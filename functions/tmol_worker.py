#!/usr/bin/env python3
"""Standalone PyTorch/tmol worker for relaxation and interface scoring.

This module intentionally avoids importing ``functions.*`` so it can execute cleanly
inside an isolated Python/tmol environment (for example ``$CONDA_PREFIX/tmol-venv``
on Linux CUDA) without pulling in JAX, OpenMM, or ColabDesign.
"""

import argparse
import gc
import json
import math
import os
import sys
import time
import traceback
import warnings

import numpy as np
import torch
import tmol
from tmol.ops import (
    get_chain_mask,
    per_residue_block_pair_energies,
    res_mask_to_coord_mask,
)
from tmol.pack import (
    FixedAAChiSampler,
    IncludeCurrentSampler,
    create_dunbrack_sampler_from_database,
)
from tmol.pack.rotamer._conjugated_groups import _add_conjugated_group_sampler
from tmol.pose._pose_stack_builder import PoseStackBuilder
from tmol.score import create_mainchain_coordinate_constraints

warnings.filterwarnings(
    "ignore",
    message=".*Sparse invariant checks are implicitly disabled.*",
    category=UserWarning,
)


def select_tmol_device(device_spec=None):
    """Select the torch device for tmol pose stacks and score functions.

    On Apple Silicon, ``tmol`` requires poses on ``cpu`` and automatically offloads
    supported terms and simulated annealing to Metal/MPS.
    """
    if device_spec is None or str(device_spec).lower() == "auto":
        if torch.cuda.is_available() and tmol.has_cuda():
            return torch.device("cuda")
        return torch.device("cpu")

    normalized = str(device_spec).lower()
    if normalized == "cuda":
        if not (torch.cuda.is_available() and tmol.has_cuda()):
            raise RuntimeError("CUDA requested for tmol, but CUDA/tmol-cu is unavailable")
        return torch.device("cuda")
    if normalized == "cpu":
        return torch.device("cpu")
    raise ValueError(f"Unsupported tmol device specification: {device_spec!r}")


def mps_offload_active():
    """Return True when tmol Apple Silicon Metal/MPS offload is enabled."""
    try:
        from tmol.score._mps_offload import mps_offload_enabled

        return bool(mps_offload_enabled())
    except Exception:
        return False


def _atom_element_symbol(atom_name):
    """Infer element symbol from a residue-template atom name."""
    stripped = atom_name.strip()
    if not stripped:
        return ""
    if stripped[0].isdigit() or stripped[0].upper() == "H":
        return "H"
    return stripped[0].upper()


def _resolve_chain_masks(pose_stack, binder_chain="B", target_chain="A"):
    """Return boolean masks [n_poses, max_n_blocks] for target and binder chains."""
    mask_b = get_chain_mask(pose_stack, binder_chain)
    if not bool(mask_b.any()):
        raise ValueError(f"Binder chain {binder_chain!r} not found in pose stack")

    mask_a = get_chain_mask(pose_stack, target_chain)
    if not bool(mask_a.any()):
        mask_a = (pose_stack.block_type_ind >= 0) & (~mask_b)
    if not bool(mask_a.any()):
        raise ValueError(
            f"Target chain {target_chain!r} (or non-binder residues) not found in pose stack"
        )
    return mask_a, mask_b


def _extract_polar_residue_geometry(pose_stack, coords_np, block_idx):
    """Return (donor_heavy_H_pairs, acceptor_indices, residue_coords) for one block."""
    bt = pose_stack.block_type(0, block_idx)
    offset = int(pose_stack.block_coord_offset[0, block_idx])
    res_coords = coords_np[offset : offset + bt.n_atoms]

    polar_heavy = []
    hydrogens = []
    for idx, atom in enumerate(bt.atoms):
        if not np.all(np.isfinite(res_coords[idx])):
            continue
        elem = _atom_element_symbol(atom.name)
        if elem in ("N", "O", "S"):
            polar_heavy.append(idx)
        elif elem == "H":
            hydrogens.append(idx)

    donor_h_pairs = []
    for h_idx in hydrogens:
        h_coord = res_coords[h_idx]
        for heavy_idx in polar_heavy:
            if float(np.linalg.norm(h_coord - res_coords[heavy_idx])) < 1.25:
                donor_h_pairs.append((heavy_idx, h_idx))
                break

    return donor_h_pairs, polar_heavy, res_coords


def _count_atomic_hbonds_in_pair(pose_stack, coords_np, block_a, block_b):
    """Count distinct cross-residue donor-H...acceptor contacts between two blocks."""
    dh_a, acc_a, coords_a = _extract_polar_residue_geometry(pose_stack, coords_np, block_a)
    dh_b, acc_b, coords_b = _extract_polar_residue_geometry(pose_stack, coords_np, block_b)

    count = 0
    for dh_list, donor_coords, acc_list, acc_coords in (
        (dh_a, coords_a, acc_b, coords_b),
        (dh_b, coords_b, acc_a, coords_a),
    ):
        for donor_idx, h_idx in dh_list:
            h_coord = donor_coords[h_idx]
            d_coord = donor_coords[donor_idx]
            vec_hd = d_coord - h_coord
            norm_hd = float(np.linalg.norm(vec_hd))
            if norm_hd < 1e-6:
                continue
            for acc_idx in acc_list:
                a_coord = acc_coords[acc_idx]
                dist_ha = float(np.linalg.norm(h_coord - a_coord))
                dist_da = float(np.linalg.norm(d_coord - a_coord))
                if dist_ha < 2.65 and dist_da < 3.6:
                    vec_ha = a_coord - h_coord
                    cos_angle = float(np.dot(vec_hd, vec_ha)) / (
                        norm_hd * float(np.linalg.norm(vec_ha)) + 1e-8
                    )
                    if cos_angle < -0.2:
                        count += 1
    return max(count, 1)


def _count_interface_hbonds(pose_stack, sfxn, pair_terms, mask_a, mask_b, cutoff=-0.1):
    """Count interface H-bonds gated by tmol block-pair hbond energy."""
    hb_idx = sfxn.all_score_types().index(tmol.ScoreType.hbond)
    hb_matrix = pair_terms[hb_idx, 0]
    a_indices = torch.where(mask_a[0])[0].tolist()
    b_indices = torch.where(mask_b[0])[0].tolist()
    if not a_indices or not b_indices:
        return 0

    hb_sym_np = (hb_matrix + hb_matrix.T).detach().cpu().numpy()
    coords_np = pose_stack.coords[0].detach().cpu().numpy()
    total_hbonds = 0
    for ai in a_indices:
        for bi in b_indices:
            if float(hb_sym_np[ai, bi]) < cutoff:
                total_hbonds += _count_atomic_hbonds_in_pair(pose_stack, coords_np, ai, bi)
    return int(total_hbonds)


def _repack_pose(pose_stack, sfxn, pack_mask=None, seed=0):
    """Repack side chains (optionally restricted to ``pack_mask``) with Dunbrack rotamers."""
    palette = tmol.PackerPalette()
    task = tmol.PackerTask(pose_stack, palette)
    task.restrict_to_repacking()
    task.add_conformer_sampler(
        create_dunbrack_sampler_from_database(sfxn._param_db, pose_stack.device)
    )
    task.add_conformer_sampler(FixedAAChiSampler())
    task.add_conformer_sampler(IncludeCurrentSampler())
    _add_conjugated_group_sampler(task, pose_stack)
    if pack_mask is not None:
        task.disable_packing_by_block_mask(~pack_mask)
    return tmol.pack_rotamers(pose_stack, sfxn, task, verbose=False, seed=seed)


def _compute_separated_dg(
    pose_stack,
    sfxn,
    e_bound,
    pair_total,
    mask_a,
    mask_b,
    pack_separated=True,
    seed=0,
):
    """Compute E_bound - E_separated (optionally repacking former interface residues)."""
    sep_pose = pose_stack.clone()
    coord_mask_b = res_mask_to_coord_mask(sep_pose, mask_b)
    shifted_coords = sep_pose.coords.clone()
    shifted_coords[:, coord_mask_b[0], 0] += 500.0
    sep_pose.coords = shifted_coords

    if pack_separated:
        pair0 = pair_total[0]
        a_mask = mask_a[0]
        b_mask = mask_b[0]
        cross_ab = pair0[a_mask][:, b_mask].abs() + pair0[b_mask][:, a_mask].abs().T
        if_a = cross_ab.sum(dim=1) > 0.05
        if_b = cross_ab.sum(dim=0) > 0.05
        if_mask = torch.zeros_like(mask_a)
        if_mask[0, a_mask] = if_a
        if_mask[0, b_mask] = if_b
        if bool(if_mask.any()):
            sep_pose = _repack_pose(sep_pose, sfxn, pack_mask=if_mask, seed=seed)

    sep_scorer = sfxn.render_whole_pose_scoring_module(sep_pose)
    e_sep = float(sep_scorer(sep_pose.coords).sum().item())
    return float(e_bound - e_sep), e_sep


def _run_fast_relax_stack(
    pose_stack,
    sfxn,
    num_repeats=1,
    constrain_coords=True,
    constraint_weight=1.0,
    verbose=False,
):
    """Run tmol FastRelax on a single- or multi-pose PoseStack and validate finite energies."""
    t0 = time.perf_counter()
    if constrain_coords:
        pose_stack = create_mainchain_coordinate_constraints(pose_stack)
        sfxn.set_weight(tmol.ScoreType.constraint, float(constraint_weight))

    try:
        initial_energies = [
            float(v)
            for v in sfxn.render_whole_pose_scoring_module(pose_stack)(pose_stack.coords).tolist()
        ]
        if not initial_energies or not all(math.isfinite(e) for e in initial_energies):
            raise RuntimeError(f"Non-finite initial tmol energy: {initial_energies}")

        move_map = tmol.MoveMap.from_pose_stack(pose_stack)
        move_map.move_all_mc = True
        move_map.move_all_sc = True
        move_map.move_all_jumps = False
        move_map.move_all_root_jumps = False
        minimizer = tmol.KinematicRelaxMinimizer(move_map=move_map)

        relaxed_pose = tmol.fast_relax(
            pose_stack,
            sfxn,
            minimizer=minimizer,
            num_repeats=int(num_repeats),
            ramp_constraints=bool(constrain_coords),
            verbose=bool(verbose),
        )
    finally:
        if constrain_coords:
            sfxn.set_weight(tmol.ScoreType.constraint, 0.0)

    final_energies = [
        float(v)
        for v in sfxn.render_whole_pose_scoring_module(relaxed_pose)(relaxed_pose.coords).tolist()
    ]
    if not final_energies or not all(math.isfinite(e) for e in final_energies):
        raise RuntimeError(f"Non-finite final tmol energy after relax: {final_energies}")

    elapsed = time.perf_counter() - t0
    return relaxed_pose, initial_energies, final_energies, elapsed


def run_tmol_relax(
    input_pdb,
    output_pdb,
    num_repeats=1,
    constrain_coords=True,
    constraint_weight=1.0,
    device=None,
    verbose=False,
):
    """Relax ``input_pdb`` via tmol FastRelax (KinematicRelaxMinimizer) and save ``output_pdb``."""
    t0 = time.perf_counter()
    torch_device = select_tmol_device(device)
    pose_stack = tmol.pose_stack_from_file(str(input_pdb), torch_device)
    sfxn = tmol.beta2016_score_function(torch_device)

    relaxed_pose, initial_energies, final_energies, _ = _run_fast_relax_stack(
        pose_stack,
        sfxn,
        num_repeats=num_repeats,
        constrain_coords=constrain_coords,
        constraint_weight=constraint_weight,
        verbose=verbose,
    )

    tmol.file_from_pose_stack(relaxed_pose, str(output_pdb), format="pdb")
    elapsed = time.perf_counter() - t0
    return {
        "device": torch_device.type,
        "mps_offload": torch_device.type == "cpu" and mps_offload_active(),
        "minimizer_type": "kinematic",
        "num_repeats": int(num_repeats),
        "initial_energy": initial_energies[0],
        "final_energy": final_energies[0],
        "relax_seconds": elapsed,
    }


def run_tmol_relax_batch(
    input_pdbs,
    output_pdbs,
    num_repeats=1,
    constrain_coords=True,
    constraint_weight=1.0,
    device=None,
    verbose=False,
):
    """Relax a batch of PDBs together in one vectorized tmol PoseStack."""
    if len(input_pdbs) != len(output_pdbs):
        raise ValueError(
            f"Mismatched batch lengths: {len(input_pdbs)} inputs vs {len(output_pdbs)} outputs"
        )
    if not input_pdbs:
        raise ValueError("input_pdbs must not be empty")
    if len(input_pdbs) == 1:
        res = run_tmol_relax(
            input_pdbs[0],
            output_pdbs[0],
            num_repeats=num_repeats,
            constrain_coords=constrain_coords,
            constraint_weight=constraint_weight,
            device=device,
            verbose=verbose,
        )
        return {
            "device": res["device"],
            "mps_offload": res["mps_offload"],
            "minimizer_type": res["minimizer_type"],
            "num_repeats": res["num_repeats"],
            "n_poses": 1,
            "initial_energies": [res["initial_energy"]],
            "final_energies": [res["final_energy"]],
            "relax_seconds": res["relax_seconds"],
        }

    t0 = time.perf_counter()
    torch_device = select_tmol_device(device)
    single_poses = [tmol.pose_stack_from_file(str(p), torch_device) for p in input_pdbs]
    batch_pose = PoseStackBuilder.from_poses(single_poses, torch_device)
    sfxn = tmol.beta2016_score_function(torch_device)

    relaxed_batch, initial_energies, final_energies, _ = _run_fast_relax_stack(
        batch_pose,
        sfxn,
        num_repeats=num_repeats,
        constrain_coords=constrain_coords,
        constraint_weight=constraint_weight,
        verbose=verbose,
    )

    for idx, out_pdb in enumerate(output_pdbs):
        tmol.file_from_pose_stack(relaxed_batch.split(idx), str(out_pdb), format="pdb")

    elapsed = time.perf_counter() - t0
    return {
        "device": torch_device.type,
        "mps_offload": torch_device.type == "cpu" and mps_offload_active(),
        "minimizer_type": "kinematic",
        "num_repeats": int(num_repeats),
        "n_poses": len(input_pdbs),
        "initial_energies": initial_energies,
        "final_energies": final_energies,
        "relax_seconds": elapsed,
    }


def run_tmol_score(
    pdb_file,
    binder_chain="B",
    target_chain="A",
    pack_separated=True,
    hbond_energy_cutoff=-0.1,
    device=None,
    seed=0,
):
    """Score binder energy, separated interface dG, and interface H-bonds via tmol."""
    t0 = time.perf_counter()
    torch_device = select_tmol_device(device)
    pose_stack = tmol.pose_stack_from_file(str(pdb_file), torch_device)
    sfxn = tmol.beta2016_score_function(torch_device)

    e_bound = float(
        sfxn.render_whole_pose_scoring_module(pose_stack)(pose_stack.coords).sum().item()
    )
    if not math.isfinite(e_bound):
        raise RuntimeError(f"Non-finite tmol complex score: {e_bound}")

    bpm = sfxn.render_block_pair_scoring_module(pose_stack)
    pair_terms = bpm(pose_stack.coords, sum_terms=False, apply_weights=True)
    pair_total = pair_terms.sum(dim=0)

    mask_a, mask_b = _resolve_chain_masks(
        pose_stack, binder_chain=binder_chain, target_chain=target_chain
    )
    per_res = per_residue_block_pair_energies(pair_total)
    binder_score = float(per_res[mask_b].sum().item())

    interface_dG, e_sep = _compute_separated_dg(
        pose_stack,
        sfxn,
        e_bound,
        pair_total,
        mask_a,
        mask_b,
        pack_separated=bool(pack_separated),
        seed=int(seed),
    )
    interface_hbonds = _count_interface_hbonds(
        pose_stack,
        sfxn,
        pair_terms,
        mask_a,
        mask_b,
        cutoff=float(hbond_energy_cutoff),
    )

    elapsed = time.perf_counter() - t0
    return {
        "device": torch_device.type,
        "mps_offload": torch_device.type == "cpu" and mps_offload_active(),
        "bound_energy": e_bound,
        "separated_energy": e_sep,
        "binder_score": binder_score,
        "interface_dG": interface_dG,
        "interface_interface_hbonds": interface_hbonds,
        "score_seconds": elapsed,
    }


def run_tmol_check(pdb_file=None, device=None):
    """Return environment details and optional finite energy smoke test for tmol."""
    torch_device = select_tmol_device(device)
    info = {
        "tmol_version": tmol.__version__,
        "torch_version": torch.__version__,
        "tmol_has_cuda": bool(tmol.has_cuda()),
        "torch_cuda_available": bool(torch.cuda.is_available()),
        "torch_mps_available": bool(torch.backends.mps.is_available()),
        "mps_offload": torch_device.type == "cpu" and mps_offload_active(),
        "device": torch_device.type,
    }
    if pdb_file is not None:
        pose_stack = tmol.pose_stack_from_file(str(pdb_file), torch_device)
        sfxn = tmol.beta2016_score_function(torch_device)
        score = float(
            sfxn.render_whole_pose_scoring_module(pose_stack)(pose_stack.coords).sum().item()
        )
        if not math.isfinite(score):
            raise RuntimeError(f"Non-finite tmol test energy on {pdb_file}: {score}")
        info["test_pdb"] = str(pdb_file)
        info["test_energy"] = round(score, 2)
    return info


def _write_json(path, payload):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def _release_device_cache():
    """Collect unreachable tensors and empty CUDA/MPS allocator caches between jobs."""
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    elif hasattr(torch, "mps") and torch.backends.mps.is_available():
        torch.mps.empty_cache()


def _handle_serve_request(req):
    cmd = req.get("cmd")
    if cmd == "ping":
        return {"ok": True, "result": {"pid": os.getpid(), "status": "ready"}}
    if cmd == "check":
        return {
            "ok": True,
            "result": run_tmol_check(pdb_file=req.get("pdb"), device=req.get("device")),
        }
    if cmd == "relax":
        try:
            res = run_tmol_relax(
                input_pdb=req["input"],
                output_pdb=req["output"],
                num_repeats=req.get("num_repeats", 1),
                constrain_coords=req.get("constrain_coords", True),
                constraint_weight=req.get("constraint_weight", 1.0),
                device=req.get("device"),
                verbose=req.get("verbose", False),
            )
            return {"ok": True, "result": res}
        finally:
            _release_device_cache()
    if cmd == "relax_batch":
        try:
            res = run_tmol_relax_batch(
                input_pdbs=req["inputs"],
                output_pdbs=req["outputs"],
                num_repeats=req.get("num_repeats", 1),
                constrain_coords=req.get("constrain_coords", True),
                constraint_weight=req.get("constraint_weight", 1.0),
                device=req.get("device"),
                verbose=req.get("verbose", False),
            )
            return {"ok": True, "result": res}
        finally:
            _release_device_cache()
    if cmd == "score":
        try:
            res = run_tmol_score(
                pdb_file=req["input"],
                binder_chain=req.get("binder_chain", "B"),
                target_chain=req.get("target_chain", "A"),
                pack_separated=req.get("pack_separated", True),
                hbond_energy_cutoff=req.get("hbond_cutoff", -0.1),
                device=req.get("device"),
                seed=req.get("seed", 0),
            )
            return {"ok": True, "result": res}
        finally:
            _release_device_cache()
    raise ValueError(f"Unknown worker command: {cmd!r}")


def serve_loop():
    """Serve JSON-line requests on stdin and emit responses on duplicated fd 1.

    Redirects fd 1 -> fd 2 so third-party C/Python prints (e.g. atomworks, RDKit,
    FastRelax verbose logs) go to stderr and never corrupt JSON-lines IPC.
    """
    ipc_fd = os.dup(1)
    ipc_out = os.fdopen(ipc_fd, "w", encoding="utf-8", buffering=1)
    os.dup2(2, 1)
    sys.stdout = sys.stderr

    try:
        while True:
            raw_line = sys.stdin.readline()
            if not raw_line:
                break
            line = raw_line.strip()
            if not line:
                continue
            try:
                req = json.loads(line)
                if req.get("cmd") == "shutdown":
                    try:
                        ipc_out.write(
                            json.dumps({"ok": True, "result": {"shutdown": True}}) + "\n"
                        )
                        ipc_out.flush()
                    except BrokenPipeError:
                        pass
                    break
                response = _handle_serve_request(req)
            except Exception as exc:
                response = {
                    "ok": False,
                    "error": f"{type(exc).__name__}: {exc}",
                    "traceback": traceback.format_exc(),
                }
            try:
                ipc_out.write(json.dumps(response) + "\n")
                ipc_out.flush()
            except BrokenPipeError:
                break
    finally:
        try:
            ipc_out.close()
        except BrokenPipeError:
            pass
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_relax = subparsers.add_parser("relax")
    p_relax.add_argument("--input", required=True)
    p_relax.add_argument("--output", required=True)
    p_relax.add_argument("--num-repeats", type=int, default=1)
    p_relax.add_argument("--no-constrain-coords", action="store_true")
    p_relax.add_argument("--constraint-weight", type=float, default=1.0)
    p_relax.add_argument("--device", default=None)
    p_relax.add_argument("--verbose", action="store_true")
    p_relax.add_argument("--result-json", required=True)

    p_relax_batch = subparsers.add_parser("relax-batch")
    p_relax_batch.add_argument("--inputs", nargs="+", required=True)
    p_relax_batch.add_argument("--outputs", nargs="+", required=True)
    p_relax_batch.add_argument("--num-repeats", type=int, default=1)
    p_relax_batch.add_argument("--no-constrain-coords", action="store_true")
    p_relax_batch.add_argument("--constraint-weight", type=float, default=1.0)
    p_relax_batch.add_argument("--device", default=None)
    p_relax_batch.add_argument("--verbose", action="store_true")
    p_relax_batch.add_argument("--result-json", required=True)

    p_score = subparsers.add_parser("score")
    p_score.add_argument("--input", required=True)
    p_score.add_argument("--binder-chain", default="B")
    p_score.add_argument("--target-chain", default="A")
    p_score.add_argument("--no-pack-separated", action="store_true")
    p_score.add_argument("--hbond-cutoff", type=float, default=-0.1)
    p_score.add_argument("--device", default=None)
    p_score.add_argument("--seed", type=int, default=0)
    p_score.add_argument("--result-json", required=True)

    p_check = subparsers.add_parser("check")
    p_check.add_argument("--pdb", default=None)
    p_check.add_argument("--device", default=None)
    p_check.add_argument("--result-json", default=None)

    subparsers.add_parser("serve")

    args = parser.parse_args(argv)
    if args.command == "serve":
        return serve_loop()
    if args.command == "relax":
        result = run_tmol_relax(
            input_pdb=args.input,
            output_pdb=args.output,
            num_repeats=args.num_repeats,
            constrain_coords=not args.no_constrain_coords,
            constraint_weight=args.constraint_weight,
            device=args.device,
            verbose=args.verbose,
        )
        _write_json(args.result_json, result)
    elif args.command == "relax-batch":
        result = run_tmol_relax_batch(
            input_pdbs=args.inputs,
            output_pdbs=args.outputs,
            num_repeats=args.num_repeats,
            constrain_coords=not args.no_constrain_coords,
            constraint_weight=args.constraint_weight,
            device=args.device,
            verbose=args.verbose,
        )
        _write_json(args.result_json, result)
    elif args.command == "score":
        result = run_tmol_score(
            pdb_file=args.input,
            binder_chain=args.binder_chain,
            target_chain=args.target_chain,
            pack_separated=not args.no_pack_separated,
            hbond_energy_cutoff=args.hbond_cutoff,
            device=args.device,
            seed=args.seed,
        )
        _write_json(args.result_json, result)
    elif args.command == "check":
        result = run_tmol_check(pdb_file=args.pdb, device=args.device)
        if args.result_json:
            _write_json(args.result_json, result)
        else:
            print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
