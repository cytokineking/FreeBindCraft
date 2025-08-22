#!/usr/bin/env python3
"""
Compute PyRosetta, FreeSASA, and Biopython interface metrics for a folder of PDBs.

Modeled on compare_pyrosetta_bypass_scores.py, this script emits a CSV with
side-by-side columns for each engine (rosetta_, freesasa_, biopython_ prefixes).

Usage:
  python compare_interface_metrics_all.py --pdb-dir /path/to/pdbs --binder-chain B \
      [--output out.csv] [--recursive] [--workers 4] [--no-pyrosetta]
"""

import os
import sys
import time
import argparse
import json
import pandas as pd

from functions.pyrosetta_utils import score_interface, PYROSETTA_AVAILABLE, pr
from functions.biopython_utils import calculate_clash_score
from functions import pr_alternative_utils as alt

_PROCESS_PYROSETTA_READY = False


def ensure_binaries_executable():
    repo_root = os.path.dirname(os.path.abspath(__file__))
    functions_dir = os.path.join(repo_root, "functions")
    binaries = ["dssp", "DAlphaBall.gcc", "sc"]
    for binary in binaries:
        binary_path = os.path.join(functions_dir, binary)
        if os.path.isfile(binary_path):
            try:
                if not os.access(binary_path, os.X_OK):
                    current_mode = os.stat(binary_path).st_mode
                    os.chmod(binary_path, current_mode | 0o755)
                    print(f"Made {binary} executable", flush=True)
            except Exception as e:
                print(f"Warning: Failed to make {binary} executable: {e}", flush=True)


def try_init_pyrosetta(dalphaball_path: str, verbose: bool = True) -> bool:
    if not PYROSETTA_AVAILABLE or pr is None:
        if verbose:
            print("PyRosetta not available in this environment; skipping Rosetta metrics.", flush=True)
        return False
    init_flags = f"-ignore_unrecognized_res -ignore_zero_occupancy -mute all -holes:dalphaball {dalphaball_path} -corrections::beta_nov16 true -relax:default_repeats 1"
    try:
        if verbose:
            print(f"Attempting to initialize PyRosetta with DAlphaBall at: {dalphaball_path}", flush=True)
        pr.init(init_flags)
        if verbose:
            print("PyRosetta initialized successfully.", flush=True)
        return True
    except Exception:
        if verbose:
            print("Failed to initialize PyRosetta; continuing without Rosetta metrics.", flush=True)
        return False


def score_all_engines(pdb_path: str, binder_chain: str, enable_pyrosetta: bool, dalphaball_path: str, verbose: bool = True) -> dict:
    global _PROCESS_PYROSETTA_READY
    row: dict = {
        "File": os.path.basename(pdb_path),
        "Path": os.path.abspath(pdb_path),
        "BinderChain": binder_chain,
    }

    # Lazy per-process PyRosetta init
    if enable_pyrosetta and not _PROCESS_PYROSETTA_READY:
        _PROCESS_PYROSETTA_READY = try_init_pyrosetta(dalphaball_path, verbose=False)

    # Biopython metrics
    try:
        t0 = time.time()
        if verbose:
            print(f"[Biopython] Scoring {row['File']}...", flush=True)
        bio_scores, bio_interface_aa, bio_interface_residues = alt.pr_alternative_score_interface(
            pdb_path, binder_chain=binder_chain, sasa_engine="biopython"
        )
        for k, v in bio_scores.items():
            row[f"biopython_{k}"] = v
        row["biopython_InterfaceAAs"] = json.dumps(bio_interface_aa)
        row["biopython_InterfaceResidues"] = bio_interface_residues
        if verbose:
            print(f"[Biopython] Done {row['File']} in {time.time()-t0:.2f}s", flush=True)
    except Exception as e:
        row["biopython_error"] = str(e)

    # FreeSASA metrics
    try:
        t0 = time.time()
        if verbose:
            print(f"[FreeSASA] Scoring {row['File']}...", flush=True)
        fs_scores, fs_interface_aa, fs_interface_residues = alt.pr_alternative_score_interface(
            pdb_path, binder_chain=binder_chain, sasa_engine="freesasa"
        )
        for k, v in fs_scores.items():
            row[f"freesasa_{k}"] = v
        row["freesasa_InterfaceAAs"] = json.dumps(fs_interface_aa)
        row["freesasa_InterfaceResidues"] = fs_interface_residues
        if verbose:
            print(f"[FreeSASA] Done {row['File']} in {time.time()-t0:.2f}s", flush=True)
    except Exception as e:
        row["freesasa_error"] = str(e)

    # PyRosetta metrics
    if enable_pyrosetta:
        try:
            t0 = time.time()
            if verbose:
                print(f"[Rosetta] Scoring {row['File']}...", flush=True)
            rosetta_scores, rosetta_interface_aa, rosetta_interface_residues = score_interface(
                pdb_path, binder_chain=binder_chain, use_pyrosetta=True
            )
            for k, v in rosetta_scores.items():
                row[f"rosetta_{k}"] = v
            row["rosetta_InterfaceAAs"] = json.dumps(rosetta_interface_aa)
            row["rosetta_InterfaceResidues"] = rosetta_interface_residues
            if verbose:
                print(f"[Rosetta] Done {row['File']} in {time.time()-t0:.2f}s", flush=True)
        except Exception as e:
            row["rosetta_error"] = str(e)
    else:
        row["rosetta_unavailable"] = True

    # Common reference metrics
    try:
        t0 = time.time()
        row["Clashes_AllAtoms"] = int(calculate_clash_score(pdb_path, threshold=2.4, only_ca=False))
        if verbose:
            print(f"[Clashes] All atoms: {row['Clashes_AllAtoms']} ({time.time()-t0:.2f}s)", flush=True)
    except Exception:
        row["Clashes_AllAtoms"] = None
    try:
        t0 = time.time()
        row["Clashes_CAOnly"] = int(calculate_clash_score(pdb_path, threshold=2.5, only_ca=True))
        if verbose:
            print(f"[Clashes] CA-only: {row['Clashes_CAOnly']} ({time.time()-t0:.2f}s)", flush=True)
    except Exception:
        row["Clashes_CAOnly"] = None

    return row


def collect_pdbs(pdb_dir: str, recursive: bool) -> list:
    pdbs = []
    if recursive:
        for root, _, files in os.walk(pdb_dir):
            for f in files:
                if f.lower().endswith(".pdb"):
                    pdbs.append(os.path.join(root, f))
    else:
        for f in os.listdir(pdb_dir):
            if f.lower().endswith(".pdb"):
                pdbs.append(os.path.join(pdb_dir, f))
    return sorted(pdbs)


def main():
    parser = argparse.ArgumentParser(description="Compute PyRosetta, FreeSASA, and Biopython interface metrics across PDBs")
    parser.add_argument("--pdb-dir", type=str, default=None, help="Folder containing PDB files")
    parser.add_argument("--output", type=str, default=None, help="Output CSV path; defaults to <pdb-dir>/interface_metrics_all.csv")
    parser.add_argument("--binder-chain", type=str, default="B", help="Binder chain ID in the complex PDB (default: B)")
    parser.add_argument("--recursive", action="store_true", help="Scan folder recursively for .pdb files")
    parser.add_argument("--no-pyrosetta", action="store_true", help="Skip attempting to use PyRosetta even if installed")
    parser.add_argument("--workers", "-w", type=int, default=1, help="Number of worker processes (default: 1)")
    args = parser.parse_args()

    # Ensure binaries
    ensure_binaries_executable()

    pdb_dir = args.pdb_dir
    if not pdb_dir:
        try:
            pdb_dir = input("Enter path to folder containing PDB files: ").strip()
        except Exception:
            print("Error: --pdb-dir must be provided in non-interactive environments.")
            sys.exit(1)

    if not os.path.isdir(pdb_dir):
        print(f"Error: directory not found: {pdb_dir}")
        sys.exit(1)

    out_csv = args.output or os.path.join(pdb_dir, "interface_metrics_all.csv")

    # Attempt PyRosetta init unless disabled
    enable_pyrosetta = False
    dalphaball_path = None
    if not args.no_pyrosetta:
        repo_root = os.path.dirname(os.path.abspath(__file__))
        dalphaball_path = os.path.join(repo_root, "functions", "DAlphaBall.gcc")
        enable_pyrosetta = try_init_pyrosetta(dalphaball_path, verbose=True)
        if not enable_pyrosetta:
            print("Warning: PyRosetta unavailable or failed to initialize; proceeding without Rosetta metrics.", flush=True)

    pdb_files = collect_pdbs(pdb_dir, recursive=args.recursive)
    if not pdb_files:
        print(f"No PDB files found in {pdb_dir} (recursive={args.recursive}).", flush=True)
        sys.exit(0)
    print(f"Found {len(pdb_files)} PDBs in {pdb_dir} (recursive={args.recursive}). Binder chain: {args.binder_chain}", flush=True)

    rows = []
    total = len(pdb_files)
    workers = max(1, int(args.workers))
    if workers == 1:
        for idx, pdb_path in enumerate(pdb_files, start=1):
            print("")
            print(f"=== [{idx}/{total}] Processing: {os.path.basename(pdb_path)} ===", flush=True)
            row = score_all_engines(pdb_path, binder_chain=args.binder_chain, enable_pyrosetta=enable_pyrosetta, dalphaball_path=dalphaball_path, verbose=True)
            rows.append(row)
    else:
        print(f"Using {workers} workers for parallel scoring...", flush=True)
        from concurrent.futures import ProcessPoolExecutor, as_completed
        submitted = 0
        with ProcessPoolExecutor(max_workers=workers) as ex:
            future_to_path = {}
            for pdb_path in pdb_files:
                fut = ex.submit(score_all_engines, pdb_path, args.binder_chain, enable_pyrosetta, dalphaball_path, False)
                future_to_path[fut] = pdb_path
                submitted += 1
            completed = 0
            for fut in as_completed(future_to_path):
                pdb_path = future_to_path[fut]
                try:
                    row = fut.result()
                    rows.append(row)
                    completed += 1
                    print(f"Completed [{completed}/{total}]: {os.path.basename(pdb_path)}", flush=True)
                except Exception as e:
                    completed += 1
                    print(f"ERROR in [{completed}/{total}] {os.path.basename(pdb_path)}: {e}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    print("")
    print(f"Wrote {len(df)} rows to {out_csv}", flush=True)


if __name__ == "__main__":
    main()


