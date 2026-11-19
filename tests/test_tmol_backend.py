"""Unit and integration tests for the PyRosetta-free tmol relaxation/scoring backend."""

import argparse
import ast
import contextlib
import importlib.util
import io
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
import zipfile

import numpy as np
from Bio.PDB import PDBIO, PDBParser, Polypeptide, Superimposer

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from functions import pyrosetta_utils, tmol_utils
from functions.biopython_utils import calculate_clash_score, is_aa

PDL1_PYROSETTA_ZIP = ROOT / "performance_data/pdl1_miniprotein/PDL1-pyrosetta.zip"
PDL1_SAMPLE_MEMBER = "3_PDL1_l66_s627125_mpnn8_model1.pdb"

_source = ast.parse((ROOT / "bindcraft.py").read_text())
_start = next(
    i
    for i, node in enumerate(_source.body)
    if isinstance(node, ast.Assign)
    and any(isinstance(target, ast.Name) and target.id == "parser" for target in node.targets)
)
_end = next(
    i
    for i, node in enumerate(_source.body[_start:], _start)
    if isinstance(node, ast.If)
    and isinstance(node.test, ast.Compare)
    and isinstance(node.test.left, ast.Attribute)
    and node.test.left.attr == "seed"
)
_CLI_PARSER = compile(
    ast.fix_missing_locations(ast.Module(body=_source.body[_start:_end], type_ignores=[])),
    str(ROOT / "bindcraft.py"),
    "exec",
)


def _extract_sample_pdb(dest_dir: Path) -> Path:
    dest = dest_dir / PDL1_SAMPLE_MEMBER
    with zipfile.ZipFile(PDL1_PYROSETTA_ZIP) as archive:
        dest.write_bytes(archive.read(PDL1_SAMPLE_MEMBER))
    return dest


def _set_chain_bfactors(pdb_path: Path, chain_bfactors: dict[str, float]) -> None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("bf", str(pdb_path))
    for chain in structure[0]:
        if chain.id in chain_bfactors:
            val = chain_bfactors[chain.id]
            for residue in chain:
                for atom in residue:
                    atom.set_bfactor(val)
    io_obj = PDBIO()
    io_obj.set_structure(structure)
    io_obj.save(str(pdb_path))


def _ca_rmsd(ref_pdb: Path, test_pdb: Path) -> float:
    parser = PDBParser(QUIET=True)
    ref_struct = parser.get_structure("ref", str(ref_pdb))
    test_struct = parser.get_structure("test", str(test_pdb))
    ref_ca = [
        res["CA"]
        for ch in ref_struct[0]
        for res in ch
        if is_aa(res, standard=True) and "CA" in res
    ]
    test_ca = [
        res["CA"]
        for ch in test_struct[0]
        for res in ch
        if is_aa(res, standard=True) and "CA" in res
    ]
    assert len(ref_ca) == len(test_ca) and len(ref_ca) > 0
    sup = Superimposer()
    sup.set_atoms(ref_ca, test_ca)
    return float(sup.rms)


class TmolBackendCliAndRoutingTests(unittest.TestCase):
    def parse_cli(self, *arguments):
        env = {
            "argparse": argparse,
            "DEFAULT_TRAJECTORY_REPEATS": 5,
            "DEFAULT_TRAJECTORY_ACCEPTED_DESIGNS": 1,
        }
        with patch.object(sys, "argv", ["bindcraft.py", *arguments]):
            exec(_CLI_PARSER, env)
        return env["args"]

    def test_bindcraft_cli_backend_flags_and_validation(self):
        self.assertIsNone(self.parse_cli().backend)
        for choice in ("pyrosetta", "openmm", "tmol"):
            with self.subTest(choice=choice):
                args = self.parse_cli("--backend", choice)
                self.assertEqual(args.backend, choice)

        # Launcher auto-prepends --no-pyrosetta before user arguments; both must parse cleanly.
        combined = self.parse_cli("--no-pyrosetta", "--backend", "tmol")
        self.assertTrue(combined.no_pyrosetta)
        self.assertEqual(combined.backend, "tmol")

        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as caught:
                self.parse_cli("--backend", "invalid-engine")
            self.assertEqual(caught.exception.code, 2)

    def test_resolve_scoring_backend_matrix(self):
        cases = [
            (True, None, True, "pyrosetta"),
            (True, None, False, "openmm"),
            (False, None, True, "openmm"),
            (False, None, False, "openmm"),
            (False, "tmol", True, "tmol"),
            (False, "tmol", False, "tmol"),
            (True, "openmm", True, "openmm"),
            (True, "pyrosetta", True, "pyrosetta"),
            (True, "pyrosetta", False, "openmm"),
        ]
        for use_pr, backend, pr_avail, expected in cases:
            with self.subTest(use_pr=use_pr, backend=backend, pr_avail=pr_avail):
                with patch.object(pyrosetta_utils, "PYROSETTA_AVAILABLE", pr_avail):
                    self.assertEqual(
                        pyrosetta_utils.resolve_scoring_backend(use_pyrosetta=use_pr, backend=backend),
                        expected,
                    )

        with self.assertRaises(ValueError):
            pyrosetta_utils.resolve_scoring_backend(backend="does-not-exist")

    def test_resolve_tmol_python_precedence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env_py = root / "env_python"
            env_py.write_text("#!/bin/sh\nexit 0\n")
            env_py.chmod(0o755)

            conda_prefix = root / "conda_env"
            venv_py = conda_prefix / "tmol-venv/bin/python"
            venv_py.parent.mkdir(parents=True)
            venv_py.write_text("#!/bin/sh\nexit 0\n")
            venv_py.chmod(0o755)

            with patch.dict(
                os.environ,
                {"TMOL_PYTHON": str(env_py), "CONDA_PREFIX": str(conda_prefix)},
                clear=True,
            ):
                self.assertEqual(tmol_utils.resolve_tmol_python(), str(env_py))

            with patch.dict(os.environ, {"CONDA_PREFIX": str(conda_prefix)}, clear=True):
                self.assertEqual(tmol_utils.resolve_tmol_python(), str(venv_py))

            with patch.dict(os.environ, {}, clear=True):
                self.assertEqual(tmol_utils.resolve_tmol_python(), sys.executable)

    def test_tmol_utils_import_does_not_load_torch_or_tmol(self):
        proc = subprocess.run(
            [
                sys.executable,
                "-c",
                (
                    "import sys, functions.tmol_utils as tu; "
                    "assert 'torch' not in sys.modules and 'tmol' not in sys.modules; "
                    "assert tu.is_tmol_available(); "
                    "assert 'torch' not in sys.modules and 'tmol' not in sys.modules"
                ),
            ],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
        )
        self.assertEqual(proc.returncode, 0, proc.stderr or proc.stdout)

    def test_fallbacks_when_tmol_worker_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            src_pdb = _extract_sample_pdb(tmp_path)
            out_pdb = tmp_path / "fallback_out.pdb"

            with patch.object(
                tmol_utils, "_invoke_tmol_worker", side_effect=RuntimeError("boom")
            ), patch.object(
                tmol_utils.alt, "openmm_relax_subprocess", return_value="MockOpenMM"
            ) as mock_openmm:
                res = tmol_utils.tmol_relax_subprocess(str(src_pdb), str(out_pdb), fallback_to_openmm=True)
                self.assertEqual(res, "MockOpenMM")
                mock_openmm.assert_called_once()

            sentinel_scores = {
                "binder_score": -1.0,
                "interface_dG": -10.0,
                "interface_packstat": 0.65,
                "interface_delta_unsat_hbonds": 1,
                "interface_delta_unsat_hbonds_percentage": 0.0,
            }
            with patch.object(
                tmol_utils.alt,
                "pr_alternative_score_interface",
                return_value=(dict(sentinel_scores), {"A": 1}, "B4"),
            ), patch.object(
                tmol_utils, "_invoke_tmol_worker", side_effect=RuntimeError("boom")
            ):
                scores, aa, res_ids = tmol_utils.tmol_score_interface_subprocess(
                    str(src_pdb), binder_chain="B", fallback_to_alt=True
                )
                self.assertEqual(scores["binder_score"], -1.0)
                self.assertEqual(scores["interface_dG"], -10.0)
                self.assertEqual(scores["interface_packstat"], 0.65)
                self.assertEqual(aa, {"A": 1})
                self.assertEqual(res_ids, "B4")

    def test_linux_activate_preserves_host_driver_dir_and_strips_cuda_toolkit(self):
        activate = ROOT / "installation/linux-activate.sh"
        deactivate = ROOT / "installation/linux-deactivate.sh"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            driver_dir = root / "lib64-nvidia"
            toolkit_dir = root / "cuda-12-both"
            plain_toolkit = root / "cuda-toolkit"
            for d in (driver_dir, toolkit_dir, plain_toolkit):
                d.mkdir()
            (driver_dir / "libcuda.so.1").write_text("")
            (toolkit_dir / "libcuda.so.1").write_text("")
            (toolkit_dir / "libcudart.so.12").write_text("")
            (plain_toolkit / "libcudart.so.12").write_text("")

            env_mixed = os.environ.copy()
            env_mixed["LD_LIBRARY_PATH"] = f"{driver_dir}:{toolkit_dir}:{plain_toolkit}"
            proc = subprocess.run(
                [
                    "bash",
                    "-uc",
                    (
                        f'source "{activate}" && source "{activate}" && '
                        'printf "ACT=%s\\n" "${LD_LIBRARY_PATH-<UNSET>}" && '
                        f'source "{deactivate}" && '
                        'printf "DEACT=%s\\n" "${LD_LIBRARY_PATH-<UNSET>}"'
                    ),
                ],
                env=env_mixed,
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertIn(f"ACT={driver_dir}\n", proc.stdout)
            self.assertIn(f"DEACT={driver_dir}:{toolkit_dir}:{plain_toolkit}\n", proc.stdout)

            env_toolkit_only = os.environ.copy()
            env_toolkit_only["LD_LIBRARY_PATH"] = f"{toolkit_dir}:{plain_toolkit}"
            proc_unset = subprocess.run(
                [
                    "bash",
                    "-uc",
                    f'source "{activate}" && printf "ACT=%s\\n" "${{LD_LIBRARY_PATH-<UNSET>}}"',
                ],
                env=env_toolkit_only,
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertEqual(proc_unset.stdout.strip(), "ACT=<UNSET>")


@unittest.skipUnless(tmol_utils.is_tmol_available(), "tmol is not installed in current Python")
class TmolPhysicsIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.addCleanup(tmol_utils.shutdown_tmol_worker)

    def test_persistent_worker_reuses_pid_restarts_on_crash_and_shuts_down(self):
        tmol_utils.shutdown_tmol_worker()
        self.assertIsNone(tmol_utils.get_tmol_worker_pid())

        with tempfile.TemporaryDirectory() as tmp:
            src_pdb = _extract_sample_pdb(Path(tmp))
            scores1, _, _ = tmol_utils.tmol_score_interface_subprocess(
                str(src_pdb), binder_chain="B", fallback_to_alt=False
            )
            pid1 = tmol_utils.get_tmol_worker_pid()
            self.assertIsNotNone(pid1)

            scores2, _, _ = tmol_utils.tmol_score_interface_subprocess(
                str(src_pdb), binder_chain="B", fallback_to_alt=False
            )
            pid2 = tmol_utils.get_tmol_worker_pid()
            self.assertEqual(pid2, pid1)
            self.assertEqual(scores2, scores1)

            # Kill the persistent daemon and confirm transparent auto-restart on next request.
            os.kill(pid1, signal.SIGKILL)
            time.sleep(0.1)
            scores3, _, _ = tmol_utils.tmol_score_interface_subprocess(
                str(src_pdb), binder_chain="B", fallback_to_alt=False
            )
            pid3 = tmol_utils.get_tmol_worker_pid()
            self.assertIsNotNone(pid3)
            self.assertNotEqual(pid3, pid1)
            self.assertEqual(scores3, scores1)

            tmol_utils.shutdown_tmol_worker()
            self.assertIsNone(tmol_utils.get_tmol_worker_pid())

            # One-shot opt-out must not spawn a lingering daemon.
            with patch.dict(os.environ, {"BINDCRAFT_TMOL_PERSISTENT": "0"}):
                scores_oneshot, _, _ = tmol_utils.tmol_score_interface_subprocess(
                    str(src_pdb), binder_chain="B", fallback_to_alt=False
                )
            self.assertIsNone(tmol_utils.get_tmol_worker_pid())
            self.assertEqual(scores_oneshot, scores1)
    def test_real_tmol_relax_and_score_subprocess_and_inprocess(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            src_pdb = _extract_sample_pdb(tmp_path)
            _set_chain_bfactors(src_pdb, {"A": 42.5, "B": 88.25})

            # 1. Subprocess relax via pyrosetta_utils.pr_relax(..., backend='tmol')
            rel_sub = tmp_path / "relaxed_sub.pdb"
            pyrosetta_utils.pr_relax(str(src_pdb), str(rel_sub), use_pyrosetta=False, backend="tmol")
            self.assertTrue(rel_sub.is_file() and rel_sub.stat().st_size > 0)
            self.assertLess(_ca_rmsd(src_pdb, rel_sub), 1.5)
            self.assertEqual(calculate_clash_score(str(rel_sub)), 0)

            parser = PDBParser(QUIET=True)
            rel_struct = parser.get_structure("rel", str(rel_sub))
            self.assertEqual([ch.id for ch in rel_struct[0]], ["A", "B"])
            b_a = next(res["CA"].get_bfactor() for res in rel_struct[0]["A"] if is_aa(res, standard=True))
            b_b = next(res["CA"].get_bfactor() for res in rel_struct[0]["B"] if is_aa(res, standard=True))
            self.assertAlmostEqual(b_a, 42.5, places=2)
            self.assertAlmostEqual(b_b, 88.25, places=2)

            # 2. Subprocess score via pyrosetta_utils.score_interface(..., backend='tmol')
            scores_sub, aa_sub, res_sub = pyrosetta_utils.score_interface(
                str(rel_sub), binder_chain="B", use_pyrosetta=False, backend="tmol"
            )
            self.assertLess(scores_sub["binder_score"], 0.0)
            self.assertLess(scores_sub["interface_dG"], 0.0)
            self.assertGreater(scores_sub["interface_dSASA"], 0.0)
            self.assertGreater(scores_sub["interface_nres"], 0)
            self.assertGreaterEqual(scores_sub["interface_interface_hbonds"], 0)
            self.assertAlmostEqual(
                scores_sub["interface_dG_SASA_ratio"],
                round((scores_sub["interface_dG"] / scores_sub["interface_dSASA"]) * 100.0, 2),
                delta=0.05,
            )
            self.assertAlmostEqual(
                scores_sub["interface_hbond_percentage"],
                round((scores_sub["interface_interface_hbonds"] / scores_sub["interface_nres"]) * 100.0, 2),
                delta=0.05,
            )
            self.assertEqual(scores_sub["interface_packstat"], 0.65)
            self.assertEqual(scores_sub["interface_delta_unsat_hbonds"], 1)
            self.assertEqual(scores_sub["interface_delta_unsat_hbonds_percentage"], 0.0)
            self.assertTrue(res_sub)
            self.assertEqual(sum(aa_sub.values()), scores_sub["interface_nres"])

            # 3-4. Direct in-process sanity checks when tmol is importable in the host venv
            if importlib.util.find_spec("tmol") is not None:
                scores_inp, _, _ = tmol_utils.tmol_score_interface(
                    str(rel_sub), binder_chain="B", use_subprocess=False, fallback_to_alt=False
                )
                self.assertEqual(scores_inp["binder_score"], scores_sub["binder_score"])
                self.assertEqual(scores_inp["interface_dG"], scores_sub["interface_dG"])
                self.assertEqual(
                    scores_inp["interface_interface_hbonds"],
                    scores_sub["interface_interface_hbonds"],
                )

                rel_inp = tmp_path / "relaxed_inp.pdb"
                perf = {}
                tmol_utils.tmol_relax(
                    str(src_pdb),
                    str(rel_inp),
                    use_subprocess=False,
                    fallback_to_openmm=False,
                    perf_report=perf,
                )
                self.assertTrue(rel_inp.is_file())
                self.assertFalse(perf.get("used_fallback"))
                self.assertLess(_ca_rmsd(src_pdb, rel_inp), 1.5)

    def test_multisegment_target_deconcatenation_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            src_pdb = _extract_sample_pdb(tmp_path)
            _set_chain_bfactors(src_pdb, {"A": 33.25, "B": 77.5})

            # Create a synthetic two-chain starting PDB matching chain A's 115 target residues (50 + 65).
            parser = PDBParser(QUIET=True)
            src_struct = parser.get_structure("src", str(src_pdb))
            a_residues = [res for res in src_struct[0]["A"] if Polypeptide.is_aa(res, standard=True)]
            self.assertEqual(len(a_residues), 115)

            starting_pdb = tmp_path / "starting_two_chain.pdb"
            # Split src_pdb chain A into A (50 aa) and C (65 aa) in starting_pdb
            from functions.biopython_utils import split_chain_into_subchains

            two_chain_raw = tmp_path / "two_chain_raw.pdb"
            split_chain_into_subchains(str(src_pdb), "A", [50, 65], ["A", "C"], str(two_chain_raw))
            two_struct = parser.get_structure("two", str(two_chain_raw))
            del two_struct[0]["B"]
            io_obj = PDBIO()
            io_obj.set_structure(two_struct)
            io_obj.save(str(starting_pdb))

            out_pdb = tmp_path / "deconcat_relaxed.pdb"
            perf = {}
            with patch.dict(
                os.environ,
                {"BINDCRAFT_STARTING_PDB": str(starting_pdb), "BINDCRAFT_TARGET_CHAINS": "A,C"},
            ):
                work_in, deconcat_tmp, reconcat_spec = tmol_utils._prepare_deconcat_pdb(str(src_pdb))
                self.assertIsNotNone(deconcat_tmp)
                self.assertEqual(reconcat_spec, (["C", "D"], "A"))
                if deconcat_tmp and os.path.isfile(deconcat_tmp):
                    os.remove(deconcat_tmp)

                tmol_utils.tmol_relax_subprocess(
                    str(src_pdb), str(out_pdb), fallback_to_openmm=False, perf_report=perf
                )

            self.assertFalse(perf.get("used_fallback"))
            out_struct = parser.get_structure("out", str(out_pdb))
            self.assertEqual([ch.id for ch in out_struct[0]], ["A", "B"])
            out_a_res = [r for r in out_struct[0]["A"] if is_aa(r, standard=True)]
            out_b_res = [r for r in out_struct[0]["B"] if is_aa(r, standard=True)]
            self.assertEqual(len(out_a_res), 115)
            self.assertEqual(len(out_b_res), 66)
            self.assertLess(_ca_rmsd(src_pdb, out_pdb), 1.5)
            self.assertAlmostEqual(out_a_res[0]["CA"].get_bfactor(), 33.25, places=2)
            self.assertAlmostEqual(out_b_res[0]["CA"].get_bfactor(), 77.5, places=2)

    def test_tmol_relax_batch_preserves_bfactors_and_falls_back_on_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            src1 = _extract_sample_pdb(tmp_path)
            src2 = tmp_path / "pose2.pdb"
            shutil_src = src1.read_text()
            src2.write_text(shutil_src)
            _set_chain_bfactors(src1, {"A": 21.25, "B": 64.5})
            _set_chain_bfactors(src2, {"A": 19.75, "B": 91.0})

            out1 = tmp_path / "batch_out1.pdb"
            out2 = tmp_path / "batch_out2.pdb"
            pyrosetta_utils.pr_relax_batch(
                [(str(src1), str(out1)), (str(src2), str(out2))],
                use_pyrosetta=False,
                backend="tmol",
            )
            self.assertTrue(out1.is_file() and out2.is_file())
            self.assertLess(_ca_rmsd(src1, out1), 1.5)
            self.assertLess(_ca_rmsd(src2, out2), 1.5)
            self.assertEqual(calculate_clash_score(str(out1)), 0)
            self.assertEqual(calculate_clash_score(str(out2)), 0)

            parser = PDBParser(QUIET=True)
            s1 = parser.get_structure("s1", str(out1))[0]
            s2 = parser.get_structure("s2", str(out2))[0]
            self.assertAlmostEqual(
                next(r["CA"].get_bfactor() for r in s1["B"] if is_aa(r, standard=True)),
                64.5,
                places=2,
            )
            self.assertAlmostEqual(
                next(r["CA"].get_bfactor() for r in s2["B"] if is_aa(r, standard=True)),
                91.0,
                places=2,
            )

            perf = {}
            out_fb1 = tmp_path / "fb1.pdb"
            out_fb2 = tmp_path / "fb2.pdb"
            with patch.object(tmol_utils, "_invoke_tmol_worker", side_effect=RuntimeError("batch-boom")) as mock_inv:
                with patch.object(tmol_utils, "tmol_relax", return_value={"ok": True}) as mock_single:
                    tmol_utils.tmol_relax_batch(
                        [(str(src1), str(out_fb1)), (str(src2), str(out_fb2))],
                        perf_report=perf,
                    )
            self.assertEqual(mock_inv.call_count, 1)
            self.assertEqual(mock_single.call_count, 2)
            self.assertTrue(perf.get("used_batch_fallback"))
            self.assertIn("batch-boom", perf.get("batch_fallback_reason", ""))


if __name__ == "__main__":
    unittest.main()
