#!/usr/bin/env python3
"""
Compare PyRosetta vs Biopython (PyRosetta-free) interface metrics for a folder of PDBs.

For each PDB in the provided folder, computes interface metrics twice using
functions.score_interface:
  - PyRosetta mode (if available and initialized)
  - PyRosetta-free Biopython mode (Shrake-Rupley SASA with Chothia radii)

Outputs a CSV with paired columns to facilitate comparison.

Usage:
  python compare_pyrosetta_bypass_scores.py --pdb-dir /path/to/pdb/folder [--output out.csv] [--binder-chain B] [--recursive]

If --pdb-dir is omitted, the script will prompt for it.
"""

import os
import sys
import time
import argparse
import json
import pandas as pd

from functions.pyrosetta_utils import score_interface, PYROSETTA_AVAILABLE, pr
from functions.biopython_utils import calculate_clash_score


def try_init_pyrosetta(dalphaball_path: str) -> bool:
    """Attempt to initialize PyRosetta with sane defaults. Returns True on success, False otherwise."""
    if not PYROSETTA_AVAILABLE or pr is None:
        print("PyRosetta not available in this environment; using Biopython-only metrics.", flush=True)
        return False
    init_flags = f"-ignore_unrecognized_res -ignore_zero_occupancy -mute all -holes:dalphaball {dalphaball_path} -corrections::beta_nov16 true -relax:default_repeats 1"
    try:
        print(f"Attempting to initialize PyRosetta with DAlphaBall at: {dalphaball_path}", flush=True)
        pr.init(init_flags)
        print("PyRosetta initialized successfully.", flush=True)
        return True
    except Exception:
        print("Failed to initialize PyRosetta; will continue with Biopython-only metrics.", flush=True)
        return False


def score_one_pdb(pdb_path: str, binder_chain: str, enable_pyrosetta: bool) -> dict:
    """Compute both Biopython and (optionally) PyRosetta interface metrics for a single PDB."""
    row: dict = {
        "File": os.path.basename(pdb_path),
        "Path": os.path.abspath(pdb_path),
        "BinderChain": binder_chain,
    }

    # Biopython (PyRosetta-free) scores
    try:
        t0 = time.time()
        print(f"[Bypass] Scoring {row['File']} (Biopython)...", flush=True)
        bypass_scores, bypass_interface_aa, bypass_interface_residues = score_interface(
            pdb_path, binder_chain=binder_chain, use_pyrosetta=False
        )
        # Prefix for CSV clarity
        for k, v in bypass_scores.items():
            row[f"bypass_{k}"] = v
        row["bypass_InterfaceAAs"] = json.dumps(bypass_interface_aa)
        row["bypass_InterfaceResidues"] = bypass_interface_residues
        print(f"[Bypass] Done {row['File']} in {time.time()-t0:.2f}s", flush=True)
    except Exception as e:
        row["bypass_error"] = str(e)
        print(f"[Bypass] ERROR {row['File']}: {e}", flush=True)

    # PyRosetta scores (if requested and available)
    if enable_pyrosetta:
        try:
            t0 = time.time()
            print(f"[Rosetta] Scoring {row['File']} (PyRosetta)...", flush=True)
            rosetta_scores, rosetta_interface_aa, rosetta_interface_residues = score_interface(
                pdb_path, binder_chain=binder_chain, use_pyrosetta=True
            )
            for k, v in rosetta_scores.items():
                row[f"rosetta_{k}"] = v
            row["rosetta_InterfaceAAs"] = json.dumps(rosetta_interface_aa)
            row["rosetta_InterfaceResidues"] = rosetta_interface_residues
            print(f"[Rosetta] Done {row['File']} in {time.time()-t0:.2f}s", flush=True)
        except Exception as e:
            row["rosetta_error"] = str(e)
            print(f"[Rosetta] ERROR {row['File']}: {e}", flush=True)
    else:
        row["rosetta_unavailable"] = True

    # Add clash counts (common reference metrics)
    try:
        t0 = time.time()
        row["Clashes_AllAtoms"] = int(calculate_clash_score(pdb_path, threshold=2.4, only_ca=False))
        print(f"[Clashes] All atoms: {row['Clashes_AllAtoms']} ({time.time()-t0:.2f}s)", flush=True)
    except Exception:
        row["Clashes_AllAtoms"] = None
    try:
        t0 = time.time()
        row["Clashes_CAOnly"] = int(calculate_clash_score(pdb_path, threshold=2.5, only_ca=True))
        print(f"[Clashes] CA-only: {row['Clashes_CAOnly']} ({time.time()-t0:.2f}s)", flush=True)
    except Exception:
        row["Clashes_CAOnly"] = None

    return row


def collect_pdbs(pdb_dir: str, recursive: bool) -> list:
    """Find PDB files in a directory."""
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
    parser = argparse.ArgumentParser(description="Compare PyRosetta vs Biopython interface metrics across PDBs")
    parser.add_argument("--pdb-dir", type=str, default=None, help="Folder containing PDB files")
    parser.add_argument("--output", type=str, default=None, help="Output CSV path; defaults to <pdb-dir>/pyrosetta_bypass_comparison.csv")
    parser.add_argument("--binder-chain", type=str, default="B", help="Binder chain ID in the complex PDB (default: B)")
    parser.add_argument("--recursive", action="store_true", help="Scan folder recursively for .pdb files")
    parser.add_argument("--no-pyrosetta", action="store_true", help="Skip attempting to use PyRosetta even if installed")
    args = parser.parse_args()

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

    out_csv = args.output or os.path.join(pdb_dir, "pyrosetta_bypass_comparison.csv")

    # Attempt to initialize PyRosetta unless explicitly disabled
    enable_pyrosetta = False
    if not args.no_pyrosetta:
        repo_root = os.path.dirname(os.path.abspath(__file__))
        dalphaball_path = os.path.join(repo_root, "functions", "DAlphaBall.gcc")
        enable_pyrosetta = try_init_pyrosetta(dalphaball_path)
        if not enable_pyrosetta:
            print("Warning: PyRosetta unavailable or failed to initialize; proceeding with Biopython-only metrics.", flush=True)

    pdb_files = collect_pdbs(pdb_dir, recursive=args.recursive)
    if not pdb_files:
        print(f"No PDB files found in {pdb_dir} (recursive={args.recursive}).", flush=True)
        sys.exit(0)
    print(f"Found {len(pdb_files)} PDBs in {pdb_dir} (recursive={args.recursive}). Binder chain: {args.binder_chain}", flush=True)

    rows = []
    total = len(pdb_files)
    for idx, pdb_path in enumerate(pdb_files, start=1):
        print("")
        print(f"=== [{idx}/{total}] Processing: {os.path.basename(pdb_path)} ===", flush=True)
        row = score_one_pdb(pdb_path, binder_chain=args.binder_chain, enable_pyrosetta=enable_pyrosetta)
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    print("")
    print(f"Wrote {len(df)} rows to {out_csv}", flush=True)


if __name__ == "__main__":
    main()


