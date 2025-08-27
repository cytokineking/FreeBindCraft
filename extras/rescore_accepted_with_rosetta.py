#!/usr/bin/env python3
"""
Rescore accepted designs (Ranked preferred) from a FreeBindCraft run that used --no-pyrosetta,
using PyRosetta to compute interface metrics, and evaluate whether PyRosetta-driven filters
would reject each design.

Outputs a CSV with side-by-side 'bypass_' metrics (copied from mpnn_design_stats.csv for the
specific model) and newly computed 'rosetta_' metrics, along with rosetta_filters_pass and
rosetta_failed_filters.

Usage:
  python extras/rescore_accepted_with_rosetta.py --design-path /path/to/run \
      [--filter-mode default|relaxed|design|custom] [--filters-path FILE] [--workers 4] [--output out.csv]

Notes:
  - Binder chain is assumed to be 'B'.
  - If PyRosetta initialization fails, the script aborts.
  - If a PDB fails to score, it records an error and continues.
  - Filters in 'design' mode are loaded based on the 'Filters' column value per row from mpnn_design_stats.csv;
    files are resolved in the current working directory (e.g., default_filters.json), consistent with
    extras/analyze_bindcraft_rejections.py.
"""

import os
import re
import sys
import json
import argparse
from typing import Dict, Any, List, Tuple, Optional, Set

import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

# Ensure repository root is on sys.path so that 'functions' package can be imported from extras/
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from functions.pyrosetta_utils import score_interface, PYROSETTA_AVAILABLE, pr


###############################################################################
# Helpers for binaries and PyRosetta init
###############################################################################

def ensure_binaries_executable() -> None:
    functions_dir = os.path.join(REPO_ROOT, "functions")
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
            print("PyRosetta not available in this environment; aborting.", flush=True)
        return False
    init_flags = f"-ignore_unrecognized_res -ignore_zero_occupancy -mute all -holes:dalphaball {dalphaball_path} -corrections::beta_nov16 true -relax:default_repeats 1"
    try:
        if verbose:
            print(f"Initializing PyRosetta with DAlphaBall at: {dalphaball_path}", flush=True)
        pr.init(init_flags)
        if verbose:
            print("PyRosetta initialized successfully.", flush=True)
        return True
    except Exception as e:
        if verbose:
            print(f"Failed to initialize PyRosetta: {e}", flush=True)
        return False


###############################################################################
# Filter spec parsing (adapted from extras/analyze_bindcraft_rejections.py)
###############################################################################

class FilterSpec:
    def __init__(self, json_key: str, base_key: str, threshold: float, higher: bool, is_interface_aas: bool = False, aa: Optional[str] = None):
        self.json_key = json_key
        self.base_key = base_key
        self.threshold = threshold
        self.higher = higher
        self.is_interface_aas = is_interface_aas
        self.aa = aa


def key_prefix_and_base(json_key: str) -> Tuple[Optional[str], Optional[str]]:
    """Return (prefix, base) for keys following BindCraft labeling.
    Recognizes 'Average_' and numeric prefixes '1_'..'5_'. Returns (None, None) otherwise.
    """
    if json_key.startswith("Average_"):
        return ("Average_", json_key.split("_", 1)[1])
    for n in ["1_", "2_", "3_", "4_", "5_"]:
        if json_key.startswith(n):
            return (n, json_key.split("_", 1)[1])
    return (None, None)


def extract_active_filters(filters_json: Dict[str, Any]) -> List[FilterSpec]:
    active: List[FilterSpec] = []
    for key, spec in filters_json.items():
        prefix, base = key_prefix_and_base(key)
        if prefix is None or base is None:
            continue
        if base == "InterfaceAAs":
            if not isinstance(spec, dict):
                continue
            for aa, aa_spec in spec.items():
                if not isinstance(aa_spec, dict):
                    continue
                threshold = aa_spec.get("threshold", None)
                if threshold is None:
                    continue
                higher = bool(aa_spec.get("higher", False))
                active.append(
                    FilterSpec(
                        json_key=key,
                        base_key=f"InterfaceAAs_{aa}",
                        threshold=float(threshold),
                        higher=higher,
                        is_interface_aas=True,
                        aa=aa,
                    )
                )
        else:
            if not isinstance(spec, dict):
                continue
            threshold = spec.get("threshold", None)
            if threshold is None:
                continue
            higher = bool(spec.get("higher", False))
            active.append(
                FilterSpec(json_key=key, base_key=base, threshold=float(threshold), higher=higher)
            )
    return active


###############################################################################
# Bypass (CSV) metrics extraction
###############################################################################

# Core labels as defined in functions/generic_utils.generate_dataframe_labels()
CORE_LABELS: List[str] = [
    'pLDDT', 'pTM', 'i_pTM', 'pAE', 'i_pAE', 'i_pLDDT', 'ss_pLDDT', 'Unrelaxed_Clashes',
    'Relaxed_Clashes', 'Binder_Energy_Score', 'Surface_Hydrophobicity', 'ShapeComplementarity',
    'PackStat', 'dG', 'dSASA', 'dG/dSASA', 'Interface_SASA_%', 'Interface_Hydrophobicity',
    'n_InterfaceResidues', 'n_InterfaceHbonds', 'InterfaceHbondsPercentage', 'n_InterfaceUnsatHbonds',
    'InterfaceUnsatHbondsPercentage', 'Interface_Helix%', 'Interface_BetaSheet%', 'Interface_Loop%',
    'Binder_Helix%', 'Binder_BetaSheet%', 'Binder_Loop%', 'InterfaceAAs', 'Hotspot_RMSD', 'Target_RMSD',
    'Binder_pLDDT', 'Binder_pTM', 'Binder_pAE', 'Binder_RMSD'
]


def get_bypass_metrics_for_design(row: pd.Series, model_index: int) -> Dict[str, Any]:
    """Extract per-model (N_) metrics and key per-design fields, prefixed with 'bypass_'.
    Drops the leading 'N_' from keys; uses base names to align with rosetta mapping.
    """
    result: Dict[str, Any] = {}
    # Basic per-design fields
    for base_key, col_name in [
        ("Sequence", "Sequence"),
        ("MPNN_score", "MPNN_score"),
        ("MPNN_seq_recovery", "MPNN_seq_recovery"),
        ("TargetSettings", "TargetSettings"),
        ("Filters", "Filters"),
        ("AdvancedSettings", "AdvancedSettings"),
    ]:
        if col_name in row:
            result[f"bypass_{base_key}"] = row.get(col_name)

    # Per-model metrics
    prefix = f"{model_index}_"
    for base in CORE_LABELS:
        col = prefix + base
        if col in row:
            result[f"bypass_{base}"] = row.get(col)
    return result


###############################################################################
# Rosetta rescoring and mapping
###############################################################################

_PROCESS_PYROSETTA_READY: bool = False


ROSETTA_BASE_MAP: Dict[str, str] = {
    # interface_scores keys (pyrosetta_utils) -> BindCraft base label names
    'binder_score': 'Binder_Energy_Score',
    'surface_hydrophobicity': 'Surface_Hydrophobicity',
    'interface_sc': 'ShapeComplementarity',
    'interface_packstat': 'PackStat',
    'interface_dG': 'dG',
    'interface_dSASA': 'dSASA',
    'interface_dG_SASA_ratio': 'dG/dSASA',
    'interface_fraction': 'Interface_SASA_%',
    'interface_hydrophobicity': 'Interface_Hydrophobicity',
    'interface_nres': 'n_InterfaceResidues',
    'interface_interface_hbonds': 'n_InterfaceHbonds',
    'interface_hbond_percentage': 'InterfaceHbondsPercentage',
    'interface_delta_unsat_hbonds': 'n_InterfaceUnsatHbonds',
    'interface_delta_unsat_hbonds_percentage': 'InterfaceUnsatHbondsPercentage',
}

ROSETTA_AVAILABLE_BASE_KEYS: Set[str] = set(ROSETTA_BASE_MAP.values()) | {"InterfaceAAs"}


def rescore_one_pdb(pdb_path: str, binder_chain: str, model_index: int,
                    mpnn_row: Optional[pd.Series],
                    filter_mode: str, global_filters_json: Optional[Dict[str, Any]],
                    filters_base_from_row: Optional[str],
                    dalphaball_path: Optional[str] = None,
                    verbose: bool = False) -> Dict[str, Any]:
    """Rescore a single PDB with PyRosetta, evaluate filters, and assemble an output row."""
    global _PROCESS_PYROSETTA_READY

    out: Dict[str, Any] = {}
    out["File"] = os.path.basename(pdb_path)
    out["Path"] = os.path.abspath(pdb_path)
    out["BinderChain"] = binder_chain

    # Parse Design and Model from filename for reporting
    base_no_ext = os.path.splitext(os.path.basename(pdb_path))[0]
    m = re.match(r"^(?:(?P<rank>\d+)_)?(?P<design>.+)_model(?P<model>\d+)$", base_no_ext)
    design_name = m.group("design") if m else base_no_ext
    out["Design"] = design_name
    out["Model"] = int(m.group("model")) if m and m.group("model") else model_index

    # Copy bypass per-model metrics if available
    if mpnn_row is not None:
        out.update(get_bypass_metrics_for_design(mpnn_row, model_index))
        filters_base_from_row = filters_base_from_row or str(mpnn_row.get("Filters", "")).strip()

    # PyRosetta scoring
    if not _PROCESS_PYROSETTA_READY and dalphaball_path:
        # Lazily initialize in this process if needed
        ready = try_init_pyrosetta(dalphaball_path, verbose=verbose)
        # Mark ready regardless of outcome to avoid repeated attempts; scoring will still raise if not ready
        globals()['_PROCESS_PYROSETTA_READY'] = ready

    rosetta_metrics_by_base: Dict[str, Any] = {}
    rosetta_interface_aas: Optional[Dict[str, int]] = None
    rosetta_interface_residues: Optional[str] = None
    try:
        interface_scores, interface_aa, interface_residues = score_interface(
            pdb_path, binder_chain=binder_chain, use_pyrosetta=True
        )
        # Map to base names and prefix for CSV
        for k_src, v in interface_scores.items():
            base = ROSETTA_BASE_MAP.get(k_src)
            if base is not None:
                rosetta_metrics_by_base[base] = v
                out[f"rosetta_{base}"] = v
        rosetta_interface_aas = interface_aa
        rosetta_interface_residues = interface_residues
        out["rosetta_InterfaceAAs"] = json.dumps(interface_aa)
        out["rosetta_InterfaceResidues"] = interface_residues
    except Exception as e:
        out["error"] = f"rosetta_score_failed: {e}"
        # Can't evaluate filters without rosetta metrics
        out["rosetta_filters_pass"] = None
        out["rosetta_failed_filters"] = ""
        return out

    # Resolve filters JSON
    active_filters: List[FilterSpec] = []
    filters_source = filter_mode
    filters_file_used = None
    try:
        if filter_mode in ("default", "relaxed", "custom"):
            filters_json = global_filters_json or {}
            active_filters = extract_active_filters(filters_json)
            filters_file_used = "(global)"
        elif filter_mode == "design":
            base = (filters_base_from_row or "").strip()
            if not base:
                raise RuntimeError("Design row missing Filters value; cannot resolve filter JSON in design mode")
            # Resolve JSON in CWD (primary), fallback to repo settings_filters
            cwd_candidate = os.path.join(os.getcwd(), f"{base}.json")
            repo_candidate = os.path.join(REPO_ROOT, "settings_filters", f"{base}.json")
            use_path = cwd_candidate if os.path.isfile(cwd_candidate) else repo_candidate
            if not os.path.isfile(use_path):
                raise FileNotFoundError(f"Filters JSON not found for design mode: {base}.json")
            with open(use_path, "r") as fh:
                filters_json = json.load(fh)
            active_filters = extract_active_filters(filters_json)
            filters_file_used = os.path.basename(use_path)
        else:
            raise RuntimeError(f"Unsupported filter-mode: {filter_mode}")
    except Exception as e:
        out["error"] = f"filters_load_failed: {e}"
        out["rosetta_filters_pass"] = None
        out["rosetta_failed_filters"] = ""
        out["FiltersSource"] = filters_source
        out["FiltersFileUsed"] = filters_file_used or ""
        return out

    out["FiltersSource"] = filters_source
    out["FiltersFileUsed"] = filters_file_used or ""

    # Evaluate only PyRosetta-available filters; respect model-specific prefixes
    failed: Set[str] = set()
    for fs in active_filters:
        prefix, base = key_prefix_and_base(fs.json_key)
        if prefix is None or base is None:
            continue
        # Ignore filters we cannot compute from rosetta metrics
        if fs.is_interface_aas:
            if base != "InterfaceAAs":
                continue
            # Apply only for Average_ or for matching numbered prefix
            if prefix != "Average_":
                try:
                    pref_model = int(prefix.rstrip("_"))
                except Exception:
                    continue
                if pref_model != model_index:
                    continue
            if not rosetta_interface_aas:
                continue
            val = rosetta_interface_aas.get(fs.aa) if fs.aa else None
            if val is None:
                continue
            try:
                val_f = float(val)
                thr = float(fs.threshold)
                if fs.higher:
                    if val_f < thr:
                        failed.add(f"InterfaceAAs_{fs.aa}")
                else:
                    if val_f > thr:
                        failed.add(f"InterfaceAAs_{fs.aa}")
            except Exception:
                continue
            continue

        # Non-InterfaceAAs keys
        if base not in ROSETTA_AVAILABLE_BASE_KEYS:
            continue

        # Only evaluate Average_ and the specific model's prefix
        if prefix != "Average_":
            try:
                pref_model = int(prefix.rstrip("_"))
            except Exception:
                continue
            if pref_model != model_index:
                continue

        val = rosetta_metrics_by_base.get(base)
        if val is None:
            continue
        try:
            val_f = float(val)
            thr = float(fs.threshold)
            if fs.higher:
                if val_f < thr:
                    failed.add(base)
            else:
                if val_f > thr:
                    failed.add(base)
        except Exception:
            continue

    out["rosetta_filters_pass"] = (len(failed) == 0)
    out["rosetta_failed_filters"] = ",".join(sorted(failed))
    return out


###############################################################################
# Discovery and main
###############################################################################


def discover_pdb_dir(design_path: str) -> Tuple[str, bool]:
    acc_ranked = os.path.join(design_path, "Accepted", "Ranked")
    acc = os.path.join(design_path, "Accepted")
    if os.path.isdir(acc_ranked):
        has_pdbs = any(f.lower().endswith(".pdb") for f in os.listdir(acc_ranked))
        if has_pdbs:
            return acc_ranked, True
    if os.path.isdir(acc):
        return acc, False
    raise FileNotFoundError("Neither Accepted/Ranked nor Accepted directory found under design-path")


def parse_filename_design_and_model(filename: str) -> Tuple[str, Optional[int]]:
    base = os.path.splitext(os.path.basename(filename))[0]
    m = re.match(r"^(?:(?P<rank>\d+)_)?(?P<design>.+)_model(?P<model>\d+)$", base)
    if not m:
        return base, None
    return m.group("design"), int(m.group("model"))


def load_filters_json_for_mode(filter_mode: str, filters_path: Optional[str]) -> Optional[Dict[str, Any]]:
    if filter_mode == "default":
        candidate = os.path.join(os.getcwd(), "default_filters.json")
        if not os.path.isfile(candidate):
            raise FileNotFoundError("default_filters.json not found in current working directory")
        with open(candidate, "r") as fh:
            return json.load(fh)
    if filter_mode == "relaxed":
        candidate = os.path.join(os.getcwd(), "relaxed_filters.json")
        if not os.path.isfile(candidate):
            raise FileNotFoundError("relaxed_filters.json not found in current working directory")
        with open(candidate, "r") as fh:
            return json.load(fh)
    if filter_mode == "custom":
        if not filters_path:
            raise ValueError("--filters-path is required when --filter-mode custom")
        if not os.path.isfile(filters_path):
            raise FileNotFoundError(f"Custom filters JSON not found: {filters_path}")
        with open(filters_path, "r") as fh:
            return json.load(fh)
    if filter_mode == "design":
        return None
    raise ValueError(f"Unsupported filter mode: {filter_mode}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Rescore accepted (Ranked preferred) PDBs with PyRosetta and evaluate Rosetta filters")
    parser.add_argument("--design-path", required=True, type=str, help="Run directory containing Accepted/ and mpnn_design_stats.csv")
    parser.add_argument("--filter-mode", required=False, type=str, default="design", choices=["default", "relaxed", "design", "custom"], help="Which filters to apply")
    parser.add_argument("--filters-path", required=False, type=str, default=None, help="Path to filters JSON if --filter-mode custom")
    parser.add_argument("--workers", "-w", required=False, type=int, default=1, help="Number of worker processes (default: 1)")
    parser.add_argument("--output", required=False, type=str, default=None, help="Output CSV path; defaults beside the PDB folder being analyzed")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose progress logs")
    args = parser.parse_args()

    design_path = os.path.abspath(args.design_path)
    if not os.path.isdir(design_path):
        print(f"Error: --design-path not found: {design_path}")
        sys.exit(1)

    ensure_binaries_executable()

    dalphaball_path = os.path.join(REPO_ROOT, "functions", "DAlphaBall.gcc")
    if not try_init_pyrosetta(dalphaball_path, verbose=True):
        print("Aborting: PyRosetta init failed.")
        sys.exit(1)

    pdb_dir, used_ranked = discover_pdb_dir(design_path)
    pdb_files = [os.path.join(pdb_dir, f) for f in sorted(os.listdir(pdb_dir)) if f.lower().endswith('.pdb')]
    if not pdb_files:
        print(f"No PDB files found in {pdb_dir}.")
        sys.exit(0)
    if args.verbose:
        print(f"Found {len(pdb_files)} PDBs in {pdb_dir} (Ranked={used_ranked}). Binder chain: B", flush=True)

    # Output default path
    out_csv = args.output or os.path.join(pdb_dir, "rosetta_rescore.csv")

    # Load MPNN CSV
    mpnn_csv = os.path.join(design_path, "mpnn_design_stats.csv")
    if not os.path.isfile(mpnn_csv):
        print(f"Error: mpnn_design_stats.csv not found at: {mpnn_csv}")
        sys.exit(1)
    try:
        df_mpnn = pd.read_csv(mpnn_csv)
    except Exception as e:
        print(f"Error: failed to read {mpnn_csv}: {e}")
        sys.exit(1)
    if df_mpnn.empty:
        print(f"Warning: {mpnn_csv} is empty. Proceeding but bypass metrics will be missing.")

    # Preload filters JSON for non-design modes
    global_filters_json: Optional[Dict[str, Any]] = None
    if args.filter_mode in ("default", "relaxed", "custom"):
        try:
            global_filters_json = load_filters_json_for_mode(args.filter_mode, args.filters_path)
        except Exception as e:
            print(f"Error: failed to load filters JSON for mode {args.filter_mode}: {e}")
            sys.exit(1)

    # Build lookup dict for designs
    rows_by_design: Dict[str, pd.Series] = {}
    filters_base_by_design: Dict[str, str] = {}
    if not df_mpnn.empty and "Design" in df_mpnn.columns:
        for _, row in df_mpnn.iterrows():
            d = str(row.get("Design"))
            if d:
                rows_by_design[d] = row
                filters_base_by_design[d] = str(row.get("Filters", "")).strip()

    binder_chain = "B"
    workers = max(1, int(args.workers))

    rows_out: List[Dict[str, Any]] = []
    if workers == 1:
        for idx, pdb_path in enumerate(pdb_files, start=1):
            if args.verbose:
                print("")
                print(f"=== [{idx}/{len(pdb_files)}] Rescoring: {os.path.basename(pdb_path)} ===", flush=True)
            design_name, model = parse_filename_design_and_model(pdb_path)
            if model is None:
                model = 1
            mpnn_row = rows_by_design.get(design_name)
            filters_base = filters_base_by_design.get(design_name)
            row = rescore_one_pdb(
                pdb_path, binder_chain, model, mpnn_row,
                args.filter_mode, global_filters_json, filters_base,
                dalphaball_path, args.verbose
            )
            rows_out.append(row)
    else:
        if args.verbose:
            print(f"Using {workers} workers for parallel rescoring...", flush=True)
        with ProcessPoolExecutor(max_workers=workers) as ex:
            future_to_path = {}
            for pdb_path in pdb_files:
                design_name, model = parse_filename_design_and_model(pdb_path)
                if model is None:
                    model = 1
                mpnn_row = rows_by_design.get(design_name)
                filters_base = filters_base_by_design.get(design_name)
                fut = ex.submit(
                    rescore_one_pdb,
                    pdb_path, binder_chain, model, mpnn_row,
                    args.filter_mode, global_filters_json, filters_base,
                    dalphaball_path, args.verbose
                )
                future_to_path[fut] = pdb_path
            completed = 0
            total = len(future_to_path)
            for fut in as_completed(future_to_path):
                pdb_path = future_to_path[fut]
                try:
                    row = fut.result()
                    rows_out.append(row)
                    completed += 1
                    if args.verbose:
                        print(f"Completed [{completed}/{total}]: {os.path.basename(pdb_path)}", flush=True)
                except Exception as e:
                    completed += 1
                    print(f"ERROR [{completed}/{total}] {os.path.basename(pdb_path)}: {e}", flush=True)

    df_out = pd.DataFrame(rows_out)
    df_out.to_csv(out_csv, index=False)
    print("")
    print(f"Wrote {len(df_out)} rows to {out_csv}", flush=True)


if __name__ == "__main__":
    main()


