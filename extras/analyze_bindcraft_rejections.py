import os
import sys
import argparse
import json
import ast
from typing import Dict, Any, List, Tuple, Set, Optional

import pandas as pd
import numpy as np


PYROSETTA_BASE_KEYS: Set[str] = {
    "Surface_Hydrophobicity",
    "ShapeComplementarity",
    "dG",
    "n_InterfaceHbonds",
    "n_InterfaceUnsatHbonds",
    "Binder_Energy_Score",
}


class FilterSpec:
    def __init__(self, json_key: str, base_key: str, threshold: float, higher: bool, is_interface_aas: bool = False, aa: Optional[str] = None):
        self.json_key = json_key
        self.base_key = base_key
        self.threshold = threshold
        self.higher = higher
        self.is_interface_aas = is_interface_aas
        self.aa = aa

    def signature(self) -> Tuple[Optional[float], Optional[bool]]:
        return (self.threshold, self.higher)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze BindCraft MPNN rejections and PyRosetta impact")
    parser.add_argument("--input-dir", required=True, type=str, help="Directory containing one or more BindCraft run folders")
    parser.add_argument("--filter-mode", required=True, type=str, choices=["default", "relaxed", "design", "custom"], help="Which filters to apply")
    parser.add_argument("--filters-path", required=False, type=str, default=None, help="Path to filters JSON if --filter-mode custom")
    parser.add_argument("--recursive", action="store_true", help="Search recursively for mpnn_design_stats.csv")
    parser.add_argument("--output-dir", required=False, type=str, default=None, help="Output directory (defaults to --input-dir)")
    return parser.parse_args()


def abort(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def ensure_file_exists(path: str, label: str) -> None:
    if not os.path.isfile(path):
        abort(f"Required {label} not found at: {path}")


def discover_mpnn_csvs(input_dir: str, recursive: bool) -> List[str]:
    results: List[str] = []
    if recursive:
        for root, _dirs, files in os.walk(input_dir):
            for f in files:
                if f == "mpnn_design_stats.csv":
                    results.append(os.path.join(root, f))
    else:
        candidate = os.path.join(input_dir, "mpnn_design_stats.csv")
        if os.path.isfile(candidate):
            results.append(candidate)
        else:
            # Non-recursive: also check immediate subdirectories only
            for name in os.listdir(input_dir):
                sub = os.path.join(input_dir, name)
                if os.path.isdir(sub):
                    candidate = os.path.join(sub, "mpnn_design_stats.csv")
                    if os.path.isfile(candidate):
                        results.append(candidate)
    return sorted(results)


def load_json(path: str) -> Dict[str, Any]:
    with open(path, "r") as fh:
        return json.load(fh)


def get_cwd_filters_path(basename_without_ext: str) -> str:
    # Filters JSON is expected in the current working directory
    return os.path.join(os.getcwd(), f"{basename_without_ext}.json")


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
        # InterfaceAAs is nested per amino acid
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


def safe_float(val: Any) -> Optional[float]:
    try:
        if pd.isna(val):
            return None
        return float(val)
    except Exception:
        return None


def parse_interface_aas_cell(cell: Any) -> Optional[Dict[str, float]]:
    # Cell is typically a string representation of a dict
    if isinstance(cell, dict):
        # Ensure float values
        out: Dict[str, float] = {}
        for k, v in cell.items():
            sv = safe_float(v)
            if sv is not None:
                out[str(k)] = sv
        return out
    if isinstance(cell, str):
        cell = cell.strip()
        if not cell:
            return None
        try:
            parsed = ast.literal_eval(cell)
            if isinstance(parsed, dict):
                out: Dict[str, float] = {}
                for k, v in parsed.items():
                    sv = safe_float(v)
                    if sv is not None:
                        out[str(k)] = sv
                return out
        except Exception:
            return None
    return None


def evaluate_row_against_filters(row: pd.Series, active_filters: List[FilterSpec]) -> Set[str]:
    failed: Set[str] = set()
    for fs in active_filters:
        if fs.is_interface_aas:
            cell = row.get(fs.json_key) if fs.json_key != "Average_InterfaceAAs" else row.get("Average_InterfaceAAs")
            # Note: for numbered keys, CSV has e.g. '1_InterfaceAAs'
            if fs.json_key != "Average_InterfaceAAs" and fs.json_key[0].isdigit():
                cell = row.get(fs.json_key)
            parsed = parse_interface_aas_cell(cell)
            if not parsed:
                continue
            value = parsed.get(fs.aa)
            sv = safe_float(value)
            if sv is None:
                continue
            if fs.higher:
                if sv < fs.threshold:
                    failed.add(fs.base_key)
            else:
                if sv > fs.threshold:
                    failed.add(fs.base_key)
        else:
            value = row.get(fs.json_key)
            sv = safe_float(value)
            if sv is None:
                continue
            if fs.higher:
                if sv < fs.threshold:
                    failed.add(fs.base_key)
            else:
                if sv > fs.threshold:
                    failed.add(fs.base_key)
    return failed


def compute_hypothetical_rank(pyro_only_row: pd.Series, final_csv: str) -> Optional[int]:
    try:
        if not os.path.isfile(final_csv):
            return None
        df_final = pd.read_csv(final_csv)
        if df_final.empty:
            return None
        if "Average_i_pTM" not in df_final.columns:
            return None
        accepted_scores = df_final["Average_i_pTM"].dropna().astype(float).tolist()
        if len(accepted_scores) == 0:
            return None
        candidate = safe_float(pyro_only_row.get("Average_i_pTM"))
        if candidate is None:
            return None
        # Rank: 1 + number of accepted with strictly greater score (descending)
        greater = sum(1 for s in accepted_scores if s > candidate)
        return greater + 1
    except Exception:
        return None


def main() -> None:
    args = parse_args()

    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir or input_dir)

    # Validate filter source(s)
    global_filters_path: Optional[str] = None
    design_mode = (args.filter_mode == "design")

    if args.filter_mode == "default":
        global_filters_path = os.path.join(os.getcwd(), "default_filters.json")
        ensure_file_exists(global_filters_path, "default filters JSON")
    elif args.filter_mode == "relaxed":
        global_filters_path = os.path.join(os.getcwd(), "relaxed_filters.json")
        ensure_file_exists(global_filters_path, "relaxed filters JSON")
    elif args.filter_mode == "custom":
        if not args.filters_path:
            abort("--filters-path is required when --filter-mode custom")
        global_filters_path = os.path.abspath(args.filters_path)
        ensure_file_exists(global_filters_path, "custom filters JSON")
    elif args.filter_mode == "design":
        # Will be validated per run by reading the CSV 'Filters' column
        pass

    # Discover mpnn files
    mpnn_csv_paths = discover_mpnn_csvs(input_dir, recursive=args.recursive)
    if not mpnn_csv_paths:
        abort("No mpnn_design_stats.csv files found. Use --recursive if needed.")

    # Pre-scan required filters JSONs in design mode
    design_filters_cache: Dict[str, Dict[str, Any]] = {}
    if design_mode:
        required_basenames: Set[str] = set()
        for mpnn_path in mpnn_csv_paths:
            try:
                df_tmp = pd.read_csv(mpnn_path, usecols=["Filters"])  # minimal scan
                if "Filters" in df_tmp.columns:
                    vals = df_tmp["Filters"].dropna().astype(str).unique().tolist()
                    required_basenames.update(vals)
            except Exception:
                # If Filters column missing, we can't resolve per-row filters
                abort(f"Filters column missing in {mpnn_path}; required for --filter-mode design")
        # Ensure each required file exists in CWD
        for base in required_basenames:
            path = get_cwd_filters_path(base)
            ensure_file_exists(path, f"filters JSON for design mode ('{base}.json')")

    # Prepare accumulators
    all_filter_specs_by_key: Dict[str, Set[Tuple[Optional[float], Optional[bool]]]] = {}
    violation_counts: Dict[str, int] = {}
    annotated_rows: List[Dict[str, Any]] = []

    total_designs = 0
    total_rejected = 0
    total_accepted = 0
    total_trajectories_all_runs = 0

    # Process each run
    for mpnn_path in mpnn_csv_paths:
        run_dir = os.path.dirname(mpnn_path)
        rel_run_dir = os.path.relpath(run_dir, input_dir)
        final_csv_path = os.path.join(run_dir, "final_design_stats.csv")

        # Load MPNN CSV
        try:
            df = pd.read_csv(mpnn_path)
        except Exception as e:
            print(f"WARNING: Skipping {mpnn_path} due to read error: {e}")
            continue
        if df.empty:
            continue

        # Count trajectories in this run (older BindCraft assumption)
        run_traj_csv = os.path.join(run_dir, "trajectory_stats.csv")
        if os.path.isfile(run_traj_csv):
            try:
                df_traj = pd.read_csv(run_traj_csv)
                total_trajectories_all_runs += len(df_traj.index)
            except Exception:
                pass

        # Determine filter JSON(s) to use for this run
        # - global if mode is default/relaxed/custom
        # - per-row if design
        global_active_filters: Optional[List[FilterSpec]] = None
        if global_filters_path:
            try:
                global_filters_json = load_json(global_filters_path)
            except Exception as e:
                abort(f"Failed to load filters JSON at {global_filters_path}: {e}")
            global_active_filters = extract_active_filters(global_filters_json)

        # Count accepted in run for context (from final CSV)
        accepted_in_run = 0
        if os.path.isfile(final_csv_path):
            try:
                df_final_for_count = pd.read_csv(final_csv_path)
                accepted_in_run = len(df_final_for_count.index)
            except Exception:
                accepted_in_run = 0

        # Track all base keys encountered for column layout
        encountered_base_keys: Set[str] = set()

        # Evaluate each row
        for _, row in df.iterrows():
            total_designs += 1

            # Resolve active filters for this row
            if design_mode:
                filters_base = str(row.get("Filters", "")).strip()
                if not filters_base:
                    abort(f"Row in {mpnn_path} lacks Filters value in design mode")
                if filters_base in design_filters_cache:
                    filters_json = design_filters_cache[filters_base]
                else:
                    filters_path = get_cwd_filters_path(filters_base)
                    ensure_file_exists(filters_path, f"filters JSON for design mode ('{filters_base}.json')")
                    filters_json = load_json(filters_path)
                    design_filters_cache[filters_base] = filters_json
                active_filters = extract_active_filters(filters_json)
            else:
                if global_active_filters is None:
                    abort("Active filters not initialized")
                active_filters = global_active_filters

            # Keep track of thresholds/higher per base key for breakdown
            for fs in active_filters:
                encountered_base_keys.add(fs.base_key)
                if fs.base_key not in all_filter_specs_by_key:
                    all_filter_specs_by_key[fs.base_key] = set()
                all_filter_specs_by_key[fs.base_key].add(fs.signature())

            failed_keys = evaluate_row_against_filters(row, active_filters)

            if failed_keys:
                total_rejected += 1
                for k in failed_keys:
                    violation_counts[k] = violation_counts.get(k, 0) + 1

                # pyrosetta flags
                pyro_involved = int(any(k in PYROSETTA_BASE_KEYS for k in failed_keys))
                pyro_only = int(pyro_involved == 1 and all(k in PYROSETTA_BASE_KEYS for k in failed_keys))

                # hypothetical rank if pyro-only and final CSV present
                hypothetical_rank: Optional[int] = None
                if pyro_only:
                    hypothetical_rank = compute_hypothetical_rank(row, final_csv_path)

                # Build annotated dict
                ann: Dict[str, Any] = {
                    "run_dir": rel_run_dir,
                    "design_name": row.get("Design"),
                    "sequence": row.get("Sequence"),
                    "Average_i_pTM": row.get("Average_i_pTM"),
                    "pyrosetta_involved": pyro_involved,
                    "pyrosetta_only": pyro_only,
                    "hypothetical_rank": hypothetical_rank,
                    "accepted_in_run": accepted_in_run,
                }
                # add flags per encountered key (1/0)
                for k in encountered_base_keys:
                    ann[k] = 1 if k in failed_keys else 0
                annotated_rows.append(ann)
            else:
                total_accepted += 1

    # If no rejections, still write empty outputs with headers
    # Build breakdown rows
    rejected_total = total_rejected
    breakdown_rows: List[Dict[str, Any]] = []
    all_keys_sorted = sorted(all_filter_specs_by_key.keys())

    for base_key in all_keys_sorted:
        sigs = all_filter_specs_by_key.get(base_key, set())
        thr: Optional[float] = None
        higher: Optional[bool] = None
        if len(sigs) == 1:
            only = next(iter(sigs))
            thr, higher = only
        num_viol = violation_counts.get(base_key, 0)
        percent = (num_viol / rejected_total) if rejected_total > 0 else 0.0
        breakdown_rows.append(
            {
                "filter_key": base_key,
                "is_pyrosetta": base_key in PYROSETTA_BASE_KEYS,
                "threshold": thr,
                "higher_is_better": higher,
                "num_violations": num_viol,
                "rejected_total": rejected_total,
                "percent_of_rejections": round(percent, 4),
                "total_trajectories": "",
                "estimated_total_mpnn_attempts": "",
            }
        )

    # Create DataFrames and write outputs
    os.makedirs(output_dir, exist_ok=True)

    estimated_total_mpnn_attempts = total_trajectories_all_runs * 20

    breakdown_df = pd.DataFrame(breakdown_rows)
    # Add totals row as metadata-like first row
    totals_row = pd.DataFrame([
        {
            "filter_key": "__TOTALS__",
            "is_pyrosetta": "",
            "threshold": "",
            "higher_is_better": "",
            "num_violations": "",
            "rejected_total": rejected_total,
            "percent_of_rejections": "",
            "total_trajectories": total_trajectories_all_runs,
            "estimated_total_mpnn_attempts": estimated_total_mpnn_attempts,
        }
    ])
    breakdown_df = pd.concat([totals_row, breakdown_df], ignore_index=True)

    breakdown_out = os.path.join(output_dir, "rejection_breakdown.csv")
    breakdown_df.to_csv(breakdown_out, index=False)

    # Annotated rejections
    annotated_out = os.path.join(output_dir, "annotated_rejections.csv")
    if annotated_rows:
        # Ensure consistent columns: fixed cols + sorted base keys
        fixed_cols = [
            "run_dir",
            "design_name",
            "sequence",
            "Average_i_pTM",
            "pyrosetta_involved",
            "pyrosetta_only",
            "hypothetical_rank",
            "accepted_in_run",
        ]
        # Collect all keys from annotated rows
        dynamic_keys: Set[str] = set()
        for r in annotated_rows:
            dynamic_keys.update(k for k in r.keys() if k not in fixed_cols)
        # Remove non-filter metadata keys from dynamic set (already in fixed)
        dynamic_keys -= set(fixed_cols)
        # Sort and place after fixed
        col_order = fixed_cols + sorted(dynamic_keys)
        ann_df = pd.DataFrame(annotated_rows)
        # Add missing columns with zeros
        for col in col_order:
            if col not in ann_df.columns:
                ann_df[col] = 0
        ann_df = ann_df[col_order]
        ann_df.to_csv(annotated_out, index=False)
    else:
        # Write empty with headers
        pd.DataFrame(
            columns=[
                "run_dir",
                "design_name",
                "sequence",
                "Average_i_pTM",
                "pyrosetta_involved",
                "pyrosetta_only",
                "hypothetical_rank",
                "accepted_in_run",
            ]
        ).to_csv(annotated_out, index=False)

    print(f"Wrote: {breakdown_out}")
    print(f"Wrote: {annotated_out}")
    print(f"Totals — designs: {total_designs}, rejected: {total_rejected}, accepted: {total_accepted}")
    print(f"Estimated total trajectories: {total_trajectories_all_runs}")
    print(f"Estimated total MPNN attempts (20x trajectories): {estimated_total_mpnn_attempts}")


if __name__ == "__main__":
    main()
