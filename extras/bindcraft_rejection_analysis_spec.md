## BindCraft Rejection Analysis Script — Product Requirements

### Goal
Analyze one or more BindCraft runs by re-applying filter criteria to MPNN design stats and reporting why designs were rejected, with a special focus on PyRosetta-dependent criteria and their impact when bypassed.

### Non-goals
- Do not read or depend on `rejected_mpnn_full_stats.csv`.
- Do not attempt to change BindCraft outputs; this is a read-only analysis tool.

### Inputs
- Required: `--input-dir` path containing one or more BindCraft run folders.
- Required: `--filter-mode` one of:
  - `default` → use `./default_filters.json`
  - `relaxed` → use `./relaxed_filters.json`
  - `design` → use the JSON named by the `Filters` column in each `mpnn_design_stats.csv` row (e.g., `default_filters` → `./default_filters.json`). If the referenced JSON is not found in the working directory, abort and alert the user.
  - `custom` → requires `--filters-path /absolute/or/relative/path.json`. If the JSON is not found, abort and alert the user.
- Optional: `--recursive` (ask user; if true, recursively search for `mpnn_design_stats.csv` under `--input-dir`).
- Optional: `--output-dir` (defaults to `--input-dir`).

All filter JSON files must be in the current working directory when using `default`, `relaxed`, or `design` modes. If any required JSON file is missing, abort with a clear error.

### Pairing rules and file discovery
- Discover `mpnn_design_stats.csv` files in `--input-dir` (recursively only if `--recursive` is set).
- For each discovered `mpnn_design_stats.csv`, pair exclusively with the `final_design_stats.csv` in the same directory. Never compare across different directories/subdirectories.
- If `final_design_stats.csv` is missing in a run directory, skip hypothetical rank computation for that run (do not abort the whole job).

### Filter evaluation scope
- Evaluate a single aggregated value per criterion (the “Average_*” form where applicable), mirroring the failure stats perspective. Ignore model-specific per-model keys (`1_*` to `5_*`) for this analysis.
- Interface amino acids are evaluated via the aggregated column (`Average_InterfaceAAs`) with per-amino-acid thresholds; report as separate keys, e.g., `InterfaceAAs_K`, `InterfaceAAs_M`.
- Skip any filter where the threshold is null or the row value is missing/NaN (do not count as pass or fail).
- Comparison direction follows the JSON:
  - `higher: true` → fail if value < threshold
  - `higher: false` → fail if value > threshold

### PyRosetta dependency classification
- Treat the following criteria as PyRosetta-dependent:
  - `Surface_Hydrophobicity`
  - `ShapeComplementarity`
  - `dG`
  - `n_InterfaceHbonds`
  - `n_InterfaceUnsatHbonds`
  - `Binder_Energy_Score`
- All other criteria are AF2-derived or replaceable via Biopython.

### Acceptance vs rejection logic
- A design is rejected if it violates at least one active criterion under the chosen filter source for that row.
- Percentages in the breakdown are calculated over the set of rejected designs only (not all designs).

### Hypothetical ranking for PyRosetta-only rejections
- For designs that would pass all non-PyRosetta criteria but fail at least one PyRosetta-dependent criterion (i.e., “pyrosetta-only” rejections):
  - Read `final_design_stats.csv` in the same run directory.
  - Extract the ranking metric used by BindCraft: `Average_i_pTM`.
  - Compute the hypothetical rank each pyrosetta-only design would have achieved by inserting its `Average_i_pTM` into the accepted designs’ `Average_i_pTM` list (sorted descending). Report the resulting 1-based rank relative to the existing accepted set for that run.
  - If `final_design_stats.csv` is missing, leave hypothetical rank fields blank for that run.

### Outputs (placed in `--output-dir`, defaulting to the root of `--input-dir`)
1) `rejection_breakdown.csv` (aggregated across all processed runs):
   - Columns:
     - `filter_key` (base name, e.g., `pLDDT`, `dG`, `InterfaceAAs_K`)
     - `is_pyrosetta` (true/false)
     - `threshold`
     - `higher_is_better` (true/false)
     - `num_violations` (count of rejected designs that violated this criterion)
     - `rejected_total` (total number of rejected designs across all processed runs)
     - `percent_of_rejections` = `num_violations / rejected_total`
   - Include a header row with overall totals: total designs, total rejected, total accepted.

2) `annotated_rejections.csv` (row-level detail for rejected designs):
   - Columns (minimum):
     - `run_dir` (relative path to the run directory containing the CSVs)
     - `design_name` (from `Design`)
     - `sequence` (from `Sequence`)
     - `Average_i_pTM` (for ranking context)
     - One column per reported filter key (1/0 flag if violated), using base names, e.g., `pLDDT`, `dG`, `InterfaceAAs_K`, etc.
     - `pyrosetta_involved` (1 if any PyRosetta-dependent filter violated)
     - `pyrosetta_only` (1 if it would pass when ignoring PyRosetta-dependent criteria)
     - `hypothetical_rank` (only if `pyrosetta_only` and `final_design_stats.csv` exists in the same directory)
     - `accepted_in_run` (count of accepted designs in that run, for context)

### CLI
- `--input-dir PATH` (required)
- `--filter-mode {default,relaxed,design,custom}` (required)
- `--filters-path PATH` (required if `--filter-mode custom`)
- `--recursive` (flag; ask user explicitly; default: false)
- `--output-dir PATH` (optional; default: `--input-dir`)

### Column usage from mpnn_design_stats.csv
- Required fields per row:
  - `Design`, `Sequence`
  - Aggregates: `Average_*` for all filters with non-null thresholds in the applicable filter JSON
  - For Interface amino acids: `Average_InterfaceAAs` (dict-like string to be parsed safely)
  - Ranking metric: `Average_i_pTM`
  - Optional (for `design` filter mode resolution): `Filters`

### Filter key normalization (reporting keys)
- For any JSON key `Average_<Metric>`, report as base key `<Metric>` (e.g., `Average_pLDDT` → `pLDDT`).
- For `Average_InterfaceAAs`, expand into `InterfaceAAs_<AA>` for each amino acid with a threshold.
- Ignore model-specific keys (`1_`–`5_`) in this analysis.

### Processing outline (per run directory)
- Load the applicable filter JSON according to `--filter-mode`:
  - `default`/`relaxed`: load the single JSON once from the working directory.
  - `custom`: load from `--filters-path`.
  - `design`: for each row, load the JSON named `<Filters>.json` from the working directory. Abort if missing.
- For each row in `mpnn_design_stats.csv`:
  - Build the active filter set: aggregated metrics with non-null thresholds from the chosen JSON.
  - For each criterion:
    - Fetch the row’s aggregated value (or parse `Average_InterfaceAAs` then get the `<AA>` count).
    - Skip if value is missing/NaN.
    - Evaluate pass/fail per threshold and direction.
  - Mark the row as rejected if any criterion failed.
  - Track which failed criteria are PyRosetta-dependent.
  - Determine `pyrosetta_involved` and `pyrosetta_only`.
- After processing the run directory:
  - If `final_design_stats.csv` exists and there are `pyrosetta_only` rows:
    - Load accepted designs’ `Average_i_pTM` from `final_design_stats.csv`.
    - Compute hypothetical rank for each `pyrosetta_only` row by inserting its `Average_i_pTM` into the accepted list.

### Aggregation and percentages
- Aggregation is across all processed runs.
- Denominator for `percent_of_rejections` is the total count of rejected designs across all runs.

### Error handling
- Abort with clear message if a required filter JSON cannot be found in the working directory (or at `--filters-path` for `custom`).
- If a run lacks `final_design_stats.csv`, continue processing and omit hypothetical ranks for that run.
- Robustly parse `Average_InterfaceAAs` with a safe parser; if parsing fails, treat as missing.

### Performance considerations
- Stream CSV processing (chunked if needed for very large files) though typical sizes should be fine in-memory.
- Cache loaded filter JSONs by filename to avoid repeated disk reads in `design` mode.

### Deliverable
- Script file (to be created after approval): `analyze_bindcraft_rejections.py` in the repository root.
- Two CSV outputs written to the chosen `--output-dir`:
  - `rejection_breakdown.csv`
  - `annotated_rejections.csv`

### Example usage
- Analyze a folder non-recursively using default filters:
```
python analyze_bindcraft_rejections.py --input-dir /path/to/runs --filter-mode default
```
- Analyze recursively using the filters used at design time (JSONs must be in working directory):
```
python analyze_bindcraft_rejections.py --input-dir /path/to/runs --filter-mode design --recursive
```
- Analyze with a custom filters JSON:
```
python analyze_bindcraft_rejections.py --input-dir /path/to/runs --filter-mode custom --filters-path ./my_filters.json
```
