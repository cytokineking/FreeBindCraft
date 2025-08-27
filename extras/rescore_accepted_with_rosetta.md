## Implementation plan: rescore_accepted_with_rosetta.py

### Goal
- Rescore accepted designs from a FreeBindCraft run executed with `--no-pyrosetta` using PyRosetta, determine whether PyRosetta-driven filters would reject each design, and report which filters fail.
- Emit a CSV with:
  - Design metadata and model index
  - A boolean `rosetta_filters_pass`
  - A list of `rosetta_failed_filters`
  - Side-by-side `bypass_` metrics copied from the run CSV (single-model, not averages) and newly calculated `rosetta_` metrics.

### Scope and assumptions
- Input folder: user-provided `--design-path`. Prefer PDBs in `<design-path>/Accepted/Ranked`; fallback to `<design-path>/Accepted`.
- Binder chain is always `B`.
- PyRosetta is available and must initialize; otherwise abort.
- Do not recompute bypass metrics; read them from `<design-path>/mpnn_design_stats.csv` for the specific model in the PDB filename.
- Filters are dynamic:
  - `--filter-mode` `default|relaxed|design|custom`.
  - For `design` mode, use the same filters recorded per row in `mpnn_design_stats.csv` (`Filters` column), following the behavior in `extras/analyze_bindcraft_rejections.py`.
- Evaluate only PyRosetta-derived filters, using rescored metrics. If a PDB fails to score, record error and continue.

### CLI
- `--design-path` (required): Path to a run folder (contains `Accepted/`, `Accepted/Ranked/`, `mpnn_design_stats.csv`).
- `--filter-mode` (optional, default: `design`): one of `default`, `relaxed`, `design`, `custom`.
- `--filters-path` (optional): Path to a JSON if `--filter-mode custom`.
- `--output` (optional): Output CSV path. Default:
  - If using Ranked: `<design-path>/Accepted/Ranked/rosetta_rescore.csv`
  - Else: `<design-path>/Accepted/rosetta_rescore.csv`
- `--workers`, `-w` (optional): Number of worker processes; default `1`.
- `--verbose` (optional): Print progress.

Example:
```bash
python extras/rescore_accepted_with_rosetta.py \
  --design-path /path/to/run \
  --filter-mode design \
  --workers 4
```

### Inputs and discovery
- PDB sources:
  - Prefer: `<design-path>/Accepted/Ranked/*.pdb`
  - Fallback: `<design-path>/Accepted/*.pdb`
- CSV source for bypass metrics: `<design-path>/mpnn_design_stats.csv`
- Filters source:
  - `default` → `CWD/default_filters.json`
  - `relaxed` → `CWD/relaxed_filters.json`
  - `custom` → provided `--filters-path`
  - `design` → for each design row, use `Filters` column value `V` and load `CWD/V.json`
- DAlphaBall and helper binaries: ensure `functions/DAlphaBall.gcc` (chmod +x), also `functions/dssp` and `functions/sc` if present; same approach as in the compare scripts.

### Name parsing and alignment
- Ranked filenames: `"<rank>_<Design>_model<N>.pdb"` → extract `<Design>` and `<N>`.
- Accepted filenames: `"<Design>_model<N>.pdb"` → extract `<Design>` and `<N>`.
- Use `Design` to locate the row in `mpnn_design_stats.csv` (exact string match).
- Use model index `N` to read per-model columns: `"N_<label>"`.
- Do not use averages; only the numbered columns for the model present on disk.

### Bypass metrics to copy
- From `mpnn_design_stats.csv` row for this Design and model `N`:
  - Per-model `N_` metrics from `functions/generic_utils.generate_dataframe_labels()` `core_labels` (e.g., `ShapeComplementarity`, `dSASA`, `InterfaceAAs`, etc.).
  - Binder metrics per-model (`Binder_pLDDT`, `Binder_pTM`, `Binder_pAE`, `Binder_RMSD`).
  - Also copy `Sequence`, `MPNN_score`, `MPNN_seq_recovery`, and provenance fields (`TargetSettings`, `Filters`, `AdvancedSettings`).
- Prefix all copied values as `bypass_<column>` in the output CSV (dropping the leading `N_` for clarity).

### PyRosetta rescoring
- Initialize PyRosetta once per process with flags:
  `-ignore_unrecognized_res -ignore_zero_occupancy -mute all -holes:dalphaball <path> -corrections::beta_nov16 true -relax:default_repeats 1`
- For each PDB, call `functions.pyrosetta_utils.score_interface(pdb, binder_chain='B', use_pyrosetta=True)` and collect:
  - `binder_score`, `surface_hydrophobicity`, `interface_sc`, `interface_packstat`, `interface_dG`, `interface_dSASA`, `interface_dG_SASA_ratio`, `interface_fraction`, `interface_hydrophobicity`, `interface_nres`, `interface_interface_hbonds`, `interface_hbond_percentage`, `interface_delta_unsat_hbonds`, `interface_delta_unsat_hbonds_percentage`, `interface_AA`.
- Prefix these as `rosetta_<metric>`.

### Filter evaluation (PyRosetta-driven)
- Load active filters using logic adapted from `extras/analyze_bindcraft_rejections.py`:
  - Parse JSON → collect active keys with thresholds and `higher` flags.
  - Recognize prefixes `Average_` and `1_`..`5_`.
- Evaluate only filters whose base metric exists in the PyRosetta outputs (e.g., `dG`, `PackStat`, `ShapeComplementarity`, `n_InterfaceHbonds`, `InterfaceHbondsPercentage`, `n_InterfaceUnsatHbonds`, `InterfaceUnsatHbondsPercentage`, `Binder_Energy_Score`, `Surface_Hydrophobicity`, `dSASA`, `dG/dSASA`, `Interface_SASA_%`, `Interface_Hydrophobicity`, `n_InterfaceResidues`, `InterfaceAAs`).
- For `Average_` keys, evaluate against the single PyRosetta value (caveat documented).
- For numbered keys `N_`, evaluate only if `N` matches the PDB’s model index.
- Results:
  - `rosetta_filters_pass`: True if no failures, False otherwise.
  - `rosetta_failed_filters`: comma-separated base keys for all failures (for `InterfaceAAs`, use `InterfaceAAs_<AA>`).

### Output CSV columns (proposed)
- `File`, `Path`, `Design`, `Model`, `BinderChain`, `FiltersSource`, `FiltersFileUsed`, `rosetta_filters_pass`, `rosetta_failed_filters`, `error`
- Bypass metrics: prefixed with `bypass_` (single-model values, no averages)
- Rosetta metrics: prefixed with `rosetta_`
- Optional: `rosetta_InterfaceAAs` (JSON), `rosetta_InterfaceResidues` (string)

### Parallelism
- Use `ProcessPoolExecutor` with `--workers`.
- Initialize PyRosetta lazily per process; cache the per-process load of `mpnn_design_stats.csv` and per-design `Filters` JSONs.
- If PyRosetta init fails at the top-level, abort immediately. If per-PDB scoring fails, capture `error` and continue.

### Errors and edge cases
- Missing `Ranked` folder → fallback to `Accepted`.
- Missing `mpnn_design_stats.csv` → abort.
- Design row not found in CSV for a PDB → mark `error` and continue.
- Model number not found → mark `error` and continue.
- Filters JSON missing (for selected mode) → abort (consistent with `analyze_bindcraft_rejections.py`).
- Robust `InterfaceAAs` serialization/parsing (dict or stringified dict).

### Validation notes
- Expect close agreement between `bypass_` and `rosetta_` for overlapping metrics (SC, SASA-derived, etc.); minor differences may occur due to engine differences and average vs single-model comparisons.
- Do not recompute bypass averages; only use the single model’s values.

### Implementation steps
1. Create `extras/rescore_accepted_with_rosetta.py`.
2. Implement CLI and discovery (Ranked fallback, CSV paths).
3. Reuse helpers: ensure binaries executable, PyRosetta init (adapt from `extras/compare_*`).
4. Name parsing → `(Design, Model)` via regex.
5. Per-process load of `mpnn_design_stats.csv`; lookup row by `Design`.
6. Copy bypass metrics (`N_` columns → `bypass_` base names) plus provenance fields.
7. PyRosetta rescoring via `functions.pyrosetta_utils.score_interface` → `rosetta_` metrics.
8. Load filters and extract active filters (adapt from `extras/analyze_bindcraft_rejections.py`).
9. Evaluate only PyRosetta-available filters; aggregate failures.
10. Assemble rows; parallelize; write CSV to output path.
11. Test on runs with and without Ranked; test filter modes and workers.
12. Document CLI usage and example commands.


