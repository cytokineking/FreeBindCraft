# Extras: Analysis and Utility Scripts

This folder contains accessory scripts and documents to support BindCraft runs, analysis, and testing.

## Contents

- `analyze_bindcraft_rejections.py`
  - Purpose: Analyze MPNN rejections across one or multiple BindCraft runs. Identifies which filters are most responsible, flags if violations are PyRosetta-related, and can estimate hypothetical rank in the accepted set when rejections are due solely to PyRosetta metrics.
  - Usage (examples):
    - Default filters in current working directory:
      ```bash
      python extras/analyze_bindcraft_rejections.py --input-dir /path/to/runs --filter-mode default --recursive
      ```
    - Custom filters JSON:
      ```bash
      python extras/analyze_bindcraft_rejections.py --input-dir /path/to/runs --filter-mode custom --filters-path ./my_filters.json
      ```

- `bindcraft_rejection_analysis_spec.md`
  - Purpose: Specification and explanation of the outputs produced by the rejection analysis script.

- `BUGFIX_FILE_DESCRIPTORS.md`
  - Purpose: Notes about file descriptor-related issues and fixes for high-throughput workflows (e.g., many parallel file handles).

- `check_ulimit.sh`
  - Purpose: Helper script to display current ulimit settings and assist in diagnosing resource limits that could impact large jobs.
  - Usage:
    ```bash
    bash extras/check_ulimit.sh
    ```

- `compare_interface_metrics_all.py`
  - Purpose: Compute interface metrics across a folder of PDBs using three engines side-by-side: PyRosetta (if available), FreeSASA (open-source), and Biopython (open-source). Produces a single CSV with prefixed columns (`rosetta_`, `freesasa_`, `biopython_`).
  - Usage:
    ```bash
    python extras/compare_interface_metrics_all.py --pdb-dir /path/to/pdbs --binder-chain B --workers 4
    # Disable PyRosetta attempts:
    python extras/compare_interface_metrics_all.py --pdb-dir /path/to/pdbs --no-pyrosetta
    ```

- `compare_pyrosetta_bypass_scores.py`
  - Purpose: Compare metrics between PyRosetta and the PyRosetta-free (Biopython) path to understand differences per PDB.
  - Usage:
    ```bash
    python extras/compare_pyrosetta_bypass_scores.py --pdb-dir /path/to/pdbs --binder-chain B
    ```

- `rescore_accepted_with_rosetta.py`
  - Purpose: Rescore accepted (Ranked preferred) PDBs from a `--no-pyrosetta` run with PyRosetta, evaluate whether PyRosetta-driven filters would reject each design, and report which filters failed. Outputs side-by-side `bypass_` (from run CSV) and `rosetta_` metrics and pass/fail columns.
  - Highlights: supports dynamic filters (`default|relaxed|design|custom`), parallel workers, and optional `--fast-relax` to run PyRosetta FastRelax before scoring (recommended for consistency).
  - Usage:
    ```bash
    python extras/rescore_accepted_with_rosetta.py \
      --design-path /path/to/run \
      --filter-mode design \
      --workers 8 \
      --fast-relax
    ```

- `test_openmm_relax.py`
  - Purpose: Test harness to run OpenMM and PyRosetta relaxation on a single PDB, writing two outputs for comparison. Supports optional FASPR side-chain repacking.
  - Usage:
    ```bash
    # basic
    python extras/test_openmm_relax.py ./example/PDL1.pdb ./tmp/output_relaxed

    # enable FASPR repacking (requires functions/FASPR and functions/dun2010bbdep.bin)
    python extras/test_openmm_relax.py ./example/PDL1.pdb ./tmp/output_relaxed --faspr

    # sweep grid (writes JSON if --json-dir provided)
    python extras/test_openmm_relax.py ./example/PDL1.pdb ./tmp/output_relaxed --sweep-basic --json-dir ./tmp/relax_json --faspr
    ```

## Notes on Environment and Binaries

- The scripts assume the repository root is importable; they update `sys.path` at runtime from within `extras/`.
- Some utilities expect executables in `functions/` (e.g., `dssp`, optional `DAlphaBall.gcc`, and optional `sc` from `sc-rs`). The scripts attempt to `chmod +x` these if needed.


## References

- `sc-rs` (Shape Complementarity, open-source): [https://github.com/cytokineking/sc-rs](https://github.com/cytokineking/sc-rs)
- FreeSASA (SASA, open-source): [https://github.com/mittinatten/freesasa](https://github.com/mittinatten/freesasa)

