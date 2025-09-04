# PDL1 Miniprotein: Head-to-Head Performance Comparison

This folder contains artifacts for a direct comparison between:

- Traditional BindCraft (with PyRosetta; CPU FastRelax and Rosetta metrics)
- FreeBindCraft (PyRosetta bypass; OpenMM relax and open-source metrics)

## Files

- pdl1_final_design_stats_pyrosetta.csv: Accepted-design statistics from the traditional BindCraft run (PyRosetta enabled).
- pdl1_final_design_stats_freebindcraft.csv: Accepted-design statistics from the FreeBindCraft run (--no-pyrosetta).
- PDL1-pyrosetta.zip: Ranked accepted designs for the PyRosetta run (PDBs).
- PDL1-freebindcraft.zip: Ranked accepted designs for the FreeBindCraft run (PDBs).
- freebindcraft_rosetta_rescore-relax.xlsx: Rescoring of FreeBindCraft-accepted designs using PyRosetta (FastRelax followed by standard interface metrics), including pass/fail against BindCraft’s traditional thresholds.

## Methods (abridged)

- Target: PDL1 example from the original BindCraft repository.
- Traditional pipeline: Traditional pipeline run with default filters and default_4stage_multimer advanced settings.
- FreeBindCraft pipeline: FreeBindCraft, run with --no-pyrosetta flag with default filters and default_4stage_multimer advanced settings.
- Hardware: Single B200-class GPU. Wall-clock times are end-to-end pipeline runtimes, including relax and filtering.
- Rescoring: FreeBindCraft accepted designs were optionally relaxed with PyRosetta FastRelax, then evaluated with Rosetta-native metrics to assess compliance with traditional filter thresholds.

## Results

| Metric | Traditional BindCraft (with PyRosetta) | FreeBindCraft (open-source) | Comparison |
|--------|----------------------------------------|------------------------------|------------|
| Accepted designs | 101 | 101 | Equivalent |
| Trajectories needed | 144 | 91 | 37% fewer |
| Runtime (B200 GPU) | 33.19 h | 12.25 h | 63% faster |
| Average ipTM | 0.785 | 0.792 | Equivalent |

## Interpretation

- Throughput: FreeBindCraft required fewer trajectories to reach the target number of accepted designs and reduced wall-clock runtime. The difference is attributable to GPU-accelerated OpenMM relaxation and trajectory count.
- Acceptance quality: Mean ipTM values are comparable, indicating similar AlphaFold confidence among accepted designs across both pipelines.
- Rosetta rescoring of FreeBindCraft designs: Approximately 58% satisfied all traditional Rosetta-based filters. Among the ~42% that did not, roughly two-thirds failed hydrogen-bond network assessments and the remainder were marginally below thresholds for surface hydrophobicity or shape complementarity. Differences primarily reflect methodological variance (e.g., area-weighted SASA and sc-rs vs Rosetta implementations) and subtle changes with Rosetta's FastRelax rather than large changes in interface geometry.

## Notes

- Metric definitions and implementation details for the PyRosetta-bypass path are documented in technical_overview/FreeBindCraft_Technical_Overview.md.
- CSV schemas reflect the pipeline variant; certain Rosetta-specific fields are absent or placeholders in FreeBindCraft outputs.
- Runtime and trajectory counts can vary with hardware, software stack, and stochastic elements (e.g., MPNN proposals).
