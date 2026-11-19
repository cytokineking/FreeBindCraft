## FreeBindCraft Technical Overview

### Purpose

FreeBindCraft provides two PyRosetta-free pathways alongside native PyRosetta so teams can run BindCraft end-to-end without a Rosetta commercial license: classic **OpenMM + FASPR + `sc-rs` + FreeSASA** (`--backend openmm` / `--no-pyrosetta`), and **PyTorch + [tmol](https://github.com/RosettaCommons/tmol)** (`--backend tmol`), which recovers Rosetta `beta2016` FastRelax, whole-binder energies, separated-chain repacked `dG`, and interface H-bonds using permissive open-source physics dependencies.

References:
- Original BindCraft: [https://github.com/martinpacesa/BindCraft](https://github.com/martinpacesa/BindCraft)
- FreeBindCraft fork: [https://github.com/cytokineking/FreeBindCraft](https://github.com/cytokineking/FreeBindCraft)
- tmol: [https://github.com/RosettaCommons/tmol](https://github.com/RosettaCommons/tmol)

### High-level engine matrix (`pyrosetta` vs. `tmol` vs. `openmm`)

| Component / Metric | `pyrosetta` (Rosetta) | `tmol` (`--backend tmol`) | Classic `openmm` (`--no-pyrosetta`) |
|---|---|---|---|
| **Relaxation** | CPU Rosetta `FastRelax` (`beta_nov16`) | **GPU `tmol.fast_relax`** (`KinematicRelaxMinimizer`, `beta2016`, batched `[model1, model2]`) | GPU **`OpenMM`** (Amber14 + OBC2 + MD shakes) + CPU `FASPR` |
| **`Binder_Energy_Score`** | Rosetta `TotalEnergyMetric` | Real `tmol` `beta2016` binder chain energy | Placeholder (`-1.0`) |
| **`dG` & `dG/dSASA`** | `InterfaceAnalyzerMover` (`pack_separated=True`) | Real `tmol` `E_bound - E_separated` (`pack_rotamers` at interface) | Placeholders (`-10.0`, `0.0`) |
| **`n_InterfaceHbonds` & `%`** | `InterfaceAnalyzerMover` | Real `tmol` `hbond` pair energy (`<-0.1`) + donor-H···acceptor geometry | Placeholders (`5`, `60.0%`) |
| **`ShapeComplementarity`** | Rosetta SC | `sc-rs` CLI | `sc-rs` CLI |
| **`dSASA`, `Interface_SASA_%`, `Surface_Hydrophobicity`** | Rosetta SASA + `LayerSelector` | FreeSASA (or Biopython SR fallback) + area-weighted SASA | FreeSASA (or Biopython SR fallback) + area-weighted SASA |
| **`n_InterfaceResidues` & `InterfaceAAs`** | Biopython 4.0 Å KD-tree | Biopython 4.0 Å KD-tree | Biopython 4.0 Å KD-tree |
| **`PackStat`** | Rosetta `InterfaceAnalyzer` | Placeholder (`0.65`) | Placeholder (`0.65`) |
| **`n_InterfaceUnsatHbonds` & `%`** | Rosetta `BuriedUnsatHbonds` (`DAlphaBall.gcc`) | Placeholders (`1`, `0.0%`) | Placeholders (`1`, `0.0%`) |
| **Alignment & RMSDs** | Rosetta `AlignChainMover` / `RMSDMetric` | Biopython Cα Superimposer / unaligned RMSD | Biopython Cα Superimposer / unaligned RMSD |

When `--backend openmm` / `--no-pyrosetta` is used, energetics and H-bond metrics are omitted in favor of placeholders; `--backend tmol` replaces those placeholders with `tmol` `beta2016` calculations while retaining `sc-rs`, FreeSASA/Biopython geometry, and placeholders only for `PackStat` and `n_InterfaceUnsatHbonds`.


### Motivation and context

BindCraft is a state-of-the-art, open pipeline for de novo miniprotein design that combines AlphaFold2/ColabDesign for structure-guided sequence generation with ProteinMPNN for sequence optimization. In the original pipeline, PyRosetta is used only in the final scoring step (relaxation plus several interface metrics).

However, PyRosetta itself requires a paid commercial license for non-academic use (see Rosetta Commons licensing FAQ: [https://rosettacommons.org/software/licensing-faq/](https://rosettacommons.org/software/licensing-faq/)). This creates friction for internal commercial drug-discovery teams. Our goal with FreeBindCraft is to remove this barrier by:
- Providing PyRosetta-free pathways through the pipeline using open-source physics/scoring dependencies (`--backend tmol` and `--backend openmm` / `--no-pyrosetta`).
- Preserving the scientific intent of key metrics with practical PyRosetta-free alternatives.
- Maintaining output quality and filter behavior where it matters most.


### Classic OpenMM bypass mode (`--no-pyrosetta` / `--backend openmm`)

When `--no-pyrosetta` (or `--backend openmm`) is selected without `--backend tmol`, FreeBindCraft switches the following components:

- Relaxation: PyRosetta FastRelax → OpenMM-based relax (GPU-accelerated)
- SASA and surface metrics: PyRosetta metrics → FreeSASA (preferred) or Biopython fallback
- Shape complementarity (SC): PyRosetta SC → `sc-rs` CLI (Rust, MIT-licensed)
- Structural alignment and RMSDs: PyRosetta → Biopython-based alignment and RMSD
- Interface energetics and some Rosetta-specific metrics: not recomputed; conservative placeholders used where needed to preserve CSV schemas and not block filters

Where possible, we validated replacements against PyRosetta-derived values. For metrics that did not materially affect accept/reject decisions in our empirical analysis, we prioritized a stable, transparent approximation over re-implementing complex internals.


## Empirical impact: what actually drives rejections

We analyzed roughly 20,000 trajectories and ~400,000 MPNN sequence variants across runs to quantify which filters determine acceptance. Highlights:

- ~88.65% of MPNN variants fail base AlphaFold metrics and are rejected before any PyRosetta scoring is involved.
- An additional ~1.63% of MPNN variants violate additional AlphaFold metrics, yielding a total of >90% of rejections being driven by AlphaFold metrics alone.
- Structural agreement checks (e.g., hotspot/binder RMSDs): ~1.66% of rejections.
- Interface amino acid composition (e.g., excessive Lys or Met residues) accounted for ~0.20 of rejections.
- Surface hydrophobicity: ~3.89% of rejections.
- Shape complementarity (SC): ~0.74% of rejections.
- Hydrogen-bond network quality (mostly “unsatisfied H-bonds”): ~3.23% of rejections.
- Interface energetic/size thresholds (e.g., dG, dSASA) were negligible drivers (only one observed rejection across these runs).

These observations guided which features we fully re-implemented, approximated, or deprioritized. The analysis workflow is codified in `extras/analyze_bindcraft_rejections.py`, which can be run on one or many result folders to reproduce breakdowns.

Example (non-recursive, default filters in current directory):
```bash
python extras/analyze_bindcraft_rejections.py \
  --input-dir /path/to/run_folder \
  --filter-mode default
```


## Implementation details

### OpenMM relaxation (replaces PyRosetta FastRelax)

Why it matters: Relaxation resolves steric clashes and poor local geometries before interface scoring. Without this step, metrics like SC and dSASA can be distorted by minor atom overlaps rather than true interface quality.

File: `functions/pr_alternative_utils.py`

- Preparation: PDBFixer removes heterogens, adds missing atoms/hydrogens.
- Force field: Amber14 with OBC2 implicit solvent.
- Backbone heavy-atom restraints: harmonic restraints ramped down across stages.
- Additional short-range repulsion: a custom nonbonded repulsive term (r^−12-like) ramped up across stages to resolve early clashes.
- Optional short MD “shakes”: applied to the first two stages to help escape poor contacts.
- Minimization: chunked L-BFGS stages with early-stopping on negligible energy improvement; accept-to-best bookkeeping across stages.
- Post-processing: align all-CA to input, copy original B-factors, and clean PDB records.
- Platform selection: prefers GPU via OpenCL; runs via a subprocess helper with retries to ensure clean driver/runtime teardown between designs.

Expected behavior: OpenMM relax is typically 2–4× faster than CPU FastRelax, and sufficient for the pipeline’s goal of removing clashes prior to interface scoring. It intentionally limits backbone drift to respect AlphaFold predictions.


### SASA, interface size, and surface hydrophobicity

Files: `functions/pr_alternative_utils.py`, `functions/biopython_utils.py`

- Why it matters: Interface size (dSASA) and the hydrophobic/hydrophilic balance at the binder surface influence stability, solubility, and specificity. Excessive surface hydrophobicity correlates with aggregation and promiscuous binding, while insufficient buried area may correlate with weak affinity.

- SASA engine: prefers FreeSASA; falls back to Biopython Shrake–Rupley. A NACCESS-style classifier config is bundled for FreeSASA.
- Interface dSASA: computed as monomer SASA minus complex SASA for binder and target, summed.
- Interface SASA fraction: total interface dSASA normalized by binder SASA within the complex.

- Surface hydrophobicity (binder): fraction of total binder SASA contributed by hydrophobic residues (A, C, F, I, L, M, P, V, W, Y). This area-weighted definition:
  - Avoids binary “exposed vs buried” thresholds (which can be brittle at boundaries).
  - Naturally down-weights partially buried residues and up-weights fully exposed residues.
  - Weighs large aromatics (W/F/Y) more than small aliphatics (A/V) in proportion to their actual surface exposure.

Comparison to original BindCraft (Rosetta layer-based exposure):
- The original approach classifies residues as “surface-exposed” or not, then computes the fraction of those residues that are hydrophobic. This binarizes exposure and treats all hydrophobics equally.
- Our area-based method produces values that correlate with the original but differ in interpretable ways: sequences rich in large aromatics or with partially buried hydrophobics may receive higher or lower scores reflecting their true SASA contributions.
- Practically, we observed close alignment between FreeSASA-derived hydrophobicity and PyRosetta-based results, with small dispersion attributable to the different exposure model (continuous area vs binary layer) and minor SASA engine differences.

Notes: FreeSASA agrees closely with PyRosetta SASA; Biopython values are slightly higher on average. The hydrophobicity definition here is intentionally area-weighted. It correlates with—though is not identical to—the original Rosetta layer-based approach and we consider it more faithful to the physical notion of “hydrophobic surface area”.

<p align="center">
  <img src="./Biopython-FreeSasa-vs-PyRosetta.png" alt="dSASA: Biopython and FreeSASA vs PyRosetta" width="900" />
</p>
<p align="center"><em>dSASA comparisons: Biopython (left) and FreeSASA (right) versus PyRosetta.</em></p>

<p align="center">
  <img src="./surface-hydrophobicity-comparison.png" alt="Surface hydrophobicity: FreeBindCraft vs PyRosetta" width="700" />
</p>
<p align="center"><em>Surface hydrophobicity: area-weighted SASA (FreeBindCraft) vs Rosetta layer‑based exposure.</em></p>

### Shape complementarity (SC)

File: `functions/pr_alternative_utils.py`

- Why it matters: SC quantifies the geometric fit of two interfaces (0–1). Higher SC often correlates with tighter, more specific binding and fewer steric voids.

- Tool: `sc-rs` CLI (MIT license) computes the Lawrence & Colman SC statistic in [0, 1].
- Integration: looks for a binary next to `functions/` (e.g., `functions/sc`), then environment variables (`SC_RS_BIN`, `SC_BIN`), then `PATH`.
- Agreement: empirical comparisons show near-identical values to PyRosetta with minor differences likely from radii/surface discretization.

<p align="center">
  <img src="./scrs-vs-pyrosetta.png" alt="Shape Complementarity: sc-rs vs PyRosetta" width="700" />
</p>
<p align="center"><em>Shape complementarity comparison: sc-rs versus PyRosetta.</em></p>


### Interface residues, composition, and structural checks

Files: `functions/biopython_utils.py`, `functions/pyrosetta_utils.py`

- Why it matters: Interface residue identity and density (and structural agreement with intended hotspots) are basic sanity checks for binding-mode fidelity. RMSDs help ensure the binder remains in the designed site.

- Interface residues: identified by KD-tree neighbor search between binder and target heavy atoms within 4.0 Å.
- Interface amino-acid composition: per-residue counts over the detected interface residues.
- Alignment and RMSDs: Biopython-based alignment (CA atoms) and unaligned RMSD for hotspot and binder comparisons.


### Metrics not recomputed without PyRosetta

Some Rosetta-native metrics are not reproduced in bypass mode due to complexity and limited impact on accept/reject decisions in our analysis:

- Energetics (e.g., `dG`, `Binder_Energy_Score`)
  - Why they matter: Interface energetics are useful for post hoc interpretation, but our analyses show they rarely act as the deciding filter thresholds within BindCraft’s default settings (~1:400,000 rejections).
  - Rationale for omission: Robust, general-purpose open re-implementations would require nontrivial modeling choices; given negligible impact on pass/fail decisions, we prefer transparency and speed.

- Buried unsatisfied hydrogen bonds (i.e., the absence of an H-bond when an H-bond donor/acceptor is present) and related network measures
  - Why they matter: In tightly-packed, hydrophobic protein cores, unsatisfied H-bond donors/acceptors are highly energetically unfavorable and often signal structural issues.
  - Why we omitted in FreeBindCraft: Protein–protein interfaces are more solvent exposed and less hydrophobic than protein cores. Implicit/explicit solvent can satisfy many potential donors/acceptors, making the “unsatisfied” call subjective and model-dependent. Burial estimation itself is noisy at interfaces.

For schema compatibility, placeholders are emitted that satisfy default thresholds but should not be interpreted as true physical values. If your workflow critically depends on these, run with PyRosetta enabled.


## Limitations and guidance

- Experimental validation directly comparing FreeBindCraft bypass mode to PyRosetta-inclusive runs has not yet been conducted; the wet‑lab impact of scoring modifications is TBD.
- SC and SASA metrics may differ slightly from Rosetta due to surface generation and radii models; these differences were small and not outcome-driving in our tests.
- The OpenMM relax focuses on clash resolution with minimal backbone drift. If your downstream workflow requires full Rosetta refinement for energetic analysis, consider enabling PyRosetta for those specific steps.


## `tmol` Engine Architecture & Parity

Files: `functions/tmol_worker.py`, `functions/tmol_utils.py`, `functions/pyrosetta_utils.py`

1. **Isolated Persistent Daemon (`tmol_worker.py` / `tmol_utils.py`):**
   - Zero top-level `torch`/`tmol` imports occur in `functions/` so JAX/Metal or JAX/CUDA initializes cleanly in the parent process.
   - On Linux CUDA, the worker runs out of `$CONDA_PREFIX/tmol-venv` (`torch==2.11.0+cu129`, `tmol==1.0.0`, `tmol-ops-cu12==1.0.0`); on Apple Silicon, it runs inside `.venv-metaljax` with PyTorch Metal/MPS offload (`device='cpu'`, `mps_offload=True`).
   - Communication uses a long-lived stdin/stdout JSON-lines daemon (`BINDCRAFT_TMOL_PERSISTENT=1`; opt out with `=0`), duplicating `fd 1` for structured IPC and redirecting `fd 1 -> fd 2` so `atomworks`/`tmol` diagnostics go to stderr. Between calls, `try/finally` resets `ScoreType.constraint` weight and runs `gc.collect()` plus `torch.cuda.empty_cache()` / `torch.mps.empty_cache()`. Unexpected worker disconnects auto-restart once transparently.
2. **Relaxation & Multi-Pose Batching (`pr_relax_batch` / `tmol_relax_batch`):**
   - Applies mainchain coordinate constraints (`constraint_weight=1.0`), runs `tmol.fast_relax(..., minimizer=tmol.KinematicRelaxMinimizer, num_repeats=1, ramp_constraints=True)`, de-concatenates/re-concatenates multi-chain targets if `BINDCRAFT_STARTING_PDB` / `BINDCRAFT_TARGET_CHAINS` are set, aligns Cα back to input, and restores per-residue B-factors.
   - In `predict_binder_complex` (`functions/colabdesign_utils.py`), AF2-passing validation models (`[model1, model2]`) are stacked via `PoseStackBuilder.from_poses` and relaxed together in a single `PoseStack` before splitting back via `PoseStack.split`. On failure, batch execution falls back cleanly to per-structure `tmol_relax`, and individual `tmol_relax` falls back to `openmm_relax`.
3. **Interface Scoring Parity on `PDL1-pyrosetta.zip` (Scoring-only, Same Relaxed Structures):**

| Complex | `pyrosetta` `E` / `tmol` `E` | `pyrosetta` `dG` / `tmol` `dG` | `pyrosetta` `dSASA` / `tmol` `dSASA` | `hbonds` (`PR` / `tmol`) | `SC` (`PR` / `tmol`) | `nres` (`PR` / `tmol`) |
|---|---:|---:|---:|---:|---:|---:|
| `1_PDL1_l123_s784666_mpnn3_model2` | -354.10 / -285.42 | -60.09 / -59.89 *(M5)* · -62.13 / -59.89 *(CUDA)* | 2360.20 / 2326.34 | 10 / 10 | 0.59 / 0.58 | 25 / 25 |
| `2_PDL1_l66_s627125_mpnn19_model2` | -179.04 / -148.71 | -55.42 / -55.87 | 1811.86 / 1782.70 | 7 / 7 | 0.73 / 0.71 | 19 / 19 |
| `3_PDL1_l66_s627125_mpnn8_model1` | -174.35 / -152.81 | -48.65 / -48.77 *(M5)* · -48.50 / -48.77 *(CUDA)* | 1742.82 / 1721.09 | 5 / 5 | 0.73 / 0.70 | 18 / 18 |
| `4_PDL1_l123_s784666_mpnn2_model2` | -337.12 / -264.14 | -69.51 / -67.58 *(M5)* · -67.74 / -67.59 *(CUDA)* | 2426.54 / 2375.21 | 12 / 12 | 0.67 / 0.67 | 29 / 29 |
| `5_PDL1_l139_s836274_mpnn19_model1` | -410.01 / -350.81 | -60.47 / -56.81 *(M5)* · -60.27 / -56.81 *(CUDA)* | 2245.44 / 2210.07 | 11 / 11 | 0.66 / 0.64 | 22 / 22 |

`dG`, `dSASA`, `n_InterfaceResidues`, `n_InterfaceHbonds`, and `SC` align closely across both platforms; `Binder_Energy_Score` is ~12–22% less negative in `tmol` due to `beta2016` reference/terminal handling (`ALA:nterm` / `SER:nterm` block resolution differences).

## Performance & Pipeline Amortization Analysis (NVIDIA RTX PRO 6000 vs. Apple M5 Max)

### 4-Way Relax + Score Benchmark (`1_PDL1..5_PDL1`, H-stripped inputs)

- **Clean cold subprocess startups:** PyRosetta `import + init` = `0.53 s` (M5 Max) / `0.80 s` (RTX PRO 6000 host) + lazy DB/rotamer load on call #1; `tmol` CUDA subprocess = `6.56 s` (`torch=0.72 s`, `tmol=0.01 s`, `sfxn_1st=2.77 s`, `io_1st=2.70 s`, `eval_1st=0.18 s`).
- **Part A steady-state `score_interface` (`calls #2–5` mean):** `tmol_persistent` **`1.04 s`** CUDA (`0.90–1.30 s`) / **`0.92 s`** M5 Max (`0.80–1.05 s`) vs. PyRosetta **`2.08 s`** (CUDA host) / **`1.74 s`** (M5 Max).

| Backend | RTX PRO 6000 Call #1 Cold (Relax / Score) | RTX PRO 6000 Mean `#2–5` (Relax / Score) | RTX PRO 6000 All-5 Mean | M5 Max Call #1 Cold (Relax / Score) | M5 Max Mean `#2–5` (Relax / Score) | M5 Max All-5 Mean |
|---|---:|---:|---:|---:|---:|---:|
| **`tmol_persistent`** | 10.56 s (9.49 / 1.07) | **3.47 s** (2.38 / 1.09) | **4.89 s** | 16.21 s (15.27 / 0.93) | **10.25 s** (9.37 / 0.88) | **11.44 s** |
| **`openmm`** (`CUDA`/`OpenCL` + `FASPR`) | 12.05 s (11.55 / 0.51) | **11.07 s** (10.65 / 0.42) | **11.27 s** | 15.18 s (14.84 / 0.34) | **13.83 s** (13.53 / 0.29) | **14.10 s** |
| **`pyrosetta`** (CPU + `DAlphaBall`) | 20.69 s (18.16 / 2.53) | **18.24 s** (16.15 / 2.09) | **18.73 s** | 14.31 s (12.22 / 2.09) | **11.70 s** (9.92 / 1.78) | **12.22 s** |
| **`tmol_oneshot`** (`BINDCRAFT_TMOL_PERSISTENT=0`) | 18.76 s (10.27 / 8.49) | **18.51 s** (10.34 / 8.17) | **18.56 s** | 20.08 s (15.20 / 4.88) | **18.66 s** (13.91 / 4.75) | **18.94 s** |

### Why AF2 Validation Stays Sequential-with-Early-Exit (and What Actually Pays Off)

Measured on `3_PDL1_l66_s627125_mpnn8_model1` (`binder_len=66`, `num_recycles=3`, `models=[0, 1]`):

| Stage | RTX PRO 6000 Blackwell (`CUDA`) | Apple M5 Max (`MetalJAX`) |
|---|---:|---:|
| **MPNN load + prep** | 0.38 s + 0.02 s | 0.06 s + 0.02 s |
| **MPNN sample (`batch=4` / `batch=20`, 1st → 2nd)** | `4.20 → 0.02 s` / `4.36 → 0.03 s` | `1.25 → 0.15 s` (`batch=4`) |
| **AF2 validation load+prep (once/trajectory)** | `3.13 s` (complex 1.83 s + monomer 1.30 s) | `4.05 s` (complex 2.19 s + monomer 1.87 s) |
| **Variant 1 (`COLD_JIT`, complex + monomer)** | **31.86 s** (complex `16.70 + 0.41 s`, monomer `14.60 + 0.16 s`) | **6.77 s** (complex `3.14 + 1.77 s`, monomer `1.42 + 0.45 s`) |
| **Variants 2–4 (`WARM_JIT`, per variant)** | **1.12 s** (complex `0.80–0.82 s`, monomer `0.32 s`) | **4.43–4.49 s** (complex `3.49–3.51 s`, monomer `0.93–0.98 s`) |
| **Cross-trajectory rebuild (`clear_mem()`, baseline before cache)** | **32.25 s** (same len 66) · **35.14 s** (new len 74) | **4.47 s** (complex rebuild `mk+prep 1.92 s + 2.55 s`) |

These measurements motivate four concrete design choices in FreeBindCraft:
1. **Do not batch-fold all 20 MPNN variants upfront:** Inside a trajectory, AF2 weights and JIT graphs are already constructed once outside the `mpnn_sequences` loop (`bindcraft.py`). Keeping AF2 sequential lets `predict_binder_complex` reject weak candidates on `model1` (`0.40 s` CUDA / `~1.75 s` M5 Max) before running `model2`, relaxation (`3.47 s` CUDA / `10.25 s` M5 Max), scoring, or monomer prediction (`0.32 s` / `0.95 s`), and preserves early exit at `max_mpnn_sequences=2`.
2. **Cache validation models across trajectories + instantiate monomers lazily (`ValidationModelRuntime`):** Without caching, recreating AF2 validation models per trajectory triggered **`~29 s`** of repeated CUDA XLA compilation (`32.25 s` total vs. `1.12 s` warm validation) even at an unchanged binder length (`len=66`). Enabling `--trajectory-cache` now caches both `complex_prediction_model` and `binder_prediction_model` per `(target, length)` block, while constructing `binder_prediction_model` lazily only when a sequence passes complex AF2 filters (even when `--trajectory-cache` is off).
3. **Batch `[model1, model2]` inside `tmol.PoseStack` (`pr_relax_batch`):** Warm `PoseStack` relaxation+repack scaling on `3_PDL1`:
   - **RTX PRO 6000 CUDA:** `n=1` `1.81 s/pose` (`207 MB` VRAM) → `n=2` **`1.31 s/pose` (`1.38×`, `290 MB`)** → `n=4` **`0.85 s/pose` (`2.14×`, `353 MB`)** → `n=8` **`0.73 s/pose` (`2.46×`, `492 MB`)**.
   - **Apple M5 Max (`MPS` offload):** `n=1` `7.93 s/pose` → `n=2` **`5.50 s/pose` (`1.44×`)** → `n=4` **`4.07 s/pose` (`1.95×`)**.
4. **CUDA VRAM Coexistence:** With AF2 resident, JAX defaulted to `73,610 / 97,887 MiB` (~75%) while persistent `tmol` occupied `824–912 MiB` (CPU RSS `3895 MB`). To avoid OOM on smaller 16–24 GB GPUs when `XLA_PYTHON_CLIENT_MEM_FRACTION` is unset, `bindcraft.py --backend tmol` sets a `0.70` default on CUDA before JAX initializes GPU memory.

## Performance and operational notes

- The relaxation runs in a fresh subprocess per design to avoid long-run GPU driver/runtime state accumulation.
- File-descriptor and DSSP caching safeguards are in place to support high-throughput runs (see `extras/BUGFIX_FILE_DESCRIPTORS.md`).
- Overall pipeline throughput improves as more MPNN candidates are accepted given H-bond network evaluation is omitted in bypass mode (exact speedups TBD).


## Reproducing the rejection analysis

Script: `extras/analyze_bindcraft_rejections.py`

Outputs:
- `rejection_breakdown.csv`: fractions of rejections attributed to each filter family.
- `annotated_rejections.csv`: row-wise flags indicating which filters were responsible.

Examples:
```bash
# Default/relaxed filters JSON expected in current working directory
python extras/analyze_bindcraft_rejections.py --input-dir /path/to/runs --filter-mode default

# Use the Filters column recorded per-run (design-mode)
python extras/analyze_bindcraft_rejections.py --input-dir /path/to/runs --filter-mode design
```
