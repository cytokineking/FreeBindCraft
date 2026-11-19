# FreeBindCraft: BindCraft with Optional PyRosetta Bypass

<p align="center">
  <img src="./free-bindcraft.png" alt="FreeBindCraft" width="720" />
</p>

This repository contains a modified version of Martin Pacesa's BindCraft (v1.52). It adds two **PyRosetta-free** scoring/relaxation engines (`--backend tmol` using PyTorch/[tmol](https://github.com/RosettaCommons/tmol) and `--backend openmm` / `--no-pyrosetta` using OpenMM + FASPR + `sc-rs`), `ipSAE` scoring, exact-length trajectory/validation caching, and native Apple Silicon (MetalJAX/OpenCL/MPS) support.

For comprehensive details on the original BindCraft pipeline, features, advanced settings, and filter explanations, please refer to the **original BindCraft repository: [https://github.com/martinpacesa/BindCraft](https://github.com/martinpacesa/BindCraft)** and the [original preprint](https://www.biorxiv.org/content/10.1101/2024.09.30.615802). This fork is hosted at: [https://github.com/cytokineking/FreeBindCraft](https://github.com/cytokineking/FreeBindCraft).

## Scoring & Relaxation Engines (`pyrosetta`, `tmol`, `openmm`)

Select a backend at runtime with `--backend {pyrosetta,openmm,tmol}` (or `--no-pyrosetta`, which selects `openmm` for backward compatibility):

1. **`tmol` (`--backend tmol`, PyRosetta-free):** Executes Rosetta's [`beta2016`](https://github.com/RosettaCommons/tmol) energy function via an isolated persistent PyTorch/`tmol` daemon (`$CONDA_PREFIX/tmol-venv` on Linux CUDA, `.venv-metaljax` with Apple MPS offload on macOS). It replaces PyRosetta FastRelax (`KinematicRelaxMinimizer`), batches multi-model passing AF2 complexes (`[model1, model2]`) in one `PoseStack`, and computes real `Binder_Energy_Score`, separated/repacked `dG` (`dG/dSASA`), and energy+geometry-gated `n_InterfaceHbonds` alongside `sc-rs` shape complementarity, FreeSASA/Biopython SASA, and Biopython interface/RMSD checks.
2. **`openmm` ("Classic" FreeBindCraft, `--backend openmm` or `--no-pyrosetta`):** Relaxes structures on GPU with Amber14 + OBC2 implicit solvent + short-range clash repulsion + MD shakes (`OpenMM`) followed by `FASPR` side-chain packing, then scores geometry via `sc-rs`, FreeSASA/Biopython SASA, and Biopython RMSDs. Rosetta-specific energetic and H-bond network metrics emit compatibility placeholders.
3. **`pyrosetta` (`--backend pyrosetta`, default when installed without `--no-pyrosetta`):** Original Rosetta FastRelax, `InterfaceAnalyzerMover`, and `BuriedUnsatHbonds` (`DAlphaBall.gcc`). Requires a Rosetta Commons license for commercial use.

### Feature & Metric Comparison Across Engines

| Metric / Capability | `pyrosetta` (Rosetta) | `tmol` (PyTorch / PyRosetta-free) | `openmm` ("Classic" `--no-pyrosetta`) |
|---|---|---|---|
| **PyRosetta Commercial License Required** | **Yes** (non-academic) | **No** (BSD-3-Clause `tmol` / MIT deps; repo is Elastic-2.0) | **No** (MIT `FASPR` / open-source deps; repo is Elastic-2.0) |
| **Structure Relaxation** | CPU Rosetta `FastRelax` (`beta_nov16`) | GPU **`tmol` `FastRelax`** (`KinematicRelaxMinimizer`, `beta2016`, batched `[model1, model2]`) | GPU **`OpenMM`** (Amber14 + OBC2 + MD shakes) + CPU `FASPR` |
| **`Binder_Energy_Score`** | Rosetta `TotalEnergyMetric` | Real `tmol` `beta2016` chain energy *(~12–22% less negative than Rosetta; see note)* | Placeholder (`-1.0`) |
| **`dG` & `dG/dSASA`** | `InterfaceAnalyzerMover` (`pack_separated=True`) | Real `tmol` `E_bound - E_separated` w/ interface `pack_rotamers` *(within ~0.1–3.7 REU on PDL1)* | Placeholders (`-10.0`, `0.0`) |
| **`n_InterfaceHbonds` & `%`** | `InterfaceAnalyzerMover` | Real `tmol` `hbond` energy + donor-H···acceptor geometry *(exact PDL1 match)* | Placeholders (`5`, `60.0%`) |
| **`ShapeComplementarity`** | Rosetta SC | `sc-rs` (Rust Lawrence–Colman SC) | `sc-rs` (Rust Lawrence–Colman SC) |
| **`dSASA`, `Interface_SASA_%`, `Surface_Hydrophobicity`** | Rosetta SASA + `LayerSelector` | FreeSASA (or Biopython SR fallback) + area-weighted hydrophobicity | FreeSASA (or Biopython SR fallback) + area-weighted hydrophobicity |
| **`n_InterfaceResidues` & `InterfaceAAs`** | Biopython 4.0 Å KD-tree | Biopython 4.0 Å KD-tree | Biopython 4.0 Å KD-tree |
| **`PackStat`** | Rosetta `InterfaceAnalyzer` | Placeholder (`0.65`) | Placeholder (`0.65`) |
| **`n_InterfaceUnsatHbonds` & `%`** | Rosetta `BuriedUnsatHbonds` (`DAlphaBall.gcc`) | Placeholders (`1`, `0.0%`) | Placeholders (`1`, `0.0%`) |
| **Alignment & RMSDs** | Rosetta `AlignChainMover` / `RMSDMetric` | Biopython Cα Superimposer / unaligned RMSD | Biopython Cα Superimposer / unaligned RMSD |

> **Calibration & Operational Notes:**
> - Across the 5 bundled PDL1 complexes ([PDL1-pyrosetta.zip](performance_data/pdl1_miniprotein/PDL1-pyrosetta.zip)), `tmol` reproduces PyRosetta `interface_nres` (`25, 19, 18, 29, 22`) and `interface_interface_hbonds` (`10, 7, 5, 12, 11`) exactly, matches `interface_dG` within **0.12–3.66 REU** and `dSASA` within **1.4–2.1%**, and scores `Binder_Energy_Score` ~12–22% less negatively due to reference/terminal energy handling.
> - A long-lived JSON-line daemon (`BINDCRAFT_TMOL_PERSISTENT=1`, opt out with `=0`) amortizes PyTorch/`tmol` initialization (~6.5s). On Linux CUDA, `--backend tmol` defaults `XLA_PYTHON_CLIENT_MEM_FRACTION=0.70` when unset so JAX leaves headroom for the ~824–912 MiB PyTorch/`tmol` worker (explicit `XLA_PYTHON_CLIENT_*` env vars always win).
> - Bare Linux containers running `pyrosetta` must have system `libgfortran5` installed (e.g., `apt-get install libgfortran5`, already in [Dockerfile](Dockerfile)) because environment activation clears `LD_LIBRARY_PATH` when calling `functions/DAlphaBall.gcc`.

### Performance Comparison (RTX PRO 6000 CUDA vs. Apple M5 Max)

Measured head-to-head on the 5 H-stripped PDL1 complexes (`1_PDL1..5_PDL1` from `PDL1-pyrosetta.zip`):

| Engine | RTX PRO 6000 Blackwell Call #1 Cold (Relax / Score) | RTX PRO 6000 Steady-State Mean `#2–5` (Relax / Score) | Apple M5 Max Call #1 Cold (Relax / Score) | Apple M5 Max Steady-State Mean `#2–5` (Relax / Score) | Worker Memory |
|---|---:|---:|---:|---:|---|
| **`tmol` (persistent)** | 10.56 s (9.49 / 1.07) | **3.47 s** (2.38 / 1.09) | 16.21 s (15.27 / 0.93) | **10.25 s** (9.37 / 0.88) | CUDA: 912 MiB VRAM / 3895 MB RSS; M5 Max: 2.7–2.9 GB RSS |
| **`openmm` (`CUDA`/`OpenCL` + `FASPR`)** | 12.05 s (11.55 / 0.51) | **11.07 s** (10.65 / 0.42) | 15.18 s (14.84 / 0.34) | **13.83 s** (13.53 / 0.29) | Subprocess per call |
| **`pyrosetta` (CPU + `DAlphaBall`)** | 20.69 s (18.16 / 2.53) | **18.24 s** (16.15 / 2.09) | 14.31 s (12.22 / 2.09) | **11.70 s** (9.92 / 1.78) | Host process |
| **`tmol` (one-shot `=0`)** | 18.76 s (10.27 / 8.49) | **18.51 s** (10.34 / 8.17) | 20.08 s (15.20 / 4.88) | **18.66 s** (13.91 / 4.75) | Re-initialized every call (~6.5s cold) |

* **Steady-state throughput:** Persistent `tmol` is **5.26× faster than PyRosetta** and **3.19× faster than OpenMM CUDA+FASPR** on RTX PRO 6000 Blackwell (and **1.14× faster than native arm64 PyRosetta** / **1.35× faster than OpenMM OpenCL** on Apple M5 Max).
* **Batched `tmol` `PoseStack` relaxation scaling:** BindCraft automatically batches `[model1, model2]` (`n=2`) relaxation via `pr_relax_batch` when using `--backend tmol`. Warm in-process `PoseStack` timings on `3_PDL1_l66_s627125_mpnn8_model1.pdb` (`relax` + `score/repack`):
  * **RTX PRO 6000 CUDA (`torch.cuda.synchronize`):** `n=1`: `1.81 s/pose` (207 MB peak VRAM) → `n=2`: **`1.31 s/pose` (1.38×, 290 MB)** → `n=4`: **`0.85 s/pose` (2.14×, 353 MB)** → `n=8`: **`0.73 s/pose` (2.46×, 492 MB)**.
  * **Apple M5 Max (`MPS` offload):** `n=1`: `7.93 s/pose` → `n=2`: **`5.50 s/pose` (1.44×)** → `n=4`: **`4.07 s/pose` (1.95×)**.

## ipSAE Scoring

FreeBindCraft includes support for **ipSAE (interface predicted Structural Alignment Error)**, a metric for evaluating protein-protein interface quality based on AlphaFold's PAE matrix. This implementation is based on [Dunbrack et al. (2025)](https://www.biorxiv.org/content/10.1101/2025.02.10.637595v1).

**What is ipSAE?**
- ipSAE transforms interface PAE values using a PTM-style scoring function with adaptive normalization
- Higher scores indicate better predicted interface quality (range 0–1)
- Unlike i_pTM which uses a global d0, ipSAE uses per-residue d0 based on the number of interface contacts
- ipSAE is computed for both trajectory statistics and final MPNN design statistics

**Ranking by ipSAE:**
By default, final designs are ranked by `i_pTM`. You can optionally rank by `ipSAE` instead:

```bash
python bindcraft.py --settings ... --rank-by ipSAE
```

The interactive CLI (both local and Docker) also prompts for ranking method selection.

### Technical Overview & Rationale

For the motivation behind the PyRosetta bypass, implementation details (OpenMM relax, FreeSASA/Biopython SASA, `sc-rs` shape complementarity), and empirical impact on filter decisions, see the in-depth technical overview:

- `technical_overview/FreeBindCraft_Technical_Overview.md`

## Installation

### NVIDIA GPUs

1.  Clone this modified repository:
    ```bash
    git clone https://github.com/cytokineking/FreeBindCraft [install_folder]
    ```
2.  Navigate into your install folder (`cd [install_folder]`) and run the installation script. Linux x86_64, Conda or Mamba, Git, curl, a C compiler, and NVIDIA driver **575.57.08 or newer** are required.

    Conda supplies Python 3.12 and native tools. The [tested pip profile](installation/linux-cuda12-requirements.txt) supplies JAX 0.11.2, NumPy 2.5.4, Flax 0.12.10, OpenMM 8.6.1, and their pinned CUDA 12.9 dependencies. By default, the installer also creates an isolated `$CONDA_PREFIX/tmol-venv` environment containing PyTorch 2.11.0+cu129 and [tmol](https://github.com/RosettaCommons/tmol) 1.0.0 ([tmol CUDA profile](installation/tmol-linux-cuda12-requirements.txt)) so `tmol` stays isolated from JAX's CUDA runtime.

    Installation script options:
    *   `--pkg_manager MANAGER`: Specify 'mamba' or 'conda' (default: 'conda').
    *   `--env-name NAME`: Choose the Conda environment name (default: `BindCraft`). Use a new name when migrating an existing Python 3.10 environment.
    *   `--no-pyrosetta`: **Use this flag to install without PyRosetta.**
    *   `--no-tmol`: Skip creating the isolated `$CONDA_PREFIX/tmol-venv` PyTorch/`tmol[cu12]` environment.
    *   `--skip-gpu-check`: Skip GPU execution checks only when building a container image without GPU access.
    *   The legacy `--cuda` and `--fix-channels` options remain accepted. CUDA libraries now come from pip, and the installer selects conda-forge without editing your channel configuration.

    **Example: PyRosetta-Free Installation (Includes `tmol` + OpenMM/FASPR)**
    ```bash
    bash install_bindcraft.sh --no-pyrosetta
    ```

    **Example: Full Installation (PyRosetta + `tmol` + OpenMM/FASPR)**
    If you choose to install with PyRosetta, a license is required for commercial use. On bare Linux containers without `apt` `libgfortran5`, install `libgfortran5` so `functions/DAlphaBall.gcc` can run (`installation/linux-activate.sh` clears `LD_LIBRARY_PATH`).
    ```bash
    bash install_bindcraft.sh
    ```

The installer checks real JAX GPU gradients, OpenMM CUDA execution, `tmol` CUDA scoring, and (unless `--no-pyrosetta` is passed) PyRosetta + `DAlphaBall.gcc`. Activation selects CUDA for JAX/OpenMM, sets `TMOL_PYTHON`, and clears `LD_LIBRARY_PATH` so an ambient toolkit cannot replace pip's CUDA libraries; deactivation restores your previous settings. Completed model weights are reused when rerunning installation. To check an installed environment again:

```bash
conda activate BindCraft
python extras/check_linux_installation.py --no-pyrosetta
```

Omit `--no-pyrosetta` from that check for a full installation (or add `--no-tmol` if you skipped `tmol`). Installing with `--no-pyrosetta` omits the optional package; it does not uninstall PyRosetta from an existing environment.

You can also run BindCraft entirely inside a Docker container without a local installation. See the Containerized usage section below for details and an interactive wrapper.

### Apple Silicon

FreeBindCraft runs natively on Apple Silicon, with MetalJAX, OpenMM/FASPR, and PyTorch/`tmol` MPS optimizations enabled by default:

```bash
bash install_bindcraft.sh
./run_apple_silicon.sh --settings settings_target/ubiquitin_apple_small.json \
  --advanced settings_advanced/apple_silicon_ubiquitin_small.json --backend tmol --max-attempts 1
```

Running `install_bindcraft.sh` installs PyRosetta, `tmol` (PyTorch MPS offload), and OpenMM/FASPR (`install_apple_silicon.sh --no-pyrosetta` omits PyRosetta; `--no-tmol` skips PyTorch/`tmol`). Use `--backend tmol`, `--backend openmm`, or `--backend pyrosetta` (or `--no-pyrosetta` for OpenMM/FASPR).
See the [Apple Silicon guide](apple_silicon/README.md) for setup and memory options,
and the [M5 Max versus NVIDIA benchmarks](apple_silicon/benchmarks/README.md).

## Running BindCraft

Activate the Conda environment:
```bash
conda activate BindCraft
```

Change to your BindCraft installation directory:
```bash
# cd /path/to/your/bindcraft/folder/
```

**To run with the PyRosetta-free `tmol` backend (`beta2016` FastRelax + real `dG` and interface H-bonds):**
```bash
python -u ./bindcraft.py --settings './settings_target/your_target.json' --filters './settings_filters/default_filters.json' --advanced './settings_advanced/default_4stage_multimer.json' --backend tmol
```

**To use classic OpenMM + FASPR PyRosetta bypass mode:**
```bash
python -u ./bindcraft.py --settings './settings_target/your_target.json' --filters './settings_filters/default_filters.json' --advanced './settings_advanced/default_4stage_multimer.json' --backend openmm
# equivalent to --no-pyrosetta
```

**To run with PyRosetta enabled (assuming it was installed):**
```bash
python -u ./bindcraft.py --settings './settings_target/your_target.json' --filters './settings_filters/default_filters.json' --advanced './settings_advanced/default_4stage_multimer.json' --backend pyrosetta
```

**Note:** Even if you installed BindCraft *with* PyRosetta, you can freely select `--backend tmol`, `--backend openmm`, or `--no-pyrosetta` at runtime.

For details on configuring target settings (`--settings`), filters (`--filters`), advanced parameters (`--advanced`), and other operational aspects of BindCraft, please consult the documentation in the [original BindCraft repository](https://github.com/martinpacesa/BindCraft).

## Interactive CLI

If you run BindCraft without `--settings` in a terminal (TTY) or pass `--interactive`, an interactive setup wizard guides you through configuration.

- Quickstart (local):
  ```bash
  python bindcraft.py               # auto-detects TTY and enters interactive mode
  # or
  python bindcraft.py --interactive
  ```
- What it asks:
  - Project/binder name, input PDB path (validated must exist), output directory (auto-created)
  - Design type: Miniprotein (≥31 aa) or Peptide (8–30 aa)
    - Peptides: default length 8–25; filters limited to `peptide_filters`, `peptide_relaxed_filters`, `no_filters`; advanced limited to peptide profiles
    - Miniproteins: default length 65–150 (min ≥31); filters limited to `default_filters`, `relaxed_filters`, `no_filters`; advanced excludes peptide profiles
  - Verbose logging, plots, animations, and scoring/relaxation backend (`pyrosetta`, `openmm`, or `tmol`)
  - Ranking method for final designs: `i_pTM` (default) or `ipSAE`
- Behavior:
  - Generates a target settings JSON inside your chosen output directory and then runs the standard pipeline

## Additional New CLI Flags 

You can further control runtime behavior with these flags:

- `--verbose`: Enable detailed timing/progress logs for BindCraft internals.
- `--backend {pyrosetta,openmm,tmol}`: Explicitly select the relaxation and interface scoring engine (`pyrosetta`, `openmm`, or `tmol`).
- `--no-pyrosetta`: Legacy shorthand selecting the `openmm` backend (cannot be combined with `--backend pyrosetta`).
- `--debug-pdbs`: Write intermediate PDBs from OpenMM relax (`.debug_deconcat.pdb`, `.debug_pdbfixer.pdb`, `.debug_post_initial_relax.pdb`, `.debug_post_faspr.pdb`)
- `--no-plots`: Disable saving design trajectory plots (overrides advanced settings).
- `--no-animations`: Disable saving trajectory animations (overrides advanced settings).
- `--rank-by {i_pTM,ipSAE}`: Metric to rank final designs by (default: `i_pTM`). Use `ipSAE` to rank by interface predicted Structural Alignment Error instead.
- `--max-attempts N`: Run at most N new trajectory attempts across all lengths, including rejected trajectories. The last attempt finishes its downstream processing, then accepted designs are ranked using `--rank-by`, even if the requested design count was not reached. This works with caching on or off, defaults to no attempt cap, and resets for each invocation when resuming.
- `--attempts-per-length N`: With the trajectory cache enabled, try up to N trajectories at each selected length per visit (default: 5), including the first attempt.
- `--accepted-per-length N`: With the trajectory cache enabled, move to the next length after a visit reaches N final-filter-accepted designs (default: 1), or its attempt limit.
- `--seed N`: Seed NumPy's trajectory and model sampling for reproducible checks within a pinned runtime. OpenMM MD and different GPU backends need not be bitwise reproducible.

Example:
```bash
python -u ./bindcraft.py \
  --settings './settings_target/your_target.json' \
  --filters './settings_filters/default_filters.json' \
  --advanced './settings_advanced/default_4stage_multimer.json' \
  --no-animations --no-plots --verbose --rank-by ipSAE
```

### Exact-length trajectory cache

Add `--trajectory-cache` to reuse compiled design calculations on CUDA or Apple Silicon, inspired by [BindCraft2](https://github.com/PacesaLab/BindCraft2/blob/8c1a858fdc1bfc94c60846ff11a0d01c9aa44136/bindcraft/af2.py#L461-L469). Caching is **off by default**. You can also enable it with `"trajectory_cache": true` in the advanced JSON; `--no-trajectory-cache` overrides that setting for a single run.

Use `--attempts-per-length N` to set the maximum number of attempts per selected length and visit. Use `--accepted-per-length N` to set how many accepted designs end that visit. For example, `--trajectory-cache --attempts-per-length 5 --accepted-per-length 3` moves on after three accepted designs or five attempts, whichever happens first. You can set both values to 5 to target five accepted designs within five attempts. Both values must be positive integers; `--max-attempts` remains a separate limit on all new trajectories in the run.

With caching enabled:

- The length range is divided into up to four bins. FreeBindCraft alternates between bins, running a block of trajectories at one exact length before moving on.
- Each visit attempts at most the configured number of trajectories at its selected length, ending when the accepted-design target is reached. Only designs that pass the final filters count; earlier confidence checks do not. The configured MPNN output limit and overall design goal still apply.
- The accepted-design target is checked after each trajectory completes. If a trajectory produces multiple accepted designs, all count and the visit can exceed its target; the existing per-trajectory MPNN limit still applies.
- Each bin starts at its midpoint, then explores unvisited lengths furthest from previous visits. Sweeps continue until the design goal or a stopping limit is reached. For 65–200 with a goal of 20 designs, the block order begins **82 → 116 → 150 → 184 → 65 → 99 → 133 → 167**.
- The first trajectory at a new length compiles; subsequent trajectories reuse both design and validation calculations (`complex_prediction_model` and lazily initialized `binder_prediction_model`) with fresh seeds and optimizer state. This eliminates `~29–32 s/trajectory` of CUDA XLA validation recompilation at repeated lengths (measured on RTX PRO 6000) in addition to the trajectory-design compile speedup below. Monomer validation models are constructed lazily only when at least one MPNN variant passes complex AF2 filters, saving `~14–15 s` of XLA compile on CUDA trajectories with zero passing candidates. Changing length releases both design and validation runtimes and clears JAX compilation caches. There is no padding or scheduling worker process.

For example, with the default settings, a visit to length 82 advances to 116 after five unsuccessful attempts, or after a trajectory produces at least one final-filter-accepted design.

Resume preserves bins, visited lengths and partial visits when campaign settings and accepted files are unchanged. If either limit changes, the current visit is evaluated against the new limits: raising a limit can continue a visit, while lowering one can end it immediately. Older checkpoints are migrated to the selected limits. With caching off, lengths are sampled independently as before, and these two options have no effect. Caching can increase peak memory during validation, and allocators may retain freed memory for reuse.

**Design-only measurements at one repeated length:**

| Backend | First call, off → on | Later calls, off → on | Five-call mean, off → on |
|---|---:|---:|---:|
| NVIDIA RTX 6000 Pro | 89.65 → 89.21 s (1.00×) | 89.15 → 11.96 s (7.46×) | 89.25 → 27.41 s (3.26×) |
| NVIDIA RTX 5090 | 89.18 → 88.85 s (1.00×) | 88.68 → 13.68 s (6.48×) | 88.78 → 28.72 s (3.09×) |
| NVIDIA L40S | 173.58 → 179.50 s (0.97×) | 174.32 → 24.59 s (7.09×) | 174.18 → 55.57 s (3.13×) |
| Apple M5 Max | 91.16 → 93.18 s (0.98×) | 92.60 → 88.13 s (1.05×) | 92.31 → 89.14 s (1.04×) |

Each result averages two opposite-order pairs of five-call blocks on ubiquitin with a 100-residue binder and one repeated seed. These are complete calls to the production three-stage design function with a **shortened 20/10/3-iteration schedule**, not the default four-stage schedule. The first call includes model setup and compilation; the later-call mean covers four calls per block. Timings include trajectory output but exclude relaxation, MPNN and validation. Outputs matched within each backend. CUDA used JAX 0.6.0 with autotuning disabled; Metal used JAX 0.11.2 and MetalJAX 0.11.9+fbc3. Persistent compilation caching was disabled. Results are workload-specific; the small Metal difference is not a general performance guarantee.

## Containerized usage (Docker)

Run FreeBindCraft in a GPU-enabled Docker container (build locally).

### Prerequisites (on the host)

- Linux x86_64 host with NVIDIA driver **575.57.08 or newer** (check with `nvidia-smi`).
- Docker CE and NVIDIA Container Toolkit:
  - Linux (Ubuntu):
    1) Install Docker CE and restart the daemon.
    2) `sudo apt-get install -y nvidia-container-toolkit`
    3) `sudo nvidia-ctk runtime configure --runtime=docker && sudo systemctl restart docker`
  - Validate GPU availability inside containers:
    ```bash
    docker run --rm --gpus all nvidia/cuda:12.9.1-base-ubuntu22.04 nvidia-smi
    ```

### Build the image from this repository

Build without PyRosetta (default):
```bash
docker build -t freebindcraft:gpu .
```

The build downloads the model weights and installs the pinned pip CUDA 12.9
profile. A GPU is needed only at runtime. Local environments, downloaded
`params/`, and generated `runs/` are excluded from the build context.

Check JAX gradients and OpenMM CUDA execution in the built image:
```bash
docker run --rm --gpus all freebindcraft:gpu \
  python extras/check_linux_installation.py --no-pyrosetta
```

Build with PyRosetta:
```bash
docker build --build-arg WITH_PYROSETTA=true -t freebindcraft:pyrosetta .
```

Sanity test (OpenMM relax):
```bash
mkdir -p ~/bindcraft_out
docker run --gpus all --rm -it \
  --ulimit nofile=65536:65536 \
  -v ~/bindcraft_out:/work/out \
  freebindcraft:gpu \
  python extras/test_openmm_relax.py example/PDL1.pdb /work/out/relax_test
```

 

### Running BindCraft in Docker (general guidance)

Decide where your inputs and outputs live. Typical layout:
- Host input PDB: `/path/on/host/your_target.pdb`
- Host output dir: `/path/on/host/run_outputs`
- Settings files: either use repo defaults inside the container or mount your own copies from the host

Examples:
1) Use repository defaults for settings (inside container), but mount outputs to host.
```bash
mkdir -p /path/on/host/run_outputs
docker run --gpus all --rm -it \
  --ulimit nofile=65536:65536 \
  -v /path/on/host/run_outputs:/root/software/pdl1 \
  freebindcraft:gpu \
  python bindcraft.py \
    --settings settings_target/PDL1.json \
    --filters settings_filters/default_filters.json \
    --advanced settings_advanced/default_4stage_multimer.json \
    --no-pyrosetta
```

2) Mount a custom target settings JSON from the host (also mount your input PDB):
```bash
mkdir -p /path/on/host/run_outputs
docker run --gpus all --rm -it \
  --ulimit nofile=65536:65536 \
  -v /path/on/host/run_outputs:/root/software/pdl1 \
  -v /path/on/host/my_settings.json:/app/settings_target/my_settings.json:ro \
  -v /path/on/host/your_target.pdb:/app/example/your_target.pdb:ro \
  freebindcraft:gpu \
  python bindcraft.py \
    --settings /app/settings_target/my_settings.json \
    --filters settings_filters/default_filters.json \
    --advanced settings_advanced/default_4stage_multimer.json \
    --no-pyrosetta
```

3) Mount custom filters/advanced settings from the host too:
```bash
docker run --gpus all --rm -it \
  --ulimit nofile=65536:65536 \
  -v /path/on/host/run_outputs:/root/software/pdl1 \
  -v /path/on/host/my_settings.json:/app/settings_target/my_settings.json:ro \
  -v /path/on/host/my_filters.json:/app/settings_filters/my_filters.json:ro \
  -v /path/on/host/my_advanced.json:/app/settings_advanced/my_advanced.json:ro \
  -v /path/on/host/your_target.pdb:/app/example/your_target.pdb:ro \
  freebindcraft:gpu \
  python bindcraft.py \
    --settings /app/settings_target/my_settings.json \
    --filters /app/settings_filters/my_filters.json \
    --advanced /app/settings_advanced/my_advanced.json \
    --no-pyrosetta
```

Notes:
- Always mount your output directory to the path set as `design_path` in your target settings (e.g., `/root/software/pdl1`).
- Always increase file descriptor limits with `--ulimit nofile=65536:65536`. The image’s entrypoint also attempts to raise the soft limit inside the container.

### Interactive wrapper for Docker

After you build an image locally, you can use the interactive Docker wrapper:

```bash
python docker_cli.py
```

This wizard:
- Lists local Docker images and lets you choose one (defaults to `freebindcraft:gpu` if present, or enter another local image tag such as `freebindcraft:pyrosetta`)
- Asks for a single GPU index to expose
- Prompts for all inputs identical to the local interactive CLI (design type, validated PDB path, output directory auto-created, chains, hotspots, peptide/miniprotein-specific filters and advanced profiles, verbose/plots/animations, PyRosetta, ranking method)
- Writes a target settings JSON into your chosen output directory
- Launches `docker run --rm -it --gpus device=<idx>` with host paths mounted to identical locations in the container so your paths work unchanged

On completion or Ctrl+C, the container exits and is auto-removed; outputs persist on the host via the mounted output directory.

## Citations & External Tools

- `tmol` (GPU FastRelax, `beta2016` scoring, and rotamer repacking): [https://github.com/RosettaCommons/tmol](https://github.com/RosettaCommons/tmol)
- Shape Complementarity (SC): `sc-rs` — [https://github.com/cytokineking/sc-rs](https://github.com/cytokineking/sc-rs)
- FreeSASA (SASA): [https://github.com/mittinatten/freesasa](https://github.com/mittinatten/freesasa)
- FASPR (Side-chain Packing): [https://github.com/tommyhuangthu/FASPR](https://github.com/tommyhuangthu/FASPR)
- Biopython: [https://biopython.org](https://biopython.org)
- ipSAE (interface predicted Structural Alignment Error): [Dunbrack et al. 2025](https://www.biorxiv.org/content/10.1101/2025.02.10.637595v1)

## Extras: Analysis and Utility Scripts

This repository includes additional scripts and documents in `extras/` to assist with analysis and testing:

- `extras/analyze_bindcraft_rejections.py`: Analyze MPNN rejections across runs and quantify which filters are most responsible; can estimate hypothetical rankings when Rosetta-only filters trigger.
- `extras/bindcraft_rejection_analysis_spec.md`: Specification describing the rejection analysis outputs and methodology.
- `extras/BUGFIX_FILE_DESCRIPTORS.md`: Notes on file descriptor fixes in high-throughput runs.
- `extras/check_ulimit.sh`: Helper script to check system ulimit and suggest adjustments for many-parallel-file workloads.
- `extras/compare_interface_metrics_all.py`: Compute interface metrics via PyRosetta (if available), FreeSASA, and Biopython across a folder of PDBs; outputs a combined CSV.
- `extras/compare_pyrosetta_bypass_scores.py`: Side-by-side comparison of PyRosetta vs Biopython-only metrics for a folder of PDBs.
- `extras/rescore_accepted_with_rosetta.py`: Rescore accepted (Ranked preferred) PDBs from a `--no-pyrosetta` run using PyRosetta, evaluate whether PyRosetta-driven filters would reject each design, and report failing filters. Supports dynamic filters (`default|relaxed|design|custom`), parallel workers, and an optional `--fast-relax` mode to run PyRosetta FastRelax prior to scoring.
- `extras/test_openmm_relax.py`: Quick test harness for OpenMM and PyRosetta relax routines.

See `extras/README.md` for detailed usage examples and options.

## Contributing

Pull requests are welcome.

## Known Issues

- OpenMM backend selection: The Linux installer and Docker image select the tested CUDA backend. The Apple Silicon profile uses OpenCL. Relaxation runs in a subprocess for isolation.
- OpenCL JIT "log spam": You may see repeated lines such as "Failed to read file: /tmp/dep-76532a.d" emitted to stdout/stderr during OpenCL kernel builds.
  - Recommend line buffering to keep logs responsive when filtering
    - Combined output (preferred when piping/teeing):
    ```bash
    python -u ./bindcraft.py \
      --settings ./settings_target/your_target.json \
      --filters ./settings_filters/default_filters.json \
      --advanced ./settings_advanced/default_4stage_multimer.json \
      --no-pyrosetta \
      |& stdbuf -oL -eL grep -v -E 'Failed to read file: .*/dep-[0-9a-fA-F]+\.d([[:space:]]*)?$' \
      | stdbuf -oL tee -a ./output/your_target.log
    ```
    - Stderr-only (no further pipes):
    ```bash
    python -u ./bindcraft.py \
      --settings ./settings_target/your_target.json \
      --filters ./settings_filters/default_filters.json \
      --advanced ./settings_advanced/default_4stage_multimer.json \
      --no-pyrosetta \
      2> >(stdbuf -oL -eL grep -v -E 'Failed to read file: .*/dep-[0-9a-fA-F]+\.d([[:space:]]*)?$' >&2)
    ```
    - Notes:
      - `stdbuf` forces line buffering for both stdout and stderr, preventing pipeline buffering delays.
      - Use POSIX character class `[[:space:]]` for portability.
      - If preserving exit code across pipes matters, consider: `set -o pipefail`.

## License

Source-available under the Elastic License 2.0 (Elastic-2.0). It permits royalty-free commercial use, including internal pharma/biotech research, development, production, and deployment on your own hardware or cloud. The main restriction: you may not provide the software (or a service exposing a substantial set of its functionality) to third parties as a hosted or managed service without an agreement from the copyright holder. See [LICENSE](LICENSE).

FreeBindCraft adapts upstream [BindCraft](https://github.com/martinpacesa/BindCraft) (`Copyright (c) 2024 Martin Pacesa`) and [ipSAE](https://www.biorxiv.org/content/10.1101/2025.02.10.637595v1) (Roland Dunbrack, Fox Chase Cancer Center) under their MIT licenses. The `--backend tmol` engine runs on [tmol](https://github.com/RosettaCommons/tmol) and PyTorch (BSD-3-Clause), and the `--backend openmm` engine uses OpenMM, FreeSASA (MIT), `sc-rs` (MIT), and `FASPR` (MIT); optional `--backend pyrosetta` (`DAlphaBall.gcc`) remains governed by the Rosetta Commons license and requires a commercial license for non-academic use. See [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md) for preserved upstream notices and third-party license details.
