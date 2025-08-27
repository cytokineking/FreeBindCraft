# FreeBindCraft: BindCraft with Optional PyRosetta Bypass

<p align="center">
  <img src="./free-bindcraft.png" alt="FreeBindCraft" width="720" />
</p>

This repository contains a modified version of Martin Pacesa's BindCraft (v1.52). The primary change is the introduction of an **optional PyRosetta bypass mechanism**.

For comprehensive details on the original BindCraft pipeline, features, advanced settings, and filter explanations, please refer to the **original BindCraft repository: [https://github.com/martinpacesa/BindCraft](https://github.com/martinpacesa/BindCraft)** and the [original preprint](https://www.biorxiv.org/content/10.1101/2024.09.30.615802). This fork is hosted at: [https://github.com/cytokineking/FreeBindCraft](https://github.com/cytokineking/FreeBindCraft).

## Key Modification: PyRosetta Bypass Functionality

The `--no-pyrosetta` flag, usable during both installation and runtime, enables the following behavior:

*   **OpenMM Relaxation:** Structural relaxation uses an OpenMM-based protocol instead of PyRosetta's FastRelax. This includes structure preparation with PDBFixer, ramped backbone restraints, OBC2 implicit solvation, an additional short-range repulsive term to mitigate clashes, and short MD "shakes" for early stages. This function is GPU-accelerated and typically is 2-4x faster than the CPU-based FastRelax, depending on the platform.
*   **Shape Complementarity (SC):** Replaced with an open-source implementation via `sc-rs` when available. See `sc-rs` project: [https://github.com/cytokineking/sc-rs](https://github.com/cytokineking/sc-rs). The results of sc-rs are nearly identical to the same calculations performed by PyRosetta.
*   **SASA Calculations:** Surface area and derived metrics are computed using [FreeSASA](https://github.com/mittinatten/freesasa) or a Biopython Shrake–Rupley fallback. These values align closely with the PyRosetta-derived values.
*   **Interface Residues and Alignment:** Uses Biopython-based routines for interface residue identification, RMSD, and PDB alignment.

**Important Note:** Rosetta-specific metrics that lack open-source equivalents are not computed; placeholder values are used where needed for compatibility with default filters. Evaluate design quality accordingly.

### Technical Overview & Rationale

For the motivation behind the PyRosetta bypass, implementation details (OpenMM relax, FreeSASA/Biopython SASA, `sc-rs` shape complementarity), and empirical impact on filter decisions, see the in-depth technical overview:

- `extras/FreeBindCraft_Technical_Overview.md`

## Installation

1.  Clone this modified repository:
    ```bash
    git clone https://github.com/cytokineking/FreeBindCraft [install_folder]
    ```
2.  Navigate into your install folder (`cd [install_folder]`) and run the installation script. A CUDA-compatible Nvidia graphics card is required.

    Installation script options:
    *   `--cuda CUDAVERSION`: Specify your CUDA version (e.g., '12.4').
    *   `--pkg_manager MANAGER`: Specify 'mamba' or 'conda' (default: 'conda').
    *   `--no-pyrosetta`: **Use this flag to install without PyRosetta.**

    **Example: PyRosetta-Free Installation**
    ```bash
    bash install_bindcraft.sh --cuda '12.4' --pkg_manager 'conda' --no-pyrosetta
    ```

    **Example: Standard Installation (Includes PyRosetta)**
    If you choose to install with PyRosetta, a license is required for commercial use.
    ```bash
    bash install_bindcraft.sh --cuda '12.4' --pkg_manager 'conda'
    ```

## Running BindCraft

Activate the Conda environment:
```bash
conda activate BindCraft
```

Change to your BindCraft installation directory:
```bash
# cd /path/to/your/bindcraft/folder/
```

**To use the PyRosetta bypass at runtime, include the `--no-pyrosetta` flag:**
```bash
python -u ./bindcraft.py --settings './settings_target/your_target.json' --filters './settings_filters/default_filters.json' --advanced './settings_advanced/default_4stage_multimer.json' --no-pyrosetta
```

**To run with PyRosetta enabled (assuming it was installed):**
```bash
python -u ./bindcraft.py --settings './settings_target/your_target.json' --filters './settings_filters/default_filters.json' --advanced './settings_advanced/default_4stage_multimer.json'
```

**Note:** Even if you installed BindCraft *with* PyRosetta, you can still run in bypass mode by adding the `--no-pyrosetta` flag at runtime.

For details on configuring target settings (`--settings`), filters (`--filters`), advanced parameters (`--advanced`), and other operational aspects of BindCraft, please consult the documentation in the [original BindCraft repository](https://github.com/martinpacesa/BindCraft).

## Additional New CLI Flags 

You can further control runtime behavior with these flags:

- `--verbose`: Enable detailed timing/progress logs for BindCraft internals.
- `--no-plots`: Disable saving design trajectory plots (overrides advanced settings).
- `--no-animations`: Disable saving trajectory animations (overrides advanced settings).

Example:
```bash
python -u ./bindcraft.py \
  --settings './settings_target/your_target.json' \
  --filters './settings_filters/default_filters.json' \
  --advanced './settings_advanced/default_4stage_multimer.json' \
  --no-animations --no-plots --verbose
```

## Containerized usage (Docker)

Run FreeBindCraft in a GPU-enabled Docker container (either build locally or pull a prebuilt image).

### Prerequisites (on the host)

- NVIDIA driver installed (check with `nvidia-smi`).
- Docker CE and NVIDIA Container Toolkit:
  - Linux (Ubuntu):
    1) Install Docker CE and restart the daemon.
    2) `sudo apt-get install -y nvidia-container-toolkit`
    3) `sudo nvidia-ctk runtime configure --runtime=docker && sudo systemctl restart docker`
  - Validate GPU availability inside containers:
    ```bash
    docker run --rm --gpus all nvidia/cuda:12.1.1-cudnn8-runtime-ubuntu22.04 nvidia-smi
    ```

### Option A: Build the image from this repository

```bash
docker build -t freebindcraft:gpu .
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

### Option B: Pull the prebuilt image from Docker Hub

- Image: `cytokineking/freebindcraft-no-pyrosetta:latest`
- Digest: `sha256:352d0707e1d87ae24e18b2dfc555645ae00dd806f43a5099aa2d9ee79981efc5`

```bash
docker pull cytokineking/freebindcraft-no-pyrosetta:latest
```

Sanity test (OpenMM relax):
```bash
mkdir -p ~/bindcraft_out
docker run --gpus all --rm -it \
  --ulimit nofile=65536:65536 \
  -v ~/bindcraft_out:/work/out \
  cytokineking/freebindcraft-no-pyrosetta:latest \
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
  cytokineking/freebindcraft-no-pyrosetta:latest \
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
  cytokineking/freebindcraft-no-pyrosetta:latest \
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
  cytokineking/freebindcraft-no-pyrosetta:latest \
  python bindcraft.py \
    --settings /app/settings_target/my_settings.json \
    --filters /app/settings_filters/my_filters.json \
    --advanced /app/settings_advanced/my_advanced.json \
    --no-pyrosetta
```

Notes:
- Always mount your output directory to the path set as `design_path` in your target settings (e.g., `/root/software/pdl1`).
- Always increase file descriptor limits with `--ulimit nofile=65536:65536`. The image’s entrypoint also attempts to raise the soft limit inside the container.

## Citations & External Tools

- Shape Complementarity (SC): `sc-rs` — [https://github.com/cytokineking/sc-rs](https://github.com/cytokineking/sc-rs)
- FreeSASA (SASA): [https://github.com/mittinatten/freesasa](https://github.com/mittinatten/freesasa)
- Biopython: [https://biopython.org](https://biopython.org)

## Extras: Analysis and Utility Scripts

This repository includes additional scripts and documents in `extras/` to assist with analysis and testing:

- `extras/analyze_bindcraft_rejections.py`: Analyze MPNN rejections across runs and quantify which filters are most responsible; can estimate hypothetical rankings when Rosetta-only filters trigger.
- `extras/bindcraft_rejection_analysis_spec.md`: Specification describing the rejection analysis outputs and methodology.
- `extras/BUGFIX_FILE_DESCRIPTORS.md`: Notes on file descriptor fixes in high-throughput runs.
- `extras/check_ulimit.sh`: Helper script to check system ulimit and suggest adjustments for many-parallel-file workloads.
- `extras/compare_interface_metrics_all.py`: Compute interface metrics via PyRosetta (if available), FreeSASA, and Biopython across a folder of PDBs; outputs a combined CSV.
- `extras/compare_pyrosetta_bypass_scores.py`: Side-by-side comparison of PyRosetta vs Biopython-only metrics for a folder of PDBs.
- `extras/test_openmm_relax.py`: Quick test harness for OpenMM and PyRosetta relax routines.

See `extras/README.md` for detailed usage examples and options.