# FreeBindCraft on Apple Silicon

FreeBindCraft runs natively on Apple Silicon Macs. AlphaFold and ProteinMPNN use the Apple GPU via MetalJAX. Relaxation and scoring support all three backends: **`tmol` (`beta2016` on PyTorch + Metal/MPS offload)**, **`openmm`/`FASPR` (`OpenCL` GPU)**, and **CPU `pyrosetta`**.

**Apple performance optimizations are enabled by default.** The installer handles the backend setup; no extra performance flags or manual patches are needed.

## Install

From the repository directory:

```bash
xcode-select --install  # If Apple Command Line Tools are not installed
bash install_bindcraft.sh
```

Finish the Command Line Tools installation before running the second command. Dependencies and model weights install inside the checkout; no Conda environment is needed.

Running `install_bindcraft.sh` installs MetalJAX, OpenMM/FASPR, PyTorch/`tmol`, and PyRosetta by default. Calling `bash install_apple_silicon.sh` directly installs PyTorch/`tmol` and OpenMM/FASPR without PyRosetta (`--with-pyrosetta` adds it; `--no-tmol` omits PyTorch/`tmol`).

```bash
bash install_apple_silicon.sh --no-pyrosetta   # PyTorch/tmol + OpenMM/FASPR only (PyRosetta-free)
```

To update an existing installation, pull the latest code and rerun the installer.

## Run

Try the small ubiquitin example (using the PyRosetta-free `--backend tmol` engine or classic OpenMM/FASPR):

```bash
./run_apple_silicon.sh \
  --settings settings_target/ubiquitin_apple_small.json \
  --advanced settings_advanced/apple_silicon_ubiquitin_small.json \
  --backend tmol \
  --max-attempts 1 --seed 42 --no-plots --no-animations
```

This designs a 70-residue binder against ubiquitin. Results go to `runs/ubiquitin-apple-example/output`. The command runs one design attempt; it may produce no designs that pass the filters.

Use your own target and advanced settings for other jobs. Change `design_path` in the target settings to choose a new output directory. No environment activation is needed.

Pass `--backend pyrosetta`, `--backend tmol`, or `--backend openmm` (or `--no-pyrosetta` for OpenMM/FASPR). If PyRosetta is omitted at install time, the launcher selects `--no-pyrosetta` unless `--backend tmol` is passed. See the [README](../README.md#scoring--relaxation-engines-pyrosetta-tmol-openmm) for the 3-engine feature matrix and M5 Max vs. RTX PRO 6000 benchmark comparison.

## Performance and memory

On the M5 Max, the small benchmark takes about four minutes and uses about 8 GiB of memory. The larger CA2 benchmark takes about 33 minutes and uses about 29 GiB. See the [NVIDIA comparison and time breakdowns](benchmarks/README.md).

For repeated seeds at interleaved exact binder lengths, the launcher accepts `--trajectory-cache`, `--attempts-per-length N` (default 5 attempts per visit), and `--accepted-per-length N` (default 1 accepted design per visit). See [usage, scheduling and design speedups](../README.md#exact-length-trajectory-cache).

The default favors speed by running without chunking. For lower memory use on targets plus binders exceeding 384 residues, set `FREEBINDCRAFT_METAL_CHUNKING=1` before your run command. This processes smaller batches and runs more slowly. Leave memory available for macOS and other applications.

<details>
<summary>Troubleshooting</summary>

To check the GPU and native tools:

```bash
JAX_PLATFORMS=metal .venv-metaljax/bin/python extras/check_apple_installation.py
```

The backend optimizations are opt-out. For diagnosing a problem:

| Optimization | How to disable it |
|---|---|
| Attention kernels | Add `--no-attention-epilogues` to the run command |
| Conversion optimization | Set `METALJAX_SLICE_BEFORE_CONVERT=0` before the run command |

Remove these overrides to restore the defaults.

</details>
