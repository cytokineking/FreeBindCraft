# FreeBindCraft on Apple Silicon

FreeBindCraft runs natively on Apple Silicon Macs. AlphaFold and ProteinMPNN use the Apple GPU. Relaxation uses PyRosetta on the CPU, or OpenMM on the GPU.

**Apple performance optimizations are enabled by default.** The installer handles the backend setup; no extra performance flags or manual patches are needed.

## Install

From the repository directory:

```bash
xcode-select --install  # If Apple Command Line Tools are not installed
bash install_bindcraft.sh
```

Finish the Command Line Tools installation before running the second command. Dependencies and model weights install inside the checkout; no Conda environment is needed.

PyRosetta is included by default. To use OpenMM/FASPR instead:

```bash
bash install_bindcraft.sh --no-pyrosetta
```

To update an existing installation, pull the latest code and rerun the installer.

## Run

Try the small ubiquitin example:

```bash
./run_apple_silicon.sh \
  --settings settings_target/ubiquitin_apple_small.json \
  --advanced settings_advanced/apple_silicon_ubiquitin_small.json \
  --max-attempts 1 --seed 42 --no-plots --no-animations
```

This designs a 70-residue binder against ubiquitin. Results go to `runs/ubiquitin-apple-example/output`. The command runs one design attempt; it may produce no designs that pass the filters.

Use your own target and advanced settings for other jobs. Change `design_path` in the target settings to choose a new output directory. No environment activation is needed.

Add `--no-pyrosetta` to a run command to select OpenMM/FASPR. The launcher selects it automatically if PyRosetta is absent.

## Performance and memory

On the M5 Max, the small benchmark takes about four minutes and uses about 8 GiB of memory. The larger CA2 benchmark takes about 33 minutes and uses about 29 GiB. See the [NVIDIA comparison and time breakdowns](benchmarks/README.md).

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
