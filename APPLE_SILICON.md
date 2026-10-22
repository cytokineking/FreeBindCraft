# Apple Silicon (experimental)

FreeBindCraft runs natively on Apple Silicon using
[metaljax](https://github.com/eterevsky/metaljax) for AlphaFold design and
ProteinMPNN on the Apple GPU. Relaxation uses native PyRosetta on the CPU by
default, or OpenMM on the GPU through OpenCL with `--no-pyrosetta`.

## Installation

Use native arm64 macOS with Apple Command Line Tools installed. If needed, run
`xcode-select --install` and complete Apple's installer. From the repository:

```bash
bash install_bindcraft.sh
```

The installer includes PyRosetta, DAlphaBall, all 15 AlphaFold parameter sets,
and the native preparation and scoring tools. To install without PyRosetta:

```bash
bash install_bindcraft.sh --no-pyrosetta
```

Dependencies install inside the checkout under `.venv-metaljax`, `.local`, and
`params`. No Homebrew or global Conda environment is required. Rerun the default
installation command to add PyRosetta later. FFmpeg is needed only for saving
animations; the supplied small examples disable them.

## Run a small example

No environment activation is needed:

```bash
./run_apple_silicon.sh \
  --settings settings_target/ubiquitin_apple_small.json \
  --advanced settings_advanced/apple_silicon_ubiquitin_small.json \
  --max-attempts 1 --seed 42 --no-plots --no-animations --verbose
```

This example uses a 76-residue ubiquitin target, hotspot 44, and a 70-residue
binder. The small profile runs the 75/45/5/15 design schedule, generates one MPNN
sequence, uses one validation recycle, and disables extra beta optimization.
Default filter thresholds are unchanged.

Add `--no-pyrosetta` to the run command to use OpenMM relaxation. The launcher
selects this mode automatically when PyRosetta is not installed.

Outputs go to `runs/ubiquitin-apple-example/output`. Set a new `design_path` in
the target settings for an independent run. `--max-attempts 1` limits the job to
one started trajectory, including early rejections; a run can finish with zero
accepted designs.

A PDL1 example is also available:

```bash
./run_apple_silicon.sh \
  --settings settings_target/PDL1_apple_small.json \
  --advanced settings_advanced/apple_silicon_small.json \
  --max-attempts 1 --seed 42 --no-plots --no-animations --verbose
```

It uses a 65-residue binder and writes to `runs/pdl1-apple-example/output`.

## Tested configuration and performance

The pinned runtime uses Python 3.12, JAX/jaxlib 0.11.2, metaljax 0.11.9, and
OpenMM 8.6.1. Exact dependencies are recorded in the
[Apple Silicon requirements](installation/apple-silicon-requirements.txt) and
[PyRosetta requirements](installation/pyrosetta-apple-requirements.txt).

Installation and pipeline execution have been tested on an M2 MacBook Pro with
24 GiB unified memory (macOS 26.7) and an M5 Max MacBook Pro with 128 GB unified
memory (macOS 27.0.1). In a full installation on the M2, the small ubiquitin
example measured:

| Relaxation mode | Wall time | Peak process memory |
| --- | ---: | ---: |
| OpenMM / `--no-pyrosetta` | 39.7 min | 5.23 GiB |
| PyRosetta | 38.6 min | 5.51 GiB |

Both runs completed with zero accepted designs under the default filters.
These are single-run measurements of process physical footprint; larger targets
and higher recycle counts can require more memory.

## Backend behavior

The launcher uses Metal for AlphaFold and ProteinMPNN and OpenCL for OpenMM,
with no CPU fallback for those calculations. PyRosetta relaxation and native
preparation/scoring tools use the CPU. The Metal-specific alignment routines
leave model weights, losses, and default filters unchanged. Numerical equivalence
to CUDA has not been established.

In OpenMM mode, Rosetta-only metrics use compatibility placeholders rather than
measured Rosetta scores, as in FreeBindCraft's Linux bypass mode.

### Optional MetalJAX loop compilation fix

A source patch for MetalJAX 0.11.9 fixes nested AlphaFold loops that otherwise
silently fall back to eager execution. After installing the pinned environment:

```bash
.venv-metaljax/bin/python installation/build_metaljax_apple.py --install
```

This builds a local `0.11.9+fbc1` wheel using the release's unchanged MLX runtime.
The initial XLA build can take several minutes and needs network access and
Apple Command Line Tools. The patch preserves model settings and gradients;
compiled execution uses more memory on the measured 128 GB M5 Max fixture.
See the [source patch, build provenance and validation](installation/metaljax/README.md)
for reproduction commands, measurements and limitations. Smaller Macs have not
been validated with this patch. Set `FREEBINDCRAFT_METALJAX_WHEEL` to the built
wheel's absolute path to explicitly select it when rerunning the Apple installer.

## Verify the installation

Run the Metal alignment tests and native OpenMM/DSSP/ProteinMPNN check:

```bash
JAX_PLATFORMS=metal .venv-metaljax/bin/python -m unittest discover -s tests -v
JAX_PLATFORMS=metal .venv-metaljax/bin/python extras/check_apple_installation.py
```

The integration check uses a bundled PDL1 benchmark and writes its report and
structures to `runs/native-validation`. If PyRosetta is installed, check native
FastRelax, interface scoring, and DAlphaBall using that extracted benchmark:

```bash
JAX_PLATFORMS=metal .venv-metaljax/bin/python extras/check_pyrosetta_apple.py \
  --pdb runs/native-validation/reference.pdb
```

The PyRosetta report is written to `runs/pyrosetta-apple/report.json`.
