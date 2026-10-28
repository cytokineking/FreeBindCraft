# Apple Silicon (experimental)

FreeBindCraft runs natively on Apple Silicon using
[FreeBindCraft's frozen MetalJAX build](installation/metaljax/README.md) for AlphaFold design and
ProteinMPNN on the Apple GPU. Relaxation uses native PyRosetta on the CPU by
default, or OpenMM on the GPU through OpenCL with `--no-pyrosetta`.

## Installation

The default is the published, checksum-pinned `0.11.9+fbc3` release. See the
[installation guide and release verification](installation/metaljax/README.md)
for its exact source, artifacts and public-download checks.

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

The default selects the prebuilt `metaljax==0.11.9+fbc3` arm64 macOS wheel.
Its versioned release URL and SHA-256 are pinned, together with JAX/jaxlib
`0.11.2`. There is no user source-build step or automatic upstream upgrade.
The installer reports the selected backend and installed versions.

For troubleshooting, install the pinned stock upstream `0.11.9` backend:

```bash
FREEBINDCRAFT_METALJAX=stock bash install_bindcraft.sh
```

Rerun `bash install_bindcraft.sh` with that variable unset to restore the frozen
build. Both commands accept `--no-pyrosetta`. Reruns retain the selected version;
switching modes explicitly replaces it. This is an installation choice, not a
runtime fallback. A download/checksum failure stops installation before Python
packages change. An explicit `FREEBINDCRAFT_METALJAX_WHEEL=/absolute/path.whl`
developer override takes precedence over both modes; unset it to return to the
frozen build. Local overrides are user-selected artifacts, outside the
frozen-wheel checksum guarantee.

The fbc3 release enables slice-before-BF16-conversion when
`METALJAX_SLICE_BEFORE_CONVERT` is unset; set it to `0` for diagnostic opt-out.
Attention epilogues also default on; unsupported shapes use the existing
lowering. To opt out, add the Apple launcher flag:

```bash
./run_apple_silicon.sh --no-attention-epilogues \
  --settings settings_target/ubiquitin_apple_small.json \
  --advanced settings_advanced/apple_silicon_ubiquitin_small.json \
  --max-attempts 1 --seed 42 --no-plots --no-animations
```

The flag sets `METALJAX_EVO_ATTN=0` before Python starts and overrides an
inherited `=1`. An explicit environment `=0` is also preserved. Scientific
settings, precision, scheduling and the memory governor are unchanged. The
[production report](installation/metaljax/production/fbc3/README.md) records the
exact artifact's prepublication checks; earlier opt-in attention measurements
remain historical.

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

The runtime uses Python 3.12, JAX/jaxlib 0.11.2, metaljax 0.11.9+fbc3, and
OpenMM 8.6.1. Exact dependencies are recorded in the
[Apple Silicon requirements](installation/apple-silicon-requirements.txt) and
[PyRosetta requirements](installation/pyrosetta-apple-requirements.txt).

The earlier stock MetalJAX 0.11.9 installation and pipeline were tested on an M2 MacBook Pro with
24 GiB unified memory (macOS 26.7) and an M5 Max MacBook Pro with 128 GB unified
memory (macOS 27.0.1). In a full installation on the M2, the small ubiquitin
example measured:

| Relaxation mode | Wall time | Peak process memory |
| --- | ---: | ---: |
| OpenMM / `--no-pyrosetta` | 39.7 min | 5.23 GiB |
| PyRosetta | 38.6 min | 5.51 GiB |

Both runs completed with zero accepted designs under the default filters.
These are single-run measurements of process physical footprint; larger targets
and higher recycle counts can require more memory. These M2 measurements use
stock upstream MetalJAX, not the new default. The frozen build's validation is on
the 128 GiB M5 Max; support on a physical smaller-memory Mac is not established.

## Backend behavior

The launcher uses Metal for AlphaFold and ProteinMPNN and OpenCL for OpenMM,
with no CPU fallback for those calculations. PyRosetta relaxation and native
preparation/scoring tools use the CPU. The Metal-specific alignment routines
leave model weights, losses, and default filters unchanged. Numerical equivalence
to CUDA has not been established.

In OpenMM mode, Rosetta-only metrics use compatibility placeholders rather than
measured Rosetta scores, as in FreeBindCraft's Linux bypass mode.

### Frozen MetalJAX backend and memory

The default build fixes nested AlphaFold loops that otherwise silently retire
their compiled variants. It keeps the existing calculation, precision,
recycling, checkpointing and optimization behavior, and uses unchanged private
MLX runtime binaries. In the original fbc1 validation, the 146-residue fixture's
warmed gradients were 1.71× faster and its full design stage 1.73× faster on the M5 Max. Full-workflow
peak physical footprint increased from 5.49 to 8.85 GiB.

The fbc3 release adds the attention and slice defaults described above. Against
public fbc1, its fixed-work M5 Max comparison measured 17.4% less design time,
13.9% less complete-harness time, and peak process footprint of 7.90 GiB.
These small-sample results are not accepted-design throughput; see the
[validation report](installation/metaljax/production/fbc3/README.md).

Existing per-loop recovery and memory safeguards remain enabled. A memory
governor refusal propagates as an error; it does not guarantee eager recovery,
swap dependencies or restart a design. Its budget is checked at admission and
sampling points, and is not a strict peak-memory cap: the original fbc1 PDL1 test
reached 9.78 GiB transiently under a 4 GiB policy before reclaiming memory.
The original fbc1 PDL1 forward/gradient probe used 10.43 GiB versus 5.14 GiB on stock.
Use the explicit stock installation choice above when investigating
memory/performance differences.

See [release provenance, bounded larger-target and memory-budget validation,
and maintainer build instructions](installation/metaljax/README.md). Future
upgrades require a new reviewed source commit, versioned artifact and checksum;
installation never follows a moving upstream branch.

### Unchunked model preparation

The Apple installer also applies a pinned ColabDesign preparation patch. Metal
uses an unchunked configuration for both forward and gradient callables,
including models above the dependency's previous 384-residue cutoff. CPU/CUDA
preparation behavior and the published MetalJAX wheel are unchanged.

A memory preflight checks the attention buffer size and an empirical model
working-set estimate before preparing new callables. It accounts for verified
resident model memory during reuse and reports uncertain ownership explicitly.
Insufficient estimated headroom raises `MemoryError`; native runtime memory
and pressure guards remain active. Admission is not a peak-memory guarantee.
The empirical calibration comes from one 406-residue BF16 multimer fixture;
other model sizes, modes and devices need application validation.

See the [preparation patch and validation](installation/colabdesign-metal/README.md)
for scope, measurements and remaining draft qualification.

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
