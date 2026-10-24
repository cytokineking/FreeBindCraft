# fbc3 combined defaults and Apple CLI opt-out

`metaljax==0.11.9+fbc3` is a distinct, **unpublished** candidate. With variables
unset, both reviewed attention epilogues and guarded slice-before-conversion
are enabled. Explicit `METALJAX_EVO_ATTN=0` and
`METALJAX_SLICE_BEFORE_CONVERT=0` retain independent diagnostic opt-outs.
Unsupported shapes/patterns, unproven layouts and unprofitable slices use the
existing lowering; kernel-preparation failures retain their existing fallback.
Scientific settings, precision, JAX/jaxlib, MLX, scheduling and governor
thresholds are unchanged.

The Apple launcher consumes `--no-attention-epilogues`, sets the attention
variable to `0` before the first Python invocation, and forwards the remaining
arguments in order with their quoting intact. This CLI flag overrides an
inherited `=1`; without it, an explicit environment `=0` is preserved. Launcher
help exposes the option. PyRosetta detection and the existing automatic
`--no-pyrosetta` path are preserved. Shared BindCraft scientific code is unchanged.

## Frozen candidate

| Item | Identity |
| --- | --- |
| Backend source | `d5962d9ba7340f8b7b893f9332851d34da31b7d7` |
| Parent, attention-off fbc2 | `bcd19ebca757b18f12b7eb963eebfacc428ba1fe` |
| Local tag | `v0.11.9-fbc3` |
| Tag object | `64be29d55e7813ecb76e4608dd8e403a133944d2` |
| Wheel SHA-256 | `0fa13c5def15e89024c624b81f5537492d269c6cb209602a5ba1e00baf25faa2` |
| Plugin SHA-256 | `25230c9a0a5b342c6a14313a06db0844f0e9cf3ce562ac5b8543a3edfd5449e4` |
| MLX runtime SHA-256 | `72dc8d9a4b1e07df2f6bc91383e0567acf3ce292398ae6ce4ed1d19fce0b1374` |

The only new backend changes relative to fbc2 are the attention default gate,
version/provenance and corresponding unset-default test. Recognition, shape,
layout, rounding, kernel preparation and the reviewed `3T + 5D < 3N` slice cost
rule are unchanged. Applying [the patch](0001-attention-default.patch) to the
exact fbc2 parent reproduces the committed source tree; the parent retains its
three-patch provenance from frozen fbc1. [source.json](source.json) records
the immutable chain; [build-manifest.json](build-manifest.json) records the
compiler/SDK, flags and frozen build/runtime inputs.

Two successful builds produce byte-identical wheel/plugin bytes with the
recorded paths/toolchain and verified reused disk cache. This is not a second
cold XLA build. The initial build attempt was guard-stopped when its detached
Bazel Java server became active; a subsequent attempt refused the prior
generated cache-path configuration. Both attempts are retained and produced
no wheel. After preserving/verifying those generated configs, foreground
`--batch` builds with an isolated output root completed under the unchanged
activity guard. No other process was stopped, and no governor or scheduling
limit was weakened. The build recipe's optional batch mode is recorded.

The exact wheel, local source archive/bundle, provenance and unchanged Apache,
MLX MIT and third-party notices remain in
`runs/production-fbc3/unpublished-assets/`.
[SHA256SUMS.candidate.txt](SHA256SUMS.candidate.txt) pins those local bytes.
No backend source/tag, release or asset has been published. Planned URLs must
not be treated as live.

## Fresh install, default and CLI checks

A new `.venv-metaljax-fbc3-clean` was populated through the ordinary Python
installer's default frozen selection with the exact verified wheel pre-seeded
in a private cache. No local-wheel override or future release URL was used.
The exact current public fbc1 wheel remains the control in its preserved fresh
environment; its plugin is
`005ab1302a7f49ec27e1447ab31db89e78b0e22793e6e14cf66bfe7773dce2fd`.
All scientific packages match except MetalJAX. JAX/jaxlib stay `0.11.2`; frozen
MLX SDK/runtime, stock troubleshooting pins and wheel license/notice bytes are
unchanged. Corrupt candidate transport bytes fail before uv and leave the
installed package inventory unchanged.

The relevant full backend/repository/installer/launcher suite passes:
**89 passed, 2 warnings in 47.49s**. Two launcher tests exercise seven real Bash/Python
scenarios for defaults, CLI consumption, inherited `=1` precedence, inherited
`=0`, argument order/paths with spaces, help and both PyRosetta paths. Only the
application fixture and site-package boundary are isolated; the launcher and
Python probe actually execute. The final launcher-only checks also pass after
hardening test-fixture path quoting and avoiding a fixed Python-invocation count.

Independent native workers explicitly remove each variable. Unset attention
matches explicit-on recognition/execution and matches opt-out finite bytes and
the existing nonfinite contract. All 24 attention cases, unsupported-shape
fallbacks, derivatives and kernel-execution checks pass. The 47 slice cases
still prove the independent unset default, cost/layout fallback, raw bits,
gradients, rematerialization, donation and governor refusals. Normal-root tests
leave both variables unset; runner-injected `=1` is not default proof.

Two separate probes also execute the **real Apple launcher and installed fbc3
backend**, recording controls before JAX/MetalJAX imports and then running the
actual AF2 benchmark in an isolated entry folder. The default leaves both
variables unset and recognizes attention kernels. With inherited attention
`=1`, the CLI opt-out reaches that pre-import point as `=0` and has no attention
admission. Slice rewriting remains admitted in both probes (158 forward /
153 gradient conversions), and all 24 saved fields match the public control.
These debug probes are excluded from performance data and read no dispatch
counters. They complement the launcher tests; they are not another campaign.

## Representative model comparisons

A is the exact public fbc1 wheel; B is the exact fresh fbc3 wheel. Default blocks
keep A's current released defaults and B's new combined defaults. Each block
uses ABBA fresh processes, two per arm and six correlated warmed calls per
process. Both optimization variables are unset for the new defaults. The
separate attention-off block explicitly sets attention `=0` in both arms,
leaving slice unset; the CLI's equivalent pre-import setting is proved above.
Every process matches all 18 model/sequence-gradient, three forward and three
previous-state fields byte for byte, with the same inputs, RNG/dropout,
previous state, precision, rematerialization and recycling.

| Configuration | Gradient seconds A → B | Change | Forward seconds A → B | Change |
| --- | ---: | ---: | ---: | ---: |
| 146 defaults | 1.235 → 1.047 | -15.2% | 0.394 → 0.296 | -24.7% |
| 185 defaults | 2.351 → 2.013 | -14.4% | 0.684 → 0.542 | -20.7% |
| 169 defaults | 1.853 → 1.589 | -14.3% | 0.562 → 0.420 | -25.2% |
| 146 attention off | 1.336 → 1.271 | -4.8% | 0.423 → 0.373 | -12.0% |
| model169-repeat | 1.679 → 1.435 | -14.5% | 0.505 → 0.375 | -25.9% |
| model185-repeat | 2.083 → 1.737 | -16.6% | 0.611 → 0.445 | -27.2% |
| model169 all valid | 1.703 → 1.495 | -12.2% | 0.517 → 0.387 | -25.1% |
| model185 all valid | 2.083 → 1.746 | -16.2% | 0.611 → 0.449 | -26.5% |

Individual warmed calls, first-call times, process medians and memory are
retained in the receipt. All valid observations remain included; any repeats
are reported separately and pooled without removing slower valid processes.
The initial 185-residue block's later B/A pair was about 30% slower than its
earlier pair; the 169 baseline processes differed by about 15%. Both blocks
were therefore repeated once in opposite-order BAAB processes. Remaining
timing variation is reported without tuning to the fastest observations.
Normal timing uses no diagnostic hooks/private dispatch reads or persistent
JAX compilation cache. Native precision, governor and 800-operation/512-MiB
scheduling with the 1-GiB loop window are left at their unchanged defaults.
An external 24-GiB process-tree RSS bound also applies.

## New combined-default fixed-work trajectories

Four fresh ABBA processes compare exact public fbc1 against exact fbc3 with
**both optimization variables unset**. The preserved scientific snapshot is
`cf136d1033790eef7b866bf07e92ed82b5572003`; the unchanged harness SHA-256 is
`e73d95f8c09f4e5dec907ac68fd2240e2220242873250e7c07eb93cf3462075a`.
The fixed 75/45/5/15 schedule uses `no-pyrosetta` OpenMM/OpenCL relaxation and
shared downstream fixtures. Configuration, native/model assets and package
identities are pinned equally.

Every actual generated trajectory PDB and all 140 metric records match byte
for byte **before fixture substitution**. Work counts agree: 125 design
backprops, 16 design forwards including semigreedy initialization, one MPNN
sample, two complex forwards, two binder forwards, three scored/DSSP structures,
and one filter/selection.

| Process | Wheel | Design seconds | Complete harness seconds | Peak process footprint GiB | Peak tree RSS GiB |
| --- | --- | ---: | ---: | ---: | ---: |
| a1 | fbc1 | 242.399 | 298.542 | 8.853 | 1.819 |
| b1 | fbc3 | 193.335 | 250.657 | 7.900 | 1.843 |
| b2 | fbc3 | 193.365 | 247.542 | 7.894 | 1.835 |
| a2 | fbc1 | 225.624 | 280.337 | 8.855 | 1.833 |

Median design time changes **234.011 → 193.350 seconds
(-17.4%)**; complete fixed-work harness time changes
**289.439 → 249.100 seconds
(-13.9%)**. Maximum lifetime process footprint
across each arm is **8.855 →
7.900 GiB**. The receipt retains
every stage, monitor wall time, sampled/lifetime footprint and RSS, system
pressure, compression, swap and available-memory observations.
The candidate design times agree within 0.02%; baseline design spread is
7.4%. No material trajectory outlier was identified, and all four processes
remain included. All four guards/monitors exited successfully, with normal
pressure and no swap growth.

- a1: pressure samples `{'1': 519}`; swap growth 0.000 GiB; minimum available RAM 73.38 GiB.
- b1: pressure samples `{'1': 437}`; swap growth 0.000 GiB; minimum available RAM 73.48 GiB.
- b2: pressure samples `{'1': 432}`; swap growth 0.000 GiB; minimum available RAM 73.53 GiB.
- a2: pressure samples `{'1': 486}`; swap growth 0.000 GiB; minimum available RAM 72.11 GiB.

These are descriptive measurements on one shared 128-GiB M5 Max, with two
trajectory processes per arm and uncontrolled thermals/desktop activity.
The guard cannot establish GPU exclusivity. Physical footprint includes
CPU/native/Metal allocations; system memory statistics include other apps.
The harness has **zero accepted designs**, a LowConfidence actual design,
forced fixed-stage work and shared downstream fixtures. This is not
accepted-design throughput or evidence of better scientific quality. MPNN
samples and relaxed outputs may vary even with equal work.

## Evidence and release boundary

[validation.json](validation.json) contains all normal and opt-out observations,
CLI pre-import/kernel evidence, actual PDB/record/work checks, every stage and
memory result, retained failures, source/build/install provenance and raw hashes.
The two failed build attempts above remain in the ledger; successful trajectory
guards and monitor results are retained per process. Raw files remain under
`runs/production-fbc3/`.

The original [attention-off fbc2 report](../README.md), its immutable source,
tag, wheel/archive/bundle, manifests and results remain unchanged. Its 666
evidence entries still verify, resolving changed tracked runner/build inputs
from snapshots frozen at `793827c` before this work. The v4/v5 ledgers also
retain all 865 entries each. Historical opt-in attention/development timings
and attention-off fbc2 results have not been relabeled as fbc3 measurements.

This is a fresh **Python dependency/cache installation**, reusing read-only
model/native assets. The real empty-cache public-download test is still
pending separately approved backend source/tag and exact-asset publication.
After that approval, publish the retained bytes without rebuilding, verify
the real public downloads/default install, record publication status, and
obtain merge approval. This PR remains draft/unmerged with no auto-merge;
no backend source/tag, release or artifact publication is part of this task.
