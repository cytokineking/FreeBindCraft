# MetalJAX static nested-loop repair

This patch fixes the compiled-loop fallbacks observed in FreeBindCraft's Apple
AlphaFold workload. It changes only MetalJAX's Metal PJRT lowering/runtime and
adds backend regression tests. FreeBindCraft's shared design functions, CUDA
installation, model, precision, recycling, checkpointing and optimizer are
unchanged. The Apple installer can explicitly select the resulting local wheel.
On the measured 146-residue M5 Max fixture, warmed gradients improved 1.71×
and the full design stage improved 1.73×, with higher memory use.

## Cause and implementation

MetalJAX admits small static nested loops into an enclosing compiled body.
`WhileTraceable` proves their trip counts from the StableHLO constants. Region
lowering can hoist the initial counter into that enclosing body's captures;
MLX replaces **all** its arguments, including those hoisted constants, with
tracers. `Program::run_while` nevertheless called `item_int` on the counter
before checking `in_trace`. Every enclosing body then failed with:

```text
[eval] Attempting to eval an array during function transformations like compile or vmap is not allowed.
```

The patch carries the existing static-trip proof in tape attribute 13. During
tracing, the runtime unrolls exactly that count without reading a device value.
It preserves the 64-iteration unroll cap, all admission/byte gates, the eager
path's actual runtime-bound reads, and the existing recovery ladder. Unknown
trip counts still refuse tracing and retain the correct fallback. This does not
raise a trace budget or change MLX kernels.

Under `METALJAX_DEBUG=1`, a compiled-body failure now reports its exception,
call/evaluation phase, repeat count, body operation/slot counts and argument
count. The reviewed 48-layer forward body, for example, failed at call time
with 2,552 operations, 3,373 slots and 112 arguments. All nine observed failures
had the same tracer-evaluation cause. Disabling MLX fusion also retained them.

## Source and build provenance

- FreeBindCraft starting commit and freshly fetched `origin/master`:
  `c6f8994bf6deb455a5b8f449f36cd65915e6b307`; clean detached worktree initially,
  then branch `codex/metaljax-compiled-loops`.
- MetalJAX 0.11.9: `7546eb93aa60c6570efdca95836276d03fe70a94`.
- JAX/jaxlib 0.11.2; XLA `91888df6ce85102c30220e41d952065925e10886`.
- ColabDesign `8fd167defd1b0f98e29f00584ff4f3987e6935dc`.
- The wheel retains the release's **byte-identical** private MLX dylibs and
  metallib (MLX 0.32.0). `build_metaljax_apple.py` verifies their SHA-256 hashes.
- The release stamp names MLX commit `661b2e38`, absent from the published fork
  inspected here. For SDK headers only, the build uses published base
  `651c39cd3d26d0a62587ca639a52ea1efa4b8186` plus the release's checked-in
  patches 0004–0007. It does not rebuild or replace the MLX runtime.
- Bazel 8.7.0 and the XLA archive are checksum pinned. The helper replaces
  maintainer-specific source paths, sets the macOS target for both host tools
  and plugin, keeps caches local, and packages `metaljax==0.11.9+fbc1`.

After the ordinary Apple installation, build and install into this checkout's
venv (first source build downloads/builds XLA and can take several minutes):

```bash
.venv-metaljax/bin/python installation/build_metaljax_apple.py --install
```

The helper prints the wheel and `build-manifest.json` paths. To preserve that
wheel when running the Apple installer explicitly:

```bash
FREEBINDCRAFT_METALJAX_WHEEL=/absolute/path/to/metaljax-0.11.9+fbc1-py3-none-macosx_14_0_arm64.whl \
  bash install_apple_silicon.sh
```

The SDK only needs the normal Apple Command Line Tools: the release's metallib
is reused, so this repair does not require the separate Metal shader toolchain.
The helper refuses an unexpected source revision, changed prepared SDK sources,
unverified runtime/archive, or a Python environment outside the pinned venv.

## Reproduce and test

The actual build in this worktree used:

```bash
.venv-metaljax/bin/python installation/build_metaljax_apple.py \
  --source-dir .metaljax-src --cache-dir .local --jobs 16
```

The diagnostic-only build retained original execution and matched the release
bit for bit. Its captured log is `runs/metaljax-loop-fix/diagnostic-built.log`;
`diagnostic-built-gradient.mlir` is the complete failing gradient StableHLO.
`stack-diagnostic.log` independently shows the release-binary stack through
`item_int -> Program::run_while -> Program::interpret -> Program::compiled`.
The first attempt at that temporary exception logger failed; all validation
and timing results use normal, uninterposed processes.

The reduced regression is included in the patch at
`tests/test_nested_counted.py`. It disables generated MSL scans only inside its
isolated child processes, so it specifically exercises this counted-loop path.
It checks forward results and both input/weight gradients against CPU, replays
with changed input values, and asserts that compiled/chunked variants are not
retired. It also tests nonzero starts, zero/negative trips and dynamic bounds.
The numerical tolerance, `rtol=1e-4, atol=1e-5`, follows the existing long-scan
coverage; the original backend passes its numerical checks but fails its
compilation-retention assertion.

Install `pytest` into the chosen test venv if needed. These paths use the
`.metaljax-src` checkout from the explicit build command above.

```bash
JAX_PLATFORMS=metal,cpu .venv-metaljax/bin/python -m pytest -q \
  .metaljax-src/tests/test_nested_counted.py \
  .metaljax-src/tests/test_control.py \
  .metaljax-src/tests/test_msl_scan.py \
  .metaljax-src/tests/test_command_buffer.py
```

The focused run with the fixed plugin passed **57 tests**. FreeBindCraft’s
existing suite passed **10 tests** with the installed wheel. The interpreter-free
C++ runtime test also passed when run directly. Its first Bazel test invocation
could not access Metal because the persistent build server inherited the build
sandbox; this was a device-access failure, not a numerical/test assertion.
The packaged wheel was also checked to load exactly its own two private MLX
libraries, rather than the build tree’s copies. Existing tests cover
conditional/dynamic control flow, generated-kernel failure recovery, wide bodies
and command-buffer ordering. The new regression failed on the original backend
because variants were retired despite matching CPU numerically.

The fixed-input harness uses the unchanged production losses, BF16, `high`
matmul precision, checkpointing, one recycle, dropout and RNG key 42, model 1
multimer v3, and the 76+70-residue ubiquitin fixture. It saves/hashes the sequence,
RNG and previous state, materializes outputs before stopping each timer, and
reports first calls separately from warmed calls. For example:

```bash
JAX_PLATFORMS=metal METALJAX_MATMUL_PRECISION=high \
METALJAX_DEBUG=1 METALJAX_TIMING=1 \
.venv-metaljax/bin/python extras/benchmark_metaljax_loops.py \
  --assets /absolute/path/to/installed/FreeBindCraft \
  --out runs/metaljax-loop-fix --name baseline --repeats 6 --dump-hlo
# Use the baseline backend for that command, then select the patched wheel:
JAX_PLATFORMS=metal METALJAX_MATMUL_PRECISION=high \
METALJAX_DEBUG=1 METALJAX_TIMING=1 \
.venv-metaljax-fixed/bin/python extras/benchmark_metaljax_loops.py \
  --assets /absolute/path/to/installed/FreeBindCraft \
  --out runs/metaljax-loop-fix --name fixed --repeats 6 \
  --previous runs/metaljax-loop-fix/baseline-previous.npz
```

## Numerical validation

With identical sequence, RNG, previous recycling tensors and precision,
`baseline-frozen.npz` and `fixed-frozen.npz` are bit-identical for **all 18 saved
fields**: initial forward loss/coordinates/confidence, final loss/coordinates/
confidence, all individual loss terms and sequence gradient. The three initial
forward recycling tensors also match bit for bit. Sequence-gradient cosine is
1 within FP64 calculation rounding; relative L2 difference and maximum absolute
difference are zero. Loss is `13.94212818145752`; gradient norm is
`1.6165108680725098`. This validates this fixture; it is not a CUDA equivalence or
binder-success claim.

Performance measurements and full-trajectory validation are recorded in the
accompanying `validation.json` and the worktree's `runs/metaljax-loop-fix/`.
Raw logs, NPZ outputs, HLO, native binaries and wheels remain local build/test
artifacts rather than repository payloads.

## Measurements and limits

Measurements were sequential on the same M5 Max (128 GiB, macOS 27.0.1).
The final steady-state comparison used fresh processes in **A-B-B-A order**, six
warmed materialized VJPs per process (12 per backend), with the same frozen
previous-state file. Persistent JAX compilation caching was disabled. No other
test workload ran concurrently. The desktop was not quiesced; clocks and thermal
state were uncontrolled, and OS/driver caches were retained. These are small,
correlated samples, not a confidence interval or a production throughput claim.

| Fixed-input probe | Release 0.11.9 | Patched 0.11.9+fbc1 |
| --- | ---: | ---: |
| Warmed VJP median | 2.192 s | 1.283 s |
| Fresh-process medians | 2.177 / 2.194 s | 1.286 / 1.278 s |
| First forward calls | 1.665–1.776 s | 1.389–1.397 s |
| First gradient calls | 4.481–4.592 s | 3.617–3.640 s |
| Warmed dispatches, median | 293,027 | 141,615 |
| Warmed command buffers, median | 897 | 310 |
| Warmed eager flushes, median | 367 | 3 |
| Backend GPU busy, median | 1,585.4 ms | 1,129.9 ms |
| Backend GPU idle, median | 596.6 ms | 141.1 ms |
| Initial caught body failures, each process | 9 | 0 |
| Lifetime peak physical footprint | 3.77–3.93 GiB | 6.75 GiB |

First calls include compilation **and execution**; they are separated from the
warmed measurements, not presented as isolated compiler timings. All warmed
calls report zero new compilations. The baseline had permanently retired its
failed bodies by then. GPU times/dispatches are backend diagnostics, not an
independently captured Instruments trace.

The unchanged archived `benchmark-fixed-budget.py` (SHA-256
`e73d95f8c09f4e5dec907ac68fd2240e2220242873250e7c07eb93cf3462075a`) was
run once per backend, in two identical clean `c6f8994` source extractions with
read-only links to the same weights/native tools. Both used `high`, BF16, the
original 75/45/5/15 schedule, one recycle, checkpointing, seed 121958, no
PyRosetta, OpenMM/OpenCL and FASPR. Debug counters and the original 1-second
physical-memory monitor were enabled in both. The design stage includes first
compilation. Each completed 125 differentiated design calls, 16 design forward
calls, one MPNN sample, two complex and two binder predictions, three complete
relaxation/scoring/DSSP checks, filtering and best-model selection.

| Full fixed-work benchmark, one run each | Release | Patched |
| --- | ---: | ---: |
| Design stage | 404.068 s | 233.421 s |
| Complete process, including exit/monitor | 463.463 s | 295.854 s |
| Lifetime peak process physical footprint | 5.49 GiB | 8.85 GiB |
| Caught compiled-body failures across workflow | 12 | 0 |

This is **1.73×** for design and **1.57×** for the complete measured workflow.
All 140 printed design iteration records, the actual generated design PDB,
both unrelaxed complex predictions and both binder predictions are bit-identical.
The filter outcome and selected model also agree. The benchmark uses common
fixtures after design; it measures fixed work, not accepted-binder throughput.

The full runs' MPNN samples and relaxed structures/scores differ. The unchanged
harness constructs MPNN with `seed=None`; ColabDesign’s `Key` then uses Python’s
unseeded `random.randint`, despite the harness resetting NumPy’s seed. OpenMM
also uses unseeded Langevin MD shakes and `setVelocitiesToTemperature`. Thus the
full runs are not a same-RNG equivalence test of those stochastic stages.
A separate explicit-seed-42 MPNN probe matched **all six output fields bit for
bit**, including sequence, logits, score and decoding order. OpenMM is outside
the patched backend, and its input structures match exactly; its stochastic
relaxations were left unchanged.

The extra memory is a material tradeoff: compiled bodies retain more
intermediates than the previously flushing eager path. No memory/trace budget
or checkpointing policy was changed. The 24 GiB M2 and larger targets have not
been validated. CUDA/shared design source is untouched; no CUDA execution or
cross-backend numerical-equivalence claim is made.

To repeat the full benchmark with the copied assets in this worktree, choose a
fresh output source extraction for each arm (the harness refuses an existing
output directory), then run the same command with the appropriate Python and
root:

```bash
JAX_PLATFORMS=metal METALJAX_MATMUL_PRECISION=high METALJAX_DEBUG=1 \
.venv-metaljax-fixed/bin/python extras/monitor_macos_job.py \
  --output-dir runs/metaljax-loop-fix/trajectory-fixed/monitor --interval 1 -- \
  .venv-metaljax-fixed/bin/python runs/metaljax-loop-fix/benchmark-fixed-budget.py \
  --root runs/metaljax-loop-fix/trajectory-fixed \
  --fixture runs/metaljax-loop-fix/fixed-fixture --mode no-pyrosetta --backend metal
```

The executed command/environment records are in each monitor’s `metadata.json`
and `paired-commands.json`. `validation.json` contains all warmed samples,
per-process first calls/counters, memory peaks, fixture/source hashes and the
full stage timings. The source rebuild check reproduced the SDK header patch
stack from a new public checkout. Two wheel builds produced identical hashes:
plugin `005ab1302a7f49ec27e1447ab31db89e78b0e22793e6e14cf66bfe7773dce2fd`,
wheel `6e4b9ca8ad6e5dbaec3611d7fd827a48bd3bbb7710ac12911644739083eac50f`.
