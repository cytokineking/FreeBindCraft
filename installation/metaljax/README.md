# MetalJAX static nested-loop repair

This patch fixes the compiled-loop fallbacks observed in FreeBindCraft's Apple
AlphaFold workload. It changes only MetalJAX's Metal PJRT lowering/runtime and
adds backend regression tests. FreeBindCraft's shared design functions, CUDA
installation, model, precision, recycling, checkpointing and optimizer are
unchanged. The Apple installer now downloads the frozen FreeBindCraft release
wheel by default; no user build is needed.
On the measured 146-residue M5 Max fixture, warmed gradients improved 1.71×
and the full design stage improved 1.73×, with higher memory use.

The separate [attention epilogue experiment](attention/README.md) requires an
explicit development build and opt-in flag. The release described here remains
the default Apple dependency.

## Default installation and troubleshooting

```bash
bash install_bindcraft.sh
# Add --no-pyrosetta to omit native Rosetta.
```

This selects `metaljax==0.11.9+fbc1`, with JAX/jaxlib `0.11.2`. The exact
versioned wheel URL and expected SHA-256 are in `release.json` and
`../apple-silicon-requirements.txt`. `install_apple_python.py` verifies the
artifact before syncing Python dependencies, and prints the installed versions.
The default requirements also contain the URL's SHA-256 fragment for direct uv
use. There is no moving branch or latest-release lookup during installation.

To troubleshoot with pinned stock upstream MetalJAX:

```bash
FREEBINDCRAFT_METALJAX=stock bash install_bindcraft.sh
```

Rerun the ordinary command with that variable unset to restore the frozen
release. The selected wheel participates in a single dependency sync, so fresh
installs, reruns and mode switches do not temporarily install another backend.
Existing explicit `FREEBINDCRAFT_METALJAX_WHEEL=/absolute/path.whl` overrides take
precedence over `stock` and the default. Local override bytes are chosen by the
user and are outside the published release's expected-checksum guarantee.
Unknown backend modes and missing overrides fail. A download/hash failure occurs
before uv changes the target venv; a damaged cache is re-downloaded and verified.

This selection happens only during installation. The backend retains its
existing per-loop recovery. Memory-governor refusals propagate as errors and
are not automatically retried eagerly. There are no dependency swaps or
whole-job restarts. CUDA/Linux dependencies and shared scientific logic are
unchanged.

## Frozen fork and release

- Fork: [cytokineking/metaljax](https://github.com/cytokineking/metaljax), verified
  under the same owner as FreeBindCraft.
- Source commit: `678c58eaaafc33c5b31bcc6b6d0176dd078f6063`.
- Tag/release: [v0.11.9-fbc1](https://github.com/cytokineking/metaljax/releases/tag/v0.11.9-fbc1).
- [Prebuilt wheel](https://github.com/cytokineking/metaljax/releases/download/v0.11.9-fbc1/metaljax-0.11.9%2Bfbc1-py3-none-macosx_14_0_arm64.whl)
  SHA-256: `8b3ec3bd4a71e57d7c77a0f3f8b1ee9fefe28b4ca81c0b325dc6d2975938d079`.
- [Canonical release checksums](https://github.com/cytokineking/metaljax/releases/download/v0.11.9-fbc1/SHA256SUMS.release.txt)
  and [build provenance](https://github.com/cytokineking/metaljax/releases/download/v0.11.9-fbc1/build-provenance.json).

The fork commits the reviewed fix and regression tests directly on the frozen
upstream base. Its default branch is `freebindcraft-0.11.9-fbc1`; installation
uses the artifact and source commit pins, not that branch. Upstream authorship,
Apache-2.0 licensing and bundled MLX MIT notices are retained in source and wheel.
FreeBindCraft carries the reviewed patch as an audit artifact; the builder does
**not** apply it again to the fork.

The release archives the exact committed MetalJAX source, checksum-pinned XLA
archive and Bazel 8.7.0 executable, staged reconstructed MLX SDK headers and
unchanged runtime binaries, recipe and provenance. The source archive contains
the four SDK-header patches. `build-provenance.json` records every compiled
source hash, SDK file hash, compiler/SDK versions, runtime binary hashes, flags,
JAX pins and the relationship to the reviewed build. The required private MLX
runtime inputs are archived; this does not claim to build MLX from its
unavailable stamped commit. XLA's transitive Bazel downloads remain pinned by
its build rules, rather than all being mirrored in this release.

The released wheel adds fork URLs and license notices to the original reviewed
wheel. All seven runtime payload files, including the plugin, Python loaders,
MLX dylibs and metallib, are byte-identical. The updated recipe, using the
frozen source and archived inputs on the recorded compiler/SDK, reproduced the
new wheel and native plugin hashes exactly. That rebuild used the existing
Bazel action cache after verifying the source and inputs; it was not a second
cold XLA build. A fresh public fork checkout also matched all 74 recorded
compiled source files. Different compilers/SDKs/paths may produce different
binaries; a local rebuild is not automatically a validated release.

All uploaded/downloaded asset bytes were checked. Use the canonical
`SHA256SUMS.release.txt` for downloaded filenames: the initial build-directory
manifest names `.gitignore`, which GitHub renames `default.gitignore`. The
canonical manifest was added under a new name; no published artifact or tag
was overwritten.

Future upgrades are deliberate: review a concrete source commit, create a new
immutable version/tag/asset identity, record and preserve inputs, run compiler,
numerical and memory checks, then update the installer URL/hash and compatible
JAX/jaxlib pins in a reviewed FreeBindCraft change. No automatic upstream
merges, dependency upgrades or runtime switching are enabled. Existing tags and
published bytes must not be replaced.

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

Normal installation uses the prebuilt wheel. Maintainers can reproduce from the
frozen fork and the release's archived build inputs after installing the pinned
Apple environment:

```bash
.venv-metaljax/bin/python installation/build_metaljax_apple.py
# Optional: --install installs the result only into the selected venv.
```

The helper checks out `678c58e` by full commit hash, checks both JAX pins,
verifies every archived SDK/XLA file on reruns, rejects changed source/build
configuration, uses a pinned hatchling, and prints the wheel and manifest paths.
The first cold native build can take several minutes. Inputs, source checkouts
and caches stay inside the checkout by default. The release also includes a
standalone recipe archive with this script and `release.json`; unpack it in a
working directory with uv on PATH, and pass
`--python /absolute/path/to/the/pinned/venv/bin/python`.

Only the normal Apple Command Line Tools are needed: the original metallib is
reused, so this build does not need the separate Metal shader toolchain.
A maintainer can explicitly install a local rebuild on a subsequent installer
run with `FREEBINDCRAFT_METALJAX_WHEEL=/absolute/path/to/the/built.whl`.

## Reproduce and test

The reviewed development build and the final frozen-recipe rebuild used the
same source/cache locations (the fork now contains the fix):

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
or checkpointing policy was changed. The frozen build has not been validated on the physical 24 GiB M2. The
additional bounded PDL1 and controlled-budget checks below extend fixture
coverage on the M5 Max; they do not establish arbitrary target-size support. CUDA/shared design source is untouched; no CUDA execution or
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
stack from a new public checkout. The original development recipe produced two byte-identical reviewed wheels:
plugin `005ab1302a7f49ec27e1447ab31db89e78b0e22793e6e14cf66bfe7773dce2fd`,
wheel `6e4b9ca8ad6e5dbaec3611d7fd827a48bd3bbb7710ac12911644739083eac50f`.

## Published distribution validation

The actual public release wheel was downloaded without GitHub authentication,
checksum-verified and installed into a fresh isolated venv with all pinned Apple
Python dependencies (without PyRosetta). The installer reported
`0.11.9+fbc1`, JAX `0.11.2`, jaxlib `0.11.2`. Rerunning audited the existing
packages. Switching to stock, rerunning stock, a fresh stock install, local
override precedence (even over stock), switching back and corrupted-cache
recovery all passed. An isolated copy of the real installer pins, with a
deliberately wrong expected hash against the real public wheel URL, failed
before uv ran and left the installed backend/versions unchanged; partial files
were removed. Five focused installer tests cover these failure/selection
contracts. The full native-tool/PyRosetta/weights bootstrap was not repeated;
its Python dependency phase is the tested distribution path.

The installed release wheel passed **72 tests**: the prior 57 compiler/control
tests, 10 FreeBindCraft tests and five new installer tests. Two existing JAX
deprecation warnings remain. The sandbox initially denied Metal device access;
the successful GPU tests ran with the normal device-access exemption.

Raw logs, installed wheel provenance, distribution checks and bounded probe
outputs are in `runs/metaljax-release/`. The committed
`release-validation.json` records the final measurements and checks.

### Bounded larger-target and memory checks

The release installed from the public URL was also compared to the actual
stock opt-out installation on bundled PDL1 chain A: 115 target residues plus a
70-residue binder, **185 residues total**, hotspot 56. The model, BF16/high
precision, checkpointing, one recycle, production losses, dropout and seed 42
were unchanged. The sequence/key hashes agree, and both gradients used the
same frozen stock first-recycle state. All **18 saved output fields** and all
three first-recycle tensors were bit-identical; gradient relative L2 and
maximum absolute difference are zero. The released wheel also matched all
18 saved fields of the earlier 146-residue stock baseline.

Both PDL1 arms ran sequentially in fresh processes, with two warmed gradients,
a 16 GiB governor process budget, a 24 GiB external process-tree RSS guard and
a six-minute external limit. Neither guard fired. The governor's other memory
settings were retained. Process physical footprint, not isolated GPU memory,
peaked at **5.14 GiB stock / 10.43 GiB release**. Warm-gradient medians were
4.710 / 2.575 seconds in this short pair; these timings are corroborative, not
a new balanced performance estimate. No system swap growth was observed.

The same release fixture under a **4 GiB** process policy and 250 ms stall
setting completed through repeated reclamation/stalls, with all 18 fields
still bit-identical. Its lifetime peak was **9.78 GiB**. The policy acts at
sampling/admission checkpoints and can reclaim/stall after graph transients;
it is **not a strict process peak-memory cap**. This is a material limitation
when interpreting budget settings as protection on smaller machines.

A deliberately tiny **256 MiB** policy with a 250 ms stall setting produced
`RESOURCE_EXHAUSTED: metaljax out of memory at transfer`, exiting with status 1
and no result file. The external monitor did not stop the job. The refusal
propagated normally; there was no eager retry, dependency change or job restart.
This exercised the transfer refusal path, not every allocation failure or
possible driver failure. Neither imposed budget simulates a physical smaller
Mac. The 24 GiB M2, arbitrary larger targets and higher recycle counts remain
unvalidated with the frozen backend.

For the larger comparison, use separate stock/release venvs installed with
`installation/install_apple_python.py --python /path/to/venv/bin/python` and
`FREEBINDCRAFT_METALJAX=stock` for the stock venv. Run each arm sequentially:

```bash
JAX_PLATFORMS=metal METALJAX_MATMUL_PRECISION=high METALJAX_MEM_BUDGET_MB=16384 \
  /path/to/stock/venv/bin/python extras/benchmark_metaljax_loops.py \
  --assets /path/to/installed/FreeBindCraft --pdb example/PDL1.pdb \
  --hotspot 56 --binder-length 70 --out runs/my-release-check \
  --name large-stock --repeats 2
JAX_PLATFORMS=metal METALJAX_MATMUL_PRECISION=high METALJAX_MEM_BUDGET_MB=16384 \
  /path/to/release/venv/bin/python extras/benchmark_metaljax_loops.py \
  --assets /path/to/installed/FreeBindCraft --pdb example/PDL1.pdb \
  --hotspot 56 --binder-length 70 --out runs/my-release-check \
  --name large-fixed --repeats 2 \
  --previous runs/my-release-check/large-stock-previous.npz
```

For an independent guard, wrap each command with
`extras/monitor_macos_job.py --output-dir <fresh-monitor-directory> --interval 1
--max-tree-rss-gib 24 --max-hours 0.1 -- <command>`. For the controlled policy
checks, use 4096 or 256 for `METALJAX_MEM_BUDGET_MB`, set
`METALJAX_MEM_STALL_MS=250` and `METALJAX_MEMDBG=1`, and choose new output paths.
These settings are test inputs, not default installation or launcher changes.
