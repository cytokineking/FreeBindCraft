# Experimental BF16 attention epilogues

This is an opt-in native MetalJAX experiment on the frozen FreeBindCraft
`0.11.9+fbc1` backend. It combines the elementwise operations and row reductions
around Evoformer attention softmax, including the matching backward chain.
It retains the original QK, attention-value and derivative matrix products.
It is not a complete fused or FlashAttention implementation.

FreeBindCraft carries a source patch and reproducible maintainer build recipe,
not a vendored copy of MetalJAX. The recipe checks out fork commit
`678c58eaaafc33c5b31bcc6b6d0176dd078f6063` and uses the archived build inputs
in the existing `../release.json`. Local source, binaries, wheels and raw runs
are ignored build artifacts. No new public release is created here, and the
production installer still selects the existing checksummed `0.11.9+fbc1`
wheel. Distributing this experiment later would require a separate versioned
MetalJAX wheel and a reviewed dependency update, or an upstream release
containing the changes.

## Implementation and admission

`0001-bf16-attention-epilogues.patch` adds:

- A StableHLO recognizer for the precise forward softmax and its generated
  backward expression, with structural checks for reductions, constants,
  broadcasts, clipping and uses outside the matched region.
- Lowering and runtime operations for the matched forward/backward regions.
- Native Metal kernels with explicit BF16 intermediate rounding and the
  frozen MLX reduction order, plus a literal-MLX differential control.
- Differential tests requiring actual recognition and execution, including
  derivative, layout, rematerialization and unsupported-pattern checks.

`METALJAX_EVO_ATTN=1` enables recognition. It is off by default. The admitted
attention-logit tensor is BF16 `[B,H,Q,K]`, with `B >= 1`, `H` 4 or 8,
`Q >= 16`, `65 <= K <= 512`, at least 32 rows and at most `2^28` elements.
An overflow-safe element-count check precedes the kernel launch. Unsupported
shapes, dtypes or expression graphs retain the existing lowering.

The extended match combines the sequential mask and pair-bias additions,
the proven `[-1e8, 1e8]` clipping expression, maximum subtraction,
exponentiation, reduction and division. The additions still round separately
to BF16. Forward exponentials and the denominator remain available for AD.
The forward sum accumulates in FP32 and casts the denominator to BF16; the
backward reduction deliberately rounds in BF16, matching the frozen runtime.
Clipping derivatives retain the original comparisons, half derivatives at
ties, and multiplication order. A mathematically equivalent FP32 rewrite is
not numerically interchangeable with this graph.

Kernels are compiled and probed outside MLX tracing. A kernel preparation
failure declines recognition; a memory-governor refusal still propagates.
The existing compiled-loop recovery remains in place. The kernel cache is
protected during preparation and lookup.

`METALJAX_EVO_ATTN_LOGITS=0` limits the match to the softmax expression.
`METALJAX_EVO_ATTN_KERNEL=0` executes the recognized expression using literal
MLX operations for differential diagnosis. These controls do not alter the
production installer or select a dependency at runtime.

The patch does not edit FreeBindCraft's model, losses, optimizer, dropout/RNG,
recycling, checkpointing, BF16 policy, shared scientific functions or CUDA
code. Its numerical evidence is against the frozen Metal baseline, not CUDA.

## Building and checking

Use an isolated native arm64 macOS venv with JAX/jaxlib `0.11.2` and the same
FreeBindCraft scientific dependencies as the frozen baseline:

```bash
.venv-metaljax-attention/bin/python installation/build_metaljax_attention.py \
  --source-dir .local/src/metaljax-attention \
  --cache-dir .local/attention-build --jobs 12 --install
```

The builder verifies the exact source commit, patches and archived inputs,
refuses unrelated or modified source, and writes a build manifest alongside
the development wheel. `--install` affects only the selected venv. It builds
the plugin against the unchanged private MLX runtime; it does not rebuild MLX
from an unavailable stamped source revision. A local build is an experimental
artifact, not a validated production release.

```bash
JAX_PLATFORMS=metal,cpu METALJAX_MATMUL_PRECISION=high \
  .venv-metaljax-attention/bin/python -m pytest -q \
  .local/src/metaljax-attention/tests/test_evo_attention.py \
  .local/src/metaljax-attention/tests/test_nested_counted.py \
  .local/src/metaljax-attention/tests/test_control.py \
  .local/src/metaljax-attention/tests/test_msl_scan.py \
  .local/src/metaljax-attention/tests/test_command_buffer.py tests
```

For a fresh-process comparison, provide a checkout with model parameters and
the frozen baseline's first-recycle state. Baseline and candidate venvs must
have identical scientific packages:

```bash
.venv-metaljax-release/bin/python extras/run_metaljax_attention_comparison.py \
  --baseline-python .venv-metaljax-release/bin/python \
  --candidate-python .venv-metaljax-attention/bin/python \
  --assets /absolute/path/to/parameter-checkout \
  --previous /absolute/path/to/frozen-previous.npz \
  --out runs/metaljax-attention/comparison --order ABBABAAB --repeats 6
```

The driver starts sequential fresh processes, records process and aggregate
GPU activity continuously, waits for a bounded CPU quiet window before each
run, and stops only its own work if it detects sustained competing compute.
These checks cannot establish GPU exclusivity or control thermal state.
The memory monitor samples process-tree RSS, physical footprint, pressure,
compression and swap. Its sampled limit is not a strict allocation cap.

The probe times materialized forward and last-recycle sequence gradients
separately, including first calls. Private MLX dispatch-counter snapshots
occur outside each timed interval. Repeated calls in a process are correlated;
the fresh process is the sampling unit. Persistent JAX compilation caching
and diagnostic tracing are disabled during performance measurements.

## Measured coverage

On the M5 Max with 128 GiB unified memory, JAX/jaxlib `0.11.2`, BF16/high
precision, rematerialization and one recycle, the reversed fresh-process
series (`BAABABBA`, four processes per arm, six warmed calls per process)
gave these medians of process medians:

| Materialized call, 146 residues | Frozen fbc1 | Candidate | Time reduction |
| --- | ---: | ---: | ---: |
| Sequence gradient | 1.2187 s | 1.1069 s | 9.17% |
| Forward | 0.3895 s | 0.3414 s | 12.34% |

First calls did not improve: forward medians were 1.3487/1.3893 s and gradient
medians 3.5143/3.5845 s, baseline/candidate respectively. These are first calls
in fresh processes with persistent JAX compilation caching disabled, not cold
machine or driver-cache measurements. The earlier entire `ABBABAAB` series
is also recorded, including its slow candidate process concurrent with
background services. It is not substituted with selected favorable samples.

Every warmed gradient in the reversed series recorded 141,615/138,441 GPU
dispatches, baseline/candidate; forward recorded 56,530/55,412. Native
completed-command-buffer GPU busy time also decreased without profiling
fences. Candidate physical-footprint high-water marks in that series were
about 7.20 GB versus 7.38–7.44 GB for the baseline. These are process memory
measurements, not isolated GPU allocation totals.

The 185-residue PDL1 fixture passed the same 18-field and three-field warmed
forward byte comparisons. Its two warmed calls per arm took 1.969/1.720 s
for gradients and 0.587/0.498 s for forward. This is a larger-shape correctness
and memory check, not a replicated performance estimate. Observed physical
footprint high-water marks were approximately 11.85/10.91 GB. Both shape
checks reported normal sampled memory pressure and no observed swap growth.

The final packaged wheel passed 76 tests. The focused differential tests
cover 24 cases with two input variants and 338 saved arrays per arm: 324
finite arrays compare byte for byte, including signed zeros; 14 nonfinite-case
arrays compare values and NaN/Inf positions. They exercise supported widths,
eight heads, bias/mask broadcasts, reversed K/V layouts, rematerialized scans,
large/zero logits, clip ties and fully masked batches, plus seven explicit
unsupported cases. All 18 saved full-model fields and all three warmed
forward fields match the frozen baseline byte for byte at both real-model
sizes. This is tested Metal numerical coverage, not a general equivalence
proof or a CUDA comparison.

The unchanged prior fixed-work harness was then run in clean source
extractions in `ABBA` order, two fresh processes per arm. Design-stage times
were 223.93 and 218.49 s for the baseline, versus 198.23 and 196.25 s for the
candidate: median 221.21/197.24 s, or 10.84% less time. Internal whole-harness
totals were median 274.33/251.47 s; startup and downstream MPNN/OpenMM work
vary, so this whole-harness difference is not an exclusive attention cost.

All four processes completed 125 design backprop calls, 16 design forward
calls, three complex scorings/DSSP analyses, one MPNN sample, two complex and
two binder predictions, one filter evaluation and one model selection. All
140 printed iteration records, confidence checkpoints, actual generated
design PDB bytes and raw validation PDB bytes match. This compares serialized
outputs, not a dump of every internal optimizer array. The real design still
receives the same `LowConfidence` rejection. The benchmark uses its common
saved structure/sequence for downstream work after comparing the actual
design; it measures fixed work and reports zero accepted designs.

The trajectory runs used the existing 75/45/5/15 schedule, one recycle,
rematerialization, BF16/high precision and seed 121958. The in-memory
confidence-branch adaptation and fixture substitution are identical to the
prior harness (SHA-256
`e73d95f8c09f4e5dec907ac68fd2240e2220242873250e7c07eb93cf3462075a`).
Physical-footprint high-water marks were 9.50/9.76 GB for baseline processes
and 9.01/9.01 GB for candidates, with no observed swap growth or abnormal
sampled pressure. Process and aggregate GPU activity were recorded throughout;
desktop and thermal state remain uncontrolled.

## Diagnostic profiling and rejected experiments

`diagnostic-fences.patch` is optional and omitted from the ordinary candidate
wheel. Build it in a separate source directory and venv with `--profile`.
Set `FBC_ATTN_PROFILE=1`, select a label prefix with
`FBC_ATTN_PROFILE_FILTER`, and pass `--profile` to the probe. The barriers use
exported functions from the frozen private MLX runtime. They split command
buffers and change scheduling: their GPU times describe instrumented
dependency regions, not exact production cost shares. Overlapping intervals
are marked and excluded from timing aggregation. Never mix these timings
with the unfenced performance series.

For the main `[146,4,146,146]` attention region, the final filtered baseline
profile recorded the following dynamic frequencies and nonoverlapping region
GPU time in one warmed call:

| Instrumented region | Forward count / GPU time | Gradient count / GPU time |
| --- | ---: | ---: |
| QK-shaped BF16 batched dot | 104 / 18.37 ms | 208 / 32.57 ms |
| FP32 row sum | 108 / 18.06 ms | 216 / 25.68 ms |
| BF16 backward row sum | 0 | 216 / 14.58 ms |

The gradient includes rematerialized forward work. These regions can include
layout work and independent unprofiled operations; do not add them into a
production cost breakdown. All 208 intervals for the dQ/dK-shaped dot filter
overlapped, so only its frequency is retained and all its times are rejected.
The four final profiling runs retained the same 18 full-model fields and
three warmed forward fields byte for byte. Repeated, measured attention row
reductions motivated the bounded epilogue fusion while leaving matrix
products intact.

The retained local evidence includes unsuccessful approaches:

- Widths at or below 64 use a different MLX reduction order. A prototype
  changed small-row gradients; these widths are excluded.
- Recognition alone was insufficient: an early rewrite incorrectly skipped
  its own root, causing whole-tape fallback. Its timings are excluded. The
  final rewrite keeps the root, and tests and dispatch counts verify use.
- Early diagnostic intervals overlapped. A later ordering-chain experiment
  retained otherwise dead work and doubled some counts. Neither provides
  valid production cost attribution. The optional profiler carries no such
  ordering chain and explicitly marks overlaps.
- A completed timing series included a slow candidate process while macOS
  media-analysis/indexing services were active. The entire series is retained
  with its activity record; it is not silently reduced to favorable samples.

Validation results and exact artifact identities are recorded separately in
`validation.json`. Raw NPZ arrays, HLO, logs, wheels and process samples live
under the ignored `runs/metaljax-attention/` directory on the validation host.
