# Experimental submission scheduling and weight-cast profiling

These tools measure scheduling on the existing, explicitly enabled
[attention experiment](../attention/README.md). There is no new runtime
optimization or default setting in this change. The production installer keeps
the frozen `0.11.9+fbc1` wheel; this experiment uses the reviewed
`0.11.9+fbc3attndev2` wheel, JAX/jaxlib 0.11.2 and high matmul precision.

Reducing `MLX_MAX_OPS_PER_BUFFER` from 800 to 100 reduced memory use on the
measured M5 Max/128 GiB machine. Short probes improved, but the complete
fixed-work design stage did not establish a speedup. Keep 100 as a documented
opt-in experiment, with the 512 MiB buffer byte threshold, 1 GiB loop window
and 16 GiB memory governor retained. Results do not establish suitable defaults
for other Apple chips, memory capacities, sequences or shapes.

## Measurements and decision

Ordinary performance runs used fresh sequential processes without tracing,
HLO dumping, inventory collection or loop diagnostics. The 146-residue series
used ABBA followed by BAAB, four processes per arm and six warmed calls per
process. Values below are medians of process medians; repeated calls within
one process are correlated.

| Measurement | 800 operations | 100 operations | Reduction |
| --- | ---: | ---: | ---: |
| Warmed gradient | 1.1246 s | 1.0928 s | 2.82% |
| Warmed forward | 0.3452 s | 0.3296 s | 4.50% |
| First gradient | 3.7503 s | 3.6995 s | 1.36% |
| First forward | 1.4410 s | 1.4056 s | 2.46% |
| Process physical-footprint high-water mark | 7.202 GB | 4.995 GB | 30.64% |
| Fixed-work design stage, two processes per arm | 208.12 s | 207.37 s | 0.36% |
| Trajectory process footprint high-water mark | 9.019 GB | 6.953 GB | 22.91% |

GB means decimal GB, and physical footprint includes the process rather than
isolated GPU allocations. The trajectory times were 207.15/209.09 s at 800
and 209.71/205.03 s at 100. Their spread exceeds the difference between arms.
Neither this small sample nor uncontrolled desktop/thermal conditions support
a full-design speedup claim. All these runs had normal sampled memory
pressure and no observed swap growth.

The larger 185-residue correctness check used one process per arm: warmed
gradient 1.7786/1.7278 s, forward 0.5187/0.4859 s, and footprint 10.950/7.894
GB, baseline/candidate. It is not a replicated performance estimate.

The bounded screen also tested 400 and 200 operations, 4 and 8 GiB loop
windows, and transient-only byte accounting. The window/accounting screen
had appreciable baseline drift and no convincing improvement; the 8 GiB
window increased the observed footprint to 10.256 GB. These were not promoted
or combined into a default. The native gradient dispatch count remained
138,441 at 146 residues; 100 operations increased command buffers from 309 to
1,442–1,446. Smaller buffers reduced some gaps while adding submission work.
Summed queue residence times are overlapping intervals and are not wall time.

The [validation receipt](validation.json) retains every process's first/warm
timings, identities, numerical checks and source-file hashes. An initial
baseline overlapped a disassembler and was excluded. An initial trajectory
series stopped during its second process when the contention guard detected
`corespotlightd`; that incomplete series is retained locally and excluded as
a whole. The subsequent complete ABBA series supplies the trajectory numbers.

## Correctness and limits

Independent review compared all 18 saved gradient/model fields, three warmed
forward fields and three previous-recycle fields byte for byte. This passed
every scheduling configuration at 146 residues, both arms at 185 residues,
and all retained diagnostic modes against their corresponding baseline.
The memory governor refused a deliberately undersized 256 MiB budget with
`RESOURCE_EXHAUSTED`; the process exited unsuccessfully without eager retry.
The reviewed wheel with the 100-operation cap passed 85 focused backend and
repository tests, including five new CPU-only tests for profiling analysis.

The unchanged fixed-work harness ran the 75/45/5/15 schedule with seed 121958,
BF16, rematerialization, one recycle and high precision. Every extraction
matched all 122 tracked files at source commit
`cf136d1033790eef7b866bf07e92ed82b5572003`. All four runs completed 125 design
backprop calls, 16 design forward calls, three complex scorings/DSSP analyses,
one MPNN sample, two complex and two binder predictions, one filter evaluation
and one model selection. All 140 printed iteration records, confidence
checkpoints, actual generated design PDB bytes and unrelaxed validation PDB
bytes matched. This does not dump every internal optimizer array.

The real design receives `LowConfidence` and zero accepted designs. After the
actual design, the common harness uses a saved structure/sequence for fixed
downstream work. MPNN sampled sequences and relaxed/selected PDB bytes vary
across processes; no equality or accepted-design throughput claim is made for
those outputs. The harness SHA-256 is
`e73d95f8c09f4e5dec907ac68fd2240e2220242873250e7c07eb93cf3462075a`.

## Next backend optimization: slice weights before conversion

The best next experiment is a narrow MetalJAX lowering rewrite that selects
the current Evoformer layer from an FP32 weight stack **before** converting
that slice to BF16. The captured graph repeatedly converts an entire
48-layer stack and then selects one layer. For admitted cases, elementwise
conversion and selection can commute without changing BF16 rounding.
Preserve slice bounds/clamping and layouts, check every consumer, and decline
unsupported graphs. Do not introduce persistent activation or weight caches
as part of this first experiment.

The provenance analyzer follows the captured JAX argument map through
unchanged loop carries and identity barriers. It proves 61 weight-cast sites
in forward execution, each repeated 48 times. The gradient has 61 once-only
sites and 59 sites repeated 48 times. Five large element-count groups agree
with the actual Metal FP32-to-BF16 copy counts:

| Elements in full stack | Forward copies | Gradient copies |
| ---: | ---: | ---: |
| 12,582,912 | 96 | 98 |
| 6,291,456 | 48 | 49 |
| 3,145,728 | 576 | 588 |
| 1,572,864 | 192 | 196 |
| 786,432 | 672 | 686 |

These five groups account for 24.91 GB of logical reads plus writes per
forward and 25.43 GB per gradient. Those are element-access counts, not
measured DRAM traffic. Full-stack casts across all proven source sites imply
25.24/25.76 GB; the five-group trace match is the directly corroborated subset.
Selecting one layer first could remove most of this specific redundant work,
but it does not imply a 48-fold improvement for the entire call.

This ranks ahead of broader layout/donation and contraction changes because
it has a concrete, repeated redundancy and a narrow equivalence condition.
Layout/dense copies remain a second lead (13.36/56.00 GB decoded logical
traffic in forward/gradient); dynamic copies account for 0.93/4.45 GB, despite
their high counts. Gradient counters already show 1,831 donated slice updates
and 680 copied updates. Removing more copies requires an alias/lifetime proof.
Contractions remain a third lead pending exclusive kernel timing; 3,084/9,183
contractions have decoded parameters, with 159/305 left unclassified.

No reliable exclusive cost share is established. Cast-containing command
buffers cover 296.6 of 302.7 ms of forward GPU span and 953.3 of 997.1 ms of
gradient GPU span, but they also contain other families. These loose bounds
overlap and cannot be added or interpreted as cast time. CPU time inside the
original Metal dispatch API is about 1.2/1.8 ms for casts, excluding other
encoding work. Loop diagnostics observed 48 waits taking 156–160 ms in forward
and 97 waits taking 566–572 ms in gradients, including readiness queries;
those waits overlap GPU execution and are not removable overhead estimates.

Consequently there is no defensible numerical Amdahl prediction yet. If an
exclusive fraction `p` is established and this particular work improves by
`s`, the ceiling is `1 / (1 - p + p/s)`. The next experiment must measure
normal unfenced calls and the same fixed-work trajectory, preserving exact
model/gradient outputs, dispatch accounting, memory refusals and donation
behavior before any default or release change.

## Reproducing the experiments

Build the attention experiment in an isolated venv first, following its
README. The scheduling runner checks the exact reviewed plugin/MLX hashes
and JAX versions; it intentionally rejects unreviewed local builds.
Supply model assets, the saved first-recycle state and an existing baseline
output prefix containing `.npz`, `-warm-forward.npz` and `-previous.npz`.

```bash
.venv-metaljax-scheduling/bin/python extras/run_metaljax_scheduling.py \
  --python .venv-metaljax-scheduling/bin/python \
  --configs installation/metaljax/scheduling/configs.json \
  --order base,ops100,ops100,base,ops100,base,base,ops100 --repeats 6 \
  --assets /absolute/path/to/parameter-checkout \
  --previous /absolute/path/to/frozen-previous.npz \
  --reference /absolute/path/to/baseline/enabled \
  --out runs/scheduling/new-balanced-run
```

Use `base,ops400,ops200,ops100` for the operation screen and
`base,flight4,flight8,base,transient,base` for the window/accounting screen.
The runner uses a 16 GiB governor and sampled 24 GiB process-tree RSS limit,
requires 32 GiB available before the 8 GiB window, waits for a bounded quiet
window, and stops only its own work on sustained competing compute. Samples
cannot prove GPU exclusivity or enforce a hard physical-footprint limit.

`installation/build_metaljax_scheduling.py` builds optional CPU loop timing
diagnostics from the same pinned source/SDK and the attention patch. Use a
separate checkout-local venv/source/cache and `--install`; its manifest is
required by the runner's `--loop-profile --profile-manifest ...` mode. The
diagnostic wheel is `0.11.9+fbc4schedprofile1` and is not a release artifact.

For passive Metal tracing, compile the hook against the exact archived SDK
and MLX dylib. For example, after adapting these two checkout-local paths:

```bash
schedule_sdk=.local/src/metaljax-attention/src/metaljax/lib/mlx
schedule_lib=.venv-metaljax-scheduling/lib/python3.12/site-packages/metaljax/lib/mlx/lib
xcrun clang++ -std=c++20 -O2 -dynamiclib -fno-objc-arc \
  -Wno-deprecated-declarations -I "$schedule_sdk/include" \
  extras/metaljax_submission_trace.mm "$schedule_lib/libmlx_metaljax.dylib" \
  -framework Foundation -framework Metal \
  -Wl,-rpath,"$PWD/$schedule_lib" -o /tmp/scheduling-trace.dylib
```

Pass `--trace-library /tmp/scheduling-trace.dylib`; add `--copy-trace
--reflection --inventory --dump-hlo` for shape, stride and provenance data.
These modes are diagnostic and excluded from ordinary performance claims.
The hook reads a private command-encoder layout only after verifying the
frozen MLX hash, adds no GPU commands or barriers, and requires all native
buffer/dispatch counts to reconcile. Reflection and host-side hooks still
perturb execution. Their overhead has not been isolated with a balanced
on/off series, so reported trace times are diagnostic observations only.

```bash
.venv-metaljax-scheduling/bin/python extras/analyze_metaljax_submission.py \
  --trace runs/profile/base-1-trace.jsonl --out runs/profile/analysis.json
.venv-metaljax-scheduling/bin/python extras/analyze_metaljax_weight_casts.py \
  --hlo runs/profile/base-1-gradient.mlir \
  --inputs runs/profile/base-1-weights-grad_fn-inputs.json \
  --out runs/profile/gradient-weight-casts.json
```

The lineage tool validates argument indices, shapes and dtypes, handles
individual results of a multi-result identity barrier, and leaves arithmetic
or modified-loop-carry provenance unknown. Counts are source-level predictions
before MLX fusion and require trace corroboration. Raw runs, model parameters,
traces, SDKs, development wheels and source checkouts remain local ignored
artifacts; the validation receipt records their identities rather than
publishing model data or binaries.
