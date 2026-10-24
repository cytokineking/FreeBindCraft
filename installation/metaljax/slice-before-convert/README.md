# Experimental slice before BF16 conversion

This opt-in MetalJAX compiler transformation replaces `slice(convert_f32_to_bf16(x))`
with `convert_f32_to_bf16(slice(x))` before native tape generation. It uses the existing
slice and conversion handlers. No new kernel, persistent weight cache, activation cache,
or scientific model change is involved. The production installer and dependency pins are
unchanged.

Admission requires static, positive, unencoded ranked shapes; a direct, shape-preserving
FP32 to BF16 conversion; proven row-contiguous source layout; and exclusively reducing
slice uses or proven slice-helper calls in the conversion's own block. Their aggregate
logical access cost must be strictly lower than the original conversion and slices,
using checked arithmetic as derived below. Both static and dynamic slices require a
row-contiguous window; static strides must be one. Dynamic slices retain their original
indices, sizes and clamp lowering. There is no hard-coded layer count, sequence length
or model name.

The source-layout proof starts at the frozen PJRT boundary: host inputs are gathered
into dense row-major arrays, and execution materializes non-row outputs before
publishing buffers. The proof follows reshapes, alias-only optimization barriers,
function arguments with every symbol caller checked, and loop carries whose body returns
the identical argument. Unknown or cyclic provenance, transposes, broadcasts, slices and
computed sources decline. Dynamic slices need the same window proof because loop
specialization can turn them into static views. Escaping full conversions, other dtypes,
cross-region uses, empty shapes and unknown patterns retain their original lowering.

These layout restrictions are necessary for numerical correctness: a wider prototype
passed the AF2 comparisons but changed native BF16 reduction rounding on a
transposed-source static slice by materializing a different layout. That case is
preserved as a negative-admission regression. Its timings are excluded from the final
results. The final suite also checks native BF16 reductions of row-contiguous static and
dynamic windows and transposed inputs through callees and optimization barriers. No
comparison tolerance was relaxed.

The actual model calls a helper around its slice. A supported helper has one block and
one return, reads the converted argument only in one slice, then returns that slice
through optional reshapes. Every other argument and operation must be scalar integer
index arithmetic from an explicit pure-op allowlist. Only admitted call sites are cloned
into the caller before the slice/conversion rewrite; the helper and other callers remain
intact. Arithmetic on the tensor, additional returned values, unknown operations,
non-scalar index arguments and aggregate full-stack selection decline the whole
conversion's rewrite. Helper names play no part in admission.

The conversion remains an ordinary BF16 conversion after each admitted slice. All
consumers are replaced together and the full conversion is erased only after its
complete use census. Unconverted source consumers remain intact. The transformed module
is serialized by PJRT and used by the first-execute recognizers, so plain and
attention-enabled tapes see the same transformation. AD/rematerialization still happen
in JAX before this backend transformation; the resulting forward and derivative graphs
are lowered independently.

`METALJAX_SLICE_BEFORE_CONVERT=1` enables the transformation in the development wheel.
It is disabled by default. Attention is independently enabled with
`METALJAX_EVO_ATTN=1`. There are no changes to the 800-operation command-buffer control,
the 512 MiB byte threshold, the 1 GiB loop window or memory governor.

## Reproduction

Start from FreeBindCraft `cf136d1033790eef7b866bf07e92ed82b5572003` plus this change.
Create two native arm64 Python 3.12.15 virtual environments inside the checkout, with
the identical pinned scientific dependencies in
`installation/apple-silicon-requirements.txt`, plus `pytest==8.3.5`. Use only verified
frozen assets from `installation/metaljax/release.json`; corrupted historical test
caches are invalid build inputs.

```bash
.venv-metaljax-baseline/bin/python installation/build_metaljax_attention.py \
  --python .venv-metaljax-baseline/bin/python \
  --source-dir .local/src/metaljax-baseline \
  --cache-dir .local/baseline-build --jobs 6 --install

.venv-metaljax-slice-v5/bin/python installation/build_metaljax_attention.py \
  --python .venv-metaljax-slice-v5/bin/python \
  --source-dir .local/src/metaljax-slice-v5 \
  --cache-dir .local/baseline-build --jobs 6 --slice-before-convert --install
```

The candidate applies the reviewed attention patch followed by
`0001-slice-before-convert.patch` to MetalJAX
`678c58eaaafc33c5b31bcc6b6d0176dd078f6063`. It builds `0.11.9+fbc4slicedev5`; the
control builds `0.11.9+fbc3attndev2`. Both use JAX/jaxlib 0.11.2, the checksummed frozen
MLX SDK, XLA revision `91888df6ce85102c30220e41d952065925e10886`, Bazel 8.7.0 and
hatchling 1.32.4. Each source directory contains its own `build-tools/attention-wheels/`
build manifest and wheel. The builder verifies source and input identities, refuses
unrelated changes, and can install only into a selected venv in this checkout.

Run tests behind the shared-Mac guard after builds have finished:

```bash
.venv-metaljax-baseline/bin/python extras/guard_metaljax_job.py \
  --enabled --out runs/slice-tests -- \
  .venv-metaljax-slice-v5/bin/python -m pytest -q \
  .local/src/metaljax-slice-v5/tests/test_slice_before_convert.py \
  .local/src/metaljax-slice-v5/tests/test_evo_attention.py \
  .local/src/metaljax-slice-v5/tests/test_nested_counted.py \
  .local/src/metaljax-slice-v5/tests/test_control.py \
  .local/src/metaljax-slice-v5/tests/test_msl_scan.py \
  .local/src/metaljax-slice-v5/tests/test_command_buffer.py \
  .local/src/metaljax-slice-v5/tests/test_donation.py tests
```

`extras/run_metaljax_slice_comparison.py` runs sequential fresh processes in ABBA order,
checking installed plugin/MLX hashes against both build manifests and requiring
identical scientific package versions. Provide `--assets`, `--previous`, `--reference`
(saved output prefix), both manifests and a fresh `--out`. The 185-residue fixture also
uses `--pdb example/PDL1.pdb --hotspot 56`. All 18 model/sequence-gradient fields, three
warmed-forward fields and three previous-state fields must match byte for byte. Normal
runs disable private dispatch reads and all hooks. First calls and warmed calls are
reported separately; repeated calls within a process are correlated.

The normal comparison also accepts `--attention 0` to apply the same disabled attention
setting to both arms. Select the new environment explicitly with `--candidate-python
.venv-metaljax-slice-v5/bin/python` and its v5 build manifest. The [v4 report at the
reviewed
commit](https://github.com/cytokineking/FreeBindCraft/blob/8c789beb8c3236281e2f62e21ad467d863a6aabf/installation/metaljax/slice-before-convert/README.md)
retains the original fixed-work trajectory recipe and receipts.

## Diagnostic boundaries

The historical 100-operation experiment is retired. Its two 138245-dispatch gradient
observations differed from the usual 138441, despite matching saved output archives. The
affected calls did not have individual output dumps. The old constant-count claim and
+/-20 runner guard must not be reused.

For separate diagnostics, build `extras/metaljax_dispatch_trace.mm` against the frozen
MLX SDK and run the comparison with `--diagnostic-trace`. For example, from the checkout
(using the built baseline SDK):

```bash
mkdir -p .local/diagnostics
xcrun clang++ -std=c++17 -O2 -dynamiclib -framework Foundation -framework Metal \
  -I .local/src/metaljax-baseline/src/metaljax/lib/mlx/include \
  extras/metaljax_dispatch_trace.mm \
  -L .local/src/metaljax-baseline/src/metaljax/lib/mlx/lib -lmlx_metaljax \
  -Wl,-rpath,"$PWD/.local/src/metaljax-baseline/src/metaljax/lib/mlx/lib" \
  -o .local/diagnostics/libdispatch_trace.dylib
```

`extras/diagnose_metaljax_dispatch.py` records immediate and delayed snapshots, saves
every warmed call's arrays, and waits for both native completion and the hook's own
logging callbacks. The callback wait corrects a limitation of the retired tracer's final
dump. It reads counters through an export in the tracer, identifying the exact linked
MLX image. Loading a byte-identical MLX copy by another pathname can create independent
counter globals; a checksum alone does not identify the active image. The trace refuses
an MLX checksum change; its encoder prefix is a private, checksum-specific ABI. The
analyzer requires direct Metal dispatch counts to match each native command-buffer count
and requires kernel totals to cover the submitted workload.
`check_metaljax_dispatch_accounting.py` additionally requires exact per-call agreement
with native committed/completed counter deltas and retained call arrays. There is no
expected constant dispatch count or tolerance window.

Diagnostic instrumentation and per-call dumps perturb host execution and cannot supply
normal performance timings. Containing-buffer intervals overlap across operation
families; logical copy bytes are neither measured DRAM bytes nor exclusive wall-time
shares. Never add overlapping host/GPU/wait intervals or infer a percentage speedup from
cast counts.

The local process/GPU guard waits for a bounded quiet window, records ongoing activity
and stops only its own work on sustained competing CPU compute. Aggregate GPU activity
cannot establish per-process exclusivity, and desktop activity and thermal state remain
uncontrolled. Memory reports distinguish process RSS/physical footprint from sampled
system pressure and swap; they do not isolate GPU allocations.

## Profitability follow-up: v5

Independent review found that the v4 element-count test admitted a near-full
runtime-index window despite the extra FP32 temporary. On the reviewed v4 wheel,
selecting 15 of 16 blocks was 35.7% slower in a fresh ABBA synthetic comparison, with
all bits equal. This was an operation regression, not an AF2 regression. The corrected
v5 source preserves all earlier layout, clamping, rounding, use-census,
rematerialization, gradient, donation and governor gates.

The frozen runtime's `kConvert` emits `mx::astype`; `kSlice` emits a static `mx::slice`
view; ordinary `kDynamicSlice` emits a materializing dynamic slice. For N full elements,
T the sum of every static window and D the sum of every dynamic window, the approximate
logical byte costs are:

- Before: `6N + 4D` (FP32 read/BF16 write plus BF16 dynamic copies).
- After: `6T + 14D` (static-window conversions plus FP32 dynamic copies and conversions).
- Strict savings: `3T + 5D < 3N`.

The same budget covers direct, helper and mixed uses together before any IR mutation. It
reserves strictly positive savings and divides before multiplying to check each spend.
Unrepresentable full budgets and aggregate costs decline. Build-time assertions exercise
the actual budget arithmetic at int64 limits without huge device allocations. Position
specialization can remove a dynamic copy, but admission conservatively charges for it.
These are logical accesses; they do not measure physical DRAM traffic or guarantee
hardware timing.

The v5 wheel and source were built in separate directories and environment. The v4
wheel, manifests, source and all 865 raw observations remain unchanged. A fresh patch
application reproduces every v5 source fingerprint. The 84-test backend/repository suite
passes, including added near-full rejection, small windows, strict equality and adjacent
boundaries, fanout, static versus dynamic, and mixed/helper admission cases. Both
admission and raw output bytes are checked; all existing conversion-edge, layout,
gradient and lifetime tests still pass.

Candidate v5 plugin SHA-256:
`fe556d0bba093ffe82cc53d0cce2b91b63c901ce4a7b6c59d144ea9dc78830ef`. Wheel SHA-256:
`01e6d7231fe41aa0bae2a8a5b6417053ce840d7164db9277c985cc0ff9713d34`.

### Fresh synthetic reproduction

The review fixture uses `(16,1024,1024)`, alternating runtime indices, three warmups and
twenty synchronized calls per fresh process. Both arms use the same v5 binary with the
rewrite off/on, attention disabled and original scheduling. Optimized HLO verifies that
the 15-block case now declines and the one-block case still rewrites. Every saved BF16
bit matches across arms. The first block's one-block control timings varied, so the full
ABBA block was repeated; every completed original observation remains included.

| Window | Block | Process medians, rewrite off → on (ms) | Change |
| --- | --- | ---: | ---: |
| 15 / 16 blocks | Initial ABBA | 0.9930 → 0.9659 | -2.7% |
| 15 / 16 blocks | Repeat ABBA | 0.8338 → 0.9705 | +16.4% |
| 15 / 16 blocks | All valid processes | 0.9023 → 0.9659 | +7.1% |
| 1 / 16 blocks | Initial ABBA | 0.6164 → 0.2622 | -57.5% |
| 1 / 16 blocks | Repeat ABBA | 0.6835 → 0.2657 | -61.1% |
| 1 / 16 blocks | All valid processes | 0.6835 → 0.2638 | -61.4% |
| 15 / 16 blocks, isolated | BAAB | 0.9496 → 0.9733 | +2.5% |

The unfavorable window has identical optimized HLO across arms in every block, with the
extra FP32 slice temporary absent. The mixed-window repeat still varied in elapsed time
despite that identical fallback, so the 15-block case was isolated before
positive-control work in an opposite-order BAAB block. Its off/on process ranges overlap
(off 0.856–1.043 ms; on 0.890–1.057 ms), with a +2.5% difference in process medians. The
small-window positive control remains eligible and faster in every candidate process.

The extra-copy regression is removed structurally, and isolation does not reproduce the
reviewed 35.7% penalty. These small shared-Mac samples do not establish exact zero
overhead or identify the cause of the remaining timing variation. Every block is
retained. Absolute times from different sessions are not directly comparable, and this
synthetic probe does not estimate AF2 throughput.

```bash
.venv-metaljax-baseline/bin/python extras/benchmark_metaljax_slice_profitability.py \
  --python .venv-metaljax-slice-v5/bin/python --out runs/profitability-check
```

### Fresh AF2 checks with v5

A is the preserved frozen attention-development control, B the v5 candidate. Every block
uses fresh ABBA processes, two per arm, six warmed calls per process, with identical
precision, governor, scheduling and scientific settings. The attention setting is the
same within each block. Values are medians of process medians; within-process calls
correlate. First-call and memory data are retained separately in the receipt.

| Attention | Residues | Warm gradient seconds A → B | Warm forward seconds A → B | Gradient / forward time change |
| --- | ---: | ---: | ---: | ---: |
| Enabled | 146 | 1.116 → 1.042 | 0.343 → 0.293 | -6.7% / -14.4% |
| Enabled | 185 | 1.745 → 1.696 | 0.500 → 0.439 | -2.8% / -12.1% |
| Enabled | 169 | 1.460 → 1.404 | 0.424 → 0.369 | -3.8% / -12.9% |
| Disabled | 146 | 1.226 → 1.168 | 0.395 → 0.341 | -4.7% / -13.7% |
| Disabled | 185 | 1.991 → 1.948 | 0.587 → 0.536 | -2.2% / -8.8% |
| Disabled | 169 | 1.648 → 1.593 | 0.504 → 0.445 | -3.4% / -11.6% |

One disabled 169-residue process was guard-stopped for background indexing. The
completed A/B prefix was retained and the remaining B/A pair restarted in fresh
processes; the interrupted process is excluded and all raw evidence retained.

All 24 normal processes match the 18 model/gradient, three warmed-forward and three
previous-state fields byte for byte. The enabled-attention cases match the preserved
146/185/169 archives. The disabled 146 case also matches the independent review's saved
control; disabled 185/169 cases share the first control reference and frozen previous
state within their block. Every comparison uses exact bytes, without relaxed tolerances.

Separate admission diagnostics verify 158 forward and 153 gradient conversions, with
identical slice/helper counts, in all six size/attention combinations. The new cost gate
is a strict subset of the v4 element gate: `3T+5D<3N` implies `T+D<N`. With the other
admission constraints unchanged, these identical censuses preserve the intended AF2
sites. The separately instrumented 146 case also matches v4 per-phase kernel-count work,
reconciles direct/native counters exactly, and matches every full warmed call tree
across arms. Diagnostic timing is excluded from normal performance data.

Fresh representative AF2 behavior and benefits were retained, so no affected fixed-work
trajectory comparison was identified for rerun. No v5 full-trajectory timing is claimed.
The fixed-work trajectory numbers below remain explicitly historical v4 results. All
samples are from one shared M5 Max with uncontrolled thermals and desktop activity;
activity checks cannot prove GPU exclusivity.

Exact v5 build identities, source fingerprints, all measurements, the rerun decision and
raw evidence checksums are in `validation.json`. The preserved [v4
receipt](https://github.com/cytokineking/FreeBindCraft/blob/8c789beb8c3236281e2f62e21ad467d863a6aabf/installation/metaljax/slice-before-convert/validation.json)
remains tied to the reviewed v4 commit. The optimization is still disabled by default;
this change does not publish a release or change installer pins.

## Historical v4 measurements

The following observations belong to `0.11.9+fbc4slicedev4` at reviewed commit
`8c789beb8c3236281e2f62e21ad467d863a6aabf`. Throughout this historical section, “final”
refers to v4.

### Normal fixed-input measurements (v4)


Measurements used a shared Apple M5 Max with 128 GiB of memory. Each initial block has
two fresh processes per arm and six warmed calls per process. Values below are medians
of fresh-process medians; the warmed calls are correlated. Normal processes disable
private counter reads, hooks, profiling fences and persistent JAX compilation caching. A
is the rebuilt enabled-attention control; B enables the transformation in
`0.11.9+fbc4slicedev4`. All earlier prototype timings are excluded.

| Initial block, total residues | Warm gradient seconds A → B | Warm forward seconds A → B | Maximum process physical footprint GiB A → B |
| ---: | ---: | ---: | ---: |
| 146 | 1.125 → 1.052 | 0.345 → 0.293 | 6.76 → 6.44 |
| 185 | 1.768 → 1.817 | 0.499 → 0.476 | 10.14 → 9.85 |
| 169 | 1.565 → 1.507 | 0.458 → 0.402 | 8.63 → 8.29 |

At 146 residues, forward time changed by -15.0% and gradient time by -6.6%. At 185
residues, forward time changed by -4.7% and gradient time by +2.8%. The gradient process
ranges overlap. At 169 residues, forward time changed by -12.2% and gradient time by
-3.7%. The gradient process ranges overlap.

The initial 185-residue candidate processes differed by about 13% in warmed forward time
and 10% in gradient time. The slower process remains a valid observation. At the user's
request, the entire block was repeated in the opposite balanced order, with every
initial and repeated observation retained:

| 185-residue block | Fresh processes per arm | Warm gradient seconds A → B | Warm forward seconds A → B |
| --- | ---: | ---: | ---: |
| Initial BAAB | 2 | 1.768 → 1.817 | 0.499 → 0.476 |
| Repeat ABBA | 2 | 1.806 → 1.761 | 0.506 → 0.453 |
| All valid processes | 4 | 1.787 → 1.761 | 0.503 → 0.453 |

Across all valid 185-residue processes, gradient time changed by -1.5% and forward time
by -10.1%. The maximum physical footprint across all 185-residue processes is 10.16 →
9.88 GiB. The repeat does not establish why the original process was slower. Shared
desktop activity and thermals remain uncontrolled; these small samples cannot establish
a universal speedup. Physical footprint is a macOS process lifetime high-water mark,
without isolating GPU allocations. Every individual timing and range is retained in [v4
validation
receipt](https://github.com/cytokineking/FreeBindCraft/blob/8c789beb8c3236281e2f62e21ad467d863a6aabf/installation/metaljax/slice-before-convert/validation.json).

First calls include compilation and materialization:

| Block | First forward seconds A → B | First gradient seconds A → B |
| --- | ---: | ---: |
| 146 | 1.424 → 1.372 | 3.724 → 3.674 |
| 185 | 1.560 → 1.524 | 4.361 → 4.329 |
| 169 | 1.496 → 1.446 | 4.098 → 4.011 |
| 185 repeat | 1.561 → 1.536 | 4.370 → 4.321 |

All sixteen final normal processes match all 18 model/sequence-gradient fields, three
warmed-forward fields and three previous-state fields byte for byte. Both established
fixtures reproduce the preserved baseline arrays despite the rebuilt control's different
plugin hash. The additional fixture uses ubiquitin with a 93-residue binder (169 total),
with its first control supplying the shared frozen previous state. No comparison
tolerance was relaxed.

The final 84-test backend/repository suite passes, including positive/negative
admission, raw BF16 edge bits, layout-sensitive BF16 reductions, clamped indexing,
rematerialization, gradients, compiled loops, retained donated outputs and governor
refusal. The trajectory resume path verifies completed process receipts, configurations,
actual designs and records before skipping them; the completed series is checked again
without relaunching GPU work.

Separate diagnostics admit 158 conversions in the forward module and 153 in the gradient
module. The five repeated large Evoformer cast groups vanish from the forward trace;
once-only gradient groups with whole-stack consumers remain. Reflected FP32-to-BF16
logical copy bytes change from 25.34 → 0.57 GB per forward and 25.91 → 1.15 GB per
gradient. These are logical accesses, not DRAM measurements or exclusive timing shares.
Dispatch counts increase slightly (55412 → 55448 forward; 138441 → 138476 gradient).

All four warmed diagnostic calls in each arm have exact agreement between direct Metal
dispatches, per-buffer native operation counts, kernel totals and settled counters.
Immediate versus delayed reads show no counter changes, and every full warmed call tree
matches byte for byte across arms. This validates the new samples. The two historical
196-dispatch deficits remain causally unresolved because they lack per-call traces and
arrays. No constant-count assumption or widened guard is used.

### Fixed-work trajectories and remaining limits (v4)

The unchanged checksum-pinned harness runs the 75/45/5/15 schedule with the original
scheduling control. The initial ABBA series spans a user-requested pause; the verified
completed A1 control was retained when resuming:

| Process | Design stage seconds | Complete fixed-work harness seconds |
| --- | ---: | ---: |
| A1 | 211.52 | 266.12 |
| B1 | 193.36 | 246.63 |
| B2 | 191.80 | 243.97 |
| A2 | 202.53 | 255.84 |

The initial two-process median design stage is 207.03 → 192.58 seconds (-7.0%). The
initial two-process median complete fixed-work harness is 260.98 → 245.30 seconds
(-6.0%).

The recorded repeat threshold is a within-arm max/min ratio above 1.10. Within-arm
max/min ratios are 1.044/1.008 for design and 1.040/1.011 for complete harness time
(A/B). No trajectory timing crossed that threshold, so no additional trajectory pair was
required.

Maximum observed process physical footprint across valid trajectories is 8.43 → 7.95
GiB. The successful processes finish within the governor and guard limits. One partial
candidate was stopped by the guard for a background indexing spike, and a partial retry
was stopped at the user's request. Both attempts are retained and excluded. Restarted
processes use fresh snapshots and caches under the unchanged guard. No completed slow
observation is discarded.

Every completed trajectory matches all 140 actual pre-substitution design metric records
and the actual generated trajectory PDB. The generated PDB SHA-256 is
`c2a36b18202d8c1eb2a56051ec31e46c0eb4bd44000702cabbbbb519633ba62c`; the record SHA-256
is `b67dc5caace0eef477bc74a4dc437fc118058751b4fbd245da66e294401f8baf`. The rebuilt
control also matches the preserved historical control's records and generated PDB. Each
process performs 125 design backpropagations, 16 design forwards, one MPNN sample, two
complex predictions, two binder predictions, three relaxations/scored complexes and one
filter/selection pass.

Accepted designs remain zero. Downstream prediction/scoring uses the unchanged shared
fixture; MPNN samples and relaxed PDB bytes may vary. This evidence measures fixed work
and cannot establish accepted-design throughput. Shared GPU/desktop activity, thermal
drift and the pause limit causal interpretation.

Discarded prototypes and validation repairs remain documented in the receipt: an initial
GUI test abort, admission-fixture/cache mistakes, the direct-only v1 pilot, a
duplicate-image diagnostic reader, the wider v2 layout bug, and the conservative v3
proof missing gradient weights behind identity barriers. All final timings use the
corrected v4 compiler.

The next narrow priority-2 opportunity is to measure a combined dynamic slice and BF16
conversion that avoids the new FP32 slice temporary. The candidate adds about 0.725 GB
of reflected logical FP32 dynamic-slice traffic per forward/VJP while removing the
larger whole-stack casts. A follow-up must retain clamp, conversion rounding, use
census, layouts and lifetime contracts and establish an actual timing benefit. Larger
remaining activation/layout copies need provenance and cost evidence. No speedup is
forecast from bytes.

[v4 validation
receipt](https://github.com/cytokineking/FreeBindCraft/blob/8c789beb8c3236281e2f62e21ad467d863a6aabf/installation/metaljax/slice-before-convert/validation.json)
records results, exact source/build identities, raw-evidence checksums, repetition
policy and limitations. Local wheels, arrays, logs and process samples remain ignored
under `.local/` and `runs/priority2/`.
