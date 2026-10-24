# Production MetalJAX fbc2 candidate: validation results

The exact `0.11.9+fbc2` wheel packages the reviewed v5 profitability guard with
slice-before-convert enabled by default. `METALJAX_SLICE_BEFORE_CONVERT=0`
provides diagnostic opt-out. Attention remains disabled by default; its existing
`METALJAX_EVO_ATTN=1` capability still passes numerical and execution checks.

This is an **unpublished candidate** in a draft FreeBindCraft integration.
The backend source commit/tag, wheel, source archive/bundle, provenance,
notices and checksums remain local. The planned fbc2 public URLs have not been
download-tested. No backend source/tag push, release upload or merge occurred.

## Exact artifacts and reproducibility

| Artifact | Identity |
| --- | --- |
| Production version | `0.11.9+fbc2` |
| Local backend source commit | `bcd19ebca757b18f12b7eb963eebfacc428ba1fe` |
| Local annotated tag | `v0.11.9-fbc2` |
| Tag object | `f11fb2d3aedcaf6002b8cceb28e933ec0e591c39` |
| Wheel SHA-256 | `3c41045c1ae5e555d32332656af21cd1b9296e2f6aa3079cc51a088d70a2933f` |
| Installed plugin SHA-256 | `d8c582111e2b7bfef8ad6e18ccbefb07ed1a22c2e0225c0b88b1bd8c0adb70aa` |
| Frozen MLX runtime SHA-256 | `72dc8d9a4b1e07df2f6bc91383e0567acf3ce292398ae6ce4ed1d19fce0b1374` |

The frozen fbc1 source plus the attention, v5 slice and production-default
patches reproduce the exact committed backend tree. [source.json](source.json)
records the patch hashes, tree/commit metadata and tag identity.
[build-manifest.json](build-manifest.json) records the compiler, SDK, flags,
build tools, runtime libraries and frozen input checksums. All committed source
file SHA-256 values are included in [validation.json](validation.json).

Two builds of the same immutable source with the recorded toolchain and paths
produced byte-identical plugins and wheels. They reused verified Bazel caches;
this does not claim a second cold XLA build or cross-toolchain reproducibility.
The first manifest retained its provisional zero expected checksums; the final
manifest records the resulting real pins and matches both bytes exactly.
JAX/jaxlib remain `0.11.2`; frozen MLX SDK/runtime, XLA, Bazel and hatchling
inputs are unchanged. Apache-2.0, MLX MIT and third-party notice contents in
the wheel are byte-identical to the released fbc1 wheel.

The retained local publication assets and their hashes are listed in
[SHA256SUMS.candidate.txt](SHA256SUMS.candidate.txt). They live under
`runs/production-fbc2/unpublished-assets/`, including the source bundle needed
to recover the local source/tag. This inventory records local bytes, not a
claim that the assets exist on the release server.

## Installation and default behavior

The candidate was installed into a new `.venv-metaljax-fbc2-clean` through the
ordinary Python dependency installer with its default `frozen` selection and
the pinned wheel pre-seeded in a private verified download cache. There was
no local-wheel override and no fetch of the future release URL. This validates
the normal selection, checksum and dependency-install path before publication.
Model parameters and native tools were reused read-only; this was a fresh
Python environment, not a full native/model-asset installation.

The current public fbc1 wheel was independently installed into a new
`.venv-metaljax-fbc1-release` using the archived released installer/pins.
All installed scientific packages match between environments except MetalJAX.
The released plugin is `005ab1302a7f49ec27e1447ab31db89e78b0e22793e6e14cf66bfe7773dce2fd`;
the earlier attention-development control is not used as this baseline.

A real isolated installer invocation with deliberately corrupt file-URL bytes
fails its checksum before invoking uv. The target package inventory remains
unchanged. Existing stock/custom override and corrupt-cache tests pass; the
stock upstream URL/hash is unchanged. The added cache test protects operation
when the unpublished URL is unavailable, mocking only the network boundary.

The backend/repository run passed 88 tests with two existing deprecation
warnings. One redundant manifest-text assertion was subsequently removed;
all 87 retained tests were covered by that passing run, and the final six
installer tests pass again. Final collection confirms 87 tests. The receipt
preserves the original run and this distinction.

Independent workers explicitly remove the slice variable and prove that the
unset default matches opt-out raw bytes and explicit-on admission across all
47 cases. Unsupported layouts and near-full dynamic windows still decline,
small profitable windows still admit, and conversion-edge, gradient,
rematerialization, donation and governor checks pass. A separate worker removes
the attention variable, verifies no forward/backward attention recognition,
and matches explicit-off bytes. Existing opt-in attention/derivative execution
tests pass. These workers prevent a runner-injected `=1` from substituting for
the default-behavior proof.

A separate 146-residue debug run on the installed fbc2 wheel, with slice and
attention unset, confirms 158 forward and 153 gradient conversions through the
intended slice/helper sites and matches all 24 fields exactly. Its compiler-log
diagnostics are excluded from normal timing and do not read dispatch counters.

Normal model and trajectory processes leave slice, attention, precision,
governor and scheduling controls unset. Launch receipts contain only Metal
selection/JIT/cache controls. Native scheduling remains 800 operations / 512 MiB,
with the 1 GiB loop window and enabled governor. On this 128 GiB machine the
native governor budget is 75% of RAM; an external 24 GiB process-tree RSS bound
also applies. Normal timing uses neither diagnostic hooks nor private dispatch
counter reads, and persistent JAX compilation caching is disabled.

## Representative AF2 model checks

A is the exact released fbc1 wheel, B the exact freshly installed fbc2 wheel.
Each default block uses ABBA order, two fresh processes per arm and six warmed
calls per process. The sampling unit is the process; warmed calls correlate.
All 18 model/gradient, three warmed-forward and three previous-state fields
match byte for byte in every process. Inputs, fixed RNG/dropout, previous state,
precision, BF16, rematerialization and recycling settings are unchanged.

| Residues / block | Warm gradient seconds A → B | Change | Warm forward seconds A → B | Change |
| --- | ---: | ---: | ---: | ---: |
| 146 | 1.262 → 1.203 | -4.7% | 0.403 → 0.358 | -11.1% |
| 185 initial ABBA | 2.253 → 2.210 | -1.9% | 0.657 → 0.613 | -6.6% |
| 185 repeat BAAB | 2.034 → 2.045 | +0.5% | 0.605 → 0.556 | -8.1% |
| 169 | 1.829 → 1.743 | -4.7% | 0.562 → 0.483 | -14.1% |
| 185 all valid processes | 2.034 → 2.045 | +0.5% | 0.605 → 0.556 | -8.1% |

The pooled 185-residue warm gradient change is +0.5%, while forward
time changes -8.1%. These samples do not establish a gradient speedup
for that case.

The original 185-residue block's later B/A processes were both slower than
the earlier pair (gradient 2.449/2.487 s versus 1.971/2.019 s). A complete
opposite-order BAAB block was therefore repeated. Every valid original and
repeat observation remains in the report/receipt; no slower process was removed.
First-call times, individual process medians, all warmed samples and memory
measurements are retained in the receipt and raw files. The shared machine's
thermal/activity state is uncontrolled, so these small timing differences
should not be generalized to every input or machine.

A separate 146-residue two-process AB check explicitly enables attention in
both launch environments and retains the slice default unset. All 24 fields
match exactly. It is capability/composition evidence, not an additional default
timing block or a full balanced attention benchmark. The exact candidate's
backend suite independently verifies that the recognized attention kernels run.

## Fresh fixed-work trajectory comparison

Four new processes run the unchanged fixed-work harness in ABBA order, two
per exact installed wheel. Both optimization flags remain unset in every
process. This production series is separate from earlier development evidence.
The scientific snapshot is `cf136d1033790eef7b866bf07e92ed82b5572003`.
The preserved harness SHA-256 is
`e73d95f8c09f4e5dec907ac68fd2240e2220242873250e7c07eb93cf3462075a`.
Configuration, fixture, native-tool, source and package identities are pinned.
The schedule is 75/45/5/15. The preserved harness uses `--mode no-pyrosetta`
with OpenMM/OpenCL relaxation and forces the fixed stage work before substituting
the shared downstream fixture.

All actual generated trajectory PDBs and all 140 actual design metric records
match byte for byte before downstream fixture substitution. Work counts match:
125 design backprop calls, 16 design forward calls (140 nominal iterations plus
the semigreedy initial evaluation), one MPNN sample, two complex forwards,
two binder forwards, three scored/DSSP structures and one filter/selection.

| Process | Wheel | Design seconds | Complete harness seconds | Peak process footprint GiB | Peak tree RSS GiB | Generated PDB / 140 records |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| a1 | fbc1 | 230.568 | 293.559 | 8.854 | 1.852 | Exact / exact |
| b1 | fbc2 | 209.716 | 265.169 | 8.691 | 1.813 | Exact / exact |
| b2 | fbc2 | 206.935 | 259.707 | 8.693 | 1.818 | Exact / exact |
| a2 | fbc1 | 225.642 | 278.558 | 8.850 | 1.816 | Exact / exact |

The median design time changes **228.105 → 208.326 s
(-8.7%)**; complete harness time changes
**286.059 → 262.438 s
(-8.3%)**. Maximum process lifetime footprint
across each arm changes **8.854 →
8.693 GiB**. Complete time includes
startup and all fixed downstream work; monitor wall time is retained separately.
Within-arm design times differ by about 2.2% for fbc1 and 1.3% for fbc2;
complete times differ by 5.4% and 2.1%. No trajectory was identified as a
material timing outlier. The first baseline startup took 15.451 s versus
8.447–8.548 s in the other processes; that observation remains included.

| Stage | Median seconds A → B | Change |
| --- | ---: | ---: |
| startup | 12.000 → 8.495 | -29.2% |
| design 140 iterations | 228.105 → 208.326 | -8.7% |
| trajectory relaxation and scoring | 12.521 → 13.009 | +3.9% |
| mpnn | 0.531 → 0.529 | -0.4% |
| validation model setup | 2.861 → 2.822 | -1.4% |
| complex prediction and relaxation | 25.369 → 24.946 | -1.7% |
| complex scoring | 3.161 → 3.078 | -2.6% |
| binder prediction and alignment | 1.362 → 1.089 | -20.1% |
| filter evaluation and selection | 0.004 → 0.004 | -2.7% |

Every trajectory exited successfully, with no guard interruption, process
memory stop or swap growth. Memory pressure stayed at level 1 (normal).
Per-process sampled/lifetime footprint, RSS, available memory, compression,
swap baseline/peak/final change and pressure samples are in the receipt.
Footprint includes CPU, native and Metal allocations; RSS alone misses some
Metal memory. System compression/swap include other running applications.

The harness has **zero accepted designs** and deliberately uses shared
downstream fixtures. The generated design's decision remains LowConfidence.
MPNN sampled sequences and relaxed structures may vary even with equal work.
The observed total is fixed-work harness completion, not accepted-design
throughput or evidence of improved scientific quality. Four processes on one
shared M5 Max cannot establish universal speedups or exclusive GPU attribution.

## Evidence and remaining release steps

[validation.json](validation.json) includes source/build/install identities,
every normal model observation, initial/repeat/pooled 185-residue statistics,
all four production trajectory reports, per-stage timing, memory/pressure,
work/equality checks and SHA-256 values for the retained raw evidence.
The first 185-repeat launcher failed before GPU work due to a missing path
variable; the corrected launcher ran fresh and the failed launch is retained.
Raw arrays, optimized IR, logs, activity/memory samples and generated PDBs stay
under `runs/production-fbc2/` in this worktree. Older v4 and v5 environments,
sources, wheels, receipts and observations are preserved; all 865 raw checksums
in each historical ledger still verify. Earlier timings retain their original
development-version labels.

FreeBindCraft scientific functions/settings, CUDA/Linux behavior, native
scheduling/governor and campaign logic are unchanged by this integration.
The profitability rule remains conservative `3T + 5D < 3N`, with checked
aggregate accounting and the reviewed layout/indexing/rounding/all-use gates.

The remaining steps are candidate review and separately approved backend
source/tag publication, upload of the exact retained assets, then an empty-cache
real public download/checksum/fresh-install/default verification. Record that
validation and publication status before approved merge. Until then, this PR
stays draft/unmerged and the planned fbc2 URLs must not be described as live.
