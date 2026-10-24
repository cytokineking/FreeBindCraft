# Frozen MetalJAX for Apple Silicon

This draft branch selects the **unpublished production candidate**
`metaljax==0.11.9+fbc3`, with JAX/jaxlib `0.11.2`. The exact source revision,
versioned wheel URL and checksums are pinned in [release.json](release.json) and
[the Apple requirements](../apple-silicon-requirements.txt). The fbc3 source
commit, tag and artifacts remain local for review. Its planned release URL is
not yet live; a real public-download installation has not been validated.

The candidate includes the existing static nested-loop repair, the reviewed
slice-before-BF16-conversion rule and the existing attention capability. Slice
rewriting is enabled when `METALJAX_SLICE_BEFORE_CONVERT` is unset. Use
`METALJAX_SLICE_BEFORE_CONVERT=0` for diagnostic opt-out. Attention epilogues
are also enabled when `METALJAX_EVO_ATTN` is unset.
Use `./run_apple_silicon.sh --no-attention-epilogues ...` to opt out; the launcher
sets `METALJAX_EVO_ATTN=0` before Python starts, overriding an inherited `=1`.
An explicit environment `=0` is preserved when the CLI flag is absent.
Unsupported shapes/patterns retain the existing lowering.

The slice rule retains the v5 layout, clamping, rounding, all-use, profitability,
donation and governor constraints. It requires strict aggregate logical savings
`3T + 5D < 3N`, charging dynamic FP32 slice materialization. Unsupported,
noncontiguous, overflowing and unprofitable cases keep the original lowering.
There is no fused slice kernel, persistent cache or new tuning parameter.

Shared scientific functions, CUDA/Linux installation, model inputs, precision,
losses, RNG/dropout, recycling, rematerialization, optimizer and campaign logic
are unchanged. Native scheduling remains 800 operations / 512 MiB, with the
1 GiB loop window and enabled memory governor. Normal validation leaves those
backend control variables unset and uses an external process memory limit.

## Installation selection and troubleshooting

After separately approved publication and public-download validation, the
ordinary Apple command selects the new frozen artifact:

```bash
bash install_bindcraft.sh
# Add --no-pyrosetta to omit native Rosetta.
```

Before publication, reviewers can seed a private download-cache directory with
the exact candidate wheel and run the same Python dependency installer. This
uses normal `frozen` selection and the pinned checksum, without a local override:

```bash
uv venv --python python3.12 .venv-metaljax-review
python3 installation/install_apple_python.py \
  --python .venv-metaljax-review/bin/python --downloads /absolute/candidate-cache
```

The cache must contain the exact filename from `release.json`. A valid cached
artifact is verified before use; the future URL need not be contacted. This
check covers the normal dependency-selection/install path, not a public
release download or a new installation of all native model assets.

Pinned stock upstream MetalJAX remains available for troubleshooting:

```bash
FREEBINDCRAFT_METALJAX=stock bash install_bindcraft.sh
```

Explicit `FREEBINDCRAFT_METALJAX_WHEEL=/absolute/path.whl` overrides still take
precedence over stock/default selection. Those developer-selected bytes are
outside the frozen artifact's expected-checksum guarantee. Unknown modes and
missing overrides fail. The installer checks release/requirements agreement and
verifies the selected frozen/stock artifact before uv changes the target venv.
A damaged cache is re-downloaded and verified; a checksum failure preserves the
environment. No automatic dependency swap or whole-job restart is introduced.

## Immutable candidate and build provenance

- Version: `0.11.9+fbc3`.
- Backend source: `d5962d9ba7340f8b7b893f9332851d34da31b7d7`.
- Local, unpushed tag: `v0.11.9-fbc3`.
- Wheel SHA-256:
  `0fa13c5def15e89024c624b81f5537492d269c6cb209602a5ba1e00baf25faa2`.
- Plugin SHA-256:
  `25230c9a0a5b342c6a14313a06db0844f0e9cf3ce562ac5b8543a3edfd5449e4`.
- MLX runtime SHA-256:
  `72dc8d9a4b1e07df2f6bc91383e0567acf3ce292398ae6ce4ed1d19fce0b1374`.

[Source provenance](production/fbc3/source.json) records the frozen parent,
reviewed patches, defaults patch, exact commit metadata and local tag object.
The four patches reproduce the committed source tree exactly. Apache-2.0,
MLX MIT and third-party license/notice bytes are unchanged from fbc1 in the
packaged wheel. Frozen MLX SDK, XLA revision
`91888df6ce85102c30220e41d952065925e10886`, Bazel 8.7.0 and hatchling 1.32.4
inputs are retained. Their existing fbc1 asset URLs/checksums remain immutable;
there is no moving-branch or latest-version lookup.

The maintainer builder consumes the exact source commit and verified inputs:

```bash
.venv-metaljax-fbc2-build/bin/python installation/build_metaljax_apple.py \
  --python .venv-metaljax-fbc2-build/bin/python \
  --source-dir .local/src/metaljax-fbc3 --cache-dir .local/fbc3-build --jobs 6 --batch
```

Before source publication, that directory must be the verified local candidate
checkout (or a checkout cloned from its local source bundle). Its commit and
tag have not been pushed to the backend repository. The builder does not
reapply patches to a committed fork. It preserves unexpected source/config
changes, verifies archived inputs, and records flags, compiler/SDK, package,
plugin, wheel and runtime identities. Repeating the recorded build reproduced
the exact candidate bytes. This reused verified build caches; it was not a
second cold XLA build. Different toolchains or paths may produce other bytes.

The earlier fbc1 release manifest is preserved in [releases/fbc1.json](releases/fbc1.json).
The superseded attention-off fbc2 candidate remains preserved in
[its manifest](releases/fbc2.json) and [original report](production/README.md),
with its own immutable source, tag, artifacts and measurements.
Development attention/slice recipes retain that frozen base rather than
applying their patches a second time to the new production source. Historical
[fbc1 validation and public release details](https://github.com/cytokineking/FreeBindCraft/blob/cedf3f254059f6a3cbaed450cec23e1a07ebb57a/installation/metaljax/README.md)
and the [v4/v5 experiment report](slice-before-convert/README.md) remain
separate from production-candidate validation.

## Production-candidate validation and approval boundary

The [production report](production/fbc3/README.md) and
[receipt](production/fbc3/validation.json) record the exact packaged-wheel checks,
fresh dependency installation, unset/default and opt-out tests, numerical and
installer integrity tests, and the balanced end-to-end comparison against the
current released fbc1 configuration. Local wheels, source bundle/archive,
provenance and checksum files remain preserved under
`runs/production-fbc3/unpublished-assets/`; they have not been uploaded.

Remaining release steps require separate review and approval:

1. Review the frozen candidate, source/tag identities, tests and measurements.
2. Push the approved backend source revision and tag, then upload the exact
   retained wheel, source archive, provenance, notices and checksum assets.
   Do not rebuild or change the approved artifact bytes during publication.
3. From an empty cache, download the real public assets, verify their pinned
   hashes and perform the fresh public-download installer/default checks.
4. Record that public-download validation, update publication status, and merge
   the FreeBindCraft integration only after the required user review.

None of those publication or merge steps has been performed by this task.
