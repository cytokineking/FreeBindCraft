# fbc3 publication verification

The reviewed `0.11.9+fbc3` artifacts were published on 2026-10-24T05:53:23Z in
the [MetalJAX fork release](https://github.com/cytokineking/metaljax/releases/tag/v0.11.9-fbc3). No artifact was rebuilt.
The source branch is [`codex/freebindcraft-0.11.9-fbc3`](https://github.com/cytokineking/metaljax/tree/codex/freebindcraft-0.11.9-fbc3),
at commit `d5962d9ba7340f8b7b893f9332851d34da31b7d7`; annotated tag `v0.11.9-fbc3`
has object `64be29d55e7813ecb76e4608dd8e403a133944d2` and resolves to that commit.
Earlier release branches/tags/assets remain intact.

All eight public assets were downloaded and their SHA-256 hashes checked against
the retained reviewed files. The wheel hash is
`0fa13c5def15e89024c624b81f5537492d269c6cb209602a5ba1e00baf25faa2`;
the installed plugin hash is
`25230c9a0a5b342c6a14313a06db0844f0e9cf3ce562ac5b8543a3edfd5449e4`.
A fresh fetch from the fork also verified all 346 committed source files.

The ordinary Python dependency installer selected `frozen` with no developer
override, starting with a new virtual environment and empty wheel and uv caches.
Its log confirms a download from the exact URL pinned in `release.json`.
The resulting package inventory and native-library hashes match the independently
reviewed environment, including JAX/jaxlib `0.11.2` and unchanged MLX bytes.
Pytest was added afterward solely for verification. Native/model assets were
reused read-only; this was not a new download/build of every native/model asset.

The public installation passed the full relevant backend/repository/installer/
launcher suite: **89 passed, 2 warnings in 58.85s**. Both optimization variables were unset at
the root. The two warnings are existing JAX `jit(backend=...)` deprecations.
Guard and memory monitor completed normally.

Two additional checks used the real Apple launcher with the publicly installed
wheel and actual AF2 benchmark, in isolated entry folders. Before backend imports,
the default had both variables unset; `--no-attention-epilogues` overrode inherited
attention `1` to `0`. Default attention admission was present, while the opt-out
had none. Slice admission remained 158 forward / 153 gradient conversions in
both cases. All 24 saved model/gradient/forward/previous-state fields matched
the public fbc1 control exactly. Both checks completed under the unchanged
activity guard and memory monitor. Their debug timings are not performance data.

The first transport attempt occurred while the release was still draft and
returned HTTP 404 before package sync. That log and its empty environment/cache
are retained. The successful installation used a new environment and empty
caches after public availability was confirmed.

[publication.json](publication.json) records release identity, download hashes,
installation, tests, CLI evidence and raw-file hashes. The original
[candidate report](README.md), [validation receipt](validation.json), source and
build provenance remain unchanged historical records; their unpublished labels
describe the prepublication review stage. Current installation pins are in
[release.json](../../release.json). The existing performance measurements and
their single-M5-Max, fixed-work limitations remain unchanged.
