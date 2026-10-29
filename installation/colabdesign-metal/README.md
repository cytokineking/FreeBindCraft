# ColabDesign Apple integration

The Apple installer applies the ColabDesign changes in this directory automatically. Metal uses unchunked execution by default, preparing forward and gradient functions consistently. CPU and CUDA behavior is unchanged.

For lower memory use, `FREEBINDCRAFT_METAL_CHUNKING=1` restores four-residue sub-batches above 384 total target-plus-binder residues. See the [Apple Silicon guide](../../apple_silicon/README.md#performance-and-memory).
