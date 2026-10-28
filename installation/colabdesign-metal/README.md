# Apple Silicon chunking policy

ColabDesign normally uses sub-batches of four for memory-intensive Evoformer
operations when the target plus binder exceeds 384 residues. This policy saves
memory on CUDA/NVIDIA, but its many small GPU operations cause severe slowdowns
with MetalJAX on Apple Silicon.

The Apple installer disables chunking by default and prepares both forward and
gradient functions with the same configuration. Removing chunking made complete
AlphaFold recycles **about 3× faster** in our 406-residue M5 Max test, at the cost
of higher memory use. CPU/CUDA behavior is unchanged.

| Test job | Total residues | Default (unchunked) | Original policy |
| --- | ---: | ---: | ---: |
| Ubiquitin (76) + 70-aa binder | 146 | ~6.9 GiB | ~6.9 GiB; already unchunked |
| CA2 (256) + 150-aa binder | 406 | ~27.5 GiB | ~18 GiB with chunking |

These are measured peak process footprints; leave additional memory for macOS,
other applications and downstream processing. Memory checks remain enabled.

When memory savings matter more than speed, restore the original chunking policy:

```bash
FREEBINDCRAFT_METAL_CHUNKING=1 ./run_apple_silicon.sh --interactive
```

The same environment variable works with your usual job arguments. It enables
four-residue sub-batches above 384 total residues; smaller jobs remain unchunked.
