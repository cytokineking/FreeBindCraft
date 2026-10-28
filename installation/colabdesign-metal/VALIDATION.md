# Draft validation: unchunked Metal preparation

Validated on 2026-10-28 against public FreeBindCraft master
`056ff605371ce46174f2cffb7ea44e7337df4ee5`, in a fresh checkout and environment.
Host: Apple M5 Max, 128 GiB unified memory, arm64 macOS 27.0.1. Dependencies:
Python 3.12.15, MetalJAX 0.11.9+fbc3, JAX/jaxlib 0.11.2, and ColabDesign 1.1.3
at the revision in `manifest.json`. Native MetalJAX and MLX hashes remained
unchanged. [Machine-readable results](validation.json) include their identities.

## Installation and tests

- CPU suite: 35 passed; five Metal-only checks skipped in that invocation.
- Metal suite: all five alignment/numerical regression checks passed.
- Four preparation integration checks also passed against a temporary copy of
  the pristine dependency. They test both callables, model reuse across
  383/384/385/406/146 residues, failure before replacing the previous model, and
  preservation of CPU/CUDA preparation using fake runners.
- The real dependency installer rerun reported the patch already applied.
  All 98 installed package files matched the tested candidate snapshot, and
  the added/modified files' distribution `RECORD` hashes and sizes matched.
- Comparing complete dependency snapshots found only `af/prep.py` and the new
  `af/metal_execution_policy.py` changed.
- Real warm preparation and a subsequent finite full recycle passed. Preparing
  the same 406-residue model for 146 residues also retained unchunked selection.
  The selected model's verified resident pages received about 355 MiB of credit.

Run the regression suites after installing this checkout:

```bash
JAX_PLATFORMS=cpu .venv-metaljax/bin/python -m unittest discover -s tests -v
JAX_PLATFORMS=metal .venv-metaljax/bin/python -m unittest discover -s tests -p test_metal_alignment.py -v
```

## Fixed-input complete-recycle probes

Both arms use the same public application and native runtime. The baseline is
the pristine pinned ColabDesign package; the candidate differs only by this
preparation patch. These are complete `_recycle(..., num_recycles=1,
backprop=True)` calls: one forward pass and one differentiated pass. Model
weights, input arrays, and PRNG keys match by fingerprint.

The fixtures are ubiquitin 1UBQ (76 target + 70 binder residues, hotspot 44) and
CA2 1CA2 (256 observed target + 150 binder residues, hotspot 40). CA2 chain A is
renumbered 1–256, with heteroatoms excluded and original coordinates preserved.
The model is `model_1_multimer_v3`, seed/key 42, BF16, rematerialization enabled,
last-recycle gradients, dropout enabled, soft 0.9, hard 0 and temperature 1.
The application loss weights and callbacks are unchanged: pLDDT 0.1, PAE 0.4,
interface PAE 0.1, contacts 1, interface contacts 1, RG 0.3, interface pTM 0.05
and helix -0.3. Contact settings are 2 contacts at 14 Å with sequence separation
9; interface contacts are 2 at 20 Å. Initial previous-state tensors are zero.
The second state adds 0.25 to sequence parameter `[0, 0, 0]`.

One warm-up and two timed calls are run for each state. Every timed call has
explicit all-stream Metal synchronization before and after execution; first
calls are recorded separately. Persistent JAX compilation caching is disabled.
The 406-residue comparison uses fresh processes in baseline/candidate/candidate/
baseline order. Calls within a process are correlated, not independent trials.

| Probe | Published default | Draft | Interpretation |
| --- | ---: | ---: | --- |
| 146 residues | 1.352 s | 1.378 s | One process per arm; same calculation; no demonstrated speed gain |
| 406 residues, quieter pair | 38.664 s | 13.279 s | 2.91×, 65.7% less time; one process per arm |
| 406 residues, all completed processes | 39.001 s | 13.300 s | Two processes per arm; supporting observations |

Headline timing excludes any process with a detected running CPU-contention
sample. The other baseline and candidate each had one brief Codex renderer CPU
spike; their complete measurements remain in the results. Admission required a
six-second quiet interval. Continuous process and aggregate GPU samples cannot
establish per-process GPU exclusivity or control desktop/thermal effects. These
are preliminary shared-host observations, not a controlled multi-host benchmark.

Peak process physical footprint at 406 residues increased from 18.00–18.04 GiB
to 27.51 GiB. At 146 residues it was 6.82 versus 6.86 GiB. Memory pressure stayed
normal and swap remained zero. Validation retained the native 104 GiB system
ceiling plus an 80 GiB owned-process guard. One incomplete candidate attempt
was stopped when its active-turn proof expired; its partial times are excluded.
Its replacement used actual-turn checks every 25 seconds and completed. All
owned validation processes exited.

## Numerical result and remaining gates

All 33 output arrays match bit-for-bit at 146 residues for both states. At 406
residues, all 32 non-gradient arrays, including losses and structures, match
bit-for-bit. The sequence gradient differs as follows:

| State | Relative gradient L2 difference | Gradient cosine | Angle |
| --- | ---: | ---: | ---: |
| Original | 3.95% | 0.999376 | 2.025° |
| Perturbed | 1.56% | 0.999879 | 0.892° |

Each arm's complete output is reproducible across its two fresh processes and
all repeated calls. Candidate gradients also exactly match the archived
unchunked implementation. This probe establishes finite, repeatable outputs;
it does not establish derivative accuracy or equivalent final design quality.

Before promoting the draft, run the agreed fixed-work workflow benchmark:
ubiquitin and CA2, OpenMM/FASPR and native PyRosetta, three fresh processes per
workload/mode/configuration, including startup, the full 75/45/5/15 schedule,
MPNN, predictions, relaxations and scoring. Preserve its fixed-work assertions
and real trajectory metrics. Its shared downstream fixtures test runtime and
wiring, not acceptance of the candidate's own designs. Any quality claim needs
actual candidate outputs evaluated with unchanged scientific filters.

The older 2.27× research result used a private LEGACY scan configuration as its
control. The measurements above are the new comparison against the published
default; the two results must not be conflated.
