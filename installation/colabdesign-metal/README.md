# ColabDesign Metal preparation

The Apple installer applies this small, checksummed patch to ColabDesign commit
`8fd167defd1b0f98e29f00584ff4f3987e6935dc`. It changes model preparation on Metal
and adds a memory preflight. The frozen MetalJAX wheel is unchanged.

On Metal, preparation selects `subbatch_size=None` before constructing both the
forward and gradient callables. The dependency's original preparation switches
to a subbatch size of four above 384 residues and replaces only the forward
callable while sharing a mutable configuration. Preparing a private configuration
first avoids that inconsistent state. CPU and CUDA retain the pinned dependency's
original preparation behavior.

The memory preflight checks a necessary attention score-buffer bound and an
empirical whole-model estimate. Its reference is one BF16 multimer fixture with
406 total residues, rematerialization enabled, and last-recycle gradients on an
M5 Max with the frozen fbc3 runtime. The estimate scales with residue count,
head count and dtype; other models, devices and workloads are not qualified by
that calibration. A refusal raises `MemoryError` before preparing new callables.

Warm preparation credits only verified resident model pages already included in
the process footprint. It does not equate logical tensor sizes, RSS, or reported
device cache sizes with reusable physical memory. An inconclusive estimate is
reported as `unknown`; the existing native allocation and pressure guards remain
responsible for live memory decisions. Neither admission nor this estimate is a
peak-memory guarantee. Inspect `model.get_execution_policy()` for the selected
policy, measurements, estimate and uncertainty.

`manifest.json` pins the original and resulting file hashes and helper payload.
The installer checks the source revision, stages and verifies the change before
replacing either file, and updates the distribution's `RECORD`. Reapplying an
identical patch is safe. A differing or locally edited dependency is rejected
without applying another patch. Installation does not modify scientific losses,
weights, precision, recycling, checkpointing, filters, or the model kernels.

The draft validation covers installation, corruption rejection, memory-accounting
cases, preservation of CPU/CUDA preparation, and both callable selections across
the old 384-residue boundary. See [validation](VALIDATION.md) for the real Metal
checks and remaining application gate.
