#!/usr/bin/env python3
"""Check the installed CUDA runtime with real JAX gradients and OpenMM energy.

The numerical check catches the non-finite GEMM gradients observed with an
older Conda CUDA stack on RTX 5090. Run in the activated BindCraft environment.
"""
import argparse
import importlib.util
import json
import math
from importlib.metadata import version

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--no-pyrosetta', action='store_true')
parser.add_argument('--skip-gpu-check', action='store_true', help='Import checks only, for container image builds')
args = parser.parse_args()
driver_versions = []
if not args.skip_gpu_check:
    from check_linux_driver import check_driver
    driver_versions = check_driver()

import colabdesign
import freesasa
import jax
import jax.numpy as jnp
import numpy as np
import openmm
import pdbfixer

report = {'packages': {name: version(name) for name in
    ['jax', 'jaxlib', 'jax-cuda12-plugin', 'numpy', 'flax', 'colabdesign', 'openmm']},
    'pyrosetta_installed': importlib.util.find_spec('pyrosetta') is not None,
    'gpu_execution_checked': not args.skip_gpu_check, 'nvidia_driver_versions': driver_versions}
if not args.no_pyrosetta:
    import pyrosetta
    pyrosetta.init('-mute all')
    energy = pyrosetta.get_fa_scorefxn()(pyrosetta.pose_from_sequence('AAAAAA'))
    if not math.isfinite(energy):
        raise RuntimeError('PyRosetta returned a non-finite energy')
    report['pyrosetta_version'] = version('pyrosetta')

if not args.skip_gpu_check:
    if jax.default_backend() != 'gpu' or not all(d.platform == 'gpu' for d in jax.devices()):
        raise RuntimeError('The CUDA profile requires a working NVIDIA JAX GPU device')
    report['jax_devices'] = [str(d) for d in jax.devices()]
    matrix = jnp.arange(256 * 256, dtype=jnp.float32).reshape(256, 256) / 65536
    loss, gradient = jax.jit(jax.value_and_grad(lambda x: jnp.mean(jnp.sin(x @ x))))(matrix)
    host = np.asarray(matrix)
    product = host @ host
    derivative = np.cos(product) / product.size
    reference_gradient = derivative @ host.T + host.T @ derivative
    if not np.isfinite(float(loss)) or not np.isfinite(np.asarray(gradient)).all():
        raise RuntimeError('JAX produced non-finite GPU values or gradients')
    np.testing.assert_allclose(float(loss), float(np.sin(product).mean()), rtol=.005, atol=1e-5)
    np.testing.assert_allclose(np.asarray(gradient), reference_gradient, rtol=.01, atol=2e-5)
    report['jit_loss'] = float(loss)
    report['gradient_max_absolute_error'] = float(np.max(np.abs(np.asarray(gradient) - reference_gradient)))

    system = openmm.System()
    system.addParticle(12)
    system.addParticle(12)
    force = openmm.HarmonicBondForce()
    force.addBond(0, 1, .15, 100)
    system.addForce(force)
    integrator = openmm.VerletIntegrator(.001)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('CUDA'))
    context.setPositions([[0, 0, 0], [.2, 0, 0]])
    energy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(openmm.unit.kilojoule_per_mole)
    if not math.isclose(energy, .125, rel_tol=1e-4, abs_tol=1e-6):
        raise RuntimeError(f'Unexpected OpenMM CUDA energy: {energy}')
    report['openmm_platform'] = context.getPlatform().getName()
    report['openmm_energy_kj_mol'] = energy
    del context, integrator

print(json.dumps(report, indent=2))
