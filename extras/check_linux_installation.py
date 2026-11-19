#!/usr/bin/env python3
"""Check the installed CUDA runtime with real JAX gradients and OpenMM energy.

The numerical check catches the non-finite GEMM gradients observed with an
older Conda CUDA stack on RTX 5090. Run in the activated BindCraft environment.
"""
import argparse
import importlib.util
import json
import math
import os
import subprocess
import sys
import tempfile
from importlib.metadata import version
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.append(str(REPO_ROOT))

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--no-pyrosetta', action='store_true')
parser.add_argument('--no-tmol', action='store_true')
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

    dalphaball_bin = Path(os.environ.get('DALPHABALL_BIN') or (REPO_ROOT / 'functions' / 'DAlphaBall.gcc'))
    if not dalphaball_bin.is_file():
        raise RuntimeError(f'Missing DAlphaBall executable required by PyRosetta scoring: {dalphaball_bin}')
    spheres = 'NPOINTS\n2\nCOORDS\n0 0 0 1\n5 0 0 1\nEND\n'
    try:
        db_proc = subprocess.run(
            [str(dalphaball_bin), 'surf_vol'],
            input=spheres,
            text=True,
            capture_output=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError) as error:
        stderr = (getattr(error, 'stderr', '') or '').strip()
        raise RuntimeError(
            f'DAlphaBall executable failed ({dalphaball_bin}): {stderr or error}. '
            'On bare Linux containers, install system libgfortran5 (e.g. apt-get install -y libgfortran5) '
            'because activation unsets LD_LIBRARY_PATH, or run with --no-pyrosetta / --backend tmol.'
        ) from error
    rows = [list(map(float, row.split())) for row in db_proc.stdout.splitlines() if row.strip()]
    if len(rows) != 2 or not all(
        math.isclose(r[1], 4 * math.pi, rel_tol=1e-9) and math.isclose(r[2], 4 * math.pi / 3, rel_tol=1e-9)
        for r in rows
    ):
        raise RuntimeError(f'Unexpected DAlphaBall sphere output: {db_proc.stdout}')
    report['dalphaball_checked'] = True

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

if not args.no_tmol:
    from functions.tmol_utils import resolve_tmol_python

    tmol_python = resolve_tmol_python()
    worker = str(REPO_ROOT / 'functions' / 'tmol_worker.py')
    test_pdb = str(REPO_ROOT / 'example' / 'PDL1.pdb')
    device_arg = 'cpu' if args.skip_gpu_check else 'cuda'
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as tf:
        result_path = tf.name
    try:
        subprocess.run(
            [tmol_python, worker, 'check', '--pdb', test_pdb, '--device', device_arg, '--result-json', result_path],
            check=True,
            cwd=str(REPO_ROOT),
        )
        tmol_report = json.loads(Path(result_path).read_text())
    finally:
        if os.path.isfile(result_path):
            os.remove(result_path)
    if not math.isfinite(tmol_report.get('test_energy', float('nan'))):
        raise RuntimeError(f'tmol returned non-finite test energy: {tmol_report}')
    if not args.skip_gpu_check and not (tmol_report.get('tmol_has_cuda') and tmol_report.get('torch_cuda_available')):
        raise RuntimeError(f'The CUDA profile requires working tmol CUDA support: {tmol_report}')
    tmol_report['python'] = tmol_python
    report['tmol'] = tmol_report

print(json.dumps(report, indent=2))
