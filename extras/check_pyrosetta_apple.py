#!/usr/bin/env python3
"""Exercise native FastRelax, Rosetta interface scoring/DAlphaBall, and Metal MPNN.

Usage: JAX_PLATFORMS=metal .venv-pyrosetta/bin/python extras/check_pyrosetta_apple.py
Optional --pdb overrides the previously generated PDL1 candidate. This deliberately
uses CPU Rosetta relaxation while retaining GPU-only JAX model execution.
"""
import argparse
import io
import json
import logging
import os
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time
from importlib.metadata import version

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)
os.environ.setdefault('JAX_PLATFORMS', 'metal')
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'.local/matplotlib'))

import numpy as np
import jax
from functions.generic_utils import check_jax_gpu
from functions.pyrosetta_utils import PYROSETTA_AVAILABLE, pr, pr_relax, score_interface
from functions.colabdesign_utils import mpnn_gen_sequence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pdb', type=Path, default=ROOT/'runs/pdl1-apple/output/MPNN/PDL1_l65_s121958_mpnn1_model1.pdb')
    parser.add_argument('--output', type=Path, default=ROOT/'runs/pyrosetta-apple')
    parser.add_argument('--dalphaball', type=Path, default=ROOT/'.local/src/DAlphaBall/DAlphaBall.gcc')
    args = parser.parse_args()
    assert platform.machine() == 'arm64' and sys.platform == 'darwin'
    assert PYROSETTA_AVAILABLE and pr is not None, 'PyRosetta is required; no relaxation fallback is permitted'
    check_jax_gpu()
    assert jax.default_backend().lower() == 'metal'
    logging.basicConfig(level=logging.WARNING)
    logging.getLogger('functions').setLevel(logging.DEBUG)
    folder = args.output.resolve()
    folder.mkdir(parents=True, exist_ok=True)
    binary = args.dalphaball.resolve()
    assert binary.is_file() and os.access(binary, os.X_OK)
    arch = subprocess.check_output(['file', str(binary)], text=True).strip()
    assert 'arm64' in arch, arch

    # Two disjoint unit spheres have known areas/volumes: checks real native code.
    spheres = 'NPOINTS\n2\nCOORDS\n0 0 0 1\n5 0 0 1\nEND\n'
    result = subprocess.run([str(binary), 'surf_vol'], input=spheres, text=True,
                            capture_output=True, check=True)
    values = np.loadtxt(io.StringIO(result.stdout))
    np.testing.assert_allclose(values[:, 1], 4*np.pi, rtol=1e-10)
    np.testing.assert_allclose(values[:, 2], 4*np.pi/3, rtol=1e-10)

    # Record calls to the real executable without changing its input or output.
    calls = folder/'dalphaball-invocations.log'
    calls.write_text('')
    wrapper = folder/'dalphaball-trace.sh'
    wrapper.write_text('#!/bin/sh\nprintf "invoked\\n" >> '+shlex.quote(str(calls))+
                       '\nexec '+shlex.quote(str(binary))+' "$@"\n')
    wrapper.chmod(0o755)
    pr.init('-ignore_unrecognized_res -ignore_zero_occupancy -mute all '
            f'-holes:dalphaball {shlex.quote(str(wrapper))} '
            '-corrections::beta_nov16 true -relax:default_repeats 1 '
            '-constant_seed -jran 42')
    source = args.pdb.resolve()
    before = pr.pose_from_pdb(str(source))
    scorefxn = pr.get_fa_scorefxn()
    energy_before = float(scorefxn(before))
    relaxed = folder/'pdl1-pyrosetta-relaxed.pdb'
    start = time.monotonic()
    pr_relax(str(source), str(relaxed), use_pyrosetta=True)
    relax_seconds = time.monotonic()-start
    after = pr.pose_from_pdb(str(relaxed))
    assert after.total_residue() == before.total_residue() and after.sequence() == before.sequence()
    assert after.num_chains() == before.num_chains() == 2
    energy_after = float(scorefxn(after))
    assert np.isfinite([energy_before, energy_after]).all()
    assert energy_after < energy_before, (energy_before, energy_after)
    start = time.monotonic()
    metrics, interface_aa, interface_residues = score_interface(str(relaxed), 'B', use_pyrosetta=True)
    score_seconds = time.monotonic()-start
    assert np.isfinite(list(metrics.values())).all(), metrics
    assert metrics['interface_nres'] > 0 and 0 < metrics['interface_sc'] <= 1, metrics
    invocations = len(calls.read_text().splitlines())
    assert invocations > 0, 'DAlphaBall was not executed by Rosetta scoring'
    print('Native Rosetta relaxation and scoring passed.', flush=True)

    # Exercise a real model after loading and using PyRosetta in this process.
    settings = json.loads((ROOT/'settings_advanced/default_4stage_multimer.json').read_text())
    settings.update(num_seqs=1, max_mpnn_sequences=1)
    start = time.monotonic()
    mpnn = mpnn_gen_sequence(str(relaxed), 'B', interface_residues, settings)
    mpnn_seconds = time.monotonic()-start
    assert len(mpnn['seq']) == 1
    assert len(mpnn['seq'][0].replace('/', '')) == after.total_residue()
    binder_length = after.chain_end(2)-after.chain_begin(2)+1
    report = dict(pyrosetta=version('pyrosetta'), python=platform.python_version(),
                  architecture=platform.machine(), jax=version('jax'), jax_backend=jax.default_backend(),
                  input=str(source), relaxed=str(relaxed), residue_count=after.total_residue(),
                  relaxation_backend='PyRosetta FastRelax (CPU)', relax_seconds=relax_seconds,
                  total_energy_before_rosetta_units=energy_before,
                  total_energy_after_rosetta_units=energy_after,
                  score_seconds=score_seconds, interface=metrics,
                  dalphaball_binary=arch, dalphaball_invocations=invocations,
                  dalphaball_reference_spheres=values.tolist(),
                  mpnn_seconds=mpnn_seconds, mpnn_sequence=mpnn['seq'][0][-binder_length:])
    (folder/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
