#!/usr/bin/env python3
"""Check a complete AF2 binder loss/gradient on GPU against finite differences.

This is an optional scientific validation, not a design profile. FP32 and fixed
RNG/dropout make the derivative check interpretable; production profiles retain
ColabDesign's normal precision and dropout settings. Run from the repository:
JAX_PLATFORMS=metal .venv-metaljax/bin/python extras/check_metal_gradient.py
"""
import argparse
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import json
import time
import numpy as np
import jax
import jax.numpy as jnp
from colabdesign import mk_afdesign_model
from colabdesign.af.alphafold.common import residue_constants
from functions.generic_utils import check_jax_gpu
from functions.colabdesign_utils import add_rg_loss, add_i_ptm_loss, add_helix_loss

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path,
                    default=root/'runs/native-validation/gradient-report.json')
args = parser.parse_args()
check_jax_gpu()
if jax.default_backend().lower() != 'metal':
    raise RuntimeError('This validation requires the Metal GPU')
model = mk_afdesign_model(protocol='binder', use_multimer=True,
                         data_dir=str(root), num_recycles=0, use_bfloat16=False,
                         model_names=['model_1_multimer_v3'])
model.prep_inputs(pdb_filename=str(root/'example/PDL1.pdb'), chain='A',
                  binder_len=65, hotspot='56', seed=17, rm_aa='C')
# Disable the intentional straight-through PSSM estimator only for this check.
model.opt.update(dropout=False, soft=1.0, hard=0.0, temp=1.0, pssm_hard=False)
# Start from a benchmark-derived sequence with smooth softmax inputs.
# Argmax amino-acid identities remain fixed under the small perturbations.
sequence = 'SMSEREHLRRMADSIIENMERVLGSVPEEVKKAVEEAIEKEEDPDEAIVKVMQILGDYLREKTGS'
indices = np.array([residue_constants.restype_order[aa] for aa in sequence])
logits = 2.0 * np.eye(20, dtype=np.float32)[indices][None]
model._params['seq'] = jnp.asarray(logits)
model.opt['weights'].update(plddt=0.1, pae=0.4, i_pae=0.1, con=1., i_con=1.)
model.opt['con'].update(num=2, cutoff=14., binary=False, seqsep=9)
model.opt['i_con'].update(num=2, cutoff=20., binary=False)
add_rg_loss(model, 0.3)
add_i_ptm_loss(model, 0.05)
add_helix_loss(model, -0.3)
length = model._inputs['residue_index'].shape[0]
model._inputs['opt'] = model.opt
model._inputs['prev'] = dict(prev_pos=np.zeros((length,37,3),np.float32),
                            prev_msa_first_row=np.zeros((length,256),np.float32),
                            prev_pair=np.zeros((length,length,128),np.float32))
key = jax.random.PRNGKey(17)
flags = [model._model_params[0], model._inputs, key]
start = time.monotonic()
# Match the production one-recycle design path. ColabDesign differentiates only
# the final pass, so keep these previous-state tensors fixed during differences.
_, warmup = model._model['fn'](model._params, *flags)
model._inputs['prev'] = warmup['prev']
(loss, aux), grads = model._model['grad_fn'](model._params, *flags)
grad = np.asarray(grads['seq'])
assert np.isfinite(grad).all() and np.linalg.norm(grad) > 0
assert np.isfinite(np.asarray(aux['atom_positions'])).all()
direction = jnp.asarray(grad / np.linalg.norm(grad))
analytic = float(jnp.sum(grads['seq']*direction))
baseline, _ = model._model['fn'](model._params, *flags)
print('Gradient/forward baseline:', float(loss), float(baseline), flush=True)
np.testing.assert_allclose(float(baseline), float(loss), rtol=1e-4, atol=1e-4)
def scalar_loss(seq):
    return model._model['fn']({**model._params, 'seq': seq}, *flags)[0]
_, forward_derivative = jax.jvp(scalar_loss, (model._params['seq'],), (direction,))
print('Reverse/forward AD:', analytic, float(forward_derivative), flush=True)
np.testing.assert_allclose(float(forward_derivative), analytic, rtol=1e-4, atol=1e-6)
results = []
for step in (0.01, 0.003, 0.001, 0.0003):
    values = []
    for sign in (1, -1):
        params = {**model._params, 'seq':model._params['seq']+sign*step*direction}
        np.testing.assert_array_equal(np.asarray(params['seq']).argmax(-1), indices[None])
        value, _ = model._model['fn'](params, *flags)
        values.append(float(value))
    difference = (values[0]-values[1])/(2*step)
    relative_error = abs(difference-analytic)/max(abs(analytic),1e-8)
    results.append(dict(step=step, values=values, finite_difference=difference, relative_error=relative_error))
    print(results[-1], flush=True)
report = dict(backend=jax.default_backend(), jax=jax.__version__, dtype='float32',
              fixed_previous_passes=1,
              loss=float(loss), gradient_norm=float(np.linalg.norm(grad)),
              analytic_directional_derivative=analytic, checks=results,
              forward_directional_derivative=float(forward_derivative),
              elapsed_seconds=time.monotonic()-start)
print(json.dumps(report, indent=2))
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(report, indent=2)+'\n')
# Finite differences through a deep nonlinear float32 network are approximate.
assert min(r['relative_error'] for r in results) < 0.02, report
