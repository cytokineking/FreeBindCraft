"""Fixed-input forward/VJP probe for the disabled attention experiment.

No optimizer, schedule, precision, recycling, or checkpointing changes. Use the
baseline's --previous file for a same-state VJP comparison. Timers include host
materialization; first calls are reported separately from warmed calls.
"""
import argparse
import ast
import ctypes
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import resource
import time

import jax
import jax.numpy as jnp
import numpy as np
from colabdesign import mk_afdesign_model
from colabdesign.af.alphafold.common import residue_constants
from colabdesign.af.loss import get_ptm, mask_loss, get_dgram_bins, _get_con_loss
from monitor_macos_job import read_footprint


class DispatchStats(ctypes.Structure):
    # POD ABI of the unchanged private MLX runtime, metal/metal.h.
    _fields_ = [(key, ctypes.c_uint64) for key in (
        'dispatches command_buffers empty_command_buffers completed_command_buffers '
        'gpu_busy_ns gpu_span_ns gpu_gap_ns queue_ns now_ns '
        'slice_update_donated slice_update_copied').split()]


def dispatch_reader():
    path = importlib.metadata.distribution('metaljax').locate_file(
        'metaljax/lib/mlx/lib/libmlx_metaljax.dylib')
    library = ctypes.CDLL(str(path))
    read = getattr(library, '_ZN3mlx4core5metal14dispatch_statsEv')
    read.argtypes, read.restype = [], DispatchStats
    def snapshot():
        end = time.perf_counter()+0.05
        value = read()
        while value.completed_command_buffers < value.command_buffers and time.perf_counter() < end:
            time.sleep(0.0001)
            value = read()
        return {key: getattr(value, key) for key, _ in DispatchStats._fields_}
    return snapshot


def digest(array):
    a = np.asarray(array)
    return hashlib.sha256(a.dtype.str.encode() + str(a.shape).encode() + a.tobytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assets', type=Path, required=True, help='Checkout containing params/ and example/')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--previous', type=Path, help='Frozen baseline first-recycle state (.npz)')
    parser.add_argument('--repeats', type=int, default=6)
    parser.add_argument('--name', required=True)
    parser.add_argument('--dump-hlo', action='store_true')
    parser.add_argument('--profile', action='store_true', help='Enable runtime profiling fences only in warmed calls')
    parser.add_argument('--no-dispatch-stats', action='store_true',
                        help='Normal execution timings without private runtime counter reads')
    parser.add_argument('--diagnostic-arrays', action='store_true',
                        help='Save every warmed result for counter diagnosis; excludes performance use')
    parser.add_argument('--pdb', type=Path, help='Target PDB; defaults to the bundled ubiquitin fixture')
    parser.add_argument('--chain', default='A')
    parser.add_argument('--binder-length', type=int, default=70)
    parser.add_argument('--hotspot', default='44')
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    if args.binder_length < 1:
        parser.error('--binder-length must be positive')
    if jax.default_backend() != 'metal':
        parser.error('This probe requires JAX_PLATFORMS=metal')
    jax.config.update('jax_enable_compilation_cache', False)
    root = Path(__file__).resolve().parents[1]
    args.out.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location('fbc_probe_metal', root/'functions/metal_utils.py')
    metal = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(metal)
    metal.enable_metal_alignment()

    # The unchanged production loss callbacks, without native scoring imports.
    source = ast.parse((root/'functions/colabdesign_utils.py').read_text())
    selected = [n for n in source.body if isinstance(n, ast.FunctionDef) and n.name in
                {'add_rg_loss', 'add_i_ptm_loss', 'add_helix_loss'}]
    callbacks = dict(globals())
    exec(compile(ast.Module(body=selected, type_ignores=[]), 'production_loss_callbacks', 'exec'), callbacks)

    start = time.perf_counter()
    model = mk_afdesign_model(protocol='binder', use_multimer=True,
                             data_dir=str(args.assets), num_recycles=1, use_remat=True,
                             model_names=['model_1_multimer_v3'])
    target = args.pdb or args.assets/'example/ubiquitin_1ubq.pdb'
    model.prep_inputs(pdb_filename=str(target), chain=args.chain,
                      binder_len=args.binder_length, hotspot=args.hotspot, seed=42, rm_aa='C')
    model.opt.update(dropout=True, soft=0.9, hard=0.0, temp=1.0)
    model.opt['weights'].update(plddt=0.1, pae=0.4, i_pae=0.1, con=1., i_con=1.)
    model.opt['con'].update(num=2, cutoff=14., binary=False, seqsep=9)
    model.opt['i_con'].update(num=2, cutoff=20., binary=False)
    for name, weight in [('add_rg_loss', 0.3), ('add_i_ptm_loss', 0.05), ('add_helix_loss', -0.3)]:
        callbacks[name](model, weight)
    length = int(model._inputs['residue_index'].shape[0])
    model._inputs['opt'] = model.opt
    model._inputs['prev'] = dict(prev_pos=np.zeros((length,37,3),np.float32),
                                prev_msa_first_row=np.zeros((length,256),np.float32),
                                prev_pair=np.zeros((length,length,128),np.float32))
    key = jax.random.PRNGKey(42)

    def call(fn):
        return fn(model._params, model._model_params[0], model._inputs, key)

    def settle(tree):
        return jax.tree_util.tree_map(np.asarray, tree)

    setup_seconds = time.perf_counter()-start
    t = time.perf_counter()
    forward_loss, forward = settle(call(model._model['fn']))
    first_forward_seconds = time.perf_counter()-t
    np.savez(args.out/f'{args.name}-previous.npz', **forward['prev'])
    if args.previous:
        with np.load(args.previous) as previous:
            model._inputs['prev'] = {k: previous[k] for k in previous.files}
    else:
        model._inputs['prev'] = forward['prev']
    previous_hashes = {k: digest(v) for k, v in model._inputs['prev'].items()}
    if args.dump_hlo:
        (args.out/f'{args.name}-forward.mlir').write_text(model._model['fn'].lower(model._params, model._model_params[0], model._inputs, key).as_text())
        lowered = model._model['grad_fn'].lower(model._params, model._model_params[0], model._inputs, key)
        (args.out/f'{args.name}-gradient.mlir').write_text(lowered.as_text())

    t = time.perf_counter()
    first = settle(call(model._model['grad_fn']))
    first_grad_seconds = time.perf_counter()-t
    times = []
    read_dispatch = (lambda: {}) if args.no_dispatch_stats else dispatch_reader()
    gradient_dispatch, forward_dispatch = [], []
    for i in range(args.repeats):
        if args.profile: os.environ['FBC_PROFILE_ACTIVE']=f'gradient-{i}'
        before = read_dispatch()
        t = time.perf_counter()
        (loss,aux),grad = settle(call(model._model['grad_fn']))
        times.append(time.perf_counter()-t)
        after = read_dispatch()
        gradient_dispatch.append({k: after[k]-before[k] for k in before})
        if args.diagnostic_arrays:
            leaves = jax.tree_util.tree_leaves(((loss, aux), grad))
            np.savez(args.out/f'{args.name}-gradient-call-{i}.npz',
                     **{f'leaf_{n}': v for n, v in enumerate(leaves)})
    forward_times=[]
    for i in range(args.repeats):
        if args.profile: os.environ['FBC_PROFILE_ACTIVE']=f'forward-{i}'
        before = read_dispatch()
        t=time.perf_counter()
        warm_loss,warm=settle(call(model._model['fn']))
        forward_times.append(time.perf_counter()-t)
        after = read_dispatch()
        forward_dispatch.append({k: after[k]-before[k] for k in before})
        if args.diagnostic_arrays:
            leaves = jax.tree_util.tree_leaves((warm_loss, warm))
            np.savez(args.out/f'{args.name}-forward-call-{i}.npz',
                     **{f'leaf_{n}': v for n, v in enumerate(leaves)})
    os.environ['FBC_PROFILE_ACTIVE']='0'
    np.savez(args.out/f'{args.name}-warm-forward.npz',loss=warm_loss,xyz=warm['atom_positions'],plddt=warm['plddt'])
    assert np.isfinite(loss) and np.isfinite(grad['seq']).all()
    assert np.isfinite(aux['atom_positions']).all()
    assert previous_hashes == {k: digest(v) for k, v in model._inputs['prev'].items()}
    np.savez(args.out/f'{args.name}.npz', loss=loss, grad=grad['seq'], xyz=aux['atom_positions'],
             plddt=aux['plddt'], initial_seq=np.asarray(model._params['seq']),
             forward_loss=forward_loss, forward_xyz=forward['atom_positions'], forward_plddt=forward['plddt'],
             **{f'loss_{k}': v for k, v in aux['losses'].items()})
    report = dict(name=args.name, backend=jax.default_backend(), jax=jax.__version__,
        metaljax=importlib.metadata.version('metaljax'),
        precision=os.environ.get('METALJAX_MATMUL_PRECISION','high'),
        plugin=os.environ.get('METALJAX_PLUGIN_PATH', 'installed'),
        use_bfloat16=model._args['use_bfloat16'], use_remat=model._args['use_remat'],
        subbatch=model._cfg.model.global_config.subbatch_size, length=length,
        target_pdb_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        chain=args.chain, binder_length=args.binder_length, hotspot=args.hotspot,
        model='model_1_multimer_v3', fixed_rng=True, dropout=True,
        sequence_sha256=digest(model._params['seq']), key_sha256=digest(key), previous_sha256=previous_hashes,
        setup_seconds=setup_seconds, first_forward_seconds=first_forward_seconds,
        first_grad_seconds=first_grad_seconds, grad_seconds=times, median_grad_seconds=float(np.median(times)),
        loss=float(loss), grad_norm=float(np.linalg.norm(grad['seq'])),
        process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        **read_footprint(os.getpid()),
        forward_seconds=forward_times, median_forward_seconds=float(np.median(forward_times)),
        gradient_dispatch=gradient_dispatch, forward_dispatch=forward_dispatch,
        attention=os.environ.get('METALJAX_EVO_ATTN','unset'), profiling=args.profile,
        dispatch_stats_enabled=not args.no_dispatch_stats,
        diagnostic_arrays=args.diagnostic_arrays,
        slice_before_convert=os.environ.get('METALJAX_SLICE_BEFORE_CONVERT','unset'),
        backend_environment={k:os.environ.get(k) for k in ['METALJAX_SLICE_BEFORE_CONVERT',
            'METALJAX_EVO_ATTN','METALJAX_MATMUL_PRECISION','METALJAX_MEM_BUDGET_MB',
            'MLX_MAX_OPS_PER_BUFFER','MLX_MAX_MB_PER_BUFFER','METALJAX_LOOP_INFLIGHT_MB']},
        persistent_jax_cache=False,
        plugin_sha256=hashlib.sha256(Path(importlib.metadata.distribution('metaljax').locate_file('metaljax/lib/libmetal_pjrt_native.dylib')).read_bytes()).hexdigest(),
        scope='Fixed-input AF2 forward and last-recycle VJP, not a trajectory or quality benchmark')
    (args.out/f'{args.name}.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
