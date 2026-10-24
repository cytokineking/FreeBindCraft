"""Balanced fresh-process comparisons of the opt-in slice-before-convert build.

Original 800-operation scheduling and the same attention setting apply to both arms.
Normal timings do not read dispatch counters or load diagnostic hooks.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from guard_metaljax_job import environment, shipping_environment, run_guarded
from check_metaljax_dispatch_accounting import check as check_dispatch


def equal_archive(left, right):
    with np.load(left) as a, np.load(right) as b:
        return dict(equal=set(a.files) == set(b.files) and all(
            a[k].shape == b[k].shape and a[k].dtype == b[k].dtype and
            a[k].tobytes() == b[k].tobytes() for k in a.files),
            fields=len(a.files), reference_sha256=hashlib.sha256(Path(left).read_bytes()).hexdigest(),
            candidate_sha256=hashlib.sha256(Path(right).read_bytes()).hexdigest())


def identity(python):
    code = """import hashlib,importlib.metadata as m,json
d=m.distribution('metaljax')
files=['metaljax/lib/libmetal_pjrt_native.dylib','metaljax/lib/mlx/lib/libmlx_metaljax.dylib']
print(json.dumps(dict(packages={d.metadata['Name'].lower():d.version for d in m.distributions()},
files={p:hashlib.sha256(d.locate_file(p).read_bytes()).hexdigest() for p in files})))"""
    return json.loads(subprocess.check_output([str(python), '-c', code], text=True))


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-python', type=Path, default=root/'.venv-metaljax-baseline/bin/python')
    parser.add_argument('--candidate-python', type=Path, default=root/'.venv-metaljax-slice/bin/python')
    parser.add_argument('--baseline-manifest', type=Path, required=True)
    parser.add_argument('--candidate-manifest', type=Path, required=True)
    parser.add_argument('--assets', type=Path, required=True)
    parser.add_argument('--previous', type=Path)
    parser.add_argument('--reference', type=Path, help='Saved prefix without .npz; omit only for an additional fixture')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--order', default='ABBA')
    parser.add_argument('--repeats', type=int, default=6)
    parser.add_argument('--pdb', type=Path)
    parser.add_argument('--hotspot', default='44')
    parser.add_argument('--binder-length', type=int, default=70)
    parser.add_argument('--diagnostic-trace', type=Path, help='Separate diagnostic run, never normal timing evidence')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--attention', choices=['0', '1'],
                        help='Use the same attention setting in both arms')
    parser.add_argument('--shipping-defaults', action='store_true',
                        help='Leave slice, attention, precision, governor and scheduling variables unset')
    args = parser.parse_args()
    if (not args.order or set(args.order)-{'A', 'B'} or args.repeats < 1 or
        ((args.reference is None or args.previous is None) and args.order[0] != 'A')):
        parser.error('Need an A/B order and positive repeats; creating a new fixture requires A first')
    out = args.out.absolute()
    out.mkdir(parents=True, exist_ok=False)
    pythons = dict(A=args.baseline_python.absolute(), B=args.candidate_python.absolute())
    identities = {arm: identity(python) for arm, python in pythons.items()}
    manifests = dict(A=json.loads(args.baseline_manifest.read_text()),
                     B=json.loads(args.candidate_manifest.read_text()))
    for arm, ident in identities.items():
        if ident['files']['metaljax/lib/libmetal_pjrt_native.dylib'] != manifests[arm]['plugin_sha256']:
            raise RuntimeError('Installed plugin differs from selected build manifest')
        if ident['files']['metaljax/lib/mlx/lib/libmlx_metaljax.dylib'] != '72dc8d9a4b1e07df2f6bc91383e0567acf3ce292398ae6ce4ed1d19fce0b1374':
            raise RuntimeError('MLX runtime differs from the frozen input')
        if any(ident['packages'][k] != '0.11.2' for k in ['jax', 'jaxlib']):
            raise RuntimeError('JAX pins differ')
    packages = [{k: v for k, v in identities[arm]['packages'].items() if k != 'metaljax'} for arm in 'AB']
    if packages[0] != packages[1]:
        raise RuntimeError('Scientific environments differ')
    (out/'identities.json').write_text(json.dumps(dict(installed=identities, builds=manifests), indent=2)+'\n')
    counts, results = dict(A=0, B=0), []
    previous = args.previous.absolute() if args.previous else None
    reference = args.reference.absolute() if args.reference else None
    for arm in args.order:
        counts[arm] += 1
        name = arm.lower()+str(counts[arm])
        python = pythons[arm]
        env = shipping_environment() if args.shipping_defaults else environment(arm == 'B')
        if args.attention is not None:
            env['METALJAX_EVO_ATTN'] = args.attention
        if args.debug:
            env['METALJAX_DEBUG'] = '1'
        probe = 'benchmark_metaljax_attention.py'
        if args.diagnostic_trace:
            probe = 'diagnose_metaljax_dispatch.py'
            env.update(FBC_MLX_TRACE_LIBRARY=str(args.diagnostic_trace.absolute()),
                       FBC_MLX_TRACE_OUTPUT=str(out/(name+'-trace.jsonl')),
                       FBC_MLX_TRACE_COPIES='1', FBC_MLX_TRACE_REFLECTION='1')
        command = [python, root/'extras/monitor_macos_job.py', '--output-dir', out/(name+'-memory'),
                   '--interval', '.25', '--max-tree-rss-gib', '24', '--max-hours', '.25', '--',
                   python, root/'extras'/probe, '--assets', args.assets.absolute(), '--out', out,
                   '--name', name, '--repeats', str(args.repeats), '--binder-length', str(args.binder_length),
                   '--hotspot', args.hotspot]
        if previous:
            command += ['--previous', previous]
        if args.pdb:
            command += ['--pdb', args.pdb.absolute()]
        if args.diagnostic_trace:
            command += ['--profile', '--diagnostic-arrays', '--dump-hlo']
        else:
            command += ['--no-dispatch-stats']
        run_guarded(command, env, out/(name+'-guard'), max_seconds=960)
        if args.diagnostic_trace:
            accounting = check_dispatch(out/(name+'-trace.jsonl'), out/(name+'.json'))
            (out/(name+'-accounting.json')).write_text(json.dumps(accounting, indent=2)+'\n')
        if reference is None:
            reference = out/name
        if previous is None:
            previous = out/(name+'-previous.npz')
        numeric = {suffix: equal_archive(Path(str(reference)+suffix+'.npz'), out/(name+suffix+'.npz'))
                   for suffix in ['', '-warm-forward', '-previous']}
        if not all(n['equal'] for n in numeric.values()):
            (out/(name+'-numerical-failure.json')).write_text(json.dumps(numeric, indent=2))
            raise RuntimeError(f'Byte comparison failed: {name}: {numeric}')
        report = json.loads((out/(name+'.json')).read_text())
        memory = json.loads((out/(name+'-memory')/'summary.json').read_text())
        results.append(dict(arm=arm, name=name, numeric=numeric, report=report,
                            memory={k: v for k, v in memory.items() if k.startswith(('peak_', 'maximum_'))}))
        summary = dict(order=args.order, results=results, diagnostic=bool(args.diagnostic_trace or args.debug),
                       shipping_defaults=args.shipping_defaults,
                       sampling_unit='fresh process; warmed calls within a process are correlated',
                       limitations='Shared Mac; guards do not establish GPU exclusivity or thermal control.')
        for key in ['first_forward_seconds', 'first_grad_seconds', 'median_forward_seconds', 'median_grad_seconds']:
            summary[key] = {a: [r['report'][key] for r in results if r['arm'] == a] for a in 'AB'}
        (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
        print(name, report['median_grad_seconds'], report['median_forward_seconds'], 'byte equal', flush=True)


if __name__ == '__main__':
    main()
