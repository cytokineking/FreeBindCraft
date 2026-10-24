"""Guarded balanced regression probe for dynamic slice-before-convert profitability.

Reproduces the review's (16,1024,1024) fixture, three warmups and twenty
synchronized calls per window in each fresh process. These synthetic operation
timings do not predict AF2 or trajectory throughput.
"""
import argparse
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import sys
import time

import numpy as np
from guard_metaljax_job import environment, run_guarded
from run_metaljax_slice_comparison import identity

ROOT = Path(__file__).resolve().parents[1]


def worker(destination, windows):
    import jax
    import jax.numpy as jnp
    from jax import lax

    destination.mkdir()
    x = jax.device_put(np.random.default_rng(117).normal(size=(16, 1024, 1024)).astype(np.float32))
    indices = [jax.device_put(np.int32(i)) for i in [0, 1]]
    report = {}
    for selected in windows:
        fn = jax.jit(lambda a, i: lax.dynamic_slice(
            a.astype(jnp.bfloat16), (i, 0, 0), (selected, 1024, 1024),
            allow_negative_indices=False)).lower(x, indices[0]).compile()
        text = fn.as_text()
        (destination/f'optimized-{selected}.txt').write_text(text)
        shapes = re.findall(r'= bf16\[([\d,]+)\].*?\bconvert\(', text)
        rewritten = shapes == [f'{selected},1024,1024']
        expected = os.environ['METALJAX_SLICE_BEFORE_CONVERT'] == '1' and selected == 1
        if rewritten != expected:
            raise RuntimeError(f'Unexpected admission: selected={selected}, shapes={shapes}')
        for i in range(3):
            fn(x, indices[i % 2]).block_until_ready()
        elapsed = []
        for i in range(20):
            start = time.perf_counter()
            result = fn(x, indices[i % 2])
            result.block_until_ready()
            elapsed.append(time.perf_counter()-start)
        np.save(destination/f'bits-{selected}.npy',
                np.asarray(lax.bitcast_convert_type(result, jnp.uint16)))
        report[str(selected)] = dict(times=elapsed, median=statistics.median(elapsed),
                                     rewritten=rewritten, converted_shapes=shapes)
    (destination/'report.json').write_text(json.dumps(report, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python', type=Path, default=ROOT/'.venv-metaljax-slice-v5/bin/python')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--windows', choices=['1,15', '15,1', '15', '1'], default='1,15',
                        help='Fixed review windows; isolate or reverse to check preceding-work effects')
    parser.add_argument('--order', choices=['ABBA', 'BAAB'], default='ABBA')
    args = parser.parse_args()
    out, python = args.out.absolute(), args.python.absolute()
    out.mkdir(parents=True, exist_ok=False)
    ident = identity(python)
    if any(ident['packages'][p] != '0.11.2' for p in ['jax', 'jaxlib']):
        raise RuntimeError('Requires the frozen JAX/jaxlib pins')
    if ident['files']['metaljax/lib/mlx/lib/libmlx_metaljax.dylib'] != '72dc8d9a4b1e07df2f6bc91383e0567acf3ce292398ae6ce4ed1d19fce0b1374':
        raise RuntimeError('Requires the frozen MLX runtime')
    (out/'identity.json').write_text(json.dumps(ident, indent=2)+'\n')
    counts, rows = dict(A=0, B=0), []
    windows = [int(n) for n in args.windows.split(',')]
    first_name = None
    for arm in args.order:
        counts[arm] += 1
        name = arm.lower()+str(counts[arm])
        env = environment(arm == 'B')
        env['METALJAX_EVO_ATTN'] = '0'
        command = [python, ROOT/'extras/monitor_macos_job.py', '--output-dir', out/(name+'-memory'),
                   '--interval', '.25', '--max-tree-rss-gib', '24', '--max-hours', '.05', '--',
                   python, Path(__file__).resolve(), '--worker', out/name, args.windows]
        run_guarded(command, env, out/(name+'-guard'), max_seconds=120)
        report = json.loads((out/name/'report.json').read_text())
        if first_name is None:
            first_name = name
        for selected in windows:
            first = np.load(out/first_name/f'bits-{selected}.npy')
            current = np.load(out/name/f'bits-{selected}.npy')
            if (first.dtype != current.dtype or first.shape != current.shape or
                first.tobytes() != current.tobytes()):
                raise RuntimeError(f'Output bits differ: {name}, selected={selected}')
            if selected == 15 and ((out/first_name/'optimized-15.txt').read_text() !=
                                   (out/name/'optimized-15.txt').read_text()):
                raise RuntimeError(f'Near-full fallback HLO differs: {name}')
        rows.append(dict(name=name, arm=arm, report=report))
        print(name, {k: v['median'] for k, v in report.items()}, 'bytes equal', flush=True)
    summary = dict(order=args.order, windows=windows, results=rows, all_bytes_equal=True,
                   units='seconds; fresh-process median of twenty calls after three warmups')
    for selected in map(str, windows):
        values = {arm: [r['report'][selected]['median'] for r in rows if r['arm'] == arm] for arm in 'AB'}
        a, b = statistics.median(values['A']), statistics.median(values['B'])
        summary[selected] = dict(**values, A_process_median=a, B_process_median=b,
                                 relative_seconds_change=b/a-1)
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    if len(sys.argv) == 4 and sys.argv[1] == '--worker':
        worker(Path(sys.argv[2]), [int(n) for n in sys.argv[3].split(',')])
    else:
        main()
