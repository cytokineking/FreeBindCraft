"""Require exact per-call agreement between direct Metal hooks and MLX counters.

No expected constant dispatch count or widened tolerance is used. Call arrays
must be present; comparison with another arm is separate from this accounting.
"""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np


def equal(left, right):
    with np.load(left) as a, np.load(right) as b:
        return set(a.files) == set(b.files) and all(
            a[k].shape == b[k].shape and a[k].dtype == b[k].dtype and
            a[k].tobytes() == b[k].tobytes() for k in a.files)


def check(trace, report_path):
    trace, report_path = Path(trace), Path(report_path)
    report = json.loads(report_path.read_text())
    snapshots = json.loads(trace.with_suffix('.snapshots.json').read_text())
    if snapshots['pending_trace_callbacks'] != 0 or not snapshots.get('mlx_image'):
        raise ValueError('Trace callbacks or counter image identity are unverified')
    direct, committed, kernels = Counter(), Counter(), Counter()
    for line in trace.open():
        row = json.loads(line)
        if row['type'] == 'buffer':
            if not row['matched'] or row['dispatches'] != row['native_operations']:
                raise ValueError('Direct dispatches and native encoder prefix differ')
            direct[row['phase']] += row['dispatches']
            committed[row['phase']] += 1
        elif row['type'] == 'kernel':
            kernels[row['phase']] += row['count']
    rows = []
    for kind, field in [('gradient', 'gradient_dispatch'), ('forward', 'forward_dispatch')]:
        for index, native in enumerate(report[field]):
            phase = f'{kind}-{index}'
            if (direct[phase] != native['dispatches'] or direct[phase] != kernels[phase] or
                committed[phase] != native['command_buffers'] or
                native['completed_command_buffers'] != native['command_buffers']):
                raise ValueError(f'{phase}: direct/native accounting mismatch: '
                                 f'{direct[phase]}, {committed[phase]}, {native}')
            call = report_path.parent/f'{report["name"]}-{kind}-call-{index}.npz'
            first = report_path.parent/f'{report["name"]}-{kind}-call-0.npz'
            same = equal(first, call)
            if not same:
                raise ValueError(f'{phase}: call arrays differ')
            rows.append(dict(phase=phase, dispatches=direct[phase], buffers=committed[phase],
                             arrays_byte_equal_first_call=same))
    delays = []
    for boundary in snapshots['boundaries']:
        first, last = boundary['snapshots'][0], boundary['snapshots'][-1]
        delays.append({k: last[k]-first[k] for k in
                       ['dispatches', 'command_buffers', 'completed_command_buffers']})
    return dict(all_direct_native_counts_equal=True, calls=rows,
                immediate_to_settled_counter_changes=delays,
                mlx_image=snapshots['mlx_image'], mlx_sha256=snapshots['mlx_sha256'],
                limitation='Fresh diagnostic calls only; historical calls lack direct traces and call arrays.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trace', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = check(args.trace, args.report)
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
