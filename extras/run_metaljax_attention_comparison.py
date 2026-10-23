"""Sequential fresh-process attention comparison with continuous contention logs.

Samples are correlated within each process. ABBA+BAAB gives four fresh processes
per arm; this does not independently control desktop/thermal/GPU activity.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import plistlib
import signal
import subprocess
import sys
import time

import numpy as np
import psutil


def snapshot():
    processes = []
    for process in psutil.process_iter(['pid', 'name', 'cpu_times', 'memory_info']):
        try:
            p = process.info
            if p['cpu_times'] is None or p['memory_info'] is None:
                continue
            processes.append(dict(pid=p['pid'], name=p['name'],
                                  cpu_seconds=p['cpu_times'].user+p['cpu_times'].system,
                                  rss_bytes=p['memory_info'].rss))
        except (psutil.AccessDenied, psutil.NoSuchProcess):
            pass
    result = subprocess.run(['ioreg', '-a', '-r', '-c', 'AGXAccelerator', '-d', '1'],
                            capture_output=True)
    try:
        gpu = [d.get('PerformanceStatistics', {}) for d in plistlib.loads(result.stdout)]
    except (ValueError, plistlib.InvalidFileException):
        gpu = []
    return dict(time=time.time(), processes=processes, gpu=gpu)


def activity_since(previous, current, own):
    old = {p['pid']: p for p in previous['processes']}
    dt = current['time']-previous['time']
    active, busy = [], []
    for process in current['processes']:
        if process['pid'] not in old or process['pid'] in own:
            continue
        cpu = (process['cpu_seconds']-old[process['pid']]['cpu_seconds'])/dt*100
        process = dict(process, cpu_percent=cpu)
        if cpu > 1:
            active.append(process)
        compute = process['name'].lower().startswith(('python', 'clang', 'java', 'ollama', 'blender'))
        if (compute and cpu > 25) or cpu > 150:
            busy.append(process)
    return active, busy, dt


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-python', type=Path, required=True)
    parser.add_argument('--candidate-python', type=Path, required=True)
    parser.add_argument('--assets', type=Path, required=True)
    parser.add_argument('--previous', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--order', default='ABBABAAB')
    parser.add_argument('--repeats', type=int, default=6)
    parser.add_argument('--pdb', type=Path)
    parser.add_argument('--hotspot', default='44')
    parser.add_argument('--binder-length', type=int, default=70)
    args = parser.parse_args()
    if set(args.order) - {'A', 'B'} or args.repeats < 1:
        parser.error('Use an A/B order and positive repeat count')
    args.out = args.out.absolute()
    if args.out.exists():
        parser.error('Choose a fresh output directory')
    args.out.mkdir(parents=True)
    env = dict(os.environ, JAX_PLATFORMS='metal', METALJAX_MATMUL_PRECISION='high',
               METALJAX_MEM_BUDGET_MB='16384', MPLBACKEND='Agg')
    for key in ['METALJAX_PLUGIN_PATH', 'METALJAX_DEBUG', 'METALJAX_TIMING',
                'FBC_ATTN_PROFILE', 'FBC_PROFILE_ACTIVE', 'METALJAX_EVO_ATTN_KERNEL',
                'METALJAX_EVO_ATTN_LOGITS']:
        env.pop(key, None)
    # Record all visible processes over time, without an age-based exemption.
    # A sustained compute process stops only our job; other workloads are kept.
    previous = snapshot()
    counts, results, baseline_fields = {'A': 0, 'B': 0}, [], None
    with (args.out/'contention.jsonl').open('w') as activity:
        activity.write(json.dumps(dict(stage='preflight', **previous))+'\n')
        for arm in args.order:
            counts[arm] += 1
            name = f'{arm.lower()}{counts[arm]}'
            # A bounded quiet-window check before every process. CPU sampling
            # cannot establish GPU exclusivity or controlled thermal state.
            deadline, quiet_since = time.monotonic()+180, time.monotonic()
            while time.monotonic()-quiet_since < 6:
                time.sleep(1)
                current = snapshot()
                active, busy, _ = activity_since(previous, current, {os.getpid()})
                activity.write(json.dumps(dict(stage=name+'-preflight', active=active, **current))+'\n')
                activity.flush()
                previous = current
                if busy:
                    quiet_since = time.monotonic()
                if time.monotonic() > deadline:
                    raise RuntimeError('No six-second CPU quiet window within 180 seconds; other workloads retained')
            python = (args.baseline_python if arm == 'A' else args.candidate_python).absolute()
            run_env = dict(env, METALJAX_EVO_ATTN='1' if arm == 'B' else '0')
            command = [str(python), str(root/'extras/monitor_macos_job.py'),
                       '--output-dir', str(args.out/(name+'-monitor')), '--interval', '1',
                       '--max-tree-rss-gib', '24', '--max-hours', '0.12', '--',
                       str(python), str(root/'extras/benchmark_metaljax_attention.py'),
                       '--assets', str(args.assets.absolute()), '--previous', str(args.previous.absolute()),
                       '--out', str(args.out), '--name', name, '--repeats', str(args.repeats),
                       '--hotspot', args.hotspot, '--binder-length', str(args.binder_length)]
            if args.pdb:
                command.extend(['--pdb', str(args.pdb.absolute())])
            sustained = {}
            reason = None
            with (args.out/(name+'-launch.log')).open('w') as log:
                child = subprocess.Popen(command, env=run_env, stdout=log,
                                         stderr=subprocess.STDOUT, start_new_session=True)
                while child.poll() is None:
                    time.sleep(1)
                    current = snapshot()
                    try:
                        own = {os.getpid(), *[p.pid for p in psutil.Process().children(recursive=True)]}
                    except psutil.NoSuchProcess:
                        own = {os.getpid(), child.pid}
                    active, busy, dt = activity_since(previous, current, own)
                    busy_ids = {p['pid'] for p in busy}
                    sustained = {pid: duration for pid, duration in sustained.items() if pid in busy_ids}
                    for p in busy:
                        sustained[p['pid']] = sustained.get(p['pid'], 0)+dt
                        if sustained[p['pid']] >= 4:
                            reason = f"Concurrent compute: {p['name']} pid {p['pid']}"
                    activity.write(json.dumps(dict(stage=name, active=active, **current))+'\n')
                    activity.flush()
                    previous = current
                    if reason:
                        # The memory monitor gives its native child a separate
                        # process group. Stop both of our groups, never a peer.
                        groups = [child.pid]
                        pid_file = args.out/(name+'-monitor')/'pid'
                        if pid_file.exists():
                            native_pid = int(pid_file.read_text())
                            if native_pid in own:
                                groups.insert(0, native_pid)
                        for group in groups:
                            try:
                                os.killpg(group, signal.SIGTERM)
                            except ProcessLookupError:
                                pass
                        try:
                            child.wait(timeout=15)
                        except subprocess.TimeoutExpired:
                            for group in groups:
                                try:
                                    os.killpg(group, signal.SIGKILL)
                                except ProcessLookupError:
                                    pass
                            child.wait()
                        break
                child.wait()
            if reason or child.returncode:
                (args.out/'failure.json').write_text(json.dumps(dict(reason=reason, returncode=child.returncode,
                                                                    command=command), indent=2)+'\n')
                raise RuntimeError(reason or f'Probe failed: {name}')
            report = json.loads((args.out/(name+'.json')).read_text())
            with np.load(args.out/(name+'.npz')) as data:
                fields = {key: data[key] for key in data.files}
            if baseline_fields is None:
                baseline_fields = fields
            equal = (fields.keys() == baseline_fields.keys() and
                     all(fields[key].shape == baseline_fields[key].shape and
                         fields[key].dtype == baseline_fields[key].dtype and
                         fields[key].tobytes() == baseline_fields[key].tobytes() for key in fields))
            if not equal:
                raise RuntimeError(f'Full-model byte comparison failed: {name}')
            with np.load(args.out/(name+'-warm-forward.npz')) as data:
                warm_fields = {key: data[key] for key in data.files}
            results.append(dict(arm=arm, name=name, command=command, report=report,
                                all_18_fields_bit_identical=equal,
                                warm_forward_hashes={k: hashlib.sha256(v.tobytes()).hexdigest()
                                                     for k, v in warm_fields.items()}))
            summary = dict(order=args.order, repeats_per_process=args.repeats, processes=results,
                           scope='Fresh processes are sampling units; repeated calls within one are correlated.',
                           contention_note='Continuous process and aggregate GPU samples are retained; desktop and thermal state are not controlled, and GPU activity is not attributed per process.')
            for key in ('grad_seconds', 'forward_seconds'):
                summary[key] = {a: [r['report'][key] for r in results if r['arm'] == a] for a in ('A', 'B')}
            (args.out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
            print(name, 'gradient', report['median_grad_seconds'],
                  'forward', report['median_forward_seconds'], 'exact', equal, flush=True)


if __name__ == '__main__':
    main()
