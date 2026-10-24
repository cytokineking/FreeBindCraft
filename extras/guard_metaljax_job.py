"""Bounded shared-Mac guard. Wait for CPU/GPU quiet, then monitor our own job.

The guard never signals another workload. Aggregate GPU samples cannot assign
GPU work to a process or prove exclusivity; retain them with every measurement.
"""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import time

import psutil
from run_metaljax_attention_comparison import snapshot, activity_since


def environment(enabled):
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(('METALJAX_', 'MLX_', 'FBC_', 'JAX_'))}
    env.update(JAX_PLATFORMS='metal,cpu', JAX_DISABLE_JIT='0',
               JAX_ENABLE_COMPILATION_CACHE='false', METALJAX_EVO_ATTN='1',
               METALJAX_SLICE_BEFORE_CONVERT=str(int(enabled)),
               METALJAX_MATMUL_PRECISION='high', METALJAX_MEM_BUDGET_MB='16384',
               MLX_MAX_OPS_PER_BUFFER='800', MLX_MAX_MB_PER_BUFFER='512',
               MLX_MAX_MB_TRANSIENT_ONLY='0', METALJAX_LOOP_INFLIGHT_MB='1024',
               METALJAX_CBUF_BYTES='inputs', MPLBACKEND='Agg')
    return env


def shipping_environment():
    """Select Metal without overriding backend optimization or resource defaults."""
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(('METALJAX_', 'MLX_', 'FBC_', 'JAX_'))}
    env.update(JAX_PLATFORMS='metal,cpu', JAX_DISABLE_JIT='0',
               JAX_ENABLE_COMPILATION_CACHE='false', MPLBACKEND='Agg')
    return env


def stop_own(child):
    # The memory monitor starts a separate process group. Resolve descendants
    # from this live child, never from an untrusted/stale pid file.
    processes = psutil.Process(child.pid).children(recursive=True)
    groups = {child.pid}
    for process in processes:
        try:
            if os.getpgid(process.pid) == process.pid:
                groups.add(process.pid)
        except ProcessLookupError:
            pass
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


def run_guarded(command, env, out, max_seconds=1200, quiet_timeout=180):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    (out/'command.json').write_text(json.dumps(dict(command=list(map(str, command)),
        environment={k: v for k, v in env.items()
                     if k.startswith(('JAX_', 'METALJAX_', 'MLX_', 'FBC_'))}), indent=2)+'\n')
    previous = snapshot()
    with (out/'activity.jsonl').open('w') as log:
        deadline = time.monotonic() + quiet_timeout
        quiet_since = time.monotonic()
        while time.monotonic() - quiet_since < 6:
            time.sleep(1)
            current = snapshot()
            if not current['gpu'] or not current['processes']:
                raise RuntimeError('Process/GPU visibility unavailable; guard cannot admit work')
            active, busy, _ = activity_since(previous, current, {os.getpid()})
            log.write(json.dumps(dict(stage='preflight', active=active, **current))+'\n')
            log.flush()
            previous = current
            if busy or any(d.get('Device Utilization %', 100) > 20 for d in current['gpu']):
                quiet_since = time.monotonic()
            if time.monotonic() > deadline:
                raise RuntimeError('No CPU/GPU quiet window; other workloads retained')
        started, sustained, reason = time.monotonic(), {}, None
        with (out/'job.log').open('w') as output:
            child = subprocess.Popen(list(map(str, command)), env=env, stdout=output,
                                     stderr=subprocess.STDOUT, start_new_session=True)
            while child.poll() is None:
                time.sleep(1)
                current = snapshot()
                own = {os.getpid(), *[p.pid for p in psutil.Process().children(recursive=True)]}
                active, busy, dt = activity_since(previous, current, own)
                log.write(json.dumps(dict(stage='running', active=active, **current))+'\n')
                log.flush()
                previous = current
                ids = {p['pid'] for p in busy}
                sustained = {p: d for p, d in sustained.items() if p in ids}
                for process in busy:
                    pid = process['pid']
                    sustained[pid] = sustained.get(pid, 0) + dt
                    if sustained[pid] >= 4:
                        reason = f'Competing compute: {process["name"]} pid {pid}'
                if time.monotonic() - started > max_seconds:
                    reason = 'Our job exceeded its bounded runtime'
                if reason:
                    stop_own(child)
                    break
            child.wait()
        result = dict(returncode=child.returncode, reason=reason,
                      elapsed_seconds=time.monotonic()-started)
        (out/'result.json').write_text(json.dumps(result, indent=2)+'\n')
        if reason or child.returncode:
            raise RuntimeError(f'Guarded job failed: {result}; see {out}/job.log')
        return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--enabled', action='store_true')
    parser.add_argument('--max-seconds', type=int, default=1200)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command:
        parser.error('Provide a command after --')
    run_guarded(command, environment(args.enabled), args.out, args.max_seconds)
