"""Sequential, fresh-process scheduling experiments on enabled attention.

Configurations are explicit JSON, e.g. [{"name":"ops200", "env":
{"MLX_MAX_OPS_PER_BUFFER":"200"}}]. The memory governor stays enabled.
Diagnostic traces and loop timing require explicit flags and are excluded from
ordinary performance comparisons. Repeated calls within a process correlate.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time

import numpy as np
import psutil
from run_metaljax_attention_comparison import snapshot, activity_since


def equal_archive(reference, candidate):
    with np.load(reference) as left, np.load(candidate) as right:
        return (set(left.files) == set(right.files) and
                all(left[k].shape == right[k].shape and left[k].dtype == right[k].dtype
                    and left[k].tobytes() == right[k].tobytes() for k in left.files))


def stop_own(child, own, pid_file):
    groups = [child.pid]
    if pid_file.exists():
        native = int(pid_file.read_text())
        if native in own:
            groups.insert(0, native)
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


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python', type=Path, default=root/'.venv-metaljax-scheduling/bin/python')
    parser.add_argument('--configs', type=Path, required=True)
    parser.add_argument('--order', help='Comma-separated configuration names; defaults to file order')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--assets', type=Path, required=True)
    parser.add_argument('--previous', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True, help='Reference output prefix without .npz')
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--pdb', type=Path)
    parser.add_argument('--hotspot', default='44')
    parser.add_argument('--binder-length', type=int, default=70)
    parser.add_argument('--loop-profile', action='store_true')
    parser.add_argument('--trace-library', type=Path)
    parser.add_argument('--trace-configs', help='Comma-separated names to trace; default traces every configuration')
    parser.add_argument('--copy-trace', action='store_true')
    parser.add_argument('--reflection', action='store_true')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--dump-hlo', action='store_true')
    parser.add_argument('--inventory', action='store_true')
    parser.add_argument('--profile-manifest', type=Path, help='Required identity receipt for the loop diagnostic build')
    args = parser.parse_args()
    configs = json.loads(args.configs.read_text())
    names = {c['name']: c for c in configs}
    order = args.order.split(',') if args.order else [c['name'] for c in configs]
    if len(names) != len(configs) or set(order)-names.keys() or args.repeats < 1:
        parser.error('Use unique configuration names, valid order and positive repeats')
    trace_names = set(args.trace_configs.split(',')) if args.trace_configs else set(names)
    if trace_names-names.keys() or ((args.copy_trace or args.reflection or args.trace_configs) and not args.trace_library):
        parser.error('Trace options require a library and valid configuration names')
    out = args.out.absolute()
    if out.exists():
        parser.error('Choose a fresh output directory')
    out.mkdir(parents=True)
    python = args.python.absolute()
    query = """import hashlib,importlib.metadata as m,json
d=m.distribution('metaljax')
print(json.dumps(dict(versions={k:m.version(k) for k in ['metaljax','jax','jaxlib']},
plugin_sha256=hashlib.sha256(d.locate_file('metaljax/lib/libmetal_pjrt_native.dylib').read_bytes()).hexdigest(),
mlx_sha256=hashlib.sha256(d.locate_file('metaljax/lib/mlx/lib/libmlx_metaljax.dylib').read_bytes()).hexdigest())))"""
    identity = json.loads(subprocess.check_output([str(python), '-c', query], text=True))
    expected = '19b28d865a0344768c75af9d25482dbc217954fa1b202e4b4258a7fb79bf9b2a'
    if args.loop_profile:
        if not args.profile_manifest:
            parser.error('--loop-profile requires its --profile-manifest')
        manifest = json.loads(args.profile_manifest.read_text())
        expected = manifest['plugin_sha256']
        if not manifest['diagnostic_profile']:
            parser.error('Selected manifest is not diagnostic')
    if (identity['plugin_sha256'] != expected or
        identity['mlx_sha256'] != '72dc8d9a4b1e07df2f6bc91383e0567acf3ce292398ae6ce4ed1d19fce0b1374' or
        any(identity['versions'][key] != '0.11.2' for key in ['jax', 'jaxlib'])):
        parser.error('Select the reviewed attention-enabled baseline or its checked diagnostic build')
    (out/'identity.json').write_text(json.dumps(identity, indent=2)+'\n')
    env = {k: v for k, v in os.environ.items() if not k.startswith(('METALJAX_', 'MLX_', 'FBC_'))}
    env.update(JAX_PLATFORMS='metal', JAX_DISABLE_JIT='0', JAX_ENABLE_COMPILATION_CACHE='false',
               METALJAX_EVO_ATTN='1', METALJAX_MATMUL_PRECISION='high', METALJAX_MEM_BUDGET_MB='16384',
               MLX_MAX_OPS_PER_BUFFER='800', MLX_MAX_MB_PER_BUFFER='512',
               MLX_MAX_MB_TRANSIENT_ONLY='0', METALJAX_LOOP_INFLIGHT_MB='1024',
               METALJAX_CBUF_BYTES='inputs', MPLBACKEND='Agg')
    allowed = {'MLX_MAX_OPS_PER_BUFFER', 'MLX_MAX_MB_PER_BUFFER', 'METALJAX_LOOP_INFLIGHT_MB',
               'METALJAX_CBUF_BYTES'}
    previous, counts, results = snapshot(), {}, []
    with (out/'contention.jsonl').open('w') as activity:
        for config_name in order:
            config = names[config_name]
            if set(config['env'])-allowed:
                parser.error('Configuration contains an unapproved experimental setting')
            counts[config_name] = counts.get(config_name, 0)+1
            name = config_name+'-'+str(counts[config_name])
            deadline, quiet_since = time.monotonic()+180, time.monotonic()
            while time.monotonic()-quiet_since < 6:
                time.sleep(1)
                current = snapshot()
                active, busy, _ = activity_since(previous, current, {os.getpid()})
                gpu_busy = any(d.get('Device Utilization %', 0) > 20 for d in current['gpu'])
                activity.write(json.dumps(dict(stage=name+'-preflight', active=active, **current))+'\n')
                activity.flush()
                previous = current
                if busy or gpu_busy:
                    quiet_since = time.monotonic()
                if time.monotonic() > deadline:
                    raise RuntimeError('No six-second CPU/GPU quiet window; other workloads retained')
            run_env = dict(env, **{k: str(v) for k, v in config['env'].items()})
            if run_env['METALJAX_CBUF_BYTES'] not in ('inputs', 'transient'):
                parser.error('Command-buffer byte accounting must be inputs or transient')
            run_env['MLX_MAX_MB_TRANSIENT_ONLY'] = '1' if run_env['METALJAX_CBUF_BYTES'] == 'transient' else '0'
            if int(run_env['MLX_MAX_MB_PER_BUFFER']) != 512:
                parser.error('This experiment keeps the byte budget at 512 MiB')
            if int(run_env['MLX_MAX_OPS_PER_BUFFER']) not in (100, 200, 400, 800):
                parser.error('Operation cap is outside this bounded sweep')
            if int(run_env['METALJAX_LOOP_INFLIGHT_MB']) not in (1024, 4096, 8192):
                parser.error('Loop window is outside this bounded sweep')
            if int(run_env['METALJAX_LOOP_INFLIGHT_MB']) >= 8192 and psutil.virtual_memory().available < 32*2**30:
                raise RuntimeError('Deferring 8 GiB window: less than 32 GiB available')
            if args.loop_profile:
                run_env['FBC_SCHED_PROFILE'] = '1'
            tracing = bool(args.trace_library and config_name in trace_names)
            if tracing:
                run_env.update(FBC_MLX_TRACE_LIBRARY=str(args.trace_library.absolute()),
                               FBC_MLX_TRACE_OUTPUT=str(out/(name+'-trace.jsonl')))
                if args.copy_trace:
                    run_env['FBC_MLX_TRACE_COPIES'] = '1'
                if args.reflection:
                    run_env['FBC_MLX_TRACE_REFLECTION'] = '1'
            if args.debug:
                run_env['METALJAX_DEBUG'] = '1'
            if args.inventory:
                run_env['FBC_MLX_MODEL_INVENTORY'] = str(out/(name+'-weights.json'))
            command = [str(python), str(root/'extras/monitor_macos_job.py'),
                       '--output-dir', str(out/(name+'-monitor')), '--interval', '1',
                       '--max-tree-rss-gib', '24', '--max-hours', '0.12', '--',
                       str(python), str(root/'extras/benchmark_metaljax_scheduling.py'),
                       '--assets', str(args.assets.absolute()), '--previous', str(args.previous.absolute()),
                       '--out', str(out), '--name', name, '--repeats', str(args.repeats),
                       '--hotspot', args.hotspot, '--binder-length', str(args.binder_length)]
            if args.pdb:
                command += ['--pdb', str(args.pdb.absolute())]
            if args.loop_profile or tracing:
                command.append('--profile')
            if args.dump_hlo:
                command.append('--dump-hlo')
            sustained, reason = {}, None
            with (out/(name+'-launch.log')).open('w') as log:
                child = subprocess.Popen(command, env=run_env, stdout=log,
                                         stderr=subprocess.STDOUT, start_new_session=True)
                while child.poll() is None:
                    time.sleep(1)
                    current = snapshot()
                    own = {os.getpid(), *[p.pid for p in psutil.Process().children(recursive=True)]}
                    active, busy, dt = activity_since(previous, current, own)
                    busy_ids = {p['pid'] for p in busy}
                    sustained = {pid: d for pid, d in sustained.items() if pid in busy_ids}
                    for p in busy:
                        sustained[p['pid']] = sustained.get(p['pid'], 0)+dt
                        if sustained[p['pid']] >= 4:
                            reason = f"Concurrent compute: {p['name']} pid {p['pid']}"
                    activity.write(json.dumps(dict(stage=name, active=active, **current))+'\n')
                    activity.flush()
                    previous = current
                    if reason:
                        stop_own(child, own, out/(name+'-monitor/pid'))
                        break
                child.wait()
            if reason or child.returncode:
                (out/'failure.json').write_text(json.dumps(dict(name=name, reason=reason, returncode=child.returncode), indent=2)+'\n')
                raise RuntimeError(reason or f'Probe failed: {name}')
            numeric = {suffix: equal_archive(Path(str(args.reference)+suffix+'.npz'), out/(name+suffix+'.npz'))
                       for suffix in ['', '-warm-forward', '-previous']}
            if not all(numeric.values()):
                raise RuntimeError(f'Byte comparison failed: {name}: {numeric}')
            report = json.loads((out/(name+'.json')).read_text())
            if report['attention'] != '1' or report['precision'] != 'high':
                raise RuntimeError('Required numerical baseline not active')
            if report['length'] in (146, 185):
                expected_dispatches = 138441 + (report['length'] == 185)
                if any(abs(v['dispatches']-expected_dispatches) > 20 for v in report['gradient_dispatch']):
                    raise RuntimeError('Unexpected dispatch count: investigate activation/fallback before comparing timing')
            results.append(dict(name=name, config=config, environment={k: run_env[k] for k in env if k.startswith(('MLX_', 'METALJAX_'))},
                                diagnostic_flags=dict(trace=tracing, copies=tracing and args.copy_trace,
                                                      reflection=tracing and args.reflection, loop_profile=args.loop_profile,
                                                      debug=args.debug, inventory=args.inventory, dump_hlo=args.dump_hlo),
                                command=command, numerical_byte_equal=numeric, report=report))
            summary = dict(identity=identity, configs=configs, order=order, repeats_per_process=args.repeats, results=results,
                           diagnostics=bool(args.loop_profile or args.trace_library or args.debug or args.inventory or args.dump_hlo),
                           configurations_sha256=hashlib.sha256(args.configs.read_bytes()).hexdigest(),
                           sampling_unit='Fresh process; repeated calls within a process are correlated.',
                           contention_limit='Continuous CPU/process and aggregate GPU samples cannot prove GPU exclusivity or controlled thermal state.')
            (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
            print(name, 'grad', report['median_grad_seconds'], 'forward', report['median_forward_seconds'],
                  'buffers', [v['command_buffers'] for v in report['gradient_dispatch']], 'exact', numeric, flush=True)


if __name__ == '__main__':
    main()
