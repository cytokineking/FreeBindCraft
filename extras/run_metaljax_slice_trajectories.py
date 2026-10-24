"""Run the unchanged, checksum-pinned fixed-work harness in isolated snapshots.

Only paths and the opt-in compiler flag differ by arm. Actual generated design
PDBs are compared separately from the downstream shared fixture outputs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tarfile

from guard_metaljax_job import environment, shipping_environment, run_guarded
from run_metaljax_slice_comparison import identity

BASE = 'cf136d1033790eef7b866bf07e92ed82b5572003'
HARNESS_SHA = 'e73d95f8c09f4e5dec907ac68fd2240e2220242873250e7c07eb93cf3462075a'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, required=True, help='Preserved priority2 handoff directory')
    parser.add_argument('--assets', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--order', default='ABBA')
    parser.add_argument('--resume', action='store_true',
                        help='Verify and retain the completed prefix; incomplete snapshots must be moved aside')
    parser.add_argument('--baseline-python', type=Path, default=root/'.venv-metaljax-baseline/bin/python')
    parser.add_argument('--candidate-python', type=Path, default=root/'.venv-metaljax-slice/bin/python')
    parser.add_argument('--baseline-manifest', type=Path, default=root/'.local/src/metaljax-baseline/build-tools/attention-wheels/build-manifest.json')
    parser.add_argument('--candidate-manifest', type=Path, default=root/'.local/src/metaljax-slice/build-tools/attention-wheels/build-manifest.json')
    parser.add_argument('--shipping-defaults', action='store_true')
    args = parser.parse_args()
    if sorted(args.order) != list('AABB'):
        parser.error('Use a balanced order with two processes per arm')
    evidence, assets, out = args.evidence.resolve(), args.assets.resolve(), args.out.absolute()
    preserved = evidence/'worker-runs/metaljax-scheduling-20261023/repeat.noindex'
    harness = preserved/'benchmark-fixed-budget.py'
    if sha(harness) != HARNESS_SHA:
        raise RuntimeError('Fixed-work harness identity differs')
    fixture = evidence/'worker-runs/metaljax-loop-fix/fixed-fixture'
    saved = json.loads((out/'summary.json').read_text())['results'] if args.resume else []
    if len(saved) > len(args.order):
        raise RuntimeError('Saved series is longer than the selected order')
    out.mkdir(parents=True, exist_ok=args.resume)
    installed = {arm: identity(python.absolute()) for arm, python in
                 [('A', args.baseline_python), ('B', args.candidate_python)]}
    for arm, path in [('A', args.baseline_manifest), ('B', args.candidate_manifest)]:
        manifest = json.loads(path.read_text())
        if installed[arm]['files']['metaljax/lib/libmetal_pjrt_native.dylib'] != manifest['plugin_sha256']:
            raise RuntimeError('Trajectory plugin differs from its build manifest')
        if installed[arm]['files']['metaljax/lib/mlx/lib/libmlx_metaljax.dylib'] != '72dc8d9a4b1e07df2f6bc91383e0567acf3ce292398ae6ce4ed1d19fce0b1374':
            raise RuntimeError('Trajectory MLX runtime differs')
    if ({k: v for k, v in installed['A']['packages'].items() if k != 'metaljax'} !=
        {k: v for k, v in installed['B']['packages'].items() if k != 'metaljax'}):
        raise RuntimeError('Trajectory scientific environments differ')
    native = {relative: sha(assets/relative) for relative in
              ['.local/src/FASPR/FASPR', '.local/src/sc-rs/target/release/sc', '.local/dssp/bin/mkdssp']}
    identities = dict(installed=installed, native_tools=native)
    if args.shipping_defaults:
        identities['shipping_defaults'] = True
    if args.resume:
        if identities != json.loads((out/'identities.json').read_text()):
            raise RuntimeError('Resume environments/native tools differ from the saved series')
        if sha(out/'benchmark-fixed-budget.py') != HARNESS_SHA:
            raise RuntimeError('Resume harness identity differs')
    else:
        (out/'identities.json').write_text(json.dumps(identities, indent=2)+'\n')
        shutil.copyfile(harness, out/'benchmark-fixed-budget.py')
    archive = out/'scientific-source.tar'
    if not args.resume:
        subprocess.run(['git', '-C', str(root), 'archive', BASE, '-o', str(archive)], check=True)
    counts, results, original, original_records = dict(A=0, B=0), [], None, None
    for arm in args.order:
        counts[arm] += 1
        name = arm.lower()+str(counts[arm])
        snapshot = out/name
        if len(results) < len(saved):
            row = saved[len(results)]
            record = json.loads((snapshot/'source-record.json').read_text())
            guard = json.loads((snapshot/'guard/result.json').read_text())
            memory = json.loads((snapshot/'memory/summary.json').read_text())
            records = (snapshot/'actual-design-records.txt').read_text().splitlines()
            generated = snapshot/'runs/fixed-budget/output'
            designs = {str(p.relative_to(generated/'Trajectory')): dict(sha256=sha(p), bytes=p.stat().st_size)
                       for p in (generated/'Trajectory').rglob('fixed_ubiquitin_benchmark.pdb')}
            if (row['name'] != name or row['arm'] != arm or
                guard['returncode'] != 0 or guard['reason'] is not None or
                record['scientific_commit'] != BASE or record['source_tar_sha256'] != sha(archive) or
                record['harness_sha256'] != HARNESS_SHA or record['native_tools_sha256'] != native or
                record['fixture_sha256'] != {p.name: sha(p) for p in fixture.iterdir() if p.is_file()} or
                record['config_sha256'] != {p.name: sha(p) for p in (preserved/'trajectory-a1/runs/benchmark/configs').iterdir()} or
                record['config_sha256'] != {p.name: sha(p) for p in (snapshot/'runs/benchmark/configs').iterdir()} or
                row['report'] != json.loads((generated/'report.json').read_text()) or
                row['memory'] != {k: v for k, v in memory.items() if k.startswith(('peak_', 'maximum_'))} or
                len(records) != 140 or sha(snapshot/'actual-design-records.txt') != row['actual_design_records_sha256'] or
                not designs or designs != row['actual_generated_designs'] or
                not row['actual_design_bytes_equal_first_baseline'] or
                not row['actual_design_records_equal_first_baseline']):
                raise RuntimeError(f'Cannot verify completed prefix process: {name}')
            if original is None:
                original, original_records = designs, records
            if designs != original or records != original_records:
                raise RuntimeError('Completed prefix actual designs differ')
            results.append(row)
            print('Verified completed prefix', name, flush=True)
            continue
        # Never reuse partial process outputs or caches. Move an interrupted
        # snapshot aside before resuming; mkdir otherwise refuses it intact.
        snapshot.mkdir()
        with tarfile.open(archive) as tar:
            tar.extractall(snapshot, filter='data')
        config = snapshot/'runs/benchmark/configs'
        shutil.copytree(preserved/'trajectory-a1/runs/benchmark/configs', config)
        (snapshot/'params').symlink_to(assets/'params', target_is_directory=True)
        for relative in ['.local/dssp', '.local/src/FASPR', '.local/src/sc-rs']:
            dest = snapshot/relative
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.symlink_to(assets/relative, target_is_directory=True)
        python = (args.baseline_python if arm == 'A' else args.candidate_python).absolute()
        env = shipping_environment() if args.shipping_defaults else environment(arm == 'B')
        command = [python, root/'extras/monitor_macos_job.py', '--output-dir', snapshot/'memory',
                   '--interval', '.5', '--max-tree-rss-gib', '24', '--max-hours', '.5', '--',
                   python, '-u', out/'benchmark-fixed-budget.py', '--root', snapshot,
                   '--fixture', fixture, '--mode', 'no-pyrosetta', '--backend', 'metal']
        (snapshot/'source-record.json').write_text(json.dumps(dict(scientific_commit=BASE,
            source_tar_sha256=sha(archive), harness_sha256=sha(harness),
            config_sha256={p.name: sha(p) for p in config.iterdir()},
            native_tools_sha256=native,
            fixture_sha256={p.name: sha(p) for p in fixture.iterdir() if p.is_file()}), indent=2)+'\n')
        run_guarded(command, env, snapshot/'guard', max_seconds=1860)
        generated = snapshot/'runs/fixed-budget/output'
        # binder_hallucination saves this before downstream fixture substitution,
        # retaining it in Trajectory/LowConfidence or Trajectory/Clashing on rejection.
        designs = {str(p.relative_to(generated/'Trajectory')): dict(sha256=sha(p), bytes=p.stat().st_size)
                   for p in (generated/'Trajectory').rglob('fixed_ubiquitin_benchmark.pdb')}
        if not designs:
            raise RuntimeError('No actual generated design was retained')
        if original is None:
            original = designs
        same = designs == original
        log = (snapshot/'memory/job.log').read_text(errors='replace')
        before_fixture = log.split('FIXED BUDGET: using the identical saved trajectory')[0]
        records = [line for line in before_fixture.splitlines() if re.match(r'^\d+ models \[', line)]
        if len(records) != 140:
            raise RuntimeError(f'Expected 140 actual design metric records, got {len(records)}')
        record_file = snapshot/'actual-design-records.txt'
        record_file.write_text('\n'.join(records)+'\n')
        if original_records is None:
            original_records = records
        same_records = records == original_records
        report = json.loads((generated/'report.json').read_text())
        memory = json.loads((snapshot/'memory/summary.json').read_text())
        results.append(dict(name=name, arm=arm, actual_generated_designs=designs,
                            actual_design_bytes_equal_first_baseline=same,
                            actual_design_metric_records=len(records),
                            actual_design_records_sha256=sha(record_file),
                            actual_design_records_equal_first_baseline=same_records,
                            report=report,
                            memory={k: v for k, v in memory.items() if k.startswith(('peak_', 'maximum_'))}))
        (out/'summary.json').write_text(json.dumps(dict(results=results,
            note='Fixed work; zero accepted designs. Downstream uses shared fixtures and is not generated-design throughput.'), indent=2)+'\n')
        if not same or not same_records:
            raise RuntimeError('Actual generated designs differ; inspect before making a performance claim')
        print(name, report['stage_seconds'], 'actual design byte equal', same, flush=True)


if __name__ == '__main__':
    main()
