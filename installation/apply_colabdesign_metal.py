#!/usr/bin/env python3
"""Apply the pinned, checksummed ColabDesign Metal design patch in this venv."""
import argparse
import base64
import csv
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def apply(package, assets):
    manifest = json.loads((assets / 'manifest.json').read_text())
    full_patch = assets / manifest['patch']
    if sha(full_patch) != manifest['patch_sha256']:
        raise RuntimeError('ColabDesign patch checksum mismatch; package unchanged')
    files = manifest['files']
    before = {name: values['before'] for name, values in files.items()}
    after = {name: values['after'] for name, values in files.items()}
    actual = {name: sha(package / name) for name in files}
    if actual == after:
        return dict(status='already applied', **manifest)
    if actual != before:
        raise RuntimeError('ColabDesign source differs from the exact pinned baseline; package unchanged')
    patch = full_patch
    for name, values in files.items():
        if 'payload' in values and sha(assets / values['payload']) != values['after']:
            raise RuntimeError('ColabDesign helper checksum mismatch; package unchanged')
    # Patch and verify a private staging copy before changing the target files.
    with tempfile.TemporaryDirectory(prefix='fbc-colabdesign-patch-') as temporary:
        stage = Path(temporary)
        for name, digest in actual.items():
            if digest is not None:
                target = stage / 'colabdesign' / name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(package / name, target)
        for name, values in files.items():
            if 'payload' in values:
                target = stage / 'colabdesign' / name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(assets / values['payload'], target)
        env = dict(os.environ, GIT_CEILING_DIRECTORIES=str(stage))
        subprocess.run(['git', 'apply', '--no-index', str(patch)], cwd=stage, env=env, check=True)
        if {name: sha(stage / 'colabdesign' / name) for name in files} != after:
            raise RuntimeError('Patched ColabDesign checksum mismatch; package unchanged')
        for name in sorted(files):
            target = package / name
            target.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as file:
                pending = Path(file.name)
                file.write((stage / 'colabdesign' / name).read_bytes())
            try:
                pending.chmod(0o644)
                pending.replace(target)
            finally:
                pending.unlink(missing_ok=True)
    return dict(status='applied', **manifest)


def update_record(package, files):
    """Let normal package uninstall/reinstall account for the added helpers."""
    distribution = importlib.metadata.distribution('colabdesign')
    records = [distribution.locate_file(p) for p in distribution.files
               if p.name == 'RECORD' and p.parent.name.endswith('.dist-info')]
    if len(records) != 1 or package.parent not in records[0].resolve().parents:
        raise RuntimeError('Cannot identify this venv ColabDesign RECORD')
    record = records[0]
    rows = {row[0]: row for row in csv.reader(record.read_text().splitlines())}
    for name in files:
        path = package / name
        digest = base64.urlsafe_b64encode(hashlib.sha256(path.read_bytes()).digest()).decode().rstrip('=')
        relative = 'colabdesign/' + name
        rows[relative] = [relative, 'sha256=' + digest, str(path.stat().st_size)]
    with tempfile.NamedTemporaryFile(mode='w', dir=record.parent, newline='', delete=False) as file:
        pending = Path(file.name)
        csv.writer(file, lineterminator='\n').writerows(rows[name] for name in sorted(rows))
    try:
        pending.chmod(0o644)
        pending.replace(record)
    finally:
        pending.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    if sys.prefix == sys.base_prefix:
        raise RuntimeError('Expected a virtual environment; system package unchanged')
    spec = importlib.util.find_spec('colabdesign')
    if spec is None or spec.origin is None:
        raise RuntimeError('ColabDesign is not installed in this venv')
    package = Path(spec.origin).resolve().parent
    if Path(sys.prefix).resolve() not in package.parents:
        raise RuntimeError('ColabDesign is outside this venv; external package unchanged')
    assets = Path(__file__).resolve().parent / 'colabdesign-metal'
    manifest = json.loads((assets / 'manifest.json').read_text())
    origin = importlib.metadata.distribution('colabdesign').read_text('direct_url.json')
    revision = json.loads(origin or '{}').get('vcs_info', {}).get('commit_id')
    if revision != manifest['source_commit']:
        raise RuntimeError('ColabDesign is not installed from the pinned source revision; package unchanged')
    result = dict(package=str(package), **apply(package, assets))
    update_record(package, result['files'])
    if args.report:
        args.report.write_text(json.dumps(result, indent=2) + '\n')
    print('ColabDesign Metal preparation:', result['status'], result['source_commit'])


if __name__ == '__main__':
    main()
