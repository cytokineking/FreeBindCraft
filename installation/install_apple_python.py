#!/usr/bin/env python3
"""Install pinned Apple Python dependencies with an explicit MetalJAX selection.

Precedence: FREEBINDCRAFT_METALJAX_WHEEL (local developer override), then
FREEBINDCRAFT_METALJAX=stock (troubleshooting), otherwise the frozen release.
The published wheel is verified before uv can modify the target environment.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import urllib.parse
import urllib.request


def sha(path):
    with path.open('rb') as file:
        return hashlib.file_digest(file, 'sha256').hexdigest()


def download_verified(asset, directory):
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory/Path(urllib.parse.urlsplit(asset['url']).path).name
    if destination.is_file() and sha(destination) == asset['sha256']:
        return destination
    # Separate temporary files also make concurrent downloads safe.
    with tempfile.NamedTemporaryFile(dir=directory, suffix='.part', delete=False) as file:
        part = Path(file.name)
    try:
        print('Downloading MetalJAX:', asset['url'], flush=True)
        urllib.request.urlretrieve(asset['url'], part)
        if sha(part) != asset['sha256']:
            raise RuntimeError(f"MetalJAX checksum mismatch for {asset['url']}; environment unchanged")
        part.replace(destination)
    finally:
        part.unlink(missing_ok=True)
    return destination


def select_wheel(release, directory):
    mode = os.environ.get('FREEBINDCRAFT_METALJAX', 'frozen')
    if mode not in ('frozen', 'stock'):
        raise ValueError('FREEBINDCRAFT_METALJAX must be frozen or stock')
    override = os.environ.get('FREEBINDCRAFT_METALJAX_WHEEL')
    if override:
        wheel = Path(override).expanduser().absolute()
        if not wheel.is_file() or wheel.suffix != '.whl':
            raise ValueError(f'MetalJAX wheel not found: {wheel}')
        return 'local override', wheel, sha(wheel)
    asset = release['stock_wheel'] if mode == 'stock' else release['wheel']
    return mode, download_verified(asset, directory), asset['sha256']


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python', type=Path, default=root/'.venv-metaljax/bin/python')
    parser.add_argument('--uv', default=shutil.which('uv') or str(root/'.local/bin/uv'))
    parser.add_argument('--downloads', type=Path, default=root/'.local/downloads/metaljax')
    parser.add_argument('--with-pyrosetta', action='store_true')
    args = parser.parse_args()
    python = args.python.absolute()  # preserve a venv interpreter's symlink
    query = 'import sys; assert sys.prefix != sys.base_prefix, "Expected a virtual environment"'
    subprocess.run([str(python), '-c', query], check=True)
    release = json.loads((root/'installation/metaljax/release.json').read_text())
    requirements = (root/'installation/apple-silicon-requirements.txt').read_text().splitlines()
    expected = f"metaljax @ {release['wheel']['url']}#sha256={release['wheel']['sha256']}"
    backend_lines = [line for line in requirements if line.startswith('metaljax')]
    if backend_lines != [expected] or any(f'{key}=={release[key]}' not in requirements for key in ('jax', 'jaxlib')):
        raise RuntimeError('Apple requirements disagree with the frozen MetalJAX release pins')
    mode, wheel, digest = select_wheel(release, args.downloads.absolute())
    selected = f'metaljax @ {wheel.as_uri()}#sha256={digest}'
    requirements = [selected if line == expected else line for line in requirements]
    command = [args.uv, 'pip', 'sync', '--python', str(python)]
    if args.with_pyrosetta:
        command += ['--find-links', 'https://west.rosettacommons.org/pyrosetta/quarterly/release/']
        requirements += [line for line in (root/'installation/pyrosetta-apple-requirements.txt').read_text().splitlines()
                         if not line.startswith('--')]
    with tempfile.TemporaryDirectory(prefix='fbc-apple-requirements-') as temporary:
        lock = Path(temporary)/'requirements.txt'
        lock.write_text('\n'.join(requirements)+'\n')
        print(f'MetalJAX selection: {mode}; SHA-256: {digest}', flush=True)
        subprocess.run([*command, str(lock)], check=True)
    report = json.loads(subprocess.check_output([str(python), '-c',
        "import importlib.metadata as m,json; print(json.dumps({n:m.version(n) for n in ['metaljax','jax','jaxlib']}))"], text=True))
    if any(report[key] != release[key] for key in ('jax', 'jaxlib')):
        raise RuntimeError(f'Unexpected JAX versions after sync: {report}')
    wanted = '0.11.9' if mode == 'stock' else release['version']
    if mode != 'local override' and report['metaljax'] != wanted:
        raise RuntimeError(f'Unexpected MetalJAX version after sync: {report}')
    print(f"Installed MetalJAX {report['metaljax']} ({mode}); JAX {report['jax']}, jaxlib {report['jaxlib']}", flush=True)


if __name__ == '__main__':
    main()
