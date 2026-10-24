#!/usr/bin/env python3
"""Maintainer build of the frozen FreeBindCraft MetalJAX release.

Normal users receive a prebuilt wheel from install_bindcraft.sh. This recipe
checks out the exact fork commit (the fix is already committed), verifies the
archived SDK/XLA/Bazel inputs, and builds only inside the selected source/cache
directories. Requires native arm64 macOS, Command Line Tools and the pinned venv.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tarfile
import urllib.request


def run(command, **kwargs):
    print('+', ' '.join(map(str, command)), flush=True)
    return subprocess.run(list(map(str, command)), check=True, **kwargs)


def sha(path):
    with Path(path).open('rb') as file:
        return hashlib.file_digest(file, 'sha256').hexdigest()


def download(asset, path):
    if path.is_file() and sha(path) == asset['sha256']:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_suffix(path.suffix + '.part')
    try:
        print('Downloading', asset['url'], flush=True)
        urllib.request.urlretrieve(asset['url'], part)
        if sha(part) != asset['sha256']:
            raise RuntimeError(f"Checksum mismatch: {asset['url']}")
        part.replace(path)
    finally:
        part.unlink(missing_ok=True)


def checkout(url, path, revision):
    if not path.exists():
        path.mkdir(parents=True)
        run(['git', 'init', path])
        run(['git', '-C', path, 'remote', 'add', 'origin', url])
        run(['git', '-C', path, 'fetch', '--depth', '1', 'origin', revision])
        run(['git', '-C', path, 'checkout', '--detach', 'FETCH_HEAD'])
    current = subprocess.check_output(['git', '-C', str(path), 'rev-parse', 'HEAD'], text=True).strip()
    if current != revision:
        raise RuntimeError(f'Preserving {path}: expected {revision}, found {current}')
    changed = subprocess.check_output(
        ['git', '-C', str(path), 'diff', '--name-only', 'HEAD'], text=True).splitlines()
    if set(changed) - {'plugin-native/.bazelrc', 'plugin-native/WORKSPACE'}:
        raise RuntimeError(f'Preserving modified source at {path}: {changed}')
    if (path/'plugin-native/.bazelrc.user').exists():
        raise RuntimeError('Preserving custom .bazelrc.user; use a fresh --source-dir')


def prepare_archive(archive, destination, prefix):
    """Extract once; on reruns verify every archived file, preserving changes."""
    with tarfile.open(archive) as tar:
        if not destination.exists():
            destination.parent.mkdir(parents=True, exist_ok=True)
            tar.extractall(destination.parent, filter='data')
            extracted = destination.parent/prefix
            if extracted != destination:
                extracted.rename(destination)
        for member in tar.getmembers():
            if not member.isfile():
                continue
            relative = Path(member.name).relative_to(prefix)
            file = destination/relative
            with tar.extractfile(member) as original:
                expected = hashlib.file_digest(original, 'sha256').hexdigest()
            if not file.is_file() or sha(file) != expected:
                raise RuntimeError(f'Preserving changed/missing build input: {file}')


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release-manifest', type=Path, default=root/'installation/metaljax/release.json')
    parser.add_argument('--python', type=Path, default=root/'.venv-metaljax/bin/python')
    parser.add_argument('--source-dir', type=Path)
    parser.add_argument('--cache-dir', type=Path, default=root/'.local/metaljax-build')
    parser.add_argument('--jobs', type=int, default=12)
    parser.add_argument('--install', action='store_true', help='Install into --python only')
    args = parser.parse_args()
    release = json.loads(args.release_manifest.read_text())
    if args.source_dir is None:
        args.source_dir = root/'.local/src'/('metaljax-'+release['version'].split('+')[-1])
    if (platform.system(), platform.machine()) != ('Darwin', 'arm64'):
        parser.error('Requires native arm64 macOS')
    # Do not resolve the interpreter symlink: it selects the global interpreter.
    python = args.python.absolute()
    source, cache = args.source_dir.absolute(), args.cache_dir.absolute()
    if not python.is_file() or args.jobs < 1:
        parser.error('Need an existing Apple venv and a positive --jobs')
    query = ("import importlib.metadata as m, json, sys; "
             "print(json.dumps(dict(venv=sys.prefix != sys.base_prefix, "
             "jax=m.version('jax'), jaxlib=m.version('jaxlib'))))")
    target = json.loads(subprocess.check_output([str(python), '-c', query], text=True))
    if not target['venv'] or any(target[k] != release[k] for k in ['jax', 'jaxlib']):
        parser.error('--python must select the pinned venv with JAX/jaxlib 0.11.2')
    checkout(release['repository']+'.git', source, release['source_commit'])
    tools = source/'build-tools'
    tools.mkdir(exist_ok=True)
    sdk_archive, xla_archive, bazel = tools/'mlx-sdk.tar.gz', tools/'xla.tar.gz', tools/'bazel-8.7.0'
    download(release['mlx_sdk'], sdk_archive)
    download(release['xla'], xla_archive)
    download(release['bazel'], bazel)
    sdk = source/'src/metaljax/lib/mlx'
    prepare_archive(sdk_archive, sdk, 'mlx')
    xla = source/'xla-91888df6'
    prepare_archive(xla_archive, xla, 'xla-'+release['xla']['revision'])
    bazel.chmod(0o755)
    for rel in ['plugin-native/WORKSPACE', 'plugin-native/.bazelrc']:
        file = source/rel
        original = subprocess.check_output(['git', '-C', str(source), 'show', 'HEAD:'+rel], text=True)
        prepared = original.replace('/Users/oleg/metaljax', str(source)).replace(
            'build --disk_cache=~/.cache/metaljax-bazel', f'build --disk_cache={cache}/bazel-disk-cache')
        if file.read_text() not in (original, prepared):
            raise RuntimeError(f'Preserving custom build configuration: {file}')
        file.write_text(prepared)
    env = dict(os.environ, METALJAX_MLX_DIR=str(sdk), UV_CACHE_DIR=str(tools/'uv-cache'),
               SDKROOT=subprocess.check_output(['xcrun', '--show-sdk-path'], text=True).strip())
    flags = ['--repositories_without_autoloads=cc_compatibility_proxy', '--macos_minimum_os=14.0',
             '--host_copt=--target=arm64-apple-macosx14.0', '--host_linkopt=--target=arm64-apple-macosx14.0',
             f'--jobs={args.jobs}', f'--disk_cache={cache}/bazel-disk-cache']
    run([bazel, f'--output_user_root={cache}/bazel', 'build', *flags,
         '//metal:libmetal_pjrt_native.dylib'], cwd=source/'plugin-native', env=env)
    dylib = source/'plugin-native/bazel-bin/metal/libmetal_pjrt_native.dylib'
    out = tools/'frozen-wheels'
    out.mkdir(exist_ok=True)
    constraints = tools/'build-constraints.txt'
    constraints.write_text(f"hatchling=={release['hatchling']}\n")
    uv = shutil.which('uv') or str(root/'.local/bin/uv')
    env.update(METALJAX_NATIVE_DYLIB=str(dylib), METALJAX_VENDORED_MLX=str(sdk))
    run([uv, 'build', '--wheel', '--build-constraint', constraints, '--out-dir', out], cwd=source, env=env)
    wheel, = out.glob(f"metaljax-{release['version']}-*.whl")
    manifest = dict(source_commit=release['source_commit'], release=release, build_flags=flags,
                    builder_sha256=sha(Path(__file__)),
                    publication_status=release.get('publication_status', 'published'),
                    sdkroot=env['SDKROOT'],
                    uv=subprocess.check_output([uv, '--version'], text=True).strip(),
                    compiler=subprocess.check_output(['xcrun', 'clang++', '--version'], text=True).strip(),
                    python=subprocess.check_output([str(python), '--version'], text=True).strip(),
                    jax=target['jax'], jaxlib=target['jaxlib'], macos=platform.mac_ver()[0],
                    plugin_sha256=sha(dylib), wheel_sha256=sha(wheel),
                    runtime_sha256={p.name: sha(p) for p in sorted((sdk/'lib').iterdir()) if p.is_file()},
                    matches_manifest_plugin=sha(dylib) == release['plugin_sha256'],
                    matches_manifest_wheel=sha(wheel) == release['wheel']['sha256'])
    (out/'build-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    if args.install:
        run([uv, 'pip', 'install', '--python', python, '--no-deps', '--reinstall', wheel], env=env)
    print(f'Built {wheel}\nManifest: {out}/build-manifest.json', flush=True)
    print('Byte-identical to manifest wheel:', manifest['matches_manifest_wheel'], flush=True)


if __name__ == '__main__':
    main()
