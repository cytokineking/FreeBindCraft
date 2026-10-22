#!/usr/bin/env python3
"""Build the Metal-only nested-loop fix against MetalJAX 0.11.9.

Requires an installed 0.11.9 Apple environment. Reuses the release's unmodified
private MLX binaries, and builds the PJRT plugin from pinned source. Downloads,
source checkouts and wheels stay beneath --source-dir; --cache-dir holds Bazel's caches. No global
packages or another checkout's environment are changed.
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

METALJAX_REV = '7546eb93aa60c6570efdca95836276d03fe70a94'
MLX_HEADERS_REV = '651c39cd3d26d0a62587ca639a52ea1efa4b8186'
XLA_REV = '91888df6ce85102c30220e41d952065925e10886'
XLA_SHA256 = '32cb8087d88318a01f9974ead6254e1ca06d263fa1847f61503379525bd46897'
BAZEL_SHA256 = '575f20fb23955e02f73519befd180df635b4ed0960c60f0e70fcc8d74014a713'
MLX_BINARIES = {
    'lib/libmlx_metaljax.dylib': '72dc8d9a4b1e07df2f6bc91383e0567acf3ce292398ae6ce4ed1d19fce0b1374',
    'lib/libjaccl_metaljax.dylib': '48f5790b66a173170527e6ba7fac3e445787550851439870019581032d5448f9',
    'lib/mlx.metallib': 'ef1aa1ef9ef5b18bcef166b75b7b8865afe0f02f5e23980f2bdf76447fe60e6b',
}
HEADER_PATCHES = ['0004-diag-dispatch-stats.patch', '0005-donate-through-stream-pins.patch',
                  '0006-skip-empty-finalize.patch', '0007-transient-byte-cadence.patch']


def run(command, **kwargs):
    print('+', ' '.join(map(str, command)), flush=True)
    return subprocess.run(list(map(str, command)), check=True, **kwargs)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def download(url, path, expected):
    if path.is_file() and sha(path) == expected:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_suffix(path.suffix + '.part')
    print('Downloading', url, flush=True)
    urllib.request.urlretrieve(url, part)
    if sha(part) != expected:
        raise RuntimeError(f'Checksum mismatch: {url}')
    part.replace(path)


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


def apply_patch(source, patch):
    # Idempotent only for the exact patch; never reset a developer's changes.
    if subprocess.run(['git', '-C', str(source), 'apply', '--reverse', '--check', str(patch)],
                      stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0:
        return
    run(['git', '-C', source, 'apply', '--check', patch])
    run(['git', '-C', source, 'apply', patch])


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python', type=Path, default=root/'.venv-metaljax/bin/python')
    parser.add_argument('--source-dir', type=Path, default=root/'.local/src/metaljax-loop-fix')
    parser.add_argument('--cache-dir', type=Path, default=root/'.local/metaljax-build')
    parser.add_argument('--jobs', type=int, default=12)
    parser.add_argument('--install', action='store_true', help='Install the built wheel into --python only')
    args = parser.parse_args()
    if (platform.system(), platform.machine()) != ('Darwin', 'arm64'):
        parser.error('Requires native arm64 macOS')
    python = args.python.absolute()
    source = args.source_dir.absolute()
    cache = args.cache_dir.absolute()
    # Do not resolve a venv interpreter symlink: doing so selects its global
    # interpreter and would install into the wrong environment.
    if not python.is_file() or args.jobs < 1:
        parser.error('Need an existing Apple venv and a positive --jobs')
    uv = shutil.which('uv') or str(root/'.local/bin/uv')
    query = ("import importlib.util, importlib.metadata, json, sys; from pathlib import Path; "
             "print(json.dumps(dict(venv=sys.prefix != sys.base_prefix, jax=importlib.metadata.version('jax'), "
             "mlx=str(Path(importlib.util.find_spec('metaljax').submodule_search_locations[0])/'lib/mlx'))))")
    target = json.loads(subprocess.check_output([str(python), '-c', query], text=True))
    if not target['venv'] or target['jax'] != '0.11.2':
        parser.error('--python must select the pinned Apple venv with JAX 0.11.2')
    installed = Path(target['mlx'])
    for rel, expected in MLX_BINARIES.items():
        if sha(installed/rel) != expected:
            raise RuntimeError(f'Expected the unmodified MetalJAX 0.11.9 private runtime: {installed/rel}')
    checkout('https://github.com/eterevsky/metaljax.git', source, METALJAX_REV)
    patch = root/'installation/metaljax/0001-static-nested-loop-trip.patch'
    apply_patch(source, patch)

    # The release stamp's later MLX commit is absent from the public fork.
    # Its public base plus the release's checked-in patches reproduce the SDK
    # headers used by the plugin. The runtime itself stays byte-for-byte the
    # release wheel's, verified above; no MLX/kernel/precision changes.
    headers = source/'mlx-src'
    checkout('https://github.com/eterevsky/mlx.git', headers, MLX_HEADERS_REV)
    header_patches = {name: sha(source/'notes/patches'/name) for name in HEADER_PATCHES}
    header_stamp = headers/'.freebindcraft-sdk-patches.json'
    def header_diff():
        return hashlib.sha256(subprocess.check_output(
            ['git', '-C', str(headers), 'diff', '--no-ext-diff', '--no-textconv',
             '--no-color', '--full-index', '--binary', 'HEAD'])).hexdigest()
    if header_stamp.is_file():
        stamp = json.loads(header_stamp.read_text())
        if stamp != dict(patches=header_patches, diff_sha256=header_diff()):
            raise RuntimeError(f'Preserving modified SDK sources at {headers}')
    else:
        run(['git', '-C', headers, 'diff', '--exit-code', 'HEAD'])
        for name in HEADER_PATCHES:
            run(['git', '-C', headers, 'apply', source/'notes/patches'/name])
        header_stamp.write_text(json.dumps(dict(patches=header_patches, diff_sha256=header_diff()))+'\n')
    sdk = source/'src/metaljax/lib/mlx'
    sdk.mkdir(parents=True, exist_ok=True)
    shutil.copytree(installed, sdk, dirs_exist_ok=True)
    shutil.copytree(headers/'mlx', sdk/'include/mlx', dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns('*.cpp', '*.cu', '*.metal', '*.txt', 'CMakeLists.txt'))

    tools = source/'build-tools'
    tools.mkdir(exist_ok=True)
    archive = tools/'xla.tar.gz'
    download(f'https://codeload.github.com/openxla/xla/tar.gz/{XLA_REV}', archive, XLA_SHA256)
    xla = source/'xla-91888df6'
    marker = xla/'.freebindcraft-source.json'
    if not xla.exists():
        with tarfile.open(archive) as tar:
            tar.extractall(source, filter='data')
        (source/f'xla-{XLA_REV}').rename(xla)
        marker.write_text(json.dumps({'revision': XLA_REV, 'archive_sha256': XLA_SHA256})+'\n')
    if not marker.is_file() or json.loads(marker.read_text())['archive_sha256'] != XLA_SHA256:
        raise RuntimeError(f'Preserving unverified XLA source at {xla}; use a fresh --source-dir')
    bazel = tools/'bazel-8.7.0'
    download('https://github.com/bazelbuild/bazel/releases/download/8.7.0/bazel-8.7.0-darwin-arm64',
             bazel, BAZEL_SHA256)
    bazel.chmod(0o755)
    for rel in ['plugin-native/WORKSPACE', 'plugin-native/.bazelrc']:
        file = source/rel
        text = file.read_text().replace('/Users/oleg/metaljax', str(source))
        text = text.replace('build --disk_cache=~/.cache/metaljax-bazel',
                            f'build --disk_cache={cache}/bazel-disk-cache')
        file.write_text(text)
    env = dict(os.environ, METALJAX_MLX_DIR=str(sdk), UV_CACHE_DIR=str(tools/'uv-cache'),
               SDKROOT=subprocess.check_output(['xcrun', '--show-sdk-path'], text=True).strip())
    flags = ['--repositories_without_autoloads=cc_compatibility_proxy', '--macos_minimum_os=14.0',
             '--host_copt=--target=arm64-apple-macosx14.0', '--host_linkopt=--target=arm64-apple-macosx14.0',
             f'--jobs={args.jobs}', f'--disk_cache={cache}/bazel-disk-cache']
    run([bazel, f'--output_user_root={cache}/bazel', 'build', *flags,
         '//metal:libmetal_pjrt_native.dylib'], cwd=source/'plugin-native', env=env)
    dylib = source/'plugin-native/bazel-bin/metal/libmetal_pjrt_native.dylib'
    version = '0.11.9+fbc1'
    for rel in ['pyproject.toml', 'src/metaljax/__init__.py']:
        p = source/rel
        p.write_text(p.read_text().replace('"0.11.9"', f'"{version}"'))
    out = tools/'wheels'
    out.mkdir(exist_ok=True)
    env.update(METALJAX_NATIVE_DYLIB=str(dylib), METALJAX_VENDORED_MLX=str(sdk))
    run([uv, 'build', '--wheel', '--out-dir', out], cwd=source, env=env)
    wheel, = out.glob(f'metaljax-{version}-*.whl')
    manifest = dict(metaljax_revision=METALJAX_REV, mlx_header_revision=MLX_HEADERS_REV,
                    mlx_header_patches=header_patches, mlx_header_diff_sha256=header_diff(),
                    mlx_binaries=MLX_BINARIES, xla_revision=XLA_REV, xla_archive_sha256=XLA_SHA256,
                    bazel_sha256=BAZEL_SHA256, patch_sha256=sha(patch), build_flags=flags,
                    compiler_version=subprocess.check_output(['xcrun', 'clang++', '--version'], text=True).strip(),
                    python_version=subprocess.check_output([str(python), '--version'], text=True).strip(),
                    jax_version=target['jax'], macos=platform.mac_ver()[0],
                    sdkroot=env['SDKROOT'], plugin_sha256=sha(dylib), wheel_sha256=sha(wheel),
                    version=version, wheel=str(wheel))
    (out/'build-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    if args.install:
        run([uv, 'pip', 'install', '--python', python, '--no-deps', '--reinstall', wheel], env=env)
    print(f'Built {wheel}\nManifest: {out}/build-manifest.json', flush=True)


if __name__ == '__main__':
    main()
