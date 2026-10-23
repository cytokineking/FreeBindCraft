#!/usr/bin/env python3
"""Build optional scheduling diagnostics on the reviewed attention-enabled source.

Uses release.json and its archived SDK/XLA/Bazel inputs. Does not change default
dependency pins or publish a release. --install targets only the selected venv.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python', type=Path, default=root/'.venv-metaljax-scheduling-profile/bin/python')
    parser.add_argument('--source-dir', type=Path, default=root/'.local/src/metaljax-scheduling-profile')
    parser.add_argument('--cache-dir', type=Path, default=root/'.local/scheduling-build')
    parser.add_argument('--jobs', type=int, default=8)
    parser.add_argument('--shared-disk-cache', type=Path, default=root/'.local/bazel-disk-cache')
    parser.add_argument('--install', action='store_true')
    args = parser.parse_args()
    if (platform.system(), platform.machine()) != ('Darwin', 'arm64') or args.jobs < 1:
        parser.error('Requires native arm64 macOS and positive --jobs')
    python = args.python.absolute()  # retain the venv interpreter symlink
    source, cache = args.source_dir.resolve(), args.cache_dir.resolve()
    disk_cache = args.shared_disk_cache.resolve()
    for path in (source, cache, disk_cache, python.parent.parent.resolve()):
        if not path.is_relative_to(root) or path == root:
            parser.error('Source, cache and selected venv must belong to this checkout')
    spec = importlib.util.spec_from_file_location('frozen_builder', root/'installation/build_metaljax_apple.py')
    base = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)
    release = json.loads((root/'installation/metaljax/release.json').read_text())
    query = ("import importlib.metadata as m,json,sys; print(json.dumps(dict("
             "venv=sys.prefix!=sys.base_prefix,jax=m.version('jax'),jaxlib=m.version('jaxlib'))))")
    target = json.loads(subprocess.check_output([str(python), '-c', query], text=True))
    if not target['venv'] or any(target[k] != release[k] for k in ('jax', 'jaxlib')):
        parser.error('Select an isolated venv with JAX/jaxlib 0.11.2')
    if not source.exists():
        base.checkout(release['repository']+'.git', source, release['source_commit'])
    current = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
    if current != release['source_commit']:
        raise RuntimeError(f'Preserving unexpected source HEAD: {current}')
    patches = [root/'installation/metaljax/attention/0001-bf16-attention-epilogues.patch']
    patches.append(root/'installation/metaljax/scheduling/diagnostic-loop-timing.patch')
    patch_hashes = {p.name: base.sha(p) for p in patches}
    targets = sorted({name for p in patches for name in
                      re.findall(r'^\+\+\+ b/(.+)$', p.read_text(), re.M)})
    stamp = source/'.freebindcraft-scheduling-source.json'

    def fingerprint():
        return {name: base.sha(source/name) for name in targets}

    def changed_files():
        return set(subprocess.check_output(
            ['git', '-C', str(source), 'diff', '--name-only', 'HEAD'], text=True).splitlines())

    build_files = {'plugin-native/.bazelrc', 'plugin-native/WORKSPACE'}
    if stamp.exists():
        saved = json.loads(stamp.read_text())
        if saved != dict(patches=patch_hashes, files=fingerprint()):
            raise RuntimeError('Preserving changed patch/source; choose a fresh --source-dir')
        if changed_files() - set(targets) - build_files:
            raise RuntimeError('Preserving unrelated source changes')
    else:
        if changed_files() - build_files:
            raise RuntimeError('Preserving modified source; choose a fresh --source-dir')
        for patch in patches:
            base.run(['git', '-C', source, 'apply', '--check', patch])
            base.run(['git', '-C', source, 'apply', patch])
        for name in ('pyproject.toml', 'src/metaljax/__init__.py'):
            file = source/name
            file.write_text(file.read_text().replace('0.11.9+fbc3attndev2', '0.11.9+fbc4schedprofile1'))
        stamp.write_text(json.dumps(dict(patches=patch_hashes, files=fingerprint()), indent=2)+'\n')
    if (source/'plugin-native/.bazelrc.user').exists():
        raise RuntimeError('Preserving custom .bazelrc.user; choose a fresh --source-dir')

    inputs = cache/'attention-inputs'
    sdk_archive, xla_archive, bazel = inputs/'mlx-sdk.tar.gz', inputs/'xla.tar.gz', inputs/'bazel-8.7.0'
    base.download(release['mlx_sdk'], sdk_archive)
    base.download(release['xla'], xla_archive)
    base.download(release['bazel'], bazel)
    sdk = source/'src/metaljax/lib/mlx'
    base.prepare_archive(sdk_archive, sdk, 'mlx')
    base.prepare_archive(xla_archive, source/'xla-91888df6', 'xla-'+release['xla']['revision'])
    bazel.chmod(0o755)
    for name in build_files:
        file = source/name
        original = subprocess.check_output(['git', '-C', str(source), 'show', 'HEAD:'+name], text=True)
        prepared = original.replace('/Users/oleg/metaljax', str(source)).replace(
            'build --disk_cache=~/.cache/metaljax-bazel', f'build --disk_cache={disk_cache}')
        if file.read_text() not in (original, prepared):
            raise RuntimeError(f'Preserving changed build configuration: {file}')
        file.write_text(prepared)
    env = dict(os.environ, METALJAX_MLX_DIR=str(sdk), UV_CACHE_DIR=str(cache/'uv-cache'),
               SDKROOT=subprocess.check_output(['xcrun', '--show-sdk-path'], text=True).strip())
    flags = ['--repositories_without_autoloads=cc_compatibility_proxy', '--macos_minimum_os=14.0',
             '--host_copt=--target=arm64-apple-macosx14.0', '--host_linkopt=--target=arm64-apple-macosx14.0',
             f'--jobs={args.jobs}', f'--disk_cache={disk_cache}']
    base.run([bazel, f'--output_user_root={cache}/bazel', 'build', *flags,
              '//metal:libmetal_pjrt_native.dylib'], cwd=source/'plugin-native', env=env)
    plugin = source/'plugin-native/bazel-bin/metal/libmetal_pjrt_native.dylib'
    output = source/'build-tools/attention-wheels'
    output.mkdir(parents=True, exist_ok=True)
    constraints = output/'build-constraints.txt'
    constraints.write_text(f"hatchling=={release['hatchling']}\n")
    uv = shutil.which('uv') or str(root/'.local/bin/uv')
    env.update(METALJAX_NATIVE_DYLIB=str(plugin), METALJAX_VENDORED_MLX=str(sdk))
    base.run([uv, 'build', '--wheel', '--build-constraint', constraints, '--out-dir', output], cwd=source, env=env)
    version = '0.11.9+fbc4schedprofile1'
    wheel, = output.glob(f'metaljax-{version}-*.whl')
    manifest = dict(builder_sha256=base.sha(Path(__file__)), version=version, frozen_source_commit=current, frozen_release=release,
                    patches=patch_hashes, modified_source=fingerprint(), diagnostic_profile=True,
                    attention_enabled_by_default=False, build_flags=flags,
                    compiler=subprocess.check_output(['xcrun', 'clang++', '--version'], text=True).strip(),
                    python=subprocess.check_output([str(python), '--version'], text=True).strip(),
                    uv=subprocess.check_output([uv, '--version'], text=True).strip(),
                    macos=platform.mac_ver()[0], sdkroot=env['SDKROOT'], jax=target['jax'], jaxlib=target['jaxlib'],
                    plugin_sha256=base.sha(plugin), wheel_sha256=base.sha(wheel),
                    mlx_runtime={p.name: base.sha(p) for p in sorted((sdk/'lib').iterdir()) if p.is_file()})
    (output/'build-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    if args.install:
        base.run([uv, 'pip', 'install', '--python', python, '--no-deps', '--reinstall', wheel], env=env)
    print(f'Built optional scheduling diagnostics: {wheel}\nManifest: {output}/build-manifest.json', flush=True)


if __name__ == '__main__':
    main()
