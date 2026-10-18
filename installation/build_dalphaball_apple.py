#!/usr/bin/env python3
"""Build the pinned native DAlphaBall source with the checkout's Fortran tools."""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
REVISION = "369e400415dd5c585f80dcf7b647147ee520e8ba"
SOURCE_PATH = "source/external/DAlpahBall"
SOURCE = ROOT / ".local/src/DAlphaBall"
TOOLCHAIN = ROOT / ".local/dalphaball-toolchain"


def fetch(url):
    request = urllib.request.Request(url, headers={"User-Agent": "FreeBindCraft-installer"})
    with urllib.request.urlopen(request, timeout=90) as response:
        return response.read()


def prepare_source():
    marker = SOURCE / "SOURCE_REVISION"
    if SOURCE.exists():
        if not marker.is_file() or marker.read_text().strip() != REVISION:
            raise RuntimeError(f"Preserving unexpected DAlphaBall source at {SOURCE}")
        return
    api = f"https://api.github.com/repos/RosettaCommons/rosetta/contents/{SOURCE_PATH}?ref={REVISION}"
    entries = json.loads(fetch(api))
    raw_prefix = f"https://raw.githubusercontent.com/RosettaCommons/rosetta/{REVISION}/{SOURCE_PATH}/"
    with tempfile.TemporaryDirectory(prefix="dalphaball-", dir=SOURCE.parent) as folder:
        staging = Path(folder) / "source"
        staging.mkdir()

        def download(entry):
            name = entry["name"]
            url = entry["download_url"]
            if entry["type"] != "file" or Path(name).name != name or not url.startswith(raw_prefix):
                raise RuntimeError(f"Unexpected source entry: {name}")
            content = fetch(url)
            blob = b"blob " + str(len(content)).encode() + b"\0" + content
            if hashlib.sha1(blob).hexdigest() != entry["sha"]:
                raise RuntimeError(f"Source checksum mismatch: {name}")
            (staging / name).write_bytes(content)

        with ThreadPoolExecutor(max_workers=6) as workers:
            list(workers.map(download, entries))
        subprocess.run(["patch", "-p1", "-i", str(ROOT / "installation/dalphaball-apple.patch")],
                       cwd=staging, check=True)
        (staging / "SOURCE_REVISION").write_text(REVISION + "\n")
        staging.rename(SOURCE)


def main():
    prepare_source()
    sdk = os.environ.get("SDKROOT") or subprocess.check_output(
        ["xcrun", "--show-sdk-path"], text=True).strip()
    clang = subprocess.check_output(["xcrun", "--find", "clang"], text=True).strip()
    fortran = TOOLCHAIN / "bin/gfortran"
    if not fortran.is_file():
        raise RuntimeError(f"Missing Fortran compiler: {fortran}")
    env = dict(os.environ, SDKROOT=sdk, CONDA_BUILD_SYSROOT=sdk,
               MACOSX_DEPLOYMENT_TARGET="14.0",
               PATH=str(TOOLCHAIN / "bin") + os.pathsep + os.environ["PATH"])
    system_flags = f"-isysroot {shlex.quote(sdk)} -mmacosx-version-min=14.0"
    subprocess.run(["make", "clean"], cwd=SOURCE, env=env, check=True)
    subprocess.run([
        "make", "-j4", f"CC={shlex.quote(clang)}",
        f"FC={shlex.quote(str(fortran))}", f"LD={shlex.quote(str(fortran))}",
        f"CFLAGS=-c -O3 -std=gnu89 {system_flags}",
        f"FFLAGS=-c -O3 -x f77-cpp-input -fsecond-underscore -fallow-argument-mismatch {system_flags}",
        f"LDFLAGS=-O3 {system_flags} -Wl,-rpath,{shlex.quote(str(TOOLCHAIN / 'lib'))}",
        f"INCS=-I{shlex.quote(str(TOOLCHAIN / 'include'))}",
        f"LIBS=-L{shlex.quote(str(TOOLCHAIN / 'lib'))}",
    ], cwd=SOURCE, env=env, check=True)
    binary = SOURCE / "DAlphaBall.gcc"
    architecture = subprocess.check_output(["file", str(binary)], text=True).strip()
    if "arm64" not in architecture:
        raise RuntimeError(architecture)
    spheres = "NPOINTS\n2\nCOORDS\n0 0 0 1\n5 0 0 1\nEND\n"
    result = subprocess.run([str(binary), "surf_vol"], input=spheres, text=True,
                            capture_output=True, check=True)
    rows = [list(map(float, row.split())) for row in result.stdout.splitlines() if row.strip()]
    if len(rows) != 2:
        raise RuntimeError(f"Unexpected sphere result: {result.stdout}")
    for row in rows:
        if not (math.isclose(row[1], 4 * math.pi, rel_tol=1e-10) and
                math.isclose(row[2], 4 * math.pi / 3, rel_tol=1e-10)):
            raise RuntimeError(f"Incorrect DAlphaBall sphere result: {row}")
    print(architecture, flush=True)
    print("DAlphaBall reference sphere check passed.", flush=True)


if __name__ == "__main__":
    main()
