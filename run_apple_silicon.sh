#!/usr/bin/env bash
# GPU-only native launcher. All paths are relative to this installation.
set -euo pipefail
cd "$(dirname "$0")"
repo_dir="$PWD"
export JAX_PLATFORMS=metal
export JAX_DISABLE_JIT=0
export OPENMM_PLATFORM_ORDER="${OPENMM_PLATFORM_ORDER:-OpenCL}"
export OPENMM_PLATFORM_ORDER="${OPENMM_PLATFORM_ORDER//[[:space:]]/}"
IFS=',' read -r -a openmm_platforms <<< "$OPENMM_PLATFORM_ORDER"
for platform in "${openmm_platforms[@]}"; do
  case "$platform" in
    OpenCL|Metal) ;;
    *) echo 'This launcher permits only the OpenCL and Metal GPU platforms.' >&2; exit 1 ;;
  esac
done
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$repo_dir/.local/cache}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$repo_dir/.local/matplotlib}"
export FASPR_BIN="$repo_dir/.local/src/FASPR/FASPR"
export SC_RS_BIN="$repo_dir/.local/src/sc-rs/target/release/sc"
export DSSP_BIN="$repo_dir/.local/dssp/bin/mkdssp"
export DALPHABALL_BIN="$repo_dir/.local/src/DAlphaBall/DAlphaBall.gcc"
export LIBCIFPP_DATA_DIR="$repo_dir/.local/dssp/share/libcifpp"
export PATH="$repo_dir/.local/dssp/bin:$PATH"
python="$repo_dir/.venv-metaljax/bin/python"
if [[ ! -x "$python" ]]; then
  echo 'Run ./install_bindcraft.sh before starting an Apple Silicon job.' >&2
  exit 1
fi
if "$python" -c 'import importlib.util, sys; sys.exit(importlib.util.find_spec("pyrosetta") is None)'; then
  exec "$python" -u bindcraft.py "$@"
else
  exec "$python" -u bindcraft.py --no-pyrosetta "$@"
fi
