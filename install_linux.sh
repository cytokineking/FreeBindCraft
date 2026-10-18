#!/usr/bin/env bash
# Conda supplies Python and native tools; pip supplies the tested CUDA runtime.
set -eo pipefail
cd "$(dirname "$0")"
repo_dir="$PWD"
pkg_manager=conda
env_name=BindCraft
install_pyrosetta=true
skip_gpu_check=false

usage() {
  cat <<'EOF'
Usage: bash install_bindcraft.sh [options]
  --no-pyrosetta, -n       Omit the optional PyRosetta package
  --pkg_manager, -p NAME   conda or mamba (default: conda)
  --env-name NAME          Conda environment name (default: BindCraft)
  --skip-gpu-check        Install without probing a GPU, for container builds
  --cuda, -c VERSION      Legacy option; pip supplies its pinned CUDA 12 runtime
  --fix-channels, -f      Legacy option; installation uses conda-forge explicitly
  --help, -h              Show this help
EOF
}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --pkg_manager=*|--env-name=*|--cuda=*)
      set -- "${1%%=*}" "${1#*=}" "${@:2}"
      continue ;;
    --no-pyrosetta|-n) install_pyrosetta=false; shift ;;
    --with-pyrosetta) install_pyrosetta=true; shift ;;
    --pkg_manager|-p|--env-name|--cuda|-c)
      [[ $# -ge 2 && -n "$2" && "$2" != -* ]] || { echo "Missing value for $1" >&2; exit 2; }
      case "$1" in
        --pkg_manager|-p) pkg_manager="$2" ;;
        --env-name) env_name="$2" ;;
        --cuda|-c)
          case "$2" in
            12|12.*|13|13.*) ;;
            *) echo 'This profile requires an NVIDIA driver compatible with CUDA 12.' >&2; exit 2 ;;
          esac
          printf 'The pip profile supplies CUDA 12.9; --cuda %s does not select a local toolkit.\n' "$2" ;;
      esac
      shift 2 ;;
    --fix-channels|-f)
      echo 'Using conda-forge explicitly; your Conda channel configuration is unchanged.'
      shift ;;
    --skip-gpu-check) skip_gpu_check=true; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown Linux installation option: $1" >&2; usage >&2; exit 2 ;;
  esac
done
[[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]] || {
  echo 'This Linux profile requires x86_64 (including its bundled native executables).' >&2; exit 1;
}
[[ "$pkg_manager" == conda || "$pkg_manager" == mamba ]] || {
  echo '--pkg_manager must be conda or mamba.' >&2; exit 2;
}
[[ "$env_name" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ && "$env_name" != base ]] || {
  echo 'Choose a named Conda environment other than base.' >&2; exit 2;
}
for tool in conda "$pkg_manager" git curl tar cc; do
  command -v "$tool" >/dev/null || { echo "Missing prerequisite: $tool" >&2; exit 1; }
done
SECONDS=0
conda_base="$(conda info --base)"
source "$conda_base/etc/profile.d/conda.sh"
if [[ "$skip_gpu_check" == false ]]; then
  "$conda_base/bin/python" extras/check_linux_driver.py
fi

# Resolve environments by name: envs_dirs need not be beneath the Conda base.
if ! conda env list --json | "$conda_base/bin/python" -c \
  'import json,pathlib,sys; sys.exit(not any(pathlib.Path(p).name == sys.argv[1] for p in json.load(sys.stdin)["envs"]))' "$env_name"; then
  "$pkg_manager" create --name "$env_name" --override-channels -c conda-forge \
    python=3.12 pip ffmpeg libgfortran5 -y
fi
conda activate "$env_name"
python - <<'PY'
import sys
if sys.version_info[:2] != (3, 12):
    raise SystemExit('The pip CUDA profile requires Python 3.12. Preserve this environment and rerun with --env-name <new-name>.')
PY
printf 'Installing the pinned pip CUDA profile into %s\n' "$CONDA_PREFIX"

# Scope CUDA library isolation to this environment and restore it on deactivation.
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d" "$CONDA_PREFIX/etc/conda/deactivate.d"
cp installation/linux-activate.sh "$CONDA_PREFIX/etc/conda/activate.d/freebindcraft-cuda.sh"
cp installation/linux-deactivate.sh "$CONDA_PREFIX/etc/conda/deactivate.d/freebindcraft-cuda.sh"
source "$CONDA_PREFIX/etc/conda/activate.d/freebindcraft-cuda.sh"
python -m pip install --requirement installation/linux-cuda12-requirements.txt
if [[ "$install_pyrosetta" == true ]]; then
  python -m pip install --no-deps --requirement installation/pyrosetta-linux-requirements.txt
fi
python -m pip check

for binary in dssp sc FASPR; do
  [[ -s "functions/$binary" ]] || { echo "Missing native helper: functions/$binary" >&2; exit 1; }
  chmod +x "functions/$binary"
done
if [[ "$install_pyrosetta" == true ]]; then
  chmod +x functions/DAlphaBall.gcc
fi

# Retain complete weights on a repeated install and resume interrupted downloads.
mkdir -p params
need_params=false
for model in 1 2 3 4 5; do
  for variant in '' _ptm _multimer_v3; do
    [[ -s "params/params_model_${model}${variant}.npz" ]] || need_params=true
  done
done
if [[ "$need_params" == true ]]; then
  archive=params/alphafold_params_2022-12-06.tar
  curl -fL --retry 3 --continue-at - -o "$archive" \
    https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar
  tar -tf "$archive" >/dev/null
  tar -xf "$archive" -C params --no-same-owner
  for model in 1 2 3 4 5; do
    for variant in '' _ptm _multimer_v3; do
      [[ -s "params/params_model_${model}${variant}.npz" ]] || {
        echo "Missing AlphaFold parameter file for model ${model}${variant}" >&2; exit 1;
      }
    done
  done
  rm "$archive"
fi

check_args=()
[[ "$install_pyrosetta" == true ]] || check_args+=(--no-pyrosetta)
[[ "$skip_gpu_check" == false ]] || check_args+=(--skip-gpu-check)
python extras/check_linux_installation.py "${check_args[@]}"
printf '\nInstalled in %s in %s seconds.\n' "$CONDA_PREFIX" "$SECONDS"
printf 'Activate with: conda activate %s\n' "$env_name"
if [[ "$skip_gpu_check" == true ]]; then
  printf 'GPU execution has not been checked; run extras/check_linux_installation.py on the GPU host.\n'
fi
