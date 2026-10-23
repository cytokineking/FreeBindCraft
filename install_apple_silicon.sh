#!/usr/bin/env bash
# Native Apple Silicon installation; tools and environments stay in this checkout.
set -euo pipefail
cd "$(dirname "$0")"
repo_dir="$PWD"
install_pyrosetta=false
for argument in "$@"; do
  case "$argument" in
    --with-pyrosetta) install_pyrosetta=true ;;
    --no-pyrosetta|-n) install_pyrosetta=false ;;
    --help|-h)
      printf 'Usage: %s [--with-pyrosetta | --no-pyrosetta]\n' "$0"
      printf 'Requires Apple Command Line Tools. Downloads local uv/Rust, Python dependencies,\n'
      printf 'native helpers, and AlphaFold weights. The default omits PyRosetta.\n'
      printf 'MetalJAX: frozen FreeBindCraft wheel by default; FREEBINDCRAFT_METALJAX=stock opts out.\n'
      printf 'FREEBINDCRAFT_METALJAX_WHEEL=/absolute/path.whl takes precedence over that selection.\n'
      exit 0 ;;
    *) echo "Unknown Apple installation option: $argument" >&2; exit 2 ;;
  esac
done
if [[ "$(uname -s)" != Darwin || "$(uname -m)" != arm64 ]]; then
  echo 'This installer requires native Apple Silicon macOS.' >&2
  exit 1
fi
if ! xcrun --find clang++ >/dev/null 2>&1; then
  echo 'Apple Command Line Tools are missing. Run xcode-select --install, finish the Apple installer, then rerun this script.' >&2
  exit 1
fi
for tool in git curl tar shasum make; do
  command -v "$tool" >/dev/null || { echo "Missing required tool: $tool" >&2; exit 1; }
done
export UV_CACHE_DIR="$repo_dir/.uv-cache"
export UV_PYTHON_INSTALL_DIR="$repo_dir/.python"
export XDG_CACHE_HOME="$repo_dir/.local/cache"
export MAMBA_ROOT_PREFIX="$repo_dir/.micromamba"
export CONDA_PKGS_DIRS="$repo_dir/.micromamba/pkgs"
export CARGO_HOME="$repo_dir/.local/cargo"
export RUSTUP_HOME="$repo_dir/.local/rustup"
mkdir -p .local/bin .local/src .local/downloads params
export PATH="$repo_dir/.local/bin:$CARGO_HOME/bin:$PATH"

# xcrun can choose a leftover SDK newer than the installed linker understands.
# Verify a real compile, then try the Command Line Tools' default SDK symlink.
sdk="${SDKROOT:-$(xcrun --show-sdk-path)}"
check_sdk() {
  printf 'int main(void) { return 0; }\n' |
    xcrun clang -isysroot "$1" -x c - -o .local/sdk-check >.local/sdk-check.log 2>&1
}
if ! check_sdk "$sdk"; then
  fallback_sdk="$(dirname "$sdk")/MacOSX.sdk"
  if [[ ! -d "$fallback_sdk" ]] || ! check_sdk "$fallback_sdk"; then
    cat .local/sdk-check.log >&2
    echo 'The Apple compiler cannot link with the installed SDK. Repair Command Line Tools, then rerun.' >&2
    exit 1
  fi
  sdk="$fallback_sdk"
fi
export SDKROOT="$sdk"
printf 'Using Apple SDK: %s\n' "$SDKROOT"
rm -f .local/sdk-check

download_verified() {
  local url="$1" destination="$2" expected="$3" actual=''
  if [[ -f "$destination" ]]; then
    actual="$(shasum -a 256 "$destination" | awk '{print $1}')"
  fi
  if [[ "$actual" != "$expected" ]]; then
    curl -fLsS --retry 3 -o "$destination.part" "$url"
    actual="$(shasum -a 256 "$destination.part" | awk '{print $1}')"
    [[ "$actual" == "$expected" ]] || { echo "Checksum mismatch for $url" >&2; exit 1; }
    mv "$destination.part" "$destination"
  fi
}

# Bootstrap only inside this checkout; do not modify shell profiles.
if [[ ! -x .local/bin/uv ]]; then
  download_verified \
    https://github.com/astral-sh/uv/releases/download/0.13.1/uv-aarch64-apple-darwin.tar.gz \
    .local/downloads/uv-0.13.1.tar.gz \
    cfe454573412ff7380d15fc7dae5997ab53ffce72fce0a78e91ddbc84001acd2
  tar -xzf .local/downloads/uv-0.13.1.tar.gz -C .local/downloads
  install -m 755 .local/downloads/uv-aarch64-apple-darwin/uv .local/bin/uv
fi
if [[ ! -x "$CARGO_HOME/bin/rustup" ]]; then
  download_verified \
    https://static.rust-lang.org/rustup/archive/1.29.1/aarch64-apple-darwin/rustup-init \
    .local/downloads/rustup-init \
    ec1b9233e7f72990ecd8e62063fa7f6c3dfc2bec8e97f88bff165f9100ac696a
  chmod +x .local/downloads/rustup-init
  PATH=/usr/bin:/bin:/usr/sbin:/sbin .local/downloads/rustup-init \
    -y --profile minimal --default-toolchain 1.89.0 --no-modify-path
fi
"$CARGO_HOME/bin/rustup" toolchain install 1.89.0 --profile minimal
export RUSTUP_TOOLCHAIN=1.89.0
.local/bin/uv --version
"$CARGO_HOME/bin/cargo" --version
if [[ ! -x .venv-metaljax/bin/python ]]; then
  .local/bin/uv venv --python 3.12.15 --managed-python .venv-metaljax
fi
python_options=(--python "$repo_dir/.venv-metaljax/bin/python" --uv "$repo_dir/.local/bin/uv")
[[ "$install_pyrosetta" == true ]] && python_options+=(--with-pyrosetta)
.venv-metaljax/bin/python installation/install_apple_python.py \
  "${python_options[@]}"

ensure_source() {
  local url="$1" destination="$2" revision="$3"
  if [[ ! -d "$destination/.git" ]]; then
    git clone --no-checkout "$url" "$destination"
    git -C "$destination" checkout --detach "$revision"
  fi
  if [[ "$(git -C "$destination" rev-parse HEAD)" != "$revision" ]]; then
    echo "Unexpected source revision in $destination; preserving it. Expected $revision." >&2
    exit 1
  fi
}
ensure_source https://github.com/tommyhuangthu/FASPR.git .local/src/FASPR 0d55732fd6307f373018c6bddd842291c355c5f7
ensure_source https://github.com/cytokineking/sc-rs.git .local/src/sc-rs befcc6cf2d8d143ffd56219f9c92aff67485ad93
(cd .local/src/FASPR && clang++ -w -O3 -ffast-math -o FASPR src/*.cpp)
install -m 644 installation/sc-rs-Cargo.lock .local/src/sc-rs/Cargo.lock
"$CARGO_HOME/bin/cargo" build --locked --release --manifest-path .local/src/sc-rs/Cargo.toml
if [[ ! -x .local/bin/micromamba ]]; then
  curl -fLsS --retry 3 https://micro.mamba.pm/api/micromamba/osx-arm64/2.9.0 |
    tar -xj -C .local bin/micromamba
fi
if [[ ! -x .local/dssp/bin/mkdssp ]]; then
  .local/bin/micromamba create --no-rc -y -p "$repo_dir/.local/dssp" \
    --file installation/dssp-osx-arm64.txt
fi
.local/dssp/bin/mkdssp --version
export LIBCIFPP_DATA_DIR="$repo_dir/.local/dssp/share/libcifpp"
.local/dssp/bin/mkdssp --output-format dssp example/PDL1.pdb .local/pdl1-install-check.dssp

export DALPHABALL_BIN="$repo_dir/.local/src/DAlphaBall/DAlphaBall.gcc"
if [[ "$install_pyrosetta" == true ]]; then
  if [[ ! -x .local/dalphaball-toolchain/bin/gfortran ]]; then
    .local/bin/micromamba create --no-rc -y -p "$repo_dir/.local/dalphaball-toolchain" \
      --file installation/dalphaball-osx-arm64.txt
  fi
  .venv-metaljax/bin/python installation/build_dalphaball_apple.py
fi

need_params=false
for model in 1 2 3 4 5; do
  for variant in '' _ptm _multimer_v3; do
    [[ -s "params/params_model_${model}${variant}.npz" ]] || need_params=true
  done
done
if [[ "$need_params" == true ]]; then
  curl -fL --retry 3 --continue-at - -o params/alphafold_params_2022-12-06.tar \
    https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar
  tar -xf params/alphafold_params_2022-12-06.tar -C params
  rm params/alphafold_params_2022-12-06.tar
fi

# Check the actual GPU and optional native Rosetta installation.
export JAX_PLATFORMS=metal
export MPLCONFIGDIR="$repo_dir/.local/matplotlib"
export FREEBINDCRAFT_INSTALL_PYROSETTA="$install_pyrosetta"
.venv-metaljax/bin/python - <<'PY'
import math
import os
from importlib.metadata import version
import jax
import openmm
import colabdesign, freesasa, pdbfixer
from functions.generic_utils import check_jax_gpu
check_jax_gpu()
print("MetalJAX version:", version("metaljax"), "JAX:", version("jax"), "jaxlib:", version("jaxlib"))
print("JAX devices:", jax.devices())
platforms = [openmm.Platform.getPlatform(i).getName()
             for i in range(openmm.Platform.getNumPlatforms())]
print("OpenMM platforms:", platforms)
assert any(name in ("OpenCL", "Metal") for name in platforms)
if os.environ["FREEBINDCRAFT_INSTALL_PYROSETTA"] == "true":
    import pyrosetta
    pyrosetta.init("-mute all -holes:dalphaball " + os.environ["DALPHABALL_BIN"])
    pose = pyrosetta.pose_from_sequence("AAAAAA")
    score = pyrosetta.get_fa_scorefxn()(pose)
    assert math.isfinite(score)
    print("Native PyRosetta:", version("pyrosetta"), "test energy:", score)
PY
.local/bin/uv pip check --python .venv-metaljax/bin/python
printf '\nInstalled. Run: ./run_apple_silicon.sh --settings <target.json> --max-attempts 1\n'
printf 'Use --no-pyrosetta at runtime to select OpenMM/FASPR instead of Rosetta.\n'
