# Sourced by Conda. CUDA user libraries come from pip, not an ambient toolkit.
if [ "${_FREEBINDCRAFT_CUDA_ACTIVE:-}" != 1 ]; then
    _FREEBINDCRAFT_LD_SET="${LD_LIBRARY_PATH+x}"
    _FREEBINDCRAFT_LD_VALUE="${LD_LIBRARY_PATH-}"
    _FREEBINDCRAFT_JAX_SET="${JAX_PLATFORMS+x}"
    _FREEBINDCRAFT_JAX_VALUE="${JAX_PLATFORMS-}"
    _FREEBINDCRAFT_OPENMM_SET="${OPENMM_PLATFORM_ORDER+x}"
    _FREEBINDCRAFT_OPENMM_VALUE="${OPENMM_PLATFORM_ORDER-}"
    _FREEBINDCRAFT_TMOL_SET="${TMOL_PYTHON+x}"
    _FREEBINDCRAFT_TMOL_VALUE="${TMOL_PYTHON-}"
    _FREEBINDCRAFT_CUDA_ACTIVE=1
fi
# Keep host-driver-only directories (e.g. /usr/lib64-nvidia on Google Colab)
# while stripping ambient CUDA toolkits so pip-installed CUDA 12.9 wins.
_fbc_driver_ld=""
_fbc_rest="${_FREEBINDCRAFT_LD_VALUE-}"
while [ -n "$_fbc_rest" ]; do
    case "$_fbc_rest" in
        *:*) _fbc_dir="${_fbc_rest%%:*}"; _fbc_rest="${_fbc_rest#*:}" ;;
        *)   _fbc_dir="$_fbc_rest"; _fbc_rest="" ;;
    esac
    case "$_fbc_dir" in
        ""|/lib|/lib64|/usr/lib|/usr/lib64|/lib/x86_64-linux-gnu|/usr/lib/x86_64-linux-gnu)
            continue ;;
    esac
    if [ -r "$_fbc_dir/libcuda.so.1" ] \
       && [ ! -r "$_fbc_dir/libcudart.so" ] \
       && [ ! -r "$_fbc_dir/libcudart.so.11.0" ] \
       && [ ! -r "$_fbc_dir/libcudart.so.12" ] \
       && [ ! -r "$_fbc_dir/libcudart.so.13" ]; then
        _fbc_driver_ld="${_fbc_driver_ld:+${_fbc_driver_ld}:}${_fbc_dir}"
    fi
done
if [ -n "$_fbc_driver_ld" ]; then
    export LD_LIBRARY_PATH="$_fbc_driver_ld"
else
    unset LD_LIBRARY_PATH
fi
unset _fbc_driver_ld _fbc_rest _fbc_dir
export JAX_PLATFORMS=cuda
export OPENMM_PLATFORM_ORDER=CUDA
if [ -n "${CONDA_PREFIX:-}" ] && [ -x "$CONDA_PREFIX/tmol-venv/bin/python" ]; then
    export TMOL_PYTHON="$CONDA_PREFIX/tmol-venv/bin/python"
fi
