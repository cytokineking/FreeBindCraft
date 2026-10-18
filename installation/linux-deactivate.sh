# Restore the calling shell's settings when leaving this Conda environment.
if [ "${_FREEBINDCRAFT_CUDA_ACTIVE:-}" = 1 ]; then
    if [ "$_FREEBINDCRAFT_LD_SET" = x ]; then
        export LD_LIBRARY_PATH="$_FREEBINDCRAFT_LD_VALUE"
    else
        unset LD_LIBRARY_PATH
    fi
    if [ "$_FREEBINDCRAFT_JAX_SET" = x ]; then
        export JAX_PLATFORMS="$_FREEBINDCRAFT_JAX_VALUE"
    else
        unset JAX_PLATFORMS
    fi
    if [ "$_FREEBINDCRAFT_OPENMM_SET" = x ]; then
        export OPENMM_PLATFORM_ORDER="$_FREEBINDCRAFT_OPENMM_VALUE"
    else
        unset OPENMM_PLATFORM_ORDER
    fi
    unset _FREEBINDCRAFT_LD_SET _FREEBINDCRAFT_LD_VALUE
    unset _FREEBINDCRAFT_JAX_SET _FREEBINDCRAFT_JAX_VALUE
    unset _FREEBINDCRAFT_OPENMM_SET _FREEBINDCRAFT_OPENMM_VALUE
    unset _FREEBINDCRAFT_CUDA_ACTIVE
fi
