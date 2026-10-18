#!/bin/bash
if [[ "$(uname -s)" == Darwin ]]; then
    # Select the native installer before invoking GNU getopt or conda/CUDA.
    apple_mode=--with-pyrosetta
    for argument in "$@"; do
        case "$argument" in
            --no-pyrosetta|-n) apple_mode=--no-pyrosetta ;;
            --with-pyrosetta) apple_mode=--with-pyrosetta ;;
            --help|-h)
                printf 'Usage: %s [--no-pyrosetta | --with-pyrosetta]\n' "$0"
                printf 'On Apple Silicon, installs native GPU dependencies and PyRosetta by default.\n'
                printf 'Requires Apple Command Line Tools; remaining tools install inside this checkout.\n'
                exit 0 ;;
            *) echo "Unsupported macOS installation option: $argument" >&2; exit 2 ;;
        esac
    done
    exec bash "$(dirname "$0")/install_apple_silicon.sh" "$apple_mode"
fi
exec bash "$(dirname "$0")/install_linux.sh" "$@"
