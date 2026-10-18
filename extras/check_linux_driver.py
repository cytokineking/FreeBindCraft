#!/usr/bin/env python3
"""Check the driver needed by the pinned CUDA 12.9 Update 1 PTX runtime."""
import subprocess

MINIMUM_DRIVER = (575, 57, 8)

def check_driver():
    try:
        output = subprocess.check_output(
            ['nvidia-smi', '--query-gpu=driver_version', '--format=csv,noheader'],
            text=True, stderr=subprocess.STDOUT)
    except (OSError, subprocess.CalledProcessError) as error:
        raise SystemExit('Cannot query the NVIDIA driver with nvidia-smi. Install a working driver, '
                         'or use --skip-gpu-check only for a container image build.') from error
    versions = sorted(set(line.strip() for line in output.splitlines() if line.strip()))
    if not versions:
        raise SystemExit('nvidia-smi did not report a GPU driver.')
    for version in versions:
        try:
            parsed = tuple(int(part) for part in version.split('.'))
        except ValueError as error:
            raise SystemExit(f'Cannot parse NVIDIA driver version: {version}') from error
        if parsed + (0,) * max(0, 3 - len(parsed)) < MINIMUM_DRIVER:
            raise SystemExit(f'NVIDIA driver {version} is too old for the pinned CUDA 12.9 profile. '
                             'Install Linux driver 575.57.08 or newer. OpenMM compiles PTX at runtime; '
                             'CUDA 12 minor-version compatibility alone is insufficient.')
    return versions

if __name__ == '__main__':
    print('NVIDIA driver check passed:', ', '.join(check_driver()))
