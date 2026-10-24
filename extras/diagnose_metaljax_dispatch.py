"""Diagnostic-only model probe: direct Metal hooks, settled counters, call arrays.

Timing from this process is excluded from performance conclusions. The hook is
restricted to the checksum-pinned private MLX ABI; it is not a portable API.
"""
import ctypes
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import benchmark_metaljax_attention as probe


library = ctypes.CDLL(os.environ['FBC_MLX_TRACE_LIBRARY'], mode=ctypes.RTLD_GLOBAL)
library.fbc_trace_wait.argtypes = []
library.fbc_trace_wait.restype = ctypes.c_uint64
library.fbc_trace_stats.argtypes = [ctypes.POINTER(probe.DispatchStats)]
library.fbc_trace_stats.restype = None
library.fbc_trace_mlx_image.argtypes = []
library.fbc_trace_mlx_image.restype = ctypes.c_char_p
mlx_image = Path(library.fbc_trace_mlx_image().decode())
mlx_sha256 = hashlib.sha256(mlx_image.read_bytes()).hexdigest()
if mlx_sha256 != '72dc8d9a4b1e07df2f6bc91383e0567acf3ce292398ae6ce4ed1d19fce0b1374':
    raise RuntimeError('Counter reader and tracer MLX identity is unrecognized')
history = []


def diagnostic_reader():
    def snapshot():
        rows = []
        start = time.monotonic()
        for delay in [0, .001, .005, .02, .1]:
            if delay:
                time.sleep(delay)
            value = probe.DispatchStats()
            library.fbc_trace_stats(ctypes.byref(value))
            rows.append(dict(elapsed_seconds=time.monotonic()-start,
                             **{k: getattr(value, k) for k, _ in probe.DispatchStats._fields_}))
        if library.fbc_trace_wait():
            raise RuntimeError('Diagnostic trace callbacks did not settle')
        last = rows[-1]
        if last['command_buffers'] != last['completed_command_buffers']:
            raise RuntimeError('Native completion accounting did not settle')
        history.append(dict(phase=os.environ.get('FBC_PROFILE_ACTIVE'), snapshots=rows))
        return {k: last[k] for k, _ in probe.DispatchStats._fields_}
    return snapshot


probe.dispatch_reader = diagnostic_reader
try:
    probe.main()
finally:
    pending = library.fbc_trace_wait()
    library.fbc_trace_dump()
    output = Path(os.environ['FBC_MLX_TRACE_OUTPUT'])
    output.with_suffix('.snapshots.json').write_text(json.dumps(dict(
        pending_trace_callbacks=pending, boundaries=history,
        mlx_image=str(mlx_image), mlx_sha256=mlx_sha256,
        note='Every warmed call has arrays; snapshots and hooks are diagnostic-only.'), indent=2)+'\n')
