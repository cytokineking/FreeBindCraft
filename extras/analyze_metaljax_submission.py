"""Decode passive Metal traces using the captured shader reflection.

Logical copy bytes are element accesses, not measured DRAM transactions.
Containing-buffer GPU spans are inclusive upper bounds and overlap between
families; they must never be added into an exclusive cost breakdown.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import struct


FORMATS = {29: 'i', 33: 'I', 80: 'q', 81: 'Q', 3: 'f'}


def merge(intervals):
    result = []
    for start, end in sorted(intervals):
        if end < start:
            raise ValueError('Reversed timestamp interval')
        if result and start <= result[-1][1]:
            result[-1][1] = max(result[-1][1], end)
        else:
            result.append([start, end])
    return result


def duration(intervals):
    return sum(b-a for a, b in merge(intervals))


def overlap(left, right):
    left, right = merge(left), merge(right)
    i = j = total = 0
    while i < len(left) and j < len(right):
        total += max(0, min(left[i][1], right[j][1])-max(left[i][0], right[j][0]))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return total


def family(name):
    if 'dynamic_copy' in name:
        return 'dynamic_copy'
    if 'copy' in name:
        return 'cast_copy' if any(pair in name for pair in
            ['float32bfloat16', 'bfloat16float32', 'float32float16', 'float16float32']) else 'layout_or_dense_copy'
    if 'steel_gemm' in name or 'gemv' in name:
        return 'contraction'
    if 'compute_dynamic_offset' in name:
        return 'dynamic_offset'
    return 'other'


def decode(record, metadata):
    sample = record['sample']
    arguments = {a['name']: a for a in metadata['arguments']}
    constants = {int(k): bytes.fromhex(v) for k, v in sample['constants_hex'].items()}
    def values(name, count=1):
        arg = arguments.get(name)
        if arg is None or arg['index'] not in constants or arg['data_type'] not in FORMATS:
            return None
        fmt = '<'+FORMATS[arg['data_type']]*count
        data = constants[arg['index']]
        return list(struct.unpack_from(fmt, data)) if len(data) >= struct.calcsize(fmt) else None
    if 'copy' in record['pipeline']:
        count, shape = None, None
        size = values('size')
        ndim = values('ndim')
        if size:
            count = size[0]
            shape = [count]
        elif ndim and 0 < ndim[0] <= 16:
            shape = values('src_shape', ndim[0])
            if shape:
                count = math.prod(shape)
        elif not sample['grid_is_threadgroups']:
            # Fixed-rank g1/g2/g3 copies dispatch one logical element per
            # thread; vectorized general copies use reflected shape/size.
            prefix = record['pipeline'].split('_dynamic_copy')[0].split('_copy')[0]
            if prefix in {'g1', 'g2', 'g3', 'gg1', 'gg2', 'gg3'}:
                shape = list(reversed(sample['grid']))
                while len(shape) > 1 and shape[0] == 1:
                    shape.pop(0)
                count = math.prod(shape)
        if count is None or count < 0 or 'src' not in arguments or 'dst' not in arguments:
            return None
        src, dst = arguments['src'], arguments['dst']
        strides = values('src_strides', len(shape)) or values('src_stride')
        dst_strides = values('dst_strides', len(shape)) or values('dst_stride')
        return dict(kind='copy', shape=shape, elements=count,
                    source_element_bytes=src['data_size'], destination_element_bytes=dst['data_size'],
                    source_mtl_dtype=src['data_type'], destination_mtl_dtype=dst['data_type'],
                    source_strides=strides, destination_strides=dst_strides,
                    logical_read_bytes=count*src['data_size'], logical_write_bytes=count*dst['data_size'])
    params = arguments.get('params')
    if params and params['index'] in constants:
        raw = constants[params['index']]
        decoded = {}
        for member in params['members']:
            fmt = FORMATS.get(member['data_type'])
            if fmt and member['offset']+struct.calcsize(fmt) <= len(raw):
                decoded[member['name']] = struct.unpack_from('<'+fmt, raw, member['offset'])[0]
        if all(key in decoded for key in ['M', 'N', 'K']):
            # Fused GEMM z indexes a batch. Split-K variants require their
            # own partition accounting and are left unclassified here.
            batch = sample['grid'][2] if 'splitk' not in record['pipeline'].lower() else None
            return dict(kind='contraction', params=decoded, batch=batch,
                        dtype_bytes=arguments.get('A', {}).get('data_size'),
                        logical_flops=2*decoded['M']*decoded['N']*decoded['K']*batch if batch is not None else None)
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trace', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    rows = [json.loads(line) for line in args.trace.read_text().splitlines()]
    metadata = {(r['phase'], r['pipeline']): r for r in rows if r['type'] == 'kernel'}
    phases = sorted({r['phase'] for r in rows})
    result = dict(trace_sha256=hashlib.sha256(args.trace.read_bytes()).hexdigest(), phases={})
    for phase in phases:
        buffers = [r for r in rows if r['type'] == 'buffer' and r['phase'] == phase]
        kernels = [r for r in rows if r['type'] == 'kernel' and r['phase'] == phase]
        operations = [r for r in rows if r['type'] in ('copy', 'operation') and r['phase'] == phase]
        if not all(b['matched'] and b['dispatches'] == b['native_operations'] for b in buffers):
            raise ValueError('Trace/native command-buffer identity or dispatch mismatch')
        if sum(b['dispatches'] for b in buffers) != sum(k['count'] for k in kernels):
            raise ValueError('Kernel counts do not cover the submitted workload')
        gpu = [(b['gpu_start_ns'], b['gpu_end_ns']) for b in buffers]
        host = [(b['host_encode_begin_ns'], b['host_commit_ns']) for b in buffers if b.get('host_encode_begin_ns')]
        grouped, unknown = defaultdict(lambda: dict(count=0, logical_read_bytes=0, logical_write_bytes=0)), Counter()
        for record in operations:
            decoded = decode(record, metadata[(phase, record['pipeline'])])
            if decoded is None:
                unknown[record['pipeline']] += 1
                continue
            key = json.dumps(dict(pipeline=record['pipeline'], **decoded), sort_keys=True)
            value = grouped[key]
            value['count'] += 1
            value['logical_read_bytes'] += decoded.get('logical_read_bytes', 0)
            value['logical_write_bytes'] += decoded.get('logical_write_bytes', 0)
        groups = [dict(operation=json.loads(key), **value) for key, value in grouped.items()]
        groups.sort(key=lambda x: -(x['logical_read_bytes']+x['logical_write_bytes']))
        bounds = {}
        for name in ['dynamic_copy', 'cast_copy', 'layout_or_dense_copy', 'contraction', 'dynamic_offset']:
            intervals = [(b['gpu_start_ns'], b['gpu_end_ns']) for b in buffers
                         if any(family(k) == name for k in b.get('kernel_counts', {}))]
            bounds[name] = dict(containing_buffer_gpu_span_ns=duration(intervals),
                                frequency=sum(k['count'] for k in kernels if family(k['pipeline']) == name),
                                exclusive_cpu_metal_dispatch_call_ns=sum(k['cpu_dispatch_ns'] for k in kernels if family(k['pipeline']) == name))
        result['phases'][phase] = dict(buffers=len(buffers), dispatches=sum(k['count'] for k in kernels),
            thresholds=dict(Counter(b['threshold'] for b in buffers)), gpu_span_ns=duration(gpu),
            host_encoding_window_union_ns=duration(host), host_encoding_windows_overlapping_gpu_ns=overlap(host, gpu),
            inclusive_family_bounds=bounds, operations=groups, unclassified=dict(unknown),
            scope='Instrumented execution only. Host windows include pauses; CPU dispatch calls cover only the original Metal dispatch APIs. Family buffer spans overlap and are loose GPU-only upper bounds. Logical element bytes are not physical DRAM traffic.')
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    for name, phase in result['phases'].items():
        print(name, 'buffers', phase['buffers'], 'dispatches', phase['dispatches'], 'groups', len(phase['operations']),
              'unclassified', sum(phase['unclassified'].values()))


if __name__ == '__main__':
    main()
