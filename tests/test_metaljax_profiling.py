"""Reject misleading provenance/count/overlap claims without using a GPU."""
import importlib.util
import json
from pathlib import Path
import struct
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('submission_analysis', ROOT/'extras/analyze_metaljax_submission.py')
submission = importlib.util.module_from_spec(spec)
spec.loader.exec_module(submission)


HLO = '''module {
  func.func @main(%weight: tensor<48x2xf32>, %sequence: tensor<48x2xf32>)
      -> (tensor<48x2xbf16>, tensor<48x2xbf16>, tensor<48x2xbf16>, tensor<48x2xbf16>) {
    %barrier:2 = stablehlo.optimization_barrier %sequence, %weight : tensor<48x2xf32>, tensor<48x2xf32>
    %activation = stablehlo.convert %barrier#0 : (tensor<48x2xf32>) -> tensor<48x2xbf16>
    %parameter = stablehlo.convert %barrier#1 : (tensor<48x2xf32>) -> tensor<48x2xbf16>
    %sum = stablehlo.add %weight, %sequence : tensor<48x2xf32>
    %unknown = stablehlo.convert %sum : (tensor<48x2xf32>) -> tensor<48x2xbf16>
    %zero = stablehlo.constant dense<0> : tensor<i32>
    %one = stablehlo.constant dense<1> : tensor<i32>
    %bound = stablehlo.constant dense<48> : tensor<i32>
    %loop:3 = stablehlo.while(%counter = %zero, %carry = %weight, %result = %parameter)
        : tensor<i32>, tensor<48x2xf32>, tensor<48x2xbf16>
     cond {
       %predicate = stablehlo.compare LT, %counter, %bound, SIGNED : (tensor<i32>, tensor<i32>) -> tensor<i1>
       stablehlo.return %predicate : tensor<i1>
     } do {
       %converted = stablehlo.convert %carry : (tensor<48x2xf32>) -> tensor<48x2xbf16>
       %next = stablehlo.add %counter, %one : tensor<i32>
       stablehlo.return %next, %carry, %converted : tensor<i32>, tensor<48x2xf32>, tensor<48x2xbf16>
     }
    return %activation, %parameter, %unknown, %loop#2 : tensor<48x2xbf16>, tensor<48x2xbf16>, tensor<48x2xbf16>, tensor<48x2xbf16>
  }
}'''


def analyze(tmp_path, dtype='float32', hlo=HLO):
    pytest.importorskip('jax')
    (tmp_path/'input.mlir').write_text(hlo)
    signature = dict(arguments=[dict(main_argument=i, shape=[48, 2], dtype=dtype,
        kind=kind, path=kind) for i, kind in enumerate(['model_weight', 'sequence'])])
    (tmp_path/'inputs.json').write_text(json.dumps(signature))
    result = subprocess.run([sys.executable, str(ROOT/'extras/analyze_metaljax_weight_casts.py'),
        '--hlo', str(tmp_path/'input.mlir'), '--inputs', str(tmp_path/'inputs.json'),
        '--out', str(tmp_path/'output.json')], capture_output=True, text=True)
    return result, tmp_path/'output.json'


def test_barrier_results_keep_separate_origins_and_static_loop_count(tmp_path):
    result, output = analyze(tmp_path)
    assert result.returncode == 0, result.stderr
    records = json.loads(output.read_text())['records']
    assert [r['model_weight_proven'] for r in records] == [False, True, False, True]
    assert [r['expanded_hlo_executions'] for r in records] == [1, 1, 1, 48]
    assert records[0]['origins'][0]['main_argument'] == 1
    assert records[1]['origins'][0]['main_argument'] == 0
    assert records[2]['origins'] == []


def test_mismatched_inventory_cannot_establish_weight_provenance(tmp_path):
    result, output = analyze(tmp_path, dtype='float16')
    assert result.returncode != 0
    assert 'Input index/shape/dtype mapping mismatch' in result.stderr
    assert not output.exists()


def test_updated_loop_carry_is_not_an_immutable_weight(tmp_path):
    hlo = HLO.replace('%next = stablehlo.add', '%changed = stablehlo.add %carry, %sequence : tensor<48x2xf32>\n       %next = stablehlo.add')
    hlo = hlo.replace('return %next, %carry, %converted', 'return %next, %changed, %converted')
    result, output = analyze(tmp_path, hlo=hlo)
    assert result.returncode == 0, result.stderr
    assert not json.loads(output.read_text())['records'][-1]['model_weight_proven']


def test_gpu_windows_use_union_instead_of_double_counting():
    assert submission.duration([(0, 10), (5, 20), (25, 30)]) == 25
    assert submission.overlap([(0, 10), (5, 20)], [(7, 12), (10, 25)]) == 13
    with pytest.raises(ValueError, match='Reversed'):
        submission.duration([(10, 0)])


def test_copy_bytes_use_reflected_element_count_and_sizes():
    record = dict(pipeline='g_copyfloat32bfloat16', sample=dict(
        constants_hex={'2': struct.pack('<i', 96).hex()}, grid=[128, 1, 1], grid_is_threadgroups=False))
    metadata = dict(arguments=[dict(name='src', index=0, data_type=3, data_size=4),
        dict(name='dst', index=1, data_type=121, data_size=2),
        dict(name='size', index=2, data_type=29, data_size=4)])
    decoded = submission.decode(record, metadata)
    assert decoded['elements'] == 96
    assert (decoded['logical_read_bytes'], decoded['logical_write_bytes']) == (384, 192)
