"""Prove model-weight provenance of whole-stack casts in captured StableHLO.

Follows only function arguments, unchanged while carries, identity barriers and
view/slice/cast operations. Unknown provenance stays unknown. Execution counts are
source-level predictions before MLX fusion; corroborate with the Metal trace.
"""
import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import re

from jax._src.interpreters import mlir
from jaxlib.mlir import ir


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hlo', type=Path, required=True)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    signature = json.loads(args.inputs.read_text())['arguments']
    with mlir.make_ir_context():
        module = ir.Module.parse(args.hlo.read_text())
        functions = {str(f.attributes['sym_name']).strip('"'): f.operation for f in module.body.operations}
        calls, casts = defaultdict(list), []
        def walk(block):
            for view in block.operations:
                op = view.operation
                if op.name == 'func.call':
                    calls[str(op.attributes['callee']).lstrip('@')].append(op)
                if op.name == 'stablehlo.convert':
                    src, dst = ir.RankedTensorType(op.operands[0].type), ir.RankedTensorType(op.results[0].type)
                    if str(src.element_type) == 'f32' and str(dst.element_type) == 'bf16':
                        casts.append(op)
                for region in op.regions:
                    for nested in region.blocks:
                        walk(nested)
        for f in functions.values():
            walk(f.regions[0].blocks[0])
        main_block = functions['main'].regions[0].blocks[0]
        if len(main_block.arguments) != len(signature):
            raise ValueError('Input map does not match main arguments')
        dtypes = {'f32': 'float32', 'bf16': 'bfloat16', 'f16': 'float16', 'f64': 'float64',
                  'i1': 'bool', 'i8': 'int8', 'ui8': 'uint8', 'i16': 'int16', 'ui16': 'uint16',
                  'i32': 'int32', 'ui32': 'uint32', 'i64': 'int64', 'ui64': 'uint64'}
        for index, (value, row) in enumerate(zip(main_block.arguments, signature)):
            t = ir.RankedTensorType(value.type)
            if (row['main_argument'] != index or list(t.shape) != row['shape'] or
                dtypes.get(str(t.element_type)) != row['dtype']):
                raise ValueError('Input index/shape/dtype mapping mismatch')
        def integer(value):
            owner = value.owner
            if hasattr(owner, 'operation'):
                owner = owner.operation
            if not isinstance(owner, ir.Operation) or owner.name != 'stablehlo.constant':
                return None
            if not re.fullmatch(r'tensor<(?:u?i)(?:32|64)>', str(value.type)):
                return None
            m = re.match(r'dense<(-?\d+)>', str(owner.attributes['value']))
            return int(m[1]) if m else None
        def trip(loop):
            cond = loop.regions[0].blocks[0]
            compares = [v.operation for v in cond.operations if v.operation.name == 'stablehlo.compare']
            if len(compares) != 1:
                return None
            cmp = compares[0]
            counter = cmp.operands[0]
            if (not isinstance(counter, ir.BlockArgument) or counter.owner != cond or
                'LT' not in str(cmp.attributes['comparison_direction']) or
                len(cond.operations[-1].operands) != 1 or
                cond.operations[-1].operands[0] != cmp.results[0]):
                return None
            index, bound = counter.arg_number, integer(cmp.operands[1])
            start = integer(loop.operands[index])
            body = loop.regions[1].blocks[0]
            update = body.operations[-1].operands[index].owner
            if hasattr(update, 'operation'):
                update = update.operation
            if not isinstance(update, ir.Operation) or update.name != 'stablehlo.add':
                return None
            left, right = update.operands
            step = integer(right) if left == body.arguments[index] else integer(left) if right == body.arguments[index] else None
            return max(0, bound-start) if step == 1 and bound is not None and start is not None and 0 <= bound-start <= 4096 else None
        def roots(value, seen=frozenset()):
            if value in seen:
                return set()
            seen = seen | {value}
            if isinstance(value, ir.BlockArgument):
                block, index = value.owner, value.arg_number
                owner = block.owner.operation
                if owner.name == 'func.func':
                    name = str(owner.attributes['sym_name']).strip('"')
                    if name == 'main':
                        row = signature[index]
                        return {(row['kind'], row['main_argument'], row['path'])}
                    sites = calls[name]
                    if not sites:
                        return set()
                    result = set()
                    for site in sites:
                        found = roots(site.operands[index], seen)
                        if not found:
                            return set()
                        result |= found
                    return result
                if owner.name == 'stablehlo.while':
                    body = owner.regions[1].blocks[0]
                    if body.operations[-1].operands[index] != body.arguments[index]:
                        return set()
                    return roots(owner.operands[index], seen)
                return set()
            owner = value.owner
            if hasattr(owner, 'operation'):
                owner = owner.operation
            # JAX rematerialization uses a multi-result identity barrier. Each
            # result comes from the operand at the same position; tracing all
            # operands together would incorrectly label activation results as
            # model weights. This is analysis only, not a compiler rewrite.
            if owner.name == 'stablehlo.optimization_barrier':
                index = value.result_number
                if len(owner.results) != len(owner.operands) or value.type != owner.operands[index].type:
                    return set()
                return roots(owner.operands[index], seen)
            if owner.name in {'stablehlo.reshape', 'stablehlo.transpose', 'stablehlo.broadcast_in_dim',
                              'stablehlo.slice', 'stablehlo.dynamic_slice', 'stablehlo.convert'}:
                return roots(owner.operands[0], seen)
            return set()
        def frequency(op, active=frozenset()):
            multiplier = 1
            child = op
            parent = op.parent
            while parent is not None:
                if parent.name == 'stablehlo.while':
                    n = trip(parent)
                    if n is None:
                        return None
                    if child.block == parent.regions[0].blocks[0]:
                        n += 1  # The condition is evaluated once after the final body.
                    multiplier *= n
                elif parent.name not in {'func.func', 'builtin.module'}:
                    # Conditional regions and unrecognized structured control flow
                    # do not establish an exact dynamic execution count.
                    return None
                if parent.name == 'func.func':
                    name = str(parent.attributes['sym_name']).strip('"')
                    if name == 'main':
                        return multiplier
                    if name in active:
                        return None
                    counts = [frequency(site, active | {name}) for site in calls[name]]
                    return multiplier*sum(counts) if counts and all(c is not None for c in counts) else None
                child, parent = parent, parent.parent
            return None
        records = []
        for op in casts:
            output = ir.RankedTensorType(op.results[0].type)
            shape = list(output.shape)
            if not shape or shape[0] != 48:
                continue
            origins = roots(op.operands[0])
            parent = op.parent
            while parent is not None and parent.name != 'func.func':
                parent = parent.parent
            records.append(dict(shape=shape, elements=math.prod(shape),
                                model_weight_proven=bool(origins) and all(o[0] == 'model_weight' for o in origins),
                                origins=[dict(kind=o[0], main_argument=o[1], path=o[2]) for o in sorted(origins)],
                                source_function=str(parent.attributes['sym_name']).strip('"') if parent else None,
                                expanded_hlo_executions=frequency(op), operation=str(op)))
        args.out.write_text(json.dumps(dict(scope=__doc__,
            hlo_sha256=hashlib.sha256(args.hlo.read_bytes()).hexdigest(),
            inputs_sha256=hashlib.sha256(args.inputs.read_bytes()).hexdigest(),
            records=records), indent=2)+'\n')
        print('whole-stack casts', len(records), 'proven weights', sum(r['model_weight_proven'] for r in records))


if __name__ == '__main__':
    main()
