"""Run the unchanged attention probe with an optional passive Metal trace."""
import ctypes
import hashlib
import json
import os
from pathlib import Path
import runpy

library = None
if path := os.environ.get('FBC_MLX_TRACE_LIBRARY'):
    library = ctypes.CDLL(path, mode=ctypes.RTLD_GLOBAL)
if path := os.environ.get('FBC_MLX_MODEL_INVENTORY'):
    import colabdesign
    import jax
    import numpy as np
    class InspectLowering:
        def __init__(self, inner, name):
            self.inner, self.name = inner, name
        def __call__(self, *args, **kwargs):
            return self.inner(*args, **kwargs)
        def __getattr__(self, name):
            return getattr(self.inner, name)
        def lower(self, *args, **kwargs):
            lowered = self.inner.lower(*args, **kwargs)
            leaves = jax.tree_util.tree_flatten_with_path((args, kwargs))[0]
            kept = sorted(lowered._lowering.compile_args['kept_var_idx'])
            avals = lowered._lowering.compile_args['global_in_avals']
            if len(avals) != len(kept):
                raise RuntimeError('Unexpected JAX input ABI: investigate before using lineage')
            rows = []
            for main_index, index in enumerate(kept):
                key, value = leaves[index]
                positional = getattr(key[1], 'idx', None) if len(key) > 1 else None
                rows.append(dict(main_argument=main_index, original_flat_index=index,
                                 path=jax.tree_util.keystr(key), shape=list(avals[main_index].shape),
                                 dtype=str(avals[main_index].dtype),
                                 original_dtype=str(getattr(value, 'dtype', type(value).__name__)),
                                 kind='model_weight' if positional == 1 else 'sequence' if positional == 0 else 'other_input'))
            dest = Path(path).with_name(Path(path).stem+'-'+self.name+'-inputs.json')
            dest.write_text(json.dumps(dict(scope='JAX 0.11.2 kept_var_idx maps original input paths to lowered arguments; validate against the dumped main signature.', arguments=rows), indent=2)+'\n')
            return lowered
    create = colabdesign.mk_afdesign_model
    def inventory_model(*args, **kwargs):
        model = create(*args, **kwargs)
        prepare = model.prep_inputs
        def inventory_inputs(*args, **kwargs):
            result = prepare(*args, **kwargs)
            for name in ['fn', 'grad_fn']:
                model._model[name] = InspectLowering(model._model[name], name)
            rows = []
            for key, value in jax.tree_util.tree_flatten_with_path(model._model_params[0])[0]:
                rows.append(dict(path=jax.tree_util.keystr(key), shape=list(value.shape),
                                 dtype=str(value.dtype), python_type=type(value).__module__+'.'+type(value).__name__,
                                 bytes=int(value.size*value.dtype.itemsize),
                                 host_numpy_sha256=hashlib.sha256(value.tobytes()).hexdigest() if isinstance(value, np.ndarray) else None))
            Path(path).write_text(json.dumps(dict(scope='Model parameter inputs only; sequence-dependent activations are not weights. No array references retained.', weights=rows), indent=2)+'\n')
            return result
        model.prep_inputs = inventory_inputs
        return model
    colabdesign.mk_afdesign_model = inventory_model
try:
    runpy.run_path(str(Path(__file__).with_name('benchmark_metaljax_attention.py')), run_name='__main__')
finally:
    if library is not None:
        library.fbc_trace_dump()
