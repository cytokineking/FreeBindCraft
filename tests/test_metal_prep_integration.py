"""Exercise a patched dependency copy with fake runners, without compiling AF2."""
import ast
import copy
import importlib.util
import os
import json
import shutil
import tempfile
from pathlib import Path
import sys
import types
import unittest
from unittest.mock import patch


class Config(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__

    def __deepcopy__(self, memo):
        return Config({key: copy.deepcopy(value, memo) for key, value in self.items()})


class PrepIntegrationTests(unittest.TestCase):
    def setUp(self):
        environment = patch.dict(os.environ, {'FREEBINDCRAFT_METAL_CHUNKING': '0'})
        environment.start()
        self.addCleanup(environment.stop)
        helper_path = Path(__file__).resolve().parents[1] / 'installation/colabdesign-metal/metal_execution_policy.py'
        spec = importlib.util.spec_from_file_location('policy_under_test', helper_path)
        self.helper = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.helper)
        override = os.environ.get('FBC_TEST_COLABDESIGN_PREP')
        if override:
            source = Path(override)
        else:
            package = importlib.util.find_spec('colabdesign')
            if package is None:
                raise unittest.SkipTest('ColabDesign is not installed')
            installed = Path(package.origin).parent
            temporary = tempfile.TemporaryDirectory()
            self.addCleanup(temporary.cleanup)
            staged = Path(temporary.name) / 'colabdesign'
            assets = helper_path.parent
            manifest = json.loads((assets / 'manifest.json').read_text())
            for name in manifest['files']:
                if (installed / name).is_file():
                    target = staged / name
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(installed / name, target)
            installer_path = assets.parent / 'apply_colabdesign_metal.py'
            installer_spec = importlib.util.spec_from_file_location('prep_test_installer', installer_path)
            installer = importlib.util.module_from_spec(installer_spec)
            installer_spec.loader.exec_module(installer)
            installer.apply(staged, assets)
            source = staged / 'af/prep.py'
        parsed = ast.parse(source.read_text())
        owner = next(node for node in parsed.body if isinstance(node, ast.ClassDef) and node.name == '_af_prep')
        methods = [node for node in owner.body if isinstance(node, ast.FunctionDef)
                   and node.name in ('_prep_model', 'get_execution_policy')]
        self.assertTrue(any(node.name == 'get_execution_policy' for node in methods),
                        'Install the checked Metal preparation patch before running this test')
        klass = ast.ClassDef(name='Prepared', bases=[], keywords=[], body=methods, decorator_list=[])
        module = ast.fix_missing_locations(ast.Module(body=[klass], type_ignores=[]))
        self.backend = types.SimpleNamespace(name='metal')
        jax = types.SimpleNamespace(default_backend=lambda: self.backend.name,
                                   devices=lambda: [types.SimpleNamespace(memory_stats=lambda: {})])
        scope = {'jax': jax, 'copy_dict': copy.deepcopy}
        exec(compile(module, str(source), 'exec'), scope)
        self.model = scope['Prepared']()
        global_config = Config(subbatch_size=4, bfloat16=True, use_remat=True)
        self.model._cfg = Config(model=Config(global_config=global_config, num_recycle=1,
            embeddings_and_evoformer=Config(evoformer=Config(
                triangle_attention_starting_node=Config(num_head=4)))))
        self.model._lengths = [256, 150]
        self.model.opt = {'dropout': True, 'num_recycles': 1}
        self.model._args = {'recycle_mode': 'last', 'use_multimer': True}
        self.model._model_params = {'retained': object()}
        self.model.restart = lambda **kwargs: None
        self.build_sizes = []

        def build(cfg):
            size = cfg.model.global_config.subbatch_size
            self.build_sizes.append(size)
            return {'runner': types.SimpleNamespace(config=cfg), 'fn': ('forward', size), 'grad_fn': ('gradient', size)}

        self.model._get_model = build
        modules = patch.dict(sys.modules, {'colabdesign.af.metal_execution_policy': self.helper})
        modules.start()
        self.addCleanup(modules.stop)

    def test_both_metal_callables_are_unchunked_across_boundary(self):
        with patch.object(self.helper, 'current_metal_memory', return_value={'status': 'admitted'}):
            for length in (383, 384, 385, 406, 146):
                self.model._lengths = [length]
                self.model._prep_model()
                self.assertEqual(self.model._model['fn'], ('forward', None))
                self.assertEqual(self.model._model['grad_fn'], ('gradient', None))
                self.assertEqual(self.model.get_execution_policy()['total_residues'], length)
            self.assertEqual(self.build_sizes, [None])

    def test_refusal_and_failed_build_preserve_existing_config(self):
        original = self.model._cfg
        old_runner = {'runner': types.SimpleNamespace(config=original)}
        self.model._model = old_runner
        self.model._execution_policy = {'previous': True}
        with patch.object(self.helper, 'current_metal_memory', side_effect=MemoryError('low headroom')):
            with self.assertRaisesRegex(MemoryError, 'low headroom'):
                self.model._prep_model()
        self.assertIs(self.model._cfg, original)
        self.assertEqual(original.model.global_config.subbatch_size, 4)
        self.assertIs(self.model._model, old_runner)
        with patch.object(self.helper, 'current_metal_memory', return_value={'status': 'admitted'}):
            with patch.object(self.model, '_get_model', side_effect=RuntimeError('failed build')):
                with self.assertRaisesRegex(RuntimeError, 'failed build'):
                    self.model._prep_model()
            self.assertIs(self.model._cfg, original)
            self.assertEqual(original.model.global_config.subbatch_size, 4)
            self.assertEqual(self.model.get_execution_policy(), {'previous': True})
            self.model._prep_model()
        self.assertIsNot(self.model._cfg, original)
        self.assertIsNone(self.model._cfg.model.global_config.subbatch_size)
        self.assertEqual(original.model.global_config.subbatch_size, 4)

    def test_chunking_opt_in_updates_both_callables_and_can_be_disabled(self):
        with patch.dict(os.environ, {'FREEBINDCRAFT_METAL_CHUNKING': '1'}), \
             patch.object(self.helper, 'current_metal_memory', side_effect=AssertionError('Unchunked estimate must not reject chunked jobs')):
            self.model._prep_model()
            self.assertEqual(self.model._model['fn'], ('forward', 4))
            self.assertEqual(self.model._model['grad_fn'], ('gradient', 4))
            self.assertEqual(self.model.get_execution_policy()['memory']['status'], 'native_guards_only')
        with patch.object(self.helper, 'current_metal_memory', return_value={'status': 'admitted'}):
            self.model._prep_model()
            self.assertEqual(self.model._model['fn'], ('forward', None))
            self.assertEqual(self.model._model['grad_fn'], ('gradient', None))
        self.assertEqual(self.build_sizes, [4, None])

    def test_chunking_opt_in_keeps_small_jobs_unchunked_on_model_reuse(self):
        with patch.dict(os.environ, {'FREEBINDCRAFT_METAL_CHUNKING': '1'}), \
             patch.object(self.helper, 'current_metal_memory', return_value={'status': 'admitted'}) as memory:
            for length, expected in ((384, None), (385, 4), (406, 4), (146, None)):
                self.model._lengths = [length]
                self.model._prep_model()
                self.assertEqual(self.model._model['fn'], ('forward', expected))
                self.assertEqual(self.model._model['grad_fn'], ('gradient', expected))
            self.assertEqual(memory.call_count, 2)
            self.assertEqual(self.build_sizes, [None, 4, None])

    def test_nonmetal_preserves_pinned_callable_selection(self):
        with patch.object(self.helper, 'current_metal_memory', side_effect=AssertionError('Metal-only preflight')):
            for backend in ('cpu', 'gpu'):
                for length in (384, 385):
                    with self.subTest(backend=backend, length=length):
                        self.backend.name = backend
                        if hasattr(self.model, '_model'):
                            del self.model._model
                        self.build_sizes.clear()
                        self.model._lengths = [length]
                        self.model._prep_model()
                        expected = 4 if length > 384 else None
                        self.assertEqual(self.model._model['fn'], ('forward', expected))
                        self.assertEqual(self.model._model['grad_fn'], ('gradient', None))
                        self.assertEqual(self.build_sizes, [None, 4] if length > 384 else [None])

    def test_warm_model_retained_bytes_are_not_counted_twice(self):
        gib = 1024 ** 3
        samples = {'warm': False, 'process': 5*gib}

        def retained(params):
            self.assertIs(params, self.model._model_params if samples['warm'] else None)
            return (20*gib if samples['warm'] else None), {'source': 'test fixture'}

        fake_psutil = types.SimpleNamespace(virtual_memory=lambda: types.SimpleNamespace(total=128*gib, available=80*gib))
        with patch.dict(sys.modules, {'psutil': fake_psutil}), \
             patch.dict(os.environ, {'METALJAX_MEM_SYS_MB': '106496', 'METALJAX_MEM_BUDGET_MB': '40960'}), \
             patch.object(self.helper, '_verified_retained_model_pages', retained), \
             patch.object(self.helper, '_system_claimed_memory', return_value=50*gib), \
             patch.object(self.helper, '_process_physical_footprint', side_effect=lambda: samples['process']):
            self.model._prep_model()
            samples.update(warm=True, process=20*gib)
            self.model._prep_model()
            details = self.model.get_execution_policy()
            self.assertEqual(details['memory']['status'], 'admitted')
            self.assertEqual(details['memory']['projected_process_footprint_bytes'], 30*gib)
            details['memory']['status'] = 'mutated outside the model'
            self.assertEqual(self.model.get_execution_policy()['memory']['status'], 'admitted')
            self.assertEqual(self.build_sizes, [None])


if __name__ == '__main__':
    unittest.main()
