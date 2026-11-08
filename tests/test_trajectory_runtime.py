"""Exercise real ColabDesign preparation without loading AF2 weights."""
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
from colabdesign import mk_afdesign_model
from functions import colabdesign_utils as design

ROOT = Path(__file__).resolve().parents[1]


class TrajectoryRuntimeTests(unittest.TestCase):
    def setUp(self):
        self.settings = json.loads((ROOT / 'settings_advanced/default_4stage_multimer.json').read_text())
        self.pdb = str(ROOT / 'example/ubiquitin_1ubq.pdb')
        self.runtime = design.TrajectoryDesignRuntime()
        factory = patch.object(design, 'mk_afdesign_model',
                               side_effect=lambda **kw: mk_afdesign_model(model_names=[], **kw))
        self.build = factory.start()
        self.addCleanup(factory.stop)

    def prepare(self, length, seed=11, helix=-0.3, hotspot='44'):
        return self.runtime.prepare(self.pdb, 'A', hotspot, length, seed, helix, self.settings)

    def test_repeated_length_has_fresh_inputs_options_rng_and_history(self):
        model = self.prepare(66)
        callbacks = list(model._callbacks['model']['loss'])
        initial_functions = model._model
        model.opt.update(temp=.01, hard=1., num_recycles=7)
        model.aux = {'log': {}}
        model._tmp['best'] = {'sentinel': 1}
        model._inputs['prev'] = object()
        self.runtime.release_outputs()
        self.assertFalse(hasattr(model, 'aux'))
        self.assertNotIn('prev', model._inputs)
        self.assertFalse(model._tmp)

        warm = self.prepare(66, seed=22, helix=.6)
        self.assertIs(warm, model)
        self.assertIs(warm._model, initial_functions)
        self.assertEqual(warm._callbacks['model']['loss'], callbacks)
        self.assertEqual(self.build.call_count, 1)
        self.assertEqual(warm.opt['weights']['helix'], .6)
        self.assertEqual(warm._tmp['best'], {})
        self.assertEqual(warm._k, 0)
        warm_key = warm.key()
        self.runtime = design.TrajectoryDesignRuntime()
        cold = self.prepare(66, seed=22, helix=.6)
        for key in ('num_recycles', 'temp', 'hard', 'dropout', 'weights', 'con', 'i_con'):
            self.assertEqual(warm.opt[key], cold.opt[key])
        np.testing.assert_array_equal(warm._params['seq'], cold._params['seq'])
        np.testing.assert_array_equal(warm._inputs['bias'], cold._inputs['bias'])
        np.testing.assert_array_equal(warm_key, cold.key())

    def test_new_length_and_loss_topology_retire_the_old_runner(self):
        first = self.prepare(64)
        second = self.prepare(65)
        self.assertIsNot(second, first)
        self.assertEqual(self.runtime.lengths, {65})
        self.assertIs(self.prepare(65), second)
        self.assertEqual(self.build.call_count, 2)
        self.settings['use_termini_distance_loss'] = True
        updated = self.prepare(65)
        self.assertIsNot(updated, second)
        self.assertIn('NC', updated.opt['weights'])
        self.assertEqual(self.runtime.lengths, {65})

    def test_retirement_clears_jax_caches_without_deleting_other_live_buffers(self):
        self.prepare(64)
        unrelated = design.jnp.arange(5)
        with patch.object(design.jax, 'clear_caches', wraps=design.jax.clear_caches) as clear_caches, \
             patch.object(design, 'clear_mem', side_effect=AssertionError('unsafe buffer deletion')):
            self.prepare(65)
        self.assertEqual(clear_caches.call_count, 1)
        np.testing.assert_array_equal(np.asarray(unrelated), np.arange(5))

    def test_new_target_at_the_same_length_rebuilds(self):
        original = self.prepare(66)
        self.pdb = (ROOT / 'example/PDL1.pdb').read_text()
        updated = self.prepare(66)
        self.assertIsNot(original, updated)
        self.assertNotEqual(original._target_len, updated._target_len)
        self.assertIs(self.prepare(66), updated)

    def test_cuda_subbatch_boundary_uses_cold_preparation(self):
        with patch.object(design.jax, 'default_backend', return_value='gpu'):
            small = self.prepare(308)  # 76 + 308 = 384
            large = self.prepare(309)
            self.assertIsNot(small, large)
            self.assertEqual(large._cfg.model.global_config.subbatch_size, 4)
            self.assertIsNone(self.prepare(308)._cfg.model.global_config.subbatch_size)

    def test_failed_warm_prep_is_discarded(self):
        first = self.prepare(66)
        with patch.object(first, 'prep_inputs', side_effect=MemoryError('headroom')):
            with self.assertRaisesRegex(MemoryError, 'headroom'):
                self.prepare(66)
        self.assertIsNone(self.runtime.model)
        self.assertFalse(self.runtime.lengths)
        self.assertIsNot(self.prepare(66), first)


if __name__ == '__main__':
    unittest.main()
