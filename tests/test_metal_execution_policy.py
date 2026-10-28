"""Policy and memory decisions without a GPU or an installed dependency."""
import importlib.util
import hashlib
import json
import os
import sys
from pathlib import Path
import unittest
import numpy as np

path = Path(__file__).resolve().parents[1] / 'installation/colabdesign-metal/metal_execution_policy.py'
spec = importlib.util.spec_from_file_location('metal_execution_policy', path)
policy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(policy)


class Config(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__


class MetalExecutionPolicyTests(unittest.TestCase):
    def test_metal_is_unchunked_at_and_above_old_boundary(self):
        for length in (128, 383, 384, 385, 406):
            with self.subTest(length=length):
                config = Config(subbatch_size=4)
                changed, route = policy.select_global_policy(config, length, 'metal')
                self.assertTrue(changed)
                self.assertIsNone(config.subbatch_size)
                self.assertEqual(route['total_residues'], length)
                self.assertIsNone(route['global_subbatch_size'])

    def test_other_backends_keep_length_cutoff(self):
        for backend in ('cpu', 'gpu', 'cuda'):
            for length, expected in ((384, None), (385, 4)):
                with self.subTest(backend=backend, length=length):
                    config = Config(subbatch_size=None)
                    _, route = policy.select_global_policy(config, length, backend)
                    self.assertEqual(config.subbatch_size, expected)
                    self.assertEqual(route['global_subbatch_size'], expected)

    def test_reuse_detects_policy_changes_even_with_an_aliased_config(self):
        config = Config(subbatch_size=None)
        changed, _ = policy.select_global_policy(config, 384, 'cpu')
        self.assertFalse(changed)
        changed, _ = policy.select_global_policy(config, 385, 'cpu')
        self.assertTrue(changed)
        changed, _ = policy.select_global_policy(config, 384, 'cpu')
        self.assertTrue(changed)
        changed, route = policy.select_global_policy(config, 406, 'metal')
        self.assertFalse(changed)
        self.assertEqual(route['total_residues'], 406)

    def test_memory_preflight_respects_headroom_and_limits(self):
        g = 1024 ** 3
        result = policy.assess_metal_memory(406, 4, 2, total_bytes=128*g,
                                           available_bytes=80*g,
                                           system_claimed_bytes=48*g,
                                           process_bytes=0,
                                           system_ceiling_bytes=104*g,
                                           rematerialization=True, recycles=1)
        self.assertEqual(result['status'], 'admitted')
        self.assertEqual(result['score_tensor_bytes'], 406**3 * 4 * 2)
        self.assertEqual(result['reserved_headroom_bytes'], 8*g)
        self.assertEqual(result['usable_bytes'], 56*g)
        self.assertEqual(result['recycles'], 1)
        self.assertEqual(result['estimated_working_set_bytes'], 30*g)
        self.assertIn('not a memory guarantee', result['measured_reference']['note'])
        with self.assertRaisesRegex(MemoryError, 'reduce the target size'):
            policy.assess_metal_memory(2048, 4, 2, total_bytes=16*g,
                                       available_bytes=5*g)
        with self.assertRaisesRegex(MemoryError, 'larger-memory host'):
            policy.assess_metal_memory(406, 4, 2, total_bytes=128*g,
                                       available_bytes=80*g,
                                       device_limit_bytes=100*1024**2,
                                       device_used_bytes=90*1024**2)

    def test_unknown_memory_does_not_invent_a_budget(self):
        result = policy.assess_metal_memory(406, 4, 2, total_bytes=None,
                                           available_bytes=None)
        self.assertEqual(result['status'], 'unknown')
        self.assertIsNone(result['device_limit_bytes'])
        self.assertNotIn('usable_bytes', result)
        self.assertIsNone(result['process_footprint_bytes'])

    def test_working_set_estimate_and_uncalibrated_modes(self):
        g = 1024 ** 3
        with self.assertRaisesRegex(MemoryError, 'estimated unchunked Metal working set'):
            policy.assess_metal_memory(406, 4, 2, total_bytes=32*g,
                                       available_bytes=20*g, process_bytes=0)
        for options in ({'rematerialization': False}, {'recycle_mode': 'backprop'}):
            with self.subTest(options=options):
                result = policy.assess_metal_memory(406, 4, 2, total_bytes=128*g,
                                                   available_bytes=80*g, **options)
                self.assertEqual(result['status'], 'unknown')
                self.assertIsNone(result['estimated_working_set_bytes'])

    def test_warm_model_counts_only_proven_overlap(self):
        g = 1024 ** 3
        common = dict(total_bytes=128*g, available_bytes=80*g,
                      system_claimed_bytes=50*g, system_ceiling_bytes=104*g,
                      process_ceiling_bytes=40*g)
        cold = policy.assess_metal_memory(406, 4, 2, process_bytes=5*g,
                                          retained_model_bytes=5*g, **common)
        warm = policy.assess_metal_memory(406, 4, 2, process_bytes=20*g,
                                          retained_model_bytes=20*g, **common)
        for result in (cold, warm):
            self.assertEqual(result['status'], 'admitted')
            self.assertEqual(result['projected_process_footprint_bytes'], 30*g)
        self.assertEqual(warm['estimated_additional_bytes'], 10*g)
        self.assertEqual(warm['projected_system_claimed_bytes'], 60*g)
        self.assertEqual(warm['unclassified_process_footprint_bytes'], 0)
        with self.assertRaisesRegex(MemoryError, 'process physical-footprint limit'):
            policy.assess_metal_memory(406, 4, 2, process_bytes=20*g,
                                       retained_model_bytes=0,
                                       nonreusable_process_bytes=20*g, **common)
        unknown = policy.assess_metal_memory(406, 4, 2, process_bytes=20*g,
                                             retained_model_bytes=0, **common)
        self.assertEqual(unknown['status'], 'unknown')
        self.assertEqual(unknown['credited_retained_model_bytes'], 0)
        self.assertEqual(unknown['projected_process_footprint_bytes'], 50*g)
        self.assertEqual(unknown['unclassified_process_footprint_bytes'], 20*g)

    def test_cache_and_device_bytes_are_not_silently_credited(self):
        g = 1024 ** 3
        common = dict(total_bytes=128*g, available_bytes=80*g,
                      process_bytes=20*g, retained_model_bytes=10*g,
                      process_ceiling_bytes=35*g, device_limit_bytes=1*g,
                      device_used_bytes=900*1024**2)
        result = policy.assess_metal_memory(406, 4, 2,
                                           reported_device_cache_bytes=10*g, **common)
        self.assertEqual(result['status'], 'unknown')
        self.assertEqual(result['projected_process_footprint_bytes'], 40*g)
        self.assertEqual(result['reported_device_cache_bytes'], 10*g)
        self.assertTrue(any('device buffers' in text for text in result['uncertainty']))
        result = policy.assess_metal_memory(
            406, 4, 2, reclaimable_process_cache_bytes=10*g,
            **{**common, 'device_used_bytes': 0})
        self.assertEqual(result['status'], 'admitted')
        self.assertEqual(result['projected_process_footprint_bytes'], 30*g)

    def test_physical_footprint_and_claimed_limit_use_comparable_metrics(self):
        g = 1024 ** 3
        with self.assertRaisesRegex(MemoryError, 'system claimed-memory limit'):
            policy.assess_metal_memory(406, 4, 2, total_bytes=128*g,
                available_bytes=80*g, process_bytes=5*g, retained_model_bytes=5*g,
                nonreusable_process_bytes=0, system_claimed_bytes=95*g,
                system_ceiling_bytes=104*g)
        unknown = policy.assess_metal_memory(406, 4, 2, total_bytes=128*g,
            available_bytes=80*g, process_bytes=5*g, retained_model_bytes=5*g,
            system_claimed_bytes=None, system_ceiling_bytes=104*g)
        self.assertEqual(unknown['status'], 'unknown')
        self.assertTrue(any('claimed memory is unavailable' in text
                            for text in unknown['uncertainty']))
        self.assertIsNone(unknown['system_claimed_bytes'])
        # A limit below the current measured footprint is not admitted even
        # when the retained model alone exceeds the empirical estimate.
        with self.assertRaisesRegex(MemoryError, 'process physical-footprint limit'):
            policy.assess_metal_memory(406, 4, 2, total_bytes=128*g,
                available_bytes=80*g, process_bytes=45*g, retained_model_bytes=45*g,
                process_ceiling_bytes=40*g)

    @unittest.skipUnless(sys.platform == 'darwin', 'macOS resident-page query')
    def test_retained_pages_deduplicate_aliases_and_ignore_views(self):
        model_array = np.ones(8 * 1024 ** 2, dtype=np.uint8)
        page = os.sysconf('SC_PAGE_SIZE')
        retained, evidence = policy._verified_retained_model_pages(
            {'weight': model_array, 'same_weight': model_array,
             'view': model_array[100:]})
        self.assertGreaterEqual(retained, model_array.nbytes - 2 * page)
        self.assertLessEqual(retained, model_array.nbytes)
        self.assertEqual(evidence['buffers_checked'], 1)
        self.assertGreater(evidence['pages_queried'], 0)
        self.assertIsInstance(policy._process_physical_footprint(), int)
        self.assertGreater(policy._system_claimed_memory(), 0)

    def test_numpy_shape_metadata_is_json_serializable(self):
        config = Config(subbatch_size=4)
        _, route = policy.select_global_policy(config, np.int64(406), 'metal')
        route['memory'] = policy.assess_metal_memory(
            np.int64(406), np.int64(4), np.int64(2), total_bytes=None,
            available_bytes=None, recycles=np.int64(1))
        self.assertEqual(json.loads(json.dumps(route))['total_residues'], 406)
        self.assertEqual(route['memory']['recycles'], 1)

    def test_policy_asset_matches_installer_manifest(self):
        assets = path.parent
        manifest = json.loads((assets / 'manifest.json').read_text())
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),
                         manifest['files']['af/metal_execution_policy.py']['after'])
        self.assertEqual(hashlib.sha256((assets / manifest['patch']).read_bytes()).hexdigest(),
                         manifest['patch_sha256'])


if __name__ == '__main__':
    unittest.main()
