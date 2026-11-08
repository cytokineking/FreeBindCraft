"""Rank real partial outputs without importing the scientific runtime."""
import ast
import copy
import os
from pathlib import Path
import shutil
import tempfile
import time
import unittest
from unittest.mock import Mock

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
source = ast.parse((ROOT / 'functions/generic_utils.py').read_text())
function = next(n for n in source.body if isinstance(n, ast.FunctionDef)
                and n.name == 'check_accepted_designs')


class PartialRankingTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.accepted = self.root / 'Accepted'
        self.ranked = self.accepted / 'Ranked'
        self.ranked.mkdir(parents=True)
        self.paths = {'Accepted': str(self.accepted), 'Accepted/Ranked': str(self.ranked),
                      'Trajectory/Animation': 'animations', 'Trajectory/Plots': 'plots'}
        self.labels = ['Design', 'Average_i_pTM', 'Average_ipSAE']
        self.final_labels = ['Rank', *self.labels]
        self.csv = self.root / 'mpnn.csv'
        self.final = self.root / 'final.csv'
        self.target = {'number_of_final_designs': 10, 'binder_name': 'UBQ'}
        self.archive = Mock()
        env = dict(os=os, pd=pd, shutil=shutil, time=time, vprint=lambda *a: None,
                   zip_and_empty_folder=self.archive)
        exec(compile(ast.Module(body=[function], type_ignores=[]),
                     str(ROOT / 'functions/generic_utils.py'), 'exec'), env)
        self.rank = env['check_accepted_designs']

    def populate(self):
        pd.DataFrame([['low', .6, .9], ['high', .9, .6], ['rejected', .99, .99]],
                     columns=self.labels).to_csv(self.csv, index=False)
        for name in ('low', 'high'):
            (self.accepted / f'{name}_model1.pdb').write_text(name)
        (self.ranked / 'stale.pdb').write_text('stale')

    def call(self, **kwargs):
        return self.rank(self.paths, self.csv, self.final_labels, self.final,
                         {'zip_animations': True, 'zip_plots': True},
                         self.target, self.labels, **kwargs)

    def test_partial_ranking_uses_only_accepted_designs_and_preserves_inputs(self):
        self.populate()
        original = copy.deepcopy(self.target)
        self.assertFalse(self.call(finalize=True))
        result = pd.read_csv(self.final)
        self.assertEqual(result['Design'].tolist(), ['high', 'low'])
        self.assertEqual(result['Rank'].tolist(), [1, 2])
        self.assertEqual(sorted(p.name for p in self.ranked.iterdir()),
                         ['1_high_model1.pdb', '2_low_model1.pdb'])
        self.assertEqual((self.ranked / '1_high_model1.pdb').read_text(), 'high')
        self.assertEqual(self.target, original)
        self.assertTrue((self.accepted / 'low_model1.pdb').exists())
        self.archive.assert_not_called()

    def test_requested_metric_and_repeat_finalization_refresh_ranks(self):
        self.populate()
        self.call(finalize=True)
        self.call(finalize=True, rank_by='Average_ipSAE')
        self.assertEqual(pd.read_csv(self.final)['Design'].tolist(), ['low', 'high'])
        self.assertEqual(sorted(p.name for p in self.ranked.iterdir()),
                         ['1_low_model1.pdb', '2_high_model1.pdb'])

    def test_no_accepted_outputs_writes_header_without_reading_missing_mpnn_csv(self):
        (self.ranked / 'stale.pdb').write_text('stale')
        self.assertFalse(self.call(finalize=True))
        self.assertEqual(pd.read_csv(self.final).columns.tolist(), self.final_labels)
        self.assertTrue(pd.read_csv(self.final).empty)
        self.assertEqual(list(self.ranked.iterdir()), [])
        self.archive.assert_not_called()

    def test_unfinished_goal_check_does_not_rank_until_finalized(self):
        self.populate()
        self.assertFalse(self.call())
        self.assertFalse(self.final.exists())
        self.assertTrue((self.ranked / 'stale.pdb').exists())

    def test_completed_goal_still_ranks_and_archives(self):
        self.populate()
        self.target['number_of_final_designs'] = 2
        self.assertTrue(self.call())
        self.assertEqual(pd.read_csv(self.final)['Rank'].tolist(), [1, 2])
        self.assertEqual(self.archive.call_count, 2)
