"""Resume, reconciliation, and crash-safety tests without scientific dependencies."""
import importlib.util
import json
from pathlib import Path
import random
import sys
import tempfile
import unittest
from unittest.mock import patch
import warnings

path = Path(__file__).resolve().parents[1] / 'functions' / 'trajectory_scheduler.py'
spec = importlib.util.spec_from_file_location('trajectory_checkpoint_test_module', path)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
Checkpoint = module.TrajectoryScheduleCheckpoint


class CheckpointTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / '.trajectory-schedule.json'
        self.accepted = set()
        self.identity = 'campaign-one'

    def open(self, goal=20, attempts=None, low=65, high=200, identity=None, repeats=5, target=1):
        return Checkpoint(self.path, self.identity if identity is None else identity,
                          low, high, goal-len(self.accepted), attempts, self.accepted,
                          attempts_per_visit=repeats, accepted_per_visit=target)

    def attempt(self, checkpoint, number, outputs=0, remaining=20, attempts=None):
        length = checkpoint.scheduler.next_length(remaining, attempts)
        name = f'UBQ_l{length}_s{number}'
        self.assertTrue(checkpoint.begin_attempt(name, self.accepted))
        self.accepted.update(f'{name}_mpnn{i}_model1.pdb' for i in range(1, outputs+1))
        self.assertTrue(checkpoint.finish_attempt(self.accepted))
        return length

    def test_resume_after_success_block_preserves_original_regions(self):
        c = self.open()
        self.assertEqual(self.attempt(c, 0, 1, 20), 82)
        uninterrupted = c.scheduler.next_length(19)
        resumed = self.open()
        self.assertEqual(resumed.scheduler.next_length(19), uninterrupted)
        self.assertEqual(uninterrupted, 116)
        self.assertEqual([(b.low, b.high) for b in resumed.scheduler.bins],
                         [(65, 98), (99, 132), (133, 166), (167, 200)])

    def test_all_failure_blocks_resume_at_the_next_region(self):
        c = self.open()
        for i in range(5):
            self.assertEqual(self.attempt(c, i, remaining=20), 82)
        self.assertTrue(c.scheduler.current['done'])
        resumed = self.open(attempts=3)
        self.assertEqual(resumed.scheduler.next_length(20, 3), c.scheduler.next_length(20, 3))
        self.assertEqual(resumed.scheduler.current['length'], 116)
        self.assertEqual(resumed.scheduler.bins[0].attempted, 5)
        self.assertEqual(resumed.scheduler.bins[0].accepted, 0)
        self.assertNotIn('limit', resumed.scheduler.current)

    def test_mid_block_resume_retains_counters_and_length(self):
        c = self.open()
        self.attempt(c, 1, 0, 20)
        self.attempt(c, 2, 0, 20)
        self.attempt(c, 3, 0, 20)
        before = c.scheduler.to_state()
        resumed = self.open(attempts=1)
        self.assertEqual(resumed.scheduler.to_state(), before)
        self.assertEqual(resumed.scheduler.next_length(20, 1), 82)
        self.assertEqual(resumed.scheduler.current['attempts'], 3)
        self.assertEqual(resumed.scheduler.current['successful'], 0)
        self.assertEqual(resumed.scheduler.current['accepted'], 0)
        # A new invocation's one-attempt limit does not rewrite this visit's
        # historical limit; the driver counts its own new attempts separately.
        self.attempt(resumed, 4, 0, 20, 1)
        self.assertFalse(resumed.scheduler.current['done'])
        self.assertEqual(self.open().scheduler.current['attempts'], 4)

    def test_custom_limit_over_eight_survives_a_mid_visit_restart(self):
        c = self.open(repeats=10)
        for i in range(9):
            self.assertEqual(self.attempt(c, i), 82)
        resumed = self.open(repeats=10)
        self.assertEqual(resumed.scheduler.current['attempts'], 9)
        self.assertFalse(resumed.scheduler.current['done'])
        self.assertEqual(self.attempt(resumed, 9), 82)
        self.assertEqual(self.open(repeats=10).scheduler.next_length(20), 116)
        self.assertEqual(json.loads(self.path.read_text())['attempts_per_visit'], 10)

    def test_lower_limit_on_resume_finishes_an_over_limit_visit(self):
        c = self.open()
        for i in range(3):
            self.attempt(c, i)
        with self.assertWarnsRegex(RuntimeWarning, 'repeat limit changed from 5 to 2'):
            resumed = self.open(repeats=2)
        self.assertTrue(resumed.scheduler.current['done'])
        self.assertEqual(resumed.scheduler.bins[0].attempted, 3)
        self.assertEqual(self.open(repeats=2).scheduler.next_length(20), 116)

    def test_raised_limit_on_resume_continues_an_unsuccessful_visit(self):
        c = self.open(repeats=2)
        self.attempt(c, 0)
        self.attempt(c, 1)
        self.assertTrue(c.scheduler.current['done'])
        with self.assertWarnsRegex(RuntimeWarning, 'repeat limit changed from 2 to 4'):
            resumed = self.open(repeats=4)
        self.assertFalse(resumed.scheduler.current['done'])
        self.assertEqual(self.attempt(resumed, 2), 82)
        self.assertEqual(self.attempt(resumed, 3), 82)
        self.assertEqual(self.open(repeats=4).scheduler.next_length(20), 116)

    def test_raised_limit_does_not_reopen_a_successful_visit(self):
        c = self.open(repeats=2)
        self.attempt(c, 0, 1)
        with self.assertWarnsRegex(RuntimeWarning, 'repeat limit changed from 2 to 10'):
            resumed = self.open(repeats=10)
        self.assertTrue(resumed.scheduler.current['done'])
        self.assertEqual(resumed.scheduler.next_length(19), 116)

    def test_changed_limit_recovers_inflight_attempt_once(self):
        c = self.open()
        for i in range(3):
            self.attempt(c, i)
        self.assertTrue(c.begin_attempt('UBQ_l82_s99', self.accepted))
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter('always')
            resumed = self.open(repeats=2)
        self.assertEqual(len(captured), 2)
        self.assertEqual(resumed.scheduler.current['attempts'], 4)
        self.assertEqual(resumed.scheduler.next_length(20), 116)
        again = self.open(repeats=2)
        self.assertEqual(again.scheduler.current['attempts'], 4)
        self.assertEqual(again.scheduler.next_length(20), 116)

    def test_v2_checkpoint_migrates_and_accepts_a_custom_limit(self):
        c = self.open()
        self.attempt(c, 0)
        self.attempt(c, 1)
        state = json.loads(self.path.read_text())
        state['version'] = 2
        state.pop('attempts_per_visit')
        state.pop('accepted_per_visit')
        self.path.write_text(json.dumps(state))
        with self.assertWarnsRegex(RuntimeWarning, 'repeat limit changed from 5 to 8'):
            resumed = self.open(repeats=8)
        self.assertEqual(resumed.scheduler.current['attempts'], 2)
        self.assertEqual(resumed.scheduler.next_length(20), 82)
        saved = json.loads(self.path.read_text())
        self.assertEqual(saved['version'], 4)
        self.assertEqual(saved['attempts_per_visit'], 8)

    def test_v3_checkpoint_preserves_repeats_and_applies_new_target(self):
        c = self.open(repeats=9)
        self.attempt(c, 0, 1)
        state = json.loads(self.path.read_text())
        state['version'] = 3
        state.pop('accepted_per_visit')
        self.path.write_text(json.dumps(state))
        with self.assertWarnsRegex(RuntimeWarning, 'accepted-design target changed from 1 to 2'):
            resumed = self.open(repeats=9, target=2)
        self.assertEqual(resumed.scheduler.current['accepted'], 1)
        self.assertEqual(resumed.scheduler.next_length(19), 82)
        saved = json.loads(self.path.read_text())
        self.assertEqual((saved['version'], saved['attempts_per_visit'], saved['accepted_per_visit']), (4, 9, 2))

    def test_resume_retains_partial_accepted_design_target(self):
        c = self.open(target=3)
        self.attempt(c, 0, 1)
        self.attempt(c, 1, 1)
        resumed = self.open(target=3)
        self.assertFalse(resumed.scheduler.current['done'])
        self.assertEqual(resumed.scheduler.current['accepted'], 2)
        self.assertEqual(self.attempt(resumed, 2, 1, remaining=18), 82)
        self.assertEqual(self.open(target=3).scheduler.next_length(17), 116)

    def test_increasing_accepted_target_can_reopen_a_visit(self):
        c = self.open()
        self.attempt(c, 0, 1)
        self.assertTrue(c.scheduler.current['done'])
        with self.assertWarnsRegex(RuntimeWarning, 'accepted-design target changed from 1 to 3'):
            resumed = self.open(target=3)
        self.assertFalse(resumed.scheduler.current['done'])
        self.assertEqual(self.attempt(resumed, 1, 1, remaining=19), 82)
        self.assertEqual(self.attempt(resumed, 2, 1, remaining=18), 82)
        self.assertEqual(self.open(target=3).scheduler.next_length(17), 116)

    def test_decreasing_accepted_target_ends_a_visit(self):
        c = self.open(target=3)
        self.attempt(c, 0, 1)
        with self.assertWarnsRegex(RuntimeWarning, 'accepted-design target changed from 3 to 1'):
            resumed = self.open(target=1)
        self.assertTrue(resumed.scheduler.current['done'])
        self.assertEqual(self.open(target=1).scheduler.next_length(19), 116)

    def test_increased_target_does_not_override_exhausted_attempts(self):
        c = self.open(repeats=2)
        self.attempt(c, 0, 0)
        self.attempt(c, 1, 0)
        with self.assertWarnsRegex(RuntimeWarning, 'accepted-design target changed from 1 to 3'):
            resumed = self.open(repeats=2, target=3)
        self.assertTrue(resumed.scheduler.current['done'])
        self.assertEqual(resumed.scheduler.next_length(20), 116)

    def test_interrupted_attempt_counts_toward_accepted_target_once(self):
        c = self.open(target=2)
        self.attempt(c, 0, 1)
        name = 'UBQ_l82_s99'
        self.assertTrue(c.begin_attempt(name, self.accepted))
        self.accepted.add(name + '_mpnn1_model1.pdb')
        with self.assertWarnsRegex(RuntimeWarning, 'Recovered an interrupted'):
            resumed = self.open(target=2)
        self.assertEqual(resumed.scheduler.current['accepted'], 2)
        self.assertTrue(resumed.scheduler.current['done'])
        self.assertEqual(self.open(target=2).scheduler.next_length(18), 116)

    def test_invalid_saved_accepted_target_discards_the_checkpoint(self):
        c = self.open()
        self.attempt(c, 0)
        state = json.loads(self.path.read_text())
        state['accepted_per_visit'] = 0
        self.path.write_text(json.dumps(state))
        with self.assertWarnsRegex(RuntimeWarning, 'invalid saved accepted-design target'):
            resumed = self.open(target=3)
        self.assertIsNone(resumed.scheduler.current)
        self.assertEqual(resumed.scheduler.accepted_per_visit, 3)

    def test_invalid_saved_repeat_limit_discards_the_checkpoint(self):
        c = self.open()
        self.attempt(c, 0)
        state = json.loads(self.path.read_text())
        state['attempts_per_visit'] = 0
        self.path.write_text(json.dumps(state))
        with self.assertWarnsRegex(RuntimeWarning, 'invalid saved trajectory repeat limit'):
            resumed = self.open(repeats=3)
        self.assertIsNone(resumed.scheduler.current)
        self.assertEqual(resumed.scheduler.attempts_per_visit, 3)

    def test_interrupted_attempt_recovers_only_its_new_accepted_files(self):
        c = self.open()
        self.assertEqual(c.scheduler.next_length(20), 82)
        name = 'UBQ_l82_s1'
        self.assertTrue(c.begin_attempt(name, self.accepted))
        self.accepted.update({name+'_mpnn1_model2.pdb', name+'_mpnn2_model1.pdb'})
        with self.assertWarnsRegex(RuntimeWarning, 'Recovered an interrupted'):
            resumed = self.open()
        self.assertIsNone(resumed.inflight)
        self.assertEqual((resumed.scheduler.current['attempts'],
                          resumed.scheduler.current['successful'],
                          resumed.scheduler.current['accepted']), (1, 1, 2))
        self.assertEqual(resumed.scheduler.next_length(18), 116)
        self.assertIsNone(json.loads(self.path.read_text())['inflight'])
        # Reopening a second time must not count the recovered attempt twice.
        self.assertEqual(self.open().scheduler.current['attempts'], 1)

    def test_interrupted_attempt_without_output_counts_as_failure(self):
        c = self.open()
        c.scheduler.next_length(20)
        c.begin_attempt('UBQ_l82_s1', self.accepted)
        with self.assertWarnsRegex(RuntimeWarning, 'Recovered an interrupted'):
            resumed = self.open()
        self.assertEqual((resumed.scheduler.current['attempts'],
                          resumed.scheduler.current['successful'],
                          resumed.scheduler.current['accepted']), (1, 0, 0))

    def test_external_addition_or_deletion_discards_stale_history(self):
        c = self.open()
        self.attempt(c, 1, 1)
        self.accepted.clear()
        self.accepted.add('external_mpnn1_model1.pdb')  # same count, different file
        with self.assertWarnsRegex(RuntimeWarning, 'accepted files changed'):
            resumed = self.open()
        self.assertIsNone(resumed.scheduler.current)
        self.assertEqual(sum(b.attempted for b in resumed.scheduler.bins), 0)

    def test_unrelated_file_during_interrupted_attempt_is_not_attributed(self):
        c = self.open()
        c.scheduler.next_length(20)
        c.begin_attempt('UBQ_l82_s1', self.accepted)
        self.accepted.add('different_mpnn1_model1.pdb')
        with self.assertWarnsRegex(RuntimeWarning, 'accepted files changed'):
            resumed = self.open()
        self.assertIsNone(resumed.scheduler.current)
        self.assertEqual(sum(b.attempted for b in resumed.scheduler.bins), 0)

    def test_similar_but_invalid_accepted_filename_is_not_attributed(self):
        c = self.open()
        c.scheduler.next_length(20)
        c.begin_attempt('UBQ_l82_s1', self.accepted)
        self.accepted.add('UBQ_l82_s1_mpnn-not-a-design.pdb')
        with self.assertWarnsRegex(RuntimeWarning, 'accepted files changed'):
            resumed = self.open()
        self.assertIsNone(resumed.scheduler.current)

    def test_changed_goal_or_range_does_not_reuse_state(self):
        c = self.open()
        self.attempt(c, 1, 1)
        for opts in ({'identity':'different-goal', 'goal':25},
                     {'identity':'different-range', 'low':70, 'high':180}):
            with self.subTest(opts=opts), self.assertWarnsRegex(RuntimeWarning, 'campaign settings changed'):
                resumed = self.open(**opts)
                self.assertIsNone(resumed.scheduler.current)
                self.assertEqual(sum(b.attempted for b in resumed.scheduler.bins), 0)

    def test_identity_changes_for_goal_range_settings_and_target_contents(self):
        pdb = Path(self.tmp.name) / 'target.pdb'
        pdb.write_text('first')
        target = {'starting_pdb':str(pdb), 'lengths':[65,200], 'number_of_final_designs':20,
                  'design_path':'unused', 'chains':'A'}
        identity = module.scheduling_identity(target, {'weights':1}, {'threshold':2})
        for updated, advanced, filters in (
                (dict(target, number_of_final_designs=21), {'weights':1}, {'threshold':2}),
                (dict(target, lengths=[66,200]), {'weights':1}, {'threshold':2}),
                (target, {'weights':2}, {'threshold':2}),
                (target, {'weights':1}, {'threshold':3})):
            self.assertNotEqual(identity, module.scheduling_identity(updated, advanced, filters))
        self.assertEqual(identity, module.scheduling_identity(dict(target, design_path='moved'),
                                                               {'weights':1}, {'threshold':2}))
        pdb.write_text('other')
        self.assertNotEqual(identity, module.scheduling_identity(target, {'weights':1}, {'threshold':2}))

    def test_legacy_or_corrupt_state_starts_fresh_conservatively(self):
        self.accepted.add('old_mpnn1_model1.pdb')
        self.assertIsNone(self.open().scheduler.current)
        self.path.write_text('{ broken')
        with self.assertWarnsRegex(RuntimeWarning, 'Invalid trajectory'):
            self.assertIsNone(self.open().scheduler.current)
        self.path.unlink()
        c = self.open()
        self.attempt(c, 1, 0, 19)
        bad = json.loads(self.path.read_text())
        bad['scheduler']['bins'][0]['accepted'] = -1
        self.path.write_text(json.dumps(bad))
        with self.assertWarnsRegex(RuntimeWarning, 'Cannot reuse'):
            self.assertIsNone(self.open().scheduler.current)

    def test_atomic_replace_failure_leaves_last_checkpoint_intact(self):
        c = self.open()
        self.attempt(c, 1, 0)
        before = self.path.read_bytes()
        c.scheduler.next_length(20)
        with patch.object(module.os, 'replace', side_effect=OSError('disk error')):
            with self.assertRaisesRegex(OSError, 'disk error'):
                c.begin_attempt('UBQ_l82_s2', self.accepted)
        self.assertEqual(self.path.read_bytes(), before)
        self.assertEqual(list(self.path.parent.glob('.trajectory-schedule-*')), [])
        self.assertEqual(self.open().scheduler.current['attempts'], 1)

    def test_repeated_restarts_match_uninterrupted_mixed_outcomes(self):
        rng = random.Random(7712)
        direct = module.TrajectoryLengthScheduler(65, 200, 30)
        checkpoint = self.open(goal=30)
        remaining = 30
        with patch.object(module.os, 'fsync'):
            for attempt in range(100):
                if remaining == 0:
                    break
                length = direct.next_length(remaining)
                self.assertEqual(checkpoint.scheduler.next_length(remaining), length)
                outputs = min(remaining, rng.choice((0, 0, 0, 1, 2)))
                name = f'UBQ_l{length}_s{attempt}'
                self.assertTrue(checkpoint.begin_attempt(name, self.accepted))
                self.accepted.update(f'{name}_mpnn{i}_model1.pdb' for i in range(1, outputs+1))
                self.assertTrue(checkpoint.finish_attempt(self.accepted))
                direct.record(outputs)
                remaining -= outputs
                if rng.random() < 0.55:
                    checkpoint = self.open(goal=30)
                self.assertEqual(checkpoint.scheduler.to_state(), direct.to_state())

    def convert_to_v1(self, *, target, limit, done):
        state = json.loads(self.path.read_text())
        state['version'] = 1
        state.pop('attempts_per_visit')
        state.pop('accepted_per_visit')
        state['scheduler']['current'].update(target=target, limit=limit, done=done)
        self.path.write_text(json.dumps(state))

    def test_v1_exhausted_two_attempt_limit_resumes_to_five(self):
        c = self.open()
        self.attempt(c, 1, 0)
        self.attempt(c, 2, 0)
        self.convert_to_v1(target=4, limit=2, done=True)
        with self.assertWarnsRegex(RuntimeWarning, 'Migrated.*five-attempt'):
            resumed = self.open()
        self.assertEqual(resumed.scheduler.next_length(20), 82)
        self.assertEqual(resumed.scheduler.current['attempts'], 2)
        for number in (3,4,5):
            self.assertEqual(self.attempt(resumed, number, 0), 82)
        self.assertEqual(resumed.scheduler.next_length(20), 116)
        saved = json.loads(self.path.read_text())
        self.assertEqual(saved['version'], 4)
        self.assertNotIn('limit', saved['scheduler']['current'])

    def test_v1_partial_success_ends_visit_under_new_policy(self):
        c = self.open()
        self.attempt(c, 1, 1)
        self.convert_to_v1(target=5, limit=8, done=False)
        with self.assertWarnsRegex(RuntimeWarning, 'Migrated.*five-attempt'):
            resumed = self.open()
        self.assertTrue(resumed.scheduler.current['done'])
        self.assertEqual(resumed.scheduler.next_length(19), 116)
        self.assertEqual(resumed.scheduler.bins[0].accepted, 1)

    def test_v1_visit_over_five_is_historical_and_moves_on(self):
        c = self.open()
        for number in range(5):
            self.attempt(c, number, 0)
        state = json.loads(self.path.read_text())
        state['version'] = 1
        state['scheduler']['bins'][0]['attempted'] = 6
        state['scheduler']['current'].update(attempts=6, target=5, limit=8, done=False)
        self.path.write_text(json.dumps(state))
        with self.assertWarnsRegex(RuntimeWarning, 'Migrated.*five-attempt'):
            resumed = self.open()
        self.assertEqual(resumed.scheduler.current['attempts'], 6)
        self.assertEqual(resumed.scheduler.next_length(20), 116)
        self.assertEqual(self.open().scheduler.current['attempts'], 6)

    def test_v1_inflight_after_old_fifth_attempt_is_recovered_once(self):
        c = self.open()
        for number in range(4):
            self.attempt(c, number, 0)
        self.assertTrue(c.begin_attempt('UBQ_l82_s99', self.accepted))
        self.convert_to_v1(target=5, limit=8, done=False)
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter('always')
            resumed = self.open()
        self.assertEqual(len(captured), 2)
        self.assertEqual(resumed.scheduler.current['attempts'], 5)
        self.assertTrue(resumed.scheduler.current['done'])
        self.assertEqual(resumed.scheduler.next_length(20), 116)
        self.assertEqual(self.open().scheduler.current['attempts'], 5)


if __name__ == '__main__':
    unittest.main()
