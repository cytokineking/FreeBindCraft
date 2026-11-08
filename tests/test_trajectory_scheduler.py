"""Behavioral checks for exact-length visits, without JAX."""
import importlib.util
from pathlib import Path
import sys
import unittest

path = Path(__file__).resolve().parents[1] / 'functions' / 'trajectory_scheduler.py'
spec = importlib.util.spec_from_file_location('trajectory_scheduler', path)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
Scheduler = module.TrajectoryLengthScheduler


class SchedulerTests(unittest.TestCase):
    def test_five_failures_then_next_bin(self):
        s = Scheduler(65, 200, 20)
        for i in range(5):
            self.assertEqual(s.next_length(20), 82)
            s.record(0)
            self.assertEqual(s.current['done'], i == 4)
        self.assertEqual(s.next_length(20), 116)
        self.assertEqual(s.bins[0].attempted, 5)

    def test_success_at_first_middle_or_fifth_attempt_ends_visit(self):
        for success_at in (1, 3, 5):
            with self.subTest(success_at=success_at):
                s = Scheduler(65, 200, 20)
                for i in range(success_at):
                    self.assertEqual(s.next_length(20), 82)
                    s.record(3 if i == success_at - 1 else 0)
                    self.assertEqual(s.current['done'], i == success_at - 1)
                self.assertEqual((s.current['successful'], s.current['accepted']), (1, 3))
                self.assertEqual(s.next_length(17), 116)

    def test_global_attempt_limit_does_not_assign_bin_shares(self):
        s = Scheduler(65, 200, 20, max_attempts=2)
        self.assertEqual(len(s.bins), 4)
        self.assertEqual(s.next_length(20, 2), 82)
        s.record(0)
        self.assertEqual(s.next_length(20, 1), 82)
        s.record(0)
        self.assertFalse(s.current['done'])  # The driver stops at the global cap.
        with self.assertRaises(ValueError):
            s.next_length(20, 0)
        self.assertEqual(s.next_length(20, 1), 82)  # A later invocation resumes.

    def test_configurable_visit_limits_and_early_success(self):
        for limit in (1, 2, 5, 9, 17):
            with self.subTest(limit=limit):
                s = Scheduler(65, 200, 20, attempts_per_visit=limit)
                for i in range(limit):
                    self.assertEqual(s.next_length(20), 82)
                    s.record(0)
                    self.assertEqual(s.current['done'], i == limit - 1)
                self.assertEqual(s.next_length(20), 116)
                s.record(1)
                self.assertTrue(s.current['done'])
                self.assertEqual(s.next_length(19), 150)

    def test_invalid_visit_limits_are_rejected(self):
        for limit in (0, -1, 1.5, True, None):
            with self.subTest(limit=limit), self.assertRaises(ValueError):
                Scheduler(65, 200, 20, attempts_per_visit=limit)
            with self.subTest(accepted=limit), self.assertRaises(ValueError):
                Scheduler(65, 200, 20, accepted_per_visit=limit)

    def test_visit_continues_until_accepted_design_target(self):
        s = Scheduler(65, 200, 20, attempts_per_visit=5, accepted_per_visit=3)
        for i, outputs in enumerate((1, 0, 1, 1)):
            self.assertEqual(s.next_length(20), 82)
            s.record(outputs)
            self.assertEqual(s.current['done'], i == 3)
        self.assertEqual(s.current['accepted'], 3)
        self.assertEqual(s.current['successful'], 3)
        self.assertEqual(s.next_length(17), 116)

    def test_repeat_limit_can_equal_target_or_end_before_it(self):
        s = Scheduler(65, 200, 20, attempts_per_visit=3, accepted_per_visit=3)
        for i in range(3):
            self.assertEqual(s.next_length(20), 82)
            s.record(1)
            self.assertEqual(s.current['done'], i == 2)
        self.assertEqual(s.next_length(17), 116)

        s = Scheduler(65, 200, 20, attempts_per_visit=2, accepted_per_visit=3)
        s.next_length(20)
        s.record(1)
        self.assertFalse(s.current['done'])
        s.next_length(19)
        s.record(0)
        self.assertTrue(s.current['done'])
        self.assertEqual(s.next_length(19), 116)

    def test_multiple_outputs_from_one_trajectory_count_as_designs(self):
        s = Scheduler(65, 200, 20, attempts_per_visit=2, accepted_per_visit=3)
        s.next_length(20)
        s.record(2)
        self.assertFalse(s.current['done'])
        s.next_length(18)
        s.record(2)
        self.assertTrue(s.current['done'])
        self.assertEqual((s.current['successful'], s.current['accepted']), (2, 4))
        self.assertEqual(s.next_length(16), 116)

    def test_later_visits_cover_new_parts_and_repeat_after_exhaustion(self):
        s = Scheduler(1, 8, 4)
        sweeps = []
        for _ in range(3):
            sweep = []
            for remaining in (4, 3, 2, 1):
                sweep.append(s.next_length(remaining))
                s.record(1)
            sweeps.append(sweep)
        self.assertEqual(sweeps, [[2,4,6,8], [1,3,5,7], [2,4,6,8]])

    def test_single_length_and_bad_calls(self):
        s = Scheduler(70, 70, 5)
        for remaining in (5,4,3,2,1):
            self.assertEqual(s.next_length(remaining), 70)
            s.record(1)
        self.assertEqual(s.bins[0].accepted, 5)
        with self.assertRaises(ValueError):
            s.record(0)
        with self.assertRaises(ValueError):
            s.next_length(0)
        for args in ((3,2,1,None), (1,2,0,None), (1,2,1,0)):
            with self.assertRaises(ValueError):
                Scheduler(*args)
        s.next_length(1)
        for invalid in (-1, 1.5, True):
            with self.assertRaises(ValueError):
                s.record(invalid)


if __name__ == '__main__':
    unittest.main()
