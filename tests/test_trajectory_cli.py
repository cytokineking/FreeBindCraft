"""Exercise the real CLI argument definitions without loading GPU libraries."""
import argparse
import ast
import contextlib
import io
from pathlib import Path
import sys
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
source = ast.parse((ROOT / 'bindcraft.py').read_text())
start = next(i for i, node in enumerate(source.body)
             if isinstance(node, ast.Assign) and any(
                 isinstance(target, ast.Name) and target.id == 'parser' for target in node.targets))
end = next(i for i, node in enumerate(source.body[start:], start)
           if isinstance(node, ast.If) and isinstance(node.test, ast.Compare)
           and isinstance(node.test.left, ast.Attribute) and node.test.left.attr == 'seed')
PARSER = compile(ast.fix_missing_locations(ast.Module(body=source.body[start:end], type_ignores=[])),
                 str(ROOT / 'bindcraft.py'), 'exec')


class TrajectoryCliTests(unittest.TestCase):
    def parse(self, *arguments):
        env = {'argparse': argparse, 'DEFAULT_TRAJECTORY_REPEATS': 5,
               'DEFAULT_TRAJECTORY_ACCEPTED_DESIGNS': 1}
        with patch.object(sys, 'argv', ['bindcraft.py', *arguments]):
            exec(PARSER, env)
        return env['args']

    def test_repeat_default_and_custom_value(self):
        self.assertEqual(self.parse().attempts_per_length, 5)
        self.assertEqual(self.parse('--trajectory-cache', '--attempts-per-length', '12').attempts_per_length, 12)

    def test_accepted_design_default_and_custom_value(self):
        self.assertEqual(self.parse().accepted_per_length, 1)
        args = self.parse('--trajectory-cache', '--attempts-per-length', '5',
                          '--accepted-per-length', '5')
        self.assertEqual((args.attempts_per_length, args.accepted_per_length), (5, 5))

    def test_repeat_rejects_zero_negative_and_noninteger_values(self):
        for value in ('0', '-2', '1.5'):
            with self.subTest(value=value), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as caught:
                    self.parse('--attempts-per-length', value)
                self.assertEqual(caught.exception.code, 2)

    def test_accepted_design_target_rejects_invalid_values(self):
        for value in ('0', '-2', '1.5'):
            with self.subTest(value=value), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as caught:
                    self.parse('--accepted-per-length', value)
                self.assertEqual(caught.exception.code, 2)


if __name__ == '__main__':
    unittest.main()
