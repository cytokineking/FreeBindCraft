"""Keep genuine DSSP diagnostics while suppressing a successful PDB probe."""
import importlib
import unittest
import warnings
from unittest.mock import patch

from functions.biopython_utils import _run_dssp

PROBE = 'parse error at line 1: This file does not seem to be an mmCIF file\n'
DSSP_SOURCE = importlib.import_module('Bio.PDB.DSSP').__file__


class DsspWarningTest(unittest.TestCase):
    def run_case(self, message, result, filename='input.pdb', fail=False):
        def calculate(*args, **kwargs):
            warnings.warn_explicit(message, UserWarning, DSSP_SOURCE, 199)
            if fail:
                raise RuntimeError('DSSP failed to produce an output')
            return result
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter('always')
            with patch('functions.biopython_utils.DSSP', side_effect=calculate):
                if fail:
                    with self.assertRaisesRegex(RuntimeError, 'failed to produce'):
                        _run_dssp(None, filename, 'mkdssp')
                else:
                    self.assertIs(_run_dssp(None, filename, 'mkdssp'), result)
        return [str(notice.message) for notice in captured]

    def test_suppresses_only_successful_pdb_probe(self):
        self.assertEqual(self.run_case(PROBE, {'residue': 'H'}), [])

    def test_preserves_failed_and_empty_results(self):
        self.assertEqual(self.run_case(PROBE, None, fail=True), [PROBE])
        self.assertEqual(self.run_case(PROBE, {}), [PROBE])

    def test_preserves_other_diagnostics_and_mmcif_errors(self):
        combined = PROBE + 'Another DSSP problem\n'
        self.assertEqual(self.run_case(combined, {'residue': 'H'}), [combined])
        self.assertEqual(self.run_case(PROBE, {'residue': 'H'}, 'input.cif'), [PROBE])


if __name__ == '__main__':
    unittest.main()
