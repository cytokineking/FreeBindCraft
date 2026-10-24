"""Exercise the actual Bash launcher through isolated Python entry points."""
import json
import os
from pathlib import Path
import shutil
import shlex
import subprocess
import sys
import tempfile
import unittest


LAUNCHER = Path(__file__).resolve().parents[1]/'run_apple_silicon.sh'


class AppleLauncherTests(unittest.TestCase):
    def launch(self, args, inherited=None, pyrosetta=False):
        with tempfile.TemporaryDirectory(prefix='apple launcher ') as temporary:
            root = Path(temporary)
            shutil.copyfile(LAUNCHER, root/'run_apple_silicon.sh')
            python = root/'.venv-metaljax/bin/python'
            python.parent.mkdir(parents=True)
            # Observe the environment at every Python entry, then execute real
            # Python without site packages. The probe still uses find_spec;
            # no installed PyRosetta/GPU/model assets are needed for this test.
            recorder = root/'record entry.py'
            recorder.write_text(
                'import json, os, sys\n'+
                'with open("python-entry.jsonl", "a") as f:\n'+
                ' f.write(json.dumps({"attention":os.getenv("METALJAX_EVO_ATTN"),"argv":sys.argv[1:]})+"\\n")\n'+
                f'os.execv({str(sys.executable)!r}, [{str(sys.executable)!r}, "-S", *sys.argv[1:]])\n')
            python.write_text('#!/bin/sh\nexec '+shlex.quote(sys.executable)+
                              ' -S '+shlex.quote(str(recorder))+' "$@"\n')
            python.chmod(0o755)
            if pyrosetta:
                (root/'pyrosetta.py').write_text('# Probe-visible fixture module.\n')
            (root/'bindcraft.py').write_text(
                'import json,os,sys\n'
                'with open("result.json","w") as f:\n'
                ' json.dump({"argv":sys.argv[1:],"attention":os.getenv("METALJAX_EVO_ATTN"),'
                '"slice":os.getenv("METALJAX_SLICE_BEFORE_CONVERT")},f)\n')
            env = {k:v for k,v in os.environ.items() if not k.startswith(
                ('METALJAX_', 'MLX_', 'JAX_', 'PYTHON'))}
            if inherited is not None:
                env['METALJAX_EVO_ATTN'] = inherited
            completed = subprocess.run(['/bin/bash', str(root/'run_apple_silicon.sh'), *args],
                                       cwd=root.parent, env=env, text=True, capture_output=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            return (json.loads((root/'result.json').read_text()),
                    [json.loads(x) for x in (root/'python-entry.jsonl').read_text().splitlines()],
                    completed.stdout)

    def test_defaults_opt_out_precedence_and_forwarding(self):
        forwarded = ['--settings', 'path with spaces/target.json', '--seed', '42']
        cases = [
            # Arguments, inherited setting, PyRosetta availability, expected setting.
            ([], None, True, None),
            (forwarded, None, False, None),
            (['--no-pyrosetta', *forwarded], '0', True, '0'),
            ([*forwarded[:2], '--no-attention-epilogues', *forwarded[2:]], '1', True, '0'),
            (forwarded, '1', False, '1'),
            (['--no-attention-epilogues'], None, False, '0'),
        ]
        for args, inherited, pyrosetta, expected in cases:
            with self.subTest(args=args, inherited=inherited, pyrosetta=pyrosetta):
                result, entries, _ = self.launch(args, inherited, pyrosetta)
                wanted = [a for a in args if a != '--no-attention-epilogues']
                if not pyrosetta:
                    wanted.insert(0, '--no-pyrosetta')
                self.assertEqual(result['argv'], wanted)
                self.assertEqual(result['attention'], expected)
                self.assertIsNone(result['slice'])
                # Includes the PyRosetta probe, before bindcraft.py can import a backend.
                self.assertTrue(entries)
                self.assertTrue(all(e['attention'] == expected for e in entries))

    def test_help_exposes_launcher_option_and_preserves_help(self):
        result, _, output = self.launch(['--help'], pyrosetta=True)
        self.assertIn('--no-attention-epilogues', output)
        self.assertEqual(result['argv'], ['--help'])


if __name__ == '__main__':
    unittest.main()
