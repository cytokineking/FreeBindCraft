"""Verify staged, repeatable patching without touching the installed package."""
import difflib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

path = Path(__file__).resolve().parents[1] / 'installation/apply_colabdesign_metal.py'
spec = importlib.util.spec_from_file_location('apply_colabdesign_metal', path)
installer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(installer)


class InstallerTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.assets = self.root / 'assets'
        self.assets.mkdir()
        self.package = self.root / 'package'
        (self.package / 'af').mkdir(parents=True)
        self.prep = self.package / 'af/prep.py'
        self.prep.write_text('original\n')
        self.helper = self.package / 'af/metal_execution_policy.py'
        payload = self.assets / 'metal_execution_policy.py'
        payload.write_text('helper\n')
        self.patch = self.assets / 'unchunked-preparation.patch'
        self.patch.write_text(''.join(difflib.unified_diff(
            ['original\n'], ['unchunked\n'],
            fromfile='a/colabdesign/af/prep.py', tofile='b/colabdesign/af/prep.py')))
        self.manifest = {
            'patch': self.patch.name, 'patch_sha256': installer.sha(self.patch),
            'files': {
                'af/prep.py': {'before': installer.sha(self.prep),
                              'after': installer.hashlib.sha256(b'unchunked\n').hexdigest()},
                'af/metal_execution_policy.py': {'before': None, 'after': installer.sha(payload),
                                                 'payload': payload.name},
            },
        }
        self.save_manifest()

    def save_manifest(self):
        (self.assets / 'manifest.json').write_text(json.dumps(self.manifest))

    def test_fresh_install_and_repeat(self):
        self.assertEqual(installer.apply(self.package, self.assets)['status'], 'applied')
        self.assertEqual(self.prep.read_text(), 'unchunked\n')
        self.assertEqual(self.helper.read_text(), 'helper\n')
        self.assertEqual(installer.apply(self.package, self.assets)['status'], 'already applied')

    def test_local_edits_are_preserved(self):
        self.prep.write_text('user edit\n')
        with self.assertRaisesRegex(RuntimeError, 'exact pinned baseline'):
            installer.apply(self.package, self.assets)
        self.assertEqual(self.prep.read_text(), 'user edit\n')
        self.assertFalse(self.helper.exists())

    def test_corrupt_patch_is_rejected_before_writing(self):
        self.patch.write_text(self.patch.read_text() + 'corrupt')
        with self.assertRaisesRegex(RuntimeError, 'patch checksum mismatch'):
            installer.apply(self.package, self.assets)
        self.assertEqual(self.prep.read_text(), 'original\n')
        self.assertFalse(self.helper.exists())

    def test_corrupt_helper_is_rejected_before_writing(self):
        (self.assets / 'metal_execution_policy.py').write_text('corrupt')
        with self.assertRaisesRegex(RuntimeError, 'helper checksum mismatch'):
            installer.apply(self.package, self.assets)
        self.assertEqual(self.prep.read_text(), 'original\n')
        self.assertFalse(self.helper.exists())

    def test_wrong_final_hash_leaves_original_package_intact(self):
        self.manifest['files']['af/prep.py']['after'] = '0' * 64
        self.save_manifest()
        with self.assertRaisesRegex(RuntimeError, 'Patched ColabDesign checksum mismatch'):
            installer.apply(self.package, self.assets)
        self.assertEqual(self.prep.read_text(), 'original\n')
        self.assertFalse(self.helper.exists())

    def test_modified_installed_helper_is_not_overwritten(self):
        installer.apply(self.package, self.assets)
        self.helper.write_text('user helper edit\n')
        with self.assertRaisesRegex(RuntimeError, 'exact pinned baseline'):
            installer.apply(self.package, self.assets)
        self.assertEqual(self.helper.read_text(), 'user helper edit\n')


if __name__ == '__main__':
    unittest.main()
