"""Dependency selection and fail-closed checks, without installing packages."""
import hashlib
import importlib.util
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


spec = importlib.util.spec_from_file_location(
    'apple_dependencies', Path(__file__).resolve().parents[1]/'installation/install_apple_python.py')
installer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(installer)


class AppleDependencyInstallTests(unittest.TestCase):
    def test_verified_default_cache_does_not_require_publication(self):
        with tempfile.TemporaryDirectory() as temporary:
            cache = Path(temporary)
            wheel = cache/'candidate.whl'
            wheel.write_bytes(b'verified unpublished candidate')
            release = {'wheel': {'url': 'https://example.invalid/future/candidate.whl',
                                'sha256': hashlib.sha256(wheel.read_bytes()).hexdigest()}}
            with patch.dict(os.environ, {}, clear=True), patch.object(installer.urllib.request, 'urlretrieve') as network:
                mode, selected, digest = installer.select_wheel(release, cache)
            network.assert_not_called()
            self.assertEqual((mode, selected, digest), ('frozen', wheel, release['wheel']['sha256']))

    def test_checksum_failure_preserves_cache_and_removes_partial(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root/'source/backend.whl'
            source.parent.mkdir()
            source.write_bytes(b'bad download')
            cache = root/'cache'
            cache.mkdir()
            (cache/source.name).write_bytes(b'previous cached bytes')
            asset = {'url': source.as_uri(), 'sha256': '0'*64}
            with self.assertRaisesRegex(RuntimeError, 'checksum mismatch'):
                installer.download_verified(asset, cache)
            self.assertEqual((cache/source.name).read_bytes(), b'previous cached bytes')
            self.assertEqual(list(cache.iterdir()), [cache/source.name])

    def test_corrupt_cache_is_replaced_by_verified_download(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root/'source/backend.whl'
            source.parent.mkdir()
            source.write_bytes(b'correct download')
            cache = root/'cache'
            cache.mkdir()
            (cache/source.name).write_bytes(b'corrupted cache')
            asset = {'url': source.as_uri(), 'sha256': hashlib.sha256(source.read_bytes()).hexdigest()}
            wheel = installer.download_verified(asset, cache)
            self.assertEqual(wheel.read_bytes(), source.read_bytes())

    def test_local_override_takes_precedence_over_stock_without_download(self):
        with tempfile.TemporaryDirectory() as temporary:
            wheel = Path(temporary)/'override.whl'
            wheel.write_bytes(b'local developer wheel')
            with patch.dict(os.environ, {'FREEBINDCRAFT_METALJAX': 'stock',
                                        'FREEBINDCRAFT_METALJAX_WHEEL': str(wheel)}, clear=True):
                with patch.object(installer, 'download_verified') as download:
                    mode, selected, digest = installer.select_wheel({}, Path(temporary))
                    download.assert_not_called()
            self.assertEqual((mode, selected), ('local override', wheel))
            self.assertEqual(digest, hashlib.sha256(wheel.read_bytes()).hexdigest())

    def test_default_and_stock_choose_their_exact_assets(self):
        release = {'wheel': {'sha256': 'frozen'}, 'stock_wheel': {'sha256': 'stock'}}
        for mode, field in [('frozen', 'wheel'), ('stock', 'stock_wheel')]:
            with self.subTest(mode=mode), patch.dict(os.environ, {} if mode == 'frozen' else
                    {'FREEBINDCRAFT_METALJAX': mode}, clear=True):
                with patch.object(installer, 'download_verified', return_value=Path('/wheel')) as download:
                    result = installer.select_wheel(release, Path('/cache'))
                    download.assert_called_once_with(release[field], Path('/cache'))
                    self.assertEqual(result, (mode, Path('/wheel'), release[field]['sha256']))

    def test_invalid_mode_is_rejected(self):
        with patch.dict(os.environ, {'FREEBINDCRAFT_METALJAX': 'typo'}, clear=True):
            with self.assertRaisesRegex(ValueError, 'must be frozen or stock'):
                installer.select_wheel({}, Path('/unused'))


if __name__ == '__main__':
    unittest.main()
