"""Exercise the pinned ColabDesign renderer with modern Matplotlib and FFmpeg."""
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from colabdesign.shared.plot import make_animation
from functions.colabdesign_utils import save_design_animation


class DesignAnimationTest(unittest.TestCase):
    def test_real_renderer_without_legacy_colormap_api(self):
        xyz = np.array([[0., 0., 0.], [1., 2., 1.], [2., 3., 0.],
                        [3., 2., 2.], [4., 0., 1.], [5., 1., 0.]])
        seq = np.eye(20)[np.arange(6)][None, ...]
        def animate(**kwargs):
            return make_animation(seq=[seq, seq], xyz=[xyz, xyz + 0.1],
                                  plddt=[np.full(6, 0.9)] * 2, length=[3, 3],
                                  align_xyz=False, **kwargs)
        original = getattr(cm, 'get_cmap', None)
        if hasattr(cm, 'get_cmap'):
            del cm.get_cmap
        try:
            with TemporaryDirectory() as folder:
                output = Path(folder) / 'animation.html'
                self.assertTrue(save_design_animation(SimpleNamespace(animate=animate), output))
                self.assertIn('data:video/mp4;base64,', output.read_text())
                np.testing.assert_allclose(cm.get_cmap('viridis', 8)(0.5),
                                           plt.get_cmap('viridis', 8)(0.5))
                self.assertEqual(plt.get_fignums(), [])
        finally:
            if original is None:
                del cm.get_cmap
            else:
                cm.get_cmap = original

    def test_render_failure_does_not_abort_design(self):
        def fail(**kwargs):
            plt.figure()
            raise RuntimeError('encoder unavailable')
        with TemporaryDirectory() as folder:
            output = Path(folder) / 'animation.html'
            with self.assertWarnsRegex(RuntimeWarning, 'encoder unavailable'):
                self.assertFalse(save_design_animation(SimpleNamespace(animate=fail), output))
            self.assertFalse(output.exists())
            self.assertEqual(plt.get_fignums(), [])


if __name__ == '__main__':
    unittest.main()
