"""GPU regression checks for the Metal alignment implementation.

These small numerical tests cover degenerate/reflected/batched structures and
derivatives that one end-to-end design cannot reliably exercise. Run with:
JAX_PLATFORMS=metal .venv-metaljax/bin/python -m unittest discover -s tests
"""
import unittest
import numpy as np
import jax
import jax.numpy as jnp
from functions.metal_utils import enable_metal_alignment, metal_rot_to_quat
from colabdesign.shared import protein
from colabdesign.af.alphafold.model import quat_affine


@unittest.skipUnless(jax.default_backend().lower() == 'metal', 'requires Metal GPU')
class MetalAlignmentTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        enable_metal_alignment()

    def test_alignment_matches_optimal_svd(self):
        rng = np.random.default_rng(42)
        a = rng.normal(size=(20, 30, 3)).astype(np.float32)
        b = rng.normal(size=a.shape).astype(np.float32)
        # Rank-deficient, reflected, and isotropic inputs exercise ambiguous axes.
        a[0] = 0
        a[1, :, 2] = 0
        a[2, :, 1:] = 0
        b[3] = a[3] * np.array([-1, 1, 1])
        a[4] = b[4] = np.tile(np.eye(3), (10, 1))
        got = np.asarray(jax.jit(protein._np_kabsch)(jnp.array(a), jnp.array(b)))
        u, _, vh = np.linalg.svd(a.swapaxes(-1, -2) @ b)
        u[:, :, -1] *= np.linalg.det(u @ vh)[:, None]
        reference = u @ vh
        np.testing.assert_allclose(got.swapaxes(-1, -2) @ got,
                                   np.broadcast_to(np.eye(3), got.shape), atol=3e-5)
        np.testing.assert_allclose(np.linalg.det(got), 1, atol=3e-5)
        np.testing.assert_allclose(np.sum((a @ got-b)**2, axis=(1,2)),
                                   np.sum((a @ reference-b)**2, axis=(1,2)), rtol=2e-5, atol=2e-5)

    def test_rotation_quaternion_roundtrip(self):
        rng = np.random.default_rng(7)
        q = rng.normal(size=(40, 4)).astype(np.float32)
        q[:4] = np.eye(4)
        q /= np.linalg.norm(q, axis=-1, keepdims=True)
        def roundtrip(q):
            rot = quat_affine.quat_to_rot(q)
            return metal_rot_to_quat(rot)
        got = np.asarray(jax.jit(roundtrip)(jnp.array(q)))
        np.testing.assert_allclose(np.abs(np.sum(got*q, axis=-1)), 1., atol=2e-5)

    def test_alignment_gradient_matches_finite_difference(self):
        rng = np.random.default_rng(3)
        a = rng.normal(size=(20,3)).astype(np.float32)
        b = rng.normal(size=a.shape).astype(np.float32)
        direction = rng.normal(size=a.shape).astype(np.float32)
        weight = rng.normal(size=(3,3)).astype(np.float32)
        def loss(x):
            return jnp.sum(protein._np_kabsch(x, jnp.array(b)) * jnp.array(weight))
        grad = np.asarray(jax.jit(jax.grad(loss))(jnp.array(a)))
        def reference(x):
            u, _, vh = np.linalg.svd(x.T @ b)
            u[:, -1] *= np.linalg.det(u @ vh)
            return np.sum((u @ vh)*weight)
        eps = 0.002
        finite_difference = (reference(a+eps*direction)-reference(a-eps*direction))/(2*eps)
        self.assertTrue(np.isfinite(grad).all())
        np.testing.assert_allclose(np.sum(grad*direction), finite_difference, atol=2e-3, rtol=2e-3)

    def test_frame_conversion_gradient(self):
        # Identity and 180-degree rotations have degenerate *lower* eigenvalues.
        # The selected eigenvector still has a unique derivative in these cases.
        rng = np.random.default_rng(31)
        rotations = np.stack([np.eye(3), np.diag([1., -1., -1.])]).astype(np.float32)
        directions = rng.normal(size=rotations.shape).astype(np.float32)
        weight = rng.normal(size=(3,3)).astype(np.float32)
        def loss(matrices):
            q = metal_rot_to_quat(matrices, unstack_inputs=True)
            rot = quat_affine.quat_to_rot(q)
            rot = jnp.stack([jnp.stack(row, axis=-1) for row in rot], axis=-2)
            return jnp.sum(rot * weight)
        grad = np.asarray(jax.jit(jax.grad(loss))(jnp.array(rotations)))
        # A central difference checks the derivative of the full frame conversion.
        eps = 0.002
        fd = (float(loss(rotations+eps*directions)) - float(loss(rotations-eps*directions)))/(2*eps)
        np.testing.assert_allclose(np.sum(grad*directions), fd, atol=2e-3, rtol=2e-3)

    def test_host_decomposition_is_rejected(self):
        with self.assertRaisesRegex(NotImplementedError, 'Host linear algebra'):
            jax.jit(jnp.linalg.eigh)(jnp.eye(3))


if __name__ == '__main__':
    unittest.main()
