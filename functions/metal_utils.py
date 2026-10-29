"""Small-matrix alignment operations for experimental JAX Metal backends.

Apple's jax-metal lacks eigh/SVD; metaljax otherwise uses host LAPACK for them.
Use cyclic Jacobi rotations for the 4x4 symmetric quaternion problem instead;
all arithmetic stays on the GPU. See apple_silicon/README.md for setup and usage.
The rest of ColabDesign, including the AlphaFold model and losses, is unchanged.
"""
import jax
import jax.numpy as jnp


def _jacobi4(matrix):
    """Diagonalize a batch of symmetric 4x4 matrices using eight sweeps."""
    eye = jnp.broadcast_to(jnp.eye(4, dtype=matrix.dtype), matrix.shape)

    def sweep(_, state):
        a, vectors = state
        for p, q in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
            off = a[..., p, q]
            theta = 0.5 * jnp.arctan2(2 * off, a[..., q, q] - a[..., p, p])
            theta = jnp.where(off == 0, 0, theta)
            c, s = jnp.cos(theta), jnp.sin(theta)
            rotation = eye.at[..., p, p].set(c).at[..., q, q].set(c)
            rotation = rotation.at[..., p, q].set(s).at[..., q, p].set(-s)
            a = rotation.swapaxes(-1, -2) @ a @ rotation
            vectors = vectors @ rotation
        return a, vectors

    diagonal, vectors = jax.lax.fori_loop(0, 8, sweep, (matrix, eye))
    return jnp.diagonal(diagonal, axis1=-2, axis2=-1), vectors


@jax.custom_jvp
def _largest_eigenvector(matrix):
    values, vectors = _jacobi4(matrix)
    index = jnp.argmax(values, axis=-1)
    return jnp.take_along_axis(vectors, index[..., None, None], axis=-1)[..., 0]


@_largest_eigenvector.defjvp
def _largest_eigenvector_jvp(primals, tangents):
    (matrix,), (tangent,) = primals, tangents
    values, vectors = _jacobi4(matrix)
    index = jnp.argmax(values, axis=-1)
    vector = jnp.take_along_axis(vectors, index[..., None, None], axis=-1)[..., 0]
    value = jnp.take_along_axis(values, index[..., None], axis=-1)
    gap = value - values
    # Analytic eigenvector derivative; do not differentiate the Jacobi sweeps
    # or truncate small, nonzero eigenvalue gaps. A repeated eigenvalue has no
    # unique derivative; take zero within its exactly degenerate eigenspace.
    valid = gap != 0
    inverse_gap = jnp.where(valid, 1 / jnp.where(valid, gap, 1), 0)
    projected = (vectors.swapaxes(-1, -2) @ tangent @ vector[..., None])[..., 0]
    derivative = (vectors @ (inverse_gap * projected)[..., None])[..., 0]
    return vector, derivative


def _quaternion_matrix(matrix):
    xx, xy, xz = matrix[..., 0, 0], matrix[..., 0, 1], matrix[..., 0, 2]
    yx, yy, yz = matrix[..., 1, 0], matrix[..., 1, 1], matrix[..., 1, 2]
    zx, zy, zz = matrix[..., 2, 0], matrix[..., 2, 1], matrix[..., 2, 2]
    rows = ((xx+yy+zz, zy-yz, xz-zx, yx-xy),
            (zy-yz, xx-yy-zz, xy+yx, xz+zx),
            (xz-zx, xy+yx, yy-xx-zz, yz+zy),
            (yx-xy, xz+zx, yz+zy, zz-xx-yy))
    return jnp.stack([jnp.stack(row, axis=-1) for row in rows], axis=-2)


def metal_rot_to_quat(rot, unstack_inputs=False):
    if not unstack_inputs:
        rot = jnp.stack([jnp.stack(row, axis=-1) for row in rot], axis=-2)
    return _largest_eigenvector(_quaternion_matrix(rot) / 3.)


def enable_metal_alignment():
    """Patch only ColabDesign's alignment helpers, only on the Metal backend."""
    from colabdesign.shared import protein
    from colabdesign.af import loss, utils
    from colabdesign.af.alphafold.model import quat_affine

    if getattr(protein._np_kabsch, '_freebindcraft_metal', False):
        return
    original = protein._np_kabsch

    def kabsch(a, b, return_v=False, use_jax=True):
        if not use_jax:
            return original(a, b, return_v=return_v, use_jax=False)
        if return_v:
            raise NotImplementedError('Metal alignment supports rotations, not SVD principal axes')
        covariance = a.swapaxes(-1, -2) @ b
        quat = _largest_eigenvector(_quaternion_matrix(covariance))
        rot = quat_affine.quat_to_rot(quat)
        return jnp.stack([jnp.stack(row, axis=-1) for row in rot], axis=-2)

    kabsch._freebindcraft_metal = True
    protein._np_kabsch = loss._np_kabsch = utils._np_kabsch = kabsch
    quat_affine.rot_to_quat = metal_rot_to_quat

    # Some Metal plugins silently implement decompositions through host LAPACK.
    # This port must fail rather than move a model operation back to the CPU.
    from jax._src.interpreters import mlir
    from jax._src.lax import linalg
    def gpu_only(ctx, *args, **kwargs):
        raise NotImplementedError('Host linear algebra is disabled for GPU-only FreeBindCraft')
    for name in ('eigh_p', 'svd_p', 'eig_p', 'qr_p', 'cholesky_p',
                 'triangular_solve_p', 'tridiagonal_solve_p', 'schur_p',
                 'hessenberg_p', 'tridiagonal_p', 'lu_p'):
        primitive = getattr(linalg, name, None)
        if primitive is not None:
            mlir.register_lowering(primitive, gpu_only, platform=jax.default_backend())
