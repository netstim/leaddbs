# -*- coding: utf-8 -*-
"""
Tensor scaling algorithms by A. Andree, parallelization by K. Butenko.
Refactored for Windows compatibility (spawn-based multiprocessing).
"""

import os
import sys
import ctypes
import itertools
from functools import partial
from multiprocessing import RawArray, Pool, cpu_count
from typing import Optional

import numpy as np
import nibabel as nib

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
_EPS = 1e-12

# Supported tensor storage orders and the corresponding (row, col) index
# mapping for the six unique components stored as a flat 6-element vector.
#
#   NIFTI         : [xx, yx, yy, zx, zy, zz]  →  indices 0–5
#   DSI_studio    : [xx, yy, zz, yx, zx, zy]  →  indices 0–5
#   FSL           : [xx, yx, zx, yy, zy, zz]  →  indices 0–5
#   Johnson_Wistar: [yy, yx, zy, xx, zx, zz]  →  indices 0–5
#
# Each entry maps (flat_index) → position in the symmetric 3×3 matrix.
_TENSOR_LAYOUTS = {
    "NIFTI": [
        (0, 0), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2)
    ],
    "DSI_studio": [
        (0, 0), (1, 1), (2, 2), (1, 0), (2, 0), (2, 1)
    ],
    "FSL": [
        (0, 0), (1, 0), (2, 0), (1, 1), (2, 1), (2, 2)
    ],
    "Johnson_Wistar": [
        (1, 1), (1, 0), (2, 1), (0, 0), (2, 0), (2, 2)
    ],
}

# Lower conductivity boundary: WM value from Gabriel et al. at 10 Hz [S/m]
_SIGMA_ISO_LOW = 0.027512
_SIGMA_ISO_LOWER_BOUNDARY = _SIGMA_ISO_LOW / 20.0


# ---------------------------------------------------------------------------
# Worker-process state (set once per process via Pool initializer)
# ---------------------------------------------------------------------------
_worker_shared_out: Optional[np.ndarray] = None  # shaped output array (mx, my, mz, 6)
_worker_dti_flat: Optional[np.ndarray] = None    # shaped read-only DTI array
_worker_dti_shape: Optional[tuple] = None


def _init_worker(
    shared_out: RawArray,
    out_shape: tuple,
    dti_flat: RawArray,
    dti_shape: tuple,
) -> None:
    """Pool initializer — runs once in each spawned worker process."""
    global _worker_shared_out, _worker_dti_flat, _worker_dti_shape
    # Both frombuffer calls produce shaped, multi-dimensional views of the
    # underlying shared memory — no copy is made.
    _worker_shared_out = np.frombuffer(shared_out, dtype=np.float64).reshape(out_shape)
    _worker_dti_flat   = np.frombuffer(dti_flat,   dtype=np.float64).reshape(dti_shape)
    _worker_dti_shape  = dti_shape


# ---------------------------------------------------------------------------
# Scaling helpers
# ---------------------------------------------------------------------------
def _theta_star(w12: float, w13: float) -> float:
    """
    Nonlinear analytic mapping for eigenvalue scaling (load preservation).
    Based on Howell & McIntyre, 2016.
    """
    v1, u1, m = 2.15, 1.21, 8e-1
    v2, u2, n = 1.85, 1.12, 8e-1

    theta = (
        (v1 / (np.power(u1 / (w12 + _EPS), m) + 1))
        * (v2 / (np.power(u2 / (w13 + _EPS), n) + 1))
    )
    return np.round(theta, -int(np.log10(_EPS)))  # suppress round-off


def _flat_to_matrix(components: np.ndarray, tensor_order: str) -> np.ndarray:
    """Reconstruct the symmetric 3×3 diffusion tensor from six stored values."""
    layout = _TENSOR_LAYOUTS[tensor_order]
    mat = np.zeros((3, 3), dtype=float)
    for value, (row, col) in zip(components, layout):
        mat[row, col] = value
        mat[col, row] = value  # enforce symmetry
    return mat


def _scale_eigenvalues(
    eigvals: np.ndarray,
    matrix: np.ndarray,
    scaling_method: str,
    diagonal_sum: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return (scaled_tensor, eigvals_of_scaled_tensor).

    For the 'Nordin' method the tensor is returned directly (no eigenvalue
    decomposition needed); for all others a new tensor is reconstructed from
    the scaled eigenvalues and the original eigenvectors.
    """
    eigvecs = np.linalg.eig(matrix)[1]  # eigvecs already computed upstream

    if scaling_method == "Tuch":
        k_dti = 0.844e3      # S·s/mm³  (converts to S/m · s/mm²)
        d_eps = 0.124e-3     # mm²/s  (extracellular diffusivity)
        scaled = k_dti * (eigvals - d_eps)

    elif scaling_method == "NormMapping":
        det_cbrt = (eigvals[0] * eigvals[1] * eigvals[2]) ** (1.0 / 3.0)
        scaled = eigvals / det_cbrt

    elif scaling_method == "LoadPreservation":
        w12 = eigvals[0] / (eigvals[1] + _EPS)
        w13 = eigvals[0] / (eigvals[2] + _EPS)   # ← was eigvals[1] (bug fix)
        theta = _theta_star(w12, w13)
        scaled = np.array([1.0, 1.0 / (w12 + _EPS), 1.0 / (w13 + _EPS)]) * theta

    elif scaling_method == "Nordin":
        # Direct tensor scaling — no eigenvalue manipulation
        tensor = matrix / (diagonal_sum / 3.0)
        return tensor, np.linalg.eigvalsh(tensor)

    else:
        raise ValueError(f"Unknown scaling method: {scaling_method!r}")

    if np.any(scaled <= 0):
        raise RuntimeError("Scaled eigenvalues must be strictly positive.")

    tensor = eigvecs @ np.diag(scaled) @ eigvecs.T
    return tensor, np.linalg.eigvalsh(tensor)


# ---------------------------------------------------------------------------
# Per-voxel processing (one (i, j) column, all k)
# ---------------------------------------------------------------------------
def _process_column(
    k_vector: np.ndarray,
    tensor_order: str,
    scaling_method: str,
    affine: np.ndarray,
    affine_inv: np.ndarray,
    ij: tuple[int, int],
) -> None:
    """
    Process all voxels in column (i, j) across z.  Writes results into the
    shared output array.  Designed to be called by a multiprocessing Pool.
    """
    i, j = ij
    out = _worker_shared_out   # already a shaped (mx, my, mz, 6) ndarray
    dti = _worker_dti_flat

    for k in k_vector:
        # ---- extract the six tensor components --------------------------------
        if dti.ndim == 5:
            components = dti[i, j, k, 0, :]
        else:
            components = dti[i, j, k, :]

        # Background / zero voxel → identity-like placeholder
        if np.all(components == 0.0):
            out[i, j, k, :] = [1.0, 0.0, 0.0, 1.0, 0.0, 1.0]
            continue

        # ---- build symmetric 3×3 matrix --------------------------------------
        mat = _flat_to_matrix(components, tensor_order)

        eigvals, _ = np.linalg.eig(mat)
        if np.any(eigvals < 0):
            print(
                "Warning: negative eigenvalue encountered; "
                "taking absolute value.  Check your DTI data!"
            )
            eigvals = np.abs(eigvals)

        # ---- scale -----------------------------------------------------------
        diag_sum = float(components[0] + components[2] + components[5])
        tensor, eigvals_tensor = _scale_eigenvalues(
            eigvals, mat, scaling_method, diag_sum
        )

        if np.any(eigvals_tensor <= 0):
            raise RuntimeError(
                f"Non-positive eigenvalue after scaling at voxel ({i},{j},{k})."
            )

        if np.any(eigvals_tensor * _SIGMA_ISO_LOW < _SIGMA_ISO_LOWER_BOUNDARY):
            raise RuntimeError(
                f"Conductivity below lower boundary at voxel ({i},{j},{k})."
            )

        # ---- re-orient to world axes -----------------------------------------
        # Build 4×4 homogeneous rotation block from the 3×3 symmetric tensor
        tc = np.zeros((4, 4), dtype=float)
        tc[:3, :3] = tensor
        tc[3, 3] = 1.0
        tensor_wa = affine @ tc @ affine_inv

        # Store as FSL order: [xx, yx, zx, yy, zy, zz]
        out[i, j, k, :] = [
            tensor_wa[0, 0],
            tensor_wa[1, 0],
            tensor_wa[2, 0],
            tensor_wa[1, 1],
            tensor_wa[2, 1],
            tensor_wa[2, 2],
        ]


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
def _run_parallel(
    dti_data: np.ndarray,
    tensor_order: str,
    scaling_method: str,
    affine: np.ndarray,
) -> np.ndarray:
    """
    Distribute voxel processing across CPU cores.

    Uses shared memory for both the read-only DTI input and the write-only
    conductivity output so that no large array copies are made per task.
    The Pool initializer makes these shared buffers available as module-level
    globals inside each worker — the only pattern that works under Windows'
    spawn-based process start method.
    """
    if np.any(np.diag(affine[:3, :3]) < 0.0):
        raise ValueError(
            "Negative diagonal in affine matrix — "
            "reverse-ordered axes are not supported.  Please flip the image."
        )

    affine_inv = np.linalg.inv(affine)
    mx, my, mz = dti_data.shape[:3]

    # ---- shared output buffer (zeros → workers write into this) --------------
    out_np = np.zeros((mx, my, mz, 6), dtype=np.float64)
    out_shape = out_np.shape
    shared_out = RawArray(ctypes.c_double, out_np.size)

    # ---- shared read-only input buffer ---------------------------------------
    dti_flat_np = np.ascontiguousarray(dti_data, dtype=np.float64)
    shared_in = RawArray(ctypes.c_double, dti_flat_np.size)
    # copy DTI data into the shared buffer once
    np.frombuffer(shared_in, dtype=np.float64)[:] = dti_flat_np.flatten()

    dti_shape = dti_flat_np.shape

    # ---- iterate over (i, j) columns, parallelised ---------------------------
    k_vector = np.arange(mz)
    ij_pairs = list(itertools.product(range(mx), range(my)))

    n_workers = max(1, cpu_count() - 1)
    with Pool(
        processes=n_workers,
        initializer=_init_worker,
        initargs=(shared_out, out_shape, shared_in, dti_shape),
    ) as pool:
        pool.map(
            partial(
                _process_column,
                k_vector,
                tensor_order,
                scaling_method,
                affine,
                affine_inv,
            ),
            ij_pairs,
        )

    return np.frombuffer(shared_out, dtype=np.float64).reshape(out_shape)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
def scale_tensor_data(
    tensor_data_name: str,
    scaling_method: str = "NormMapping",
    tensor_order: str = "NIFTI",
) -> None:
    """
    Load a NIfTI diffusion tensor image, scale it to a conductivity tensor,
    and save the result alongside the original file.

    Parameters
    ----------
    tensor_data_name:
        Path to the input NIfTI file.
    scaling_method:
        One of ``'Tuch'``, ``'NormMapping'``, ``'LoadPreservation'``,
        ``'Nordin'``.
    tensor_order:
        Storage convention of the six tensor components:
        ``'NIFTI'``, ``'DSI_studio'``, ``'FSL'``, or ``'Johnson_Wistar'``.
    """
    if tensor_order not in _TENSOR_LAYOUTS:
        raise ValueError(
            f"Unknown tensor_order {tensor_order!r}. "
            f"Choose from {list(_TENSOR_LAYOUTS)}."
        )

    filepath = os.path.realpath(tensor_data_name)
    img = nib.load(filepath)
    dti_data = img.get_fdata()

    if np.any(np.isnan(dti_data)):
        raise ValueError("NaN values detected in the DTI data — please remove them.")

    normalized = _run_parallel(dti_data, tensor_order, scaling_method, img.affine)

    stem = nib.filename_parser.splitext_addext(filepath)[0]
    out_path = f"{stem}_{scaling_method}.nii.gz"
    nib.save(nib.Nifti1Image(normalized, img.affine), out_path)
    print(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # The `if __name__ == "__main__"` guard is required on Windows so that
    # spawned worker processes do not re-execute the top-level script.
    scale_tensor_data(*sys.argv[1:])
