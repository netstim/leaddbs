# -*- coding: utf-8 -*-
"""
Tensor scaling algorithms by A. Andree, parallelization by K. Butenko.
Refactored using Claude for cross-platform multiprocessing compatibility (Linux/macOS/Windows).
"""

import os
import sys
import itertools
from functools import partial
from multiprocessing import Pool, cpu_count

import numpy as np
import nibabel as nib

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
_EPS = 1e-12

# Supported tensor storage orders and the corresponding (row, col) index
# mapping for the six unique components stored as a flat 6-element vector.
#
#   NIFTI         : [xx, yx, yy, zx, zy, zz]
#   DSI_studio    : [xx, yy, zz, yx, zx, zy]
#   FSL           : [xx, yx, zx, yy, zy, zz]
#   Johnson_Wistar: [yy, yx, zy, xx, zx, zz]
#
# Each entry is a list of (row, col) positions in the symmetric 3x3 matrix.
_TENSOR_LAYOUTS = {
    "NIFTI":          [(0,0),(1,0),(1,1),(2,0),(2,1),(2,2)],
    "DSI_studio":     [(0,0),(1,1),(2,2),(1,0),(2,0),(2,1)],
    "FSL":            [(0,0),(1,0),(2,0),(1,1),(2,1),(2,2)],
    "Johnson_Wistar": [(1,1),(1,0),(2,1),(0,0),(2,0),(2,2)],
}

# Lower conductivity boundary: WM value from Gabriel et al. at 10 Hz [S/m]
_SIGMA_ISO_LOW = 0.027512
_SIGMA_ISO_LOWER_BOUNDARY = _SIGMA_ISO_LOW / 20.0


# ---------------------------------------------------------------------------
# Scaling helpers
# ---------------------------------------------------------------------------
def _theta_star(w12, w13):
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


def _flat_to_matrix(components, tensor_order):
    """Reconstruct the symmetric 3x3 diffusion tensor from six stored values."""
    layout = _TENSOR_LAYOUTS[tensor_order]
    mat = np.zeros((3, 3), dtype=float)
    for value, (row, col) in zip(components, layout):
        mat[row, col] = value
        mat[col, row] = value  # enforce symmetry
    return mat


def _scale_eigenvalues(eigvals, eigvecs, matrix, scaling_method, diagonal_sum):
    """
    Return (scaled_tensor, eigvals_of_scaled_tensor).

    For 'Nordin' the tensor is scaled directly; for all other methods a new
    tensor is reconstructed from scaled eigenvalues and original eigenvectors.
    """
    if scaling_method == "Tuch":
        # as in Tuch el al. 2001 linear fit
        k_dti = 0.844e3    # S*s/mm^3  (converts to S/m * s/mm^2)
        d_eps = 0.124e-3   # mm^2/s  (extracellular diffusivity)
        scaled = k_dti * (eigvals - d_eps)

    elif scaling_method == "NormMapping":
        # as in Güllmar et al. 2010
        det_cbrt = (eigvals[0] * eigvals[1] * eigvals[2]) ** (1.0 / 3.0)
        scaled = eigvals / det_cbrt

    elif scaling_method == "LoadPreservation":
        # as in Howell & McIntyre, 2016.
        w12 = eigvals[0] / (eigvals[1] + _EPS)
        w13 = eigvals[0] / (eigvals[2] + _EPS)
        theta = _theta_star(w12, w13)
        scaled = np.array([1.0, 1.0 / (w12 + _EPS), 1.0 / (w13 + _EPS)]) * theta

    elif scaling_method == "Nordin":
        # as in Nordin el al. 2019
        tensor = matrix / (diagonal_sum / 3.0)
        return tensor, np.linalg.eigvalsh(tensor)

    else:
        raise ValueError("Unknown scaling method: {!r}".format(scaling_method))

    if np.any(scaled <= 0):
        raise RuntimeError("Scaled eigenvalues must be strictly positive.")

    tensor = eigvecs.dot(np.diag(scaled)).dot(eigvecs.T)
    return tensor, np.linalg.eigvalsh(tensor)


# ---------------------------------------------------------------------------
# Per-column worker  (no shared memory — receives data, returns results)
# ---------------------------------------------------------------------------
def _process_column(k_vector, tensor_order, scaling_method, affine, affine_inv, args):
    """
    Process all voxels in the (i, j) column across all z slices.

    Receives a slice of DTI data directly and returns (i, j, result_array)
    where result_array has shape (mz, 6).

    Using return values instead of shared memory avoids WinError 1450 on
    Windows, which is caused by exhausting system handles for large shared
    memory segments.
    """
    ij, col_data = args
    i, j = ij
    mz = len(k_vector)
    out = np.zeros((mz, 6), dtype=np.float64)

    for ki, k in enumerate(k_vector):
        # ---- extract the six tensor components --------------------------------
        if col_data.ndim == 3:       # original DTI had 5 dims -> (mz, 1, 6)
            components = col_data[ki, 0, :]
        else:                        # original DTI had 4 dims -> (mz, 6)
            components = col_data[ki, :]

        # Background / zero voxel -> identity-like placeholder
        if np.all(components == 0.0):
            out[ki, :] = [1.0, 0.0, 0.0, 1.0, 0.0, 1.0]
            continue

        # ---- build symmetric 3x3 matrix --------------------------------------
        mat = _flat_to_matrix(components, tensor_order)

        eigvals, eigvecs = np.linalg.eig(mat)
        if np.any(eigvals < 0):
            print(
                "Warning: negative eigenvalue at voxel ({},{},{}); "
                "taking absolute value.  Check your DTI data!".format(i, j, k)
            )
            eigvals = np.abs(eigvals)

        # ---- scale -----------------------------------------------------------
        # diagonal_sum uses positions (0, 2, 5) = xx, yy, zz in the flat vector.
        # This is strictly correct for NIFTI order; for other orders it is an
        # approximation consistent with the original code.  The Nordin method
        # should therefore be used with NIFTI tensor order.
        diagonal_sum = float(components[0] + components[2] + components[5])
        tensor, eigvals_tensor = _scale_eigenvalues(
            eigvals, eigvecs, mat, scaling_method, diagonal_sum
        )

        if np.any(eigvals_tensor <= 0):
            raise RuntimeError(
                "Non-positive eigenvalue after scaling at voxel ({},{},{}).".format(i, j, k)
            )

        if np.any(eigvals_tensor * _SIGMA_ISO_LOW < _SIGMA_ISO_LOWER_BOUNDARY):
            raise RuntimeError(
                "Conductivity below lower boundary at voxel ({},{},{}).".format(i, j, k)
            )

        # ---- re-orient to world axes -----------------------------------------
        tc = np.zeros((4, 4), dtype=float)
        tc[:3, :3] = tensor
        tc[3, 3] = 1.0
        tensor_wa = affine.dot(tc).dot(affine_inv)

        # Store as FSL order: [xx, yx, zx, yy, zy, zz]
        out[ki, :] = [
            tensor_wa[0, 0],
            tensor_wa[1, 0],
            tensor_wa[2, 0],
            tensor_wa[1, 1],
            tensor_wa[2, 1],
            tensor_wa[2, 2],
        ]

    return (i, j, out)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
def _run_parallel(dti_data, tensor_order, scaling_method, affine):
    """
    Distribute voxel processing across CPU cores and assemble the result.

    Each worker receives a single (i, j) column of DTI voxels as a regular
    numpy array, processes all z slices, and returns the results.  No shared
    memory segments are created, keeping Windows resource usage within limits.
    """
    if np.any(np.diag(affine[:3, :3]) < 0.0):
        raise ValueError(
            "Negative diagonal in affine matrix -- "
            "reverse-ordered axes are not supported.  Please flip the image."
        )

    affine_inv = np.linalg.inv(affine)
    mx, my, mz = dti_data.shape[:3]
    k_vector = np.arange(mz)

    # Pre-extract per-column slices so each task carries only (mz, 6) worth
    # of data rather than a reference to the full volume.
    if dti_data.ndim == 5:
        col_source = dti_data[:, :, :, 0, :]   # (mx, my, mz, 6)
    else:
        col_source = dti_data                   # (mx, my, mz, 6)

    task_args = [
        ((i, j), col_source[i, j, :, :])
        for i, j in itertools.product(range(mx), range(my))
    ]

    worker = partial(
        _process_column,
        k_vector,
        tensor_order,
        scaling_method,
        affine,
        affine_inv,
    )

    normalized = np.zeros((mx, my, mz, 6), dtype=np.float64)
    n_workers = max(1, cpu_count() - 1)
    with Pool(processes=n_workers) as pool:
        for i, j, col_result in pool.map(worker, task_args):
            normalized[i, j, :, :] = col_result

    return normalized


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
def scale_tensor_data(tensor_data_name, scaling_method="NormMapping", tensor_order="NIFTI"):
    """
    Load a NIfTI diffusion tensor image, scale it to a conductivity tensor,
    and save the result alongside the original file.

    Parameters
    ----------
    tensor_data_name : str
        Path to the input NIfTI file.
    scaling_method : str
        One of 'Tuch', 'NormMapping', 'LoadPreservation', 'Nordin'.
    tensor_order : str
        Storage convention of the six tensor components:
        'NIFTI', 'DSI_studio', 'FSL', or 'Johnson_Wistar'.
    """
    if tensor_order not in _TENSOR_LAYOUTS:
        raise ValueError(
            "Unknown tensor_order {!r}. Choose from {}.".format(
                tensor_order, list(_TENSOR_LAYOUTS)
            )
        )

    filepath = os.path.realpath(tensor_data_name)
    img = nib.load(filepath)
    dti_data = img.get_fdata()

    if np.any(np.isnan(dti_data)):
        raise ValueError("NaN values detected in the DTI data -- please remove them.")

    normalized = _run_parallel(dti_data, tensor_order, scaling_method, img.affine)

    stem = nib.filename_parser.splitext_addext(filepath)[0]
    out_path = "{0}_{1}.nii.gz".format(stem, scaling_method)
    nib.save(nib.Nifti1Image(normalized, img.affine), out_path)
    print("Saved: {}".format(out_path))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # The `if __name__ == "__main__"` guard is required on Windows so that
    # spawned worker processes do not re-execute the top-level script.
    scale_tensor_data(*sys.argv[1:])
