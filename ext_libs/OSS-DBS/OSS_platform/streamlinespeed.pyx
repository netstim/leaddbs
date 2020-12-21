# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 15:29:15 2020

@author: konstantin
"""

#taken from Dipy (https://dipy.org/)
#E. Garyfallidis et al., “Dipy, a library for the analysis of diffusion MRI data,” Frontiers in Neuroinformatics, vol. 8, 2014. 

import cython
import numpy as np

from libc.math cimport sqrt
from libc.stdlib cimport malloc, free

cimport numpy as np

        
def set_number_of_points(streamlines, nb_points=3):
    
    only_one_streamlines = True
    streamlines = [streamlines]
    
    if len(streamlines) == 0:
        return []

    if nb_points < 2:
        raise ValueError("nb_points must be at least 2")
        
    dtype = streamlines[0].dtype
    for streamline in streamlines:
        if streamline.dtype != dtype:
            dtype = None
            
        if len(streamline) < 2:
            raise ValueError("All streamlines must have at least 2 points.")
            
    # Allocate memory for each modified streamline
    new_streamlines = []
    cdef np.npy_intp i


    if dtype == np.float32:
        # All streamlines have composed of float32 points
        for i in range(len(streamlines)):
            streamline = streamlines[i].astype(dtype)
            modified_streamline = np.empty((nb_points, streamline.shape[1]),
                                           dtype=streamline.dtype)
            c_set_number_of_points[float2d](streamline, modified_streamline)
            # HACK: To avoid memleaks we have to recast with astype(dtype).
            new_streamlines.append(modified_streamline.astype(dtype))

    elif dtype == np.float64:
        for i in range(len(streamlines)):
            streamline = streamlines[i].astype(dtype)
            modified_streamline = np.empty((nb_points, streamline.shape[1]),
                                           dtype=streamline.dtype)
            c_set_number_of_points[double2d](streamline, modified_streamline)
            # HACK: To avoid memleaks we have to recast with astype(dtype).
            new_streamlines.append(modified_streamline.astype(dtype))
            
    else:
        print("Wrong data type, check out Dipy github to expand the functionality")
        
    if only_one_streamlines:
        return new_streamlines[0]
    else:
        return new_streamlines
        
cdef void c_arclengths(Streamline streamline, double* out) nogil:
    cdef np.npy_intp i = 0
    cdef double dn

    out[0] = 0.0
    for i in range(1, streamline.shape[0]):
        out[i] = 0.0
        for j in range(streamline.shape[1]):
            dn = streamline[i, j] - streamline[i-1, j]
            out[i] += dn*dn

        out[i] = out[i-1] + sqrt(out[i])
        
cdef void c_set_number_of_points(Streamline streamline, Streamline out) nogil:
    cdef:
        np.npy_intp N = streamline.shape[0]
        np.npy_intp D = streamline.shape[1]
        np.npy_intp new_N = out.shape[0]
        double ratio, step, next_point, delta
        np.npy_intp i, j, k, dim

    # Get arclength at each point.
    arclengths = <double*> malloc(streamline.shape[0] * sizeof(double))
    c_arclengths(streamline, arclengths)

    step = arclengths[N-1] / (new_N-1)

    next_point = 0.0
    i = 0
    j = 0
    k = 0

    while next_point < arclengths[N-1]:
        if next_point == arclengths[k]:
            for dim in range(D):
                out[i, dim] = streamline[j, dim]

            next_point += step
            i += 1
            j += 1
            k += 1
        elif next_point < arclengths[k]:
            ratio = 1 - ((arclengths[k]-next_point) /
                         (arclengths[k]-arclengths[k-1]))

            for dim in range(D):
                delta = streamline[j, dim] - streamline[j-1, dim]
                out[i, dim] = streamline[j-1, dim] + ratio * delta

            next_point += step
            i += 1
        else:
            j += 1
            k += 1

    # Last resampled point always the one from original streamline.
    for dim in range(D):
        out[new_N-1, dim] = streamline[N-1, dim]

    free(arclengths)