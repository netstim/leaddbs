#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 14:09:42 2020

@author: butenko
"""

import h5py
import os
import numpy as np
import logging
logging.getLogger('UFL').setLevel(logging.WARNING)
logging.getLogger('FFC').setLevel(logging.WARNING)

from dolfin import Mesh, Point
import matplotlib.pyplot as plt
from pandas import read_csv
import time as time_lib

import nibabel as nib

#This script allows to use VTA arrays of points instead of axons in OSS-DBS

def create_VTA_array(seed_coords, rodent_electrode, MRI_param):    #in mm, in the MRI space

    # set VTA resolution
    VTA_res = max(MRI_param.voxel_dims)
    if VTA_res > 0.5:
        VTA_res = 0.5

    if rodent_electrode == True:
        VTA_box_length = 2.0
    else:
        VTA_box_length = 20.0  # can be adjusted

    # seed probe points in a cube centered at the tip (or the first contact) of the electrode
    x_vector = np.arange(seed_coords[0] - VTA_box_length / 2.0, seed_coords[0] + VTA_box_length / 2.0 + VTA_res, VTA_res)
    y_vector = np.arange(seed_coords[1] - VTA_box_length / 2.0, seed_coords[1] + VTA_box_length / 2.0 + VTA_res, VTA_res)
    z_vector = np.arange(seed_coords[2] - VTA_box_length / 2.0, seed_coords[2] + VTA_box_length / 2.0 + VTA_res, VTA_res)

    VTA_array = np.zeros((x_vector.shape[0]*y_vector.shape[0]*z_vector.shape[0], 3), float)

    total_counter = 0
    for i in range(x_vector.shape[0]):
       for j in range(y_vector.shape[0]):
           for k in range(z_vector.shape[0]):
               VTA_array[total_counter,:] = (x_vector[i], y_vector[j], z_vector[k])
               total_counter += 1

    np.savetxt(os.environ['PATIENTDIR']+'/VTA_default_array.csv', VTA_array, delimiter=" ")

    # check if probe points are inside the computational domain
    mesh = Mesh(os.environ['PATIENTDIR']+"/Meshes/Mesh_unref.xml")
    for j in range(VTA_array.shape[0]):
        VTA_array[j,:] =  VTA_array[j,:] - MRI_param.first_vox_coords[:] # shift to the positive octant space
        pnt = Point(VTA_array[j, 0], VTA_array[j, 1], VTA_array[j, 2])
        # mark the point if it is not inside the computational domain (e.g. intersects with the electrode)
        if not (mesh.bounding_box_tree().compute_first_entity_collision(pnt) < mesh.num_cells() * 10):
            VTA_array[j, :] = -100000000.0
    # remove marked
    VTA_array = VTA_array[~np.all(VTA_array == -100000000.0, axis=1)]

    np.savetxt(os.environ['PATIENTDIR']+'/Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', VTA_array, delimiter=" ")
    return(x_vector.shape[0], VTA_array.shape[0])



#here we need to add function that will check that these vertices do not intersect with CSF, encap and so on

# do not use parallelization here



def get_VTA(d,vox_along_axis,Max_signal_for_point,shift_to_MRI_space):
    import os

    example_filename = os.path.join(os.environ['PATIENTDIR']+'/'+d['MRI_data_name'])
    img = nib.load(example_filename)

    VTA_res = max([x for x in img.header.get_zooms()[:]])
    if VTA_res > 0.5:
        VTA_res = 0.5

    VTA_Vertices_get = read_csv(os.environ['PATIENTDIR']+'/Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', delimiter=' ', header=None)    # get only physiologically correct neuron models
    VTA_Vertices = VTA_Vertices_get.values

    Array_full_coord_get = read_csv(os.environ['PATIENTDIR']+'/'+'VTA_default_array.csv', delimiter=' ', header=None)    # get only physiologically correct neuron models
    Array_full_coord = Array_full_coord_get.values

    VTA_affected = np.zeros((VTA_Vertices.shape[0],4),float)
    VTA_affected[:,:3] = VTA_Vertices

    E_field_MRI_space = np.zeros((VTA_Vertices.shape[0], 4), float)
    E_field_MRI_space[:,:3] = VTA_affected[:, :3] + shift_to_MRI_space

    VTA_size = 0.0

    for i in range(VTA_Vertices.shape[0]):
        if abs(Max_signal_for_point[i]) >= d["Activation_threshold_VTA"] / 1000.0:  # activation threshold is in V/m
            VTA_affected[i,3] = 1.0
            VTA_size += VTA_res**3

        E_field_MRI_space[i,3] = Max_signal_for_point[i] * 1000.0  # switch to V/m


    logging.critical("VTA_size: {}".format(VTA_size))

    VTA_affected_MRI_space=np.zeros((VTA_Vertices.shape[0],4),float)
    VTA_affected_MRI_space[:,:3] = VTA_affected[:,:3] + shift_to_MRI_space
    VTA_affected_MRI_space[:,3] = VTA_affected[:,3]

    np.savetxt(os.environ['PATIENTDIR']+'/Field_solutions/VTA_affected.csv', VTA_affected, delimiter=" ")

    VTA_nifti = np.zeros((vox_along_axis,vox_along_axis,vox_along_axis),int)
    E_field_nifti = np.zeros((vox_along_axis,vox_along_axis,vox_along_axis),float)

    # will throw an error, because we need to have the same number of points (no extractions)
    counter_truncated = 0

    logging.critical("vox_along_axis : {}".format(vox_along_axis))

    for i in range(vox_along_axis):  #go over all voxels
        for j in range(vox_along_axis):  #go over all voxels
            for k in range(vox_along_axis):  #go over all voxels
                #total_counter=i+j*vox_along_axis+k*vox_along_axis*vox_along_axis
                total_counter=k+j*vox_along_axis+i*vox_along_axis*vox_along_axis

                if np.all(np.round(VTA_affected_MRI_space[counter_truncated,:3],6)==np.round(Array_full_coord[total_counter,:],6)):            # if coordinates match, then
                    VTA_nifti[i,j,k]=int(VTA_affected_MRI_space[counter_truncated,3])
                    E_field_nifti[i,j,k] = Max_signal_for_point[counter_truncated] * 1000.0  # switch to V/m
                    counter_truncated += 1
                else:
                    VTA_nifti[i,j,k] = 0
                    E_field_nifti[i,j,k] = 0.0

    if counter_truncated != VTA_affected_MRI_space.shape[0]:
        logging.critical("Hasn't iterated over whole VTA_affected_MRI_space, check the algorithm")
        raise SystemExit


    # IMPORTANT:
    # If we computed in MNI, then we just use nifti directly
    # otherwise easier to operate with a 2-D spatial array
    # Lead-DBS will choose what to use

    affine_info = np.eye(4)
    affine_info[0,0] = VTA_res
    affine_info[1,1] = VTA_res
    affine_info[2,2] = VTA_res
    affine_info[0:3,3] = VTA_affected_MRI_space[0,0:3]

    #img.header.set_data_dtype(np.double)

    img.header.set_zooms((VTA_res, VTA_res, VTA_res))

    # obviously, the header is wrong if computed in native
    img3 = nib.Nifti1Image(VTA_nifti, affine_info, img.header)
    if d['Stim_side']==0:
        nib.save(img3, os.environ['PATIENTDIR']+'/Results_rh/VTA_solution.nii')
    else:
        nib.save(img3, os.environ['PATIENTDIR']+'/Results_lh/VTA_solution.nii')

    img4 = nib.Nifti1Image(E_field_nifti, affine_info,img.header)
    if d['Stim_side']==0:
        np.savetxt(os.environ['PATIENTDIR'] + '/Results_rh/VTA_affected_MRI_space.csv', VTA_affected_MRI_space, delimiter=" ")
        np.savetxt(os.environ['PATIENTDIR'] + '/Results_rh/E_field_MRI_space.csv', E_field_MRI_space,delimiter=" ")
        nib.save(img4, os.environ['PATIENTDIR'] + '/Results_rh/E_field_solution.nii')
    else:
        np.savetxt(os.environ['PATIENTDIR'] + '/Results_lh/VTA_affected_MRI_space.csv', VTA_affected_MRI_space, delimiter=" ")
        np.savetxt(os.environ['PATIENTDIR'] + '/Results_lh/E_field_MRI_space.csv', E_field_MRI_space, delimiter=" ")
        nib.save(img4, os.environ['PATIENTDIR'] + '/Results_lh/E_field_solution.nii')

    return True