# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 14:56:19 2018
@author: Konstantin Butenko
Routines to convert segmented MRI/DTI data into the format supported by OSS-DBS
Either .nii or .txt (grid voxel arrays, COMSOL format) files can be imported
With minor modifications, it should be applicable for segmented CT data

.../MRI_DTI_derived_data/Tissue_array_MRI stores the segmentation as a 1-D (flat) array
.../MRI_DTI_derived_data/Tensor_array_DTI stores the tensors as flat arrays for each tensor component

"""

import logging
import numpy as np
import time
import pickle
import os


class MRI_segm_data:

    """ Instance of MRI_segm_data stores necessary meta data (voxel dimensions, affine, number of voxels per axes, etc).
        Besides the segmented MRI stored in PATIENTDIR, you need to supply a dictionary that describes the data
        (see entries with inp_dict) """

    def __init__(self, inp_dict):
        self.name = inp_dict["MRI_data_name"]  # name of the external segmented (MRI) file, string
        self.obtain_MRI_class(inp_dict)

    def process_segm_MRI(self, inp_dict):
        start_MRI_processing = time.time()
        import os

        #  PREFERRED: extract everything from the nifti file, meta data is in the header
        if self.name[-3:] == 'nii' or self.name[-6:] == 'nii.gz':
            # load the nifti file
            import nibabel as nib
            example_filename = os.path.join(os.environ['PATIENTDIR'] + '/' + self.name)
            img = nib.load(example_filename)
            tissue_array = img.get_fdata()
            voxel_array_temp = tissue_array.flatten('F')  #  OSS-DBS operates 1-D arrays
            self.affine_MRI = img.affine

            self.N_voxels = list(tissue_array.shape[:])  # store number of voxels along axes
            del tissue_array

            if inp_dict["MRI_in_m"] == 1:  # switch to mm if the MRI data is in m
                self.voxel_dims = [x * 1000 for x in img.header.get_zooms()[:]]  # dimensions of voxels along axes, list
                self.first_vox_coords = [x * 1000 for x in img.affine[:3,3]]     #  coordinates of the first voxel, list
            else:
                self.voxel_dims = [x for x in img.header.get_zooms()[:]]  # dimensions of voxels along axes, list
                self.first_vox_coords = [x for x in img.affine[:3, 3]]  # coordinates of the first voxel, list

            self.last_vox_coords = []  # coordinates of the last voxel, only important for brain approx.
            self.last_vox_coords.append(self.first_vox_coords[0] + self.voxel_dims[0] * self.N_voxels[0])
            self.last_vox_coords.append(self.first_vox_coords[1] + self.voxel_dims[1] * self.N_voxels[1])
            self.last_vox_coords.append(self.first_vox_coords[2] + self.voxel_dims[2] * self.N_voxels[2])

        else:
            # processes a .txt segmentation file (COMSOL format)
            # For an example of such a file, see https://github.com/SFB-ELAINE/OSS-DBS/blob/master/OSS_platform/Example_files/Dummy_MRI.txt
            infile = open(os.environ['PATIENTDIR'] + '/' + self.name, 'r').readlines() #  kicks out comment strings
            with open(os.environ['PATIENTDIR'] + '/MRI_DTI_derived_data/Filtered_' + self.name, 'w') as outfile:
                for index, line in enumerate(infile):
                    # if index != 0 and index != 4 and index != 5:      #Waxholm atlas space
                    if index != 0 and index != 4:  # here, 1st and 5th line are comments (tissue_full.txt)
                        outfile.write(line)

            spc = ' '    #  might be different

            # extracts vectors and values from segmented MRI slices
            with open(os.environ['PATIENTDIR'] + '/MRI_DTI_derived_data/Filtered_' + self.name, 'r') as f:
                for index, line in enumerate(f):
                    if index == 0:
                        line = line.rstrip()
                        x_vector = line
                        x_arr = np.fromstring(x_vector, dtype=np.float, sep=' ')
                    if index == 1:
                        line = line.rstrip()
                        y_vector = line
                        y_arr = np.fromstring(y_vector, dtype=np.float, sep=' ')
                    if index == 2:
                        line = line.rstrip()
                        z_vector = line
                        z_arr = np.fromstring(z_vector, dtype=np.float, sep=' ')
                    if index >= 3:
                        line = line.rstrip()
                        if index == 3:
                            voxel_val = line + spc
                        else:
                            voxel_val = voxel_val_old + line + spc  # if voxel data is divided into slices, this will collect all of them
                        voxel_val_old = voxel_val

                voxel_array_temp = np.fromstring(voxel_val, dtype=np.float, sep=' ')

            self.N_voxels = (x_arr.shape[0], y_arr.shape[0], z_arr.shape[0])  # number of voxels along axes

            if self.N_voxels[0] * self.N_voxels[0] * self.N_voxels[0] != voxel_array_temp.shape[0]:
                logging.critical("Error while processing segmentation, maybe not all slices were extracted")
                raise SystemExit

            if inp_dict["MRI_in_m"] == 1:  # switch to mm if the MRI data is in m
                x_arr[:] = [x_i * 1000.0 for x_i in x_arr]
                y_arr[:] = [y_i * 1000.0 for y_i in y_arr]
                z_arr[:] = [z_i * 1000.0 for z_i in z_arr]

            self.voxel_dims = [abs(round(x_arr[1] - x_arr[0], 6)), abs(round(y_arr[1] - y_arr[0], 6)), abs(round(z_arr[1] - z_arr[0], 6))]
            self.first_vox_coords = [round(min(x_arr), 6), round(min(y_arr), 6), round(min(z_arr), 6)]
            self.last_vox_coords = [round(max(x_arr), 6), round(max(y_arr), 6), round(max(z_arr), 6)]

            # we assume that in .txt the spatial vectors are in the ascending order and orthogonal
            self.affine_MRI = np.array([[self.voxel_dims[0], 0.0, 0.0, min(x_arr)],
                                        [0.0, self.voxel_dims[1], 0, min(y_arr)],
                                        [0.0, 0.0, self.voxel_dims[2], min(z_arr)],
                                        [0.0, 0.0, 0.0, 1.0]])

        # Explicitly defined for clarity. The shift is required to have the segmented MRI data in the positive octant starting in (0,0,0)
        self.MRI_shift = [x * -1.0 for x in self.first_vox_coords]

        voxel_arr = np.zeros(voxel_array_temp.shape[0], int)
        voxel_arr[voxel_array_temp == inp_dict["CSF_index"]] = 1  # changes indices to the internal notation
        voxel_arr[voxel_array_temp == inp_dict["WM_index"]] = 2
        voxel_arr[voxel_array_temp == inp_dict["GM_index"]] = 3

        ## For Lead-DBS we assume that the MRI segmentation is well-defined, no interpolation is needed
        #  Alternatively, we should have here an interpolation function

        np.save(os.environ['PATIENTDIR'] + '/MRI_DTI_derived_data/Tissue_array_MRI', voxel_arr.astype('b'),
                allow_pickle=False, fix_imports=False)
        del voxel_arr, voxel_array_temp

        # we actually store it in the class as well, so this is redundant
        np.save(os.environ['PATIENTDIR'] + '/MRI_DTI_derived_data/affine_MRI', self.affine_MRI)

        logging.critical("----- Processing of segmented (MRI) data took {} seconds -----\n".format(time.time() - start_MRI_processing))

    def obtain_MRI_class(self, inp_dict):
        if self.name == 0:
            logging.critical("MRI_data_name was not provided, exiting")
            raise Exception('exit')

        if inp_dict[ "Segm_MRI_processed"] == 0:  # 1 if MRI data were already processed by the platform and corresp. meta data were created
            logging.critical("--- processing provided segmented MRI data")
            print("\n Processing the segmented MRI...")
            self.process_segm_MRI(inp_dict)  #  creates Tissue_array_MRI and stores meta data in an instance of class MRI_segm_data

            '''Save meta data for the future simulations with the current segmented MRI data set'''
            with open(os.environ['PATIENTDIR'] + '/MRI_DTI_derived_data/Segm_MRI_meta.file', "wb") as f:
                pickle.dump(self.__dict__, f, pickle.HIGHEST_PROTOCOL)
            MRI_segm_param = self
        else:
            with open(os.environ['PATIENTDIR'] + '/MRI_DTI_derived_data/Segm_MRI_meta.file', "rb") as f:
                tmp_dict = pickle.load(f)
            self.__dict__.update(tmp_dict)
            MRI_segm_param = self
            logging.critical("--- Meta data for segmented MRI were loaded\n")

        return MRI_segm_param



class DTI_meta_data:

    """ Instance of DTI_meta_data stores necessary meta data (voxel dimensions, affine, number of voxels per axes, etc).
        Besides the tensor data stored in PATIENTDIR, you need to supply a dictionary that describes the data
        (see entries with inp_dict) """

    def __init__(self, inp_dict, Segm_param):
        self.name = inp_dict["DTI_data_name"]
        self.obtain_DTI_class(inp_dict, Segm_param)

    def process_DTI(self, inp_dict):

        start_DTI_processing = time.time()
        import os

        #  PREFERRED: extract everything from the nifti file, meta data is in the header
        if self.name[-3:] == 'nii' or self.name[-6:] == 'nii.gz':
            # load the nifti file
            import nibabel as nib
            import os
            example_filename = os.path.join(os.environ['PATIENTDIR'] + '/' + self.name)
            img = nib.load(example_filename)
            self.affine_DTI = img.affine

            if inp_dict["DTI_in_m"] == 1:   # switch to mm if the DTI data is in m
                self.voxel_dims = list(img.header.get_zooms()[:3] * 1000.0)
                self.first_vox_coords = [x * 1000 for x in img.affine[:3, 3]]  # coordinates of the first voxel, list
            else:
                self.voxel_dims = list(img.header.get_zooms()[:3])
                self.first_vox_coords = [x for x in img.affine[:3, 3]]  # coordinates of the first voxel, list

            #  here we extract a subset of DTI data just to cover only the brain approximation
            #  at the moment, only works for the orthogonal basis (otherwise takes the whole tensor)
            if inp_dict["Brain_shape_name"] == 0 and img.affine[0, 1] == 0.0 and img.affine[0, 2] == 0.0 and img.affine[1, 2] == 0.0:
                logging.critical("Extracting a subset of the tensor data for the approx. volume")

                if inp_dict['Approximating_Dimensions'][0] == 0:
                    tensor_array = img.get_fdata()
                else:
                    start_voxel = []
                    vox_window=[]
                    impl_coords = [inp_dict['Implantation_coordinate_X'],inp_dict['Implantation_coordinate_Y'],inp_dict['Implantation_coordinate_Z']]
                    for j in range(3):
                        start_voxel.append(int((impl_coords[j] - self.first_vox_coords[j]) / self.voxel_dims[j]) - int(inp_dict['Approximating_Dimensions'][j] / (2.0 * self.voxel_dims[j])))
                        vox_window.append(int(inp_dict['Approximating_Dimensions'][j] / (self.voxel_dims[j])))
                        if start_voxel[j] < 0:
                            logging.critical(
                                'Warning, the DTI data does not cover the whole computational domain (isotropic values will be assigned)')
                            start_voxel[j] = 0

                        img.affine[j, 3] = self.first_vox_coords[j] + self.voxel_dims[j] * start_voxel[j]
                        self.first_vox_coords[j] = img.affine[j, 3]  # re-assignment

                    tensor_array = img.dataobj[start_voxel[0]:start_voxel[0] + vox_window[0],
                                   start_voxel[1]:start_voxel[1] + vox_window[1], start_voxel[2]:start_voxel[2] + vox_window[2], ...]
            else:
                tensor_array = img.get_fdata()

            #  For clarity, here we explicitly declare all components of the symmetric DTI tensor
            voxel_arr_c11 = tensor_array[:, :, :, 0].flatten('F')
            voxel_arr_c21 = tensor_array[:, :, :, 1].flatten('F')
            voxel_arr_c31 = tensor_array[:, :, :, 2].flatten('F')
            voxel_arr_c22 = tensor_array[:, :, :, 3].flatten('F')
            voxel_arr_c32 = tensor_array[:, :, :, 4].flatten('F')
            voxel_arr_c33 = tensor_array[:, :, :, 5].flatten('F')

            self.N_voxels = list(tensor_array.shape[:])  # number of voxels along axes

        else:
            # processes a .txt tensor file (COMSOL format)
            # For an example of such a file, see https://github.com/SFB-ELAINE/OSS-DBS/blob/master/OSS_platform/Example_files/Dummy_DTI.txt
            infile = open(os.environ['PATIENTDIR'] + '/' + self.name, 'r').readlines()   #  kicks out comment strings
            with open(os.environ['PATIENTDIR'] + '/MRI_DTI_derived_data/Filtered_' + self.name, 'w') as outfile:

                # first, we check the length of z and y vectors (one line of DTI data corresponds to the same y and z coordinate, directions are separated by a comment)
                for index, line in enumerate(infile):
                    if index == 2:
                        line = line.rstrip()
                        y_vector = line
                        y_stps = np.fromstring(y_vector, dtype=np.float, sep=' ')
                        y_length = y_stps.shape[0]
                    if index == 3:
                        line = line.rstrip()
                        z_vector = line
                        z_stps = np.fromstring(z_vector, dtype=np.float, sep=' ')
                        z_length = z_stps.shape[0]

                # jumps over lines with comments and writes only vectors
                for index, line in enumerate(infile):
                    if index != 0 and index != 4 and index != 4 + (y_length * z_length + 1) and index != 4 + 2 * (
                            y_length * z_length + 1) and (index < 4 + 3 * (y_length * z_length + 1) or index > 4 + 4 * (
                            y_length * z_length + 1)) and index != 4 + 5 * (y_length * z_length + 1) and (
                            index < 4 + 6 * (y_length * z_length + 1) or index > 4 + 8 * (y_length * z_length + 1)):
                        outfile.write(line)

            spc = ' '  # delimiter might be different

            with open(os.environ['PATIENTDIR'] + '/MRI_DTI_derived_data/Filtered_' + self.name, 'r') as f:
                for index, line in enumerate(f):
                    if index == 0:
                        line = line.rstrip()
                        x_vector = line
                        x_arr = np.fromstring(x_vector, dtype=np.float, sep=' ')
                    if index == 1:
                        line = line.rstrip()
                        y_vector = line
                        y_arr = np.fromstring(y_vector, dtype=np.float, sep=' ')
                    if index == 2:
                        line = line.rstrip()
                        z_vector = line
                        z_arr = np.fromstring(z_vector, dtype=np.float, sep=' ')

                    '''Here we assumed the order c11,c21,c31,c22,c32,c33'''
                    if index >= 3 and index < 3 + y_length * z_length:  # corresponds to collection of slices (separated or combined) for c11
                        line = line.rstrip()
                        if index == 3:
                            voxel_val_c11 = line + spc
                        else:
                            voxel_val_c11 = voxel_val_old_c11 + line + spc
                        voxel_val_old_c11 = voxel_val_c11

                    if index >= 3 + y_length * z_length and index < 3 + 2 * y_length * z_length:
                        line = line.rstrip()
                        if index == 3 + y_length * z_length:
                            voxel_val_c21 = line + spc
                        else:
                            voxel_val_c21 = voxel_val_old_c21 + line + spc
                        voxel_val_old_c21 = voxel_val_c21

                    if index >= 3 + 2 * y_length * z_length and index < 3 + 3 * y_length * z_length:
                        line = line.rstrip()
                        if index == 3 + 2 * y_length * z_length:
                            voxel_val_c31 = line + spc
                        else:
                            voxel_val_c31 = voxel_val_old_c31 + line + spc
                        voxel_val_old_c31 = voxel_val_c31

                    if index >= 3 + 3 * y_length * z_length and index < 3 + 4 * y_length * z_length:
                        line = line.rstrip()
                        if index == 3 + 3 * y_length * z_length:
                            voxel_val_c22 = line + spc
                        else:
                            voxel_val_c22 = voxel_val_old_c22 + line + spc
                        voxel_val_old_c22 = voxel_val_c22

                    if index >= 3 + 4 * y_length * z_length and index < 3 + 5 * y_length * z_length:
                        line = line.rstrip()
                        if index == 3 + 4 * y_length * z_length:
                            voxel_val_c32 = line + spc
                        else:
                            voxel_val_c32 = voxel_val_old_c32 + line + spc
                        voxel_val_old_c32 = voxel_val_c32

                    if index >= 3 + 5 * y_length * z_length and index < 3 + 6 * y_length * z_length:
                        line = line.rstrip()
                        if index == 3 + 5 * y_length * z_length:
                            voxel_val_c33 = line + spc
                        else:
                            voxel_val_c33 = voxel_val_old_c33 + line + spc
                        voxel_val_old_c33 = voxel_val_c33

                voxel_arr_c11 = np.fromstring(voxel_val_c11, dtype=np.float, sep=' ')
                voxel_arr_c21 = np.fromstring(voxel_val_c21, dtype=np.float, sep=' ')
                voxel_arr_c31 = np.fromstring(voxel_val_c31, dtype=np.float, sep=' ')
                voxel_arr_c22 = np.fromstring(voxel_val_c22, dtype=np.float, sep=' ')
                voxel_arr_c32 = np.fromstring(voxel_val_c32, dtype=np.float, sep=' ')
                voxel_arr_c33 = np.fromstring(voxel_val_c33, dtype=np.float, sep=' ')

            self.N_voxels = [x_arr.shape[0], y_arr.shape[0], z_arr.shape[0]]  # number of voxels along axes

            if inp_dict["DTI_in_m"] == 1:
                self.voxel_dims = [abs(1000 * round(x_arr[1] - x_arr[0], 6)), abs(1000 * round(y_arr[1] - y_arr[0], 6)),
                                   abs(1000 * round(z_arr[1] - z_arr[0], 6))]  # only positive direction
                self.first_vox_coords = [1000 * round(min(x_arr), 6), 1000 * round(min(y_arr), 6), 1000 * round(min(z_arr), 6)]
            else:
                self.voxel_dims = [abs(round(x_arr[1] - x_arr[0], 6)), abs(round(y_arr[1] - y_arr[0], 6)), abs(round(z_arr[1] - z_arr[0], 6))]
                self.first_vox_coords = [round(min(x_arr), 6), round(min(y_arr), 6), round(min(z_arr), 6)]

            # we assume that in .txt the spatial vectors are in the ascending order and orthogonal
            self.affine_DTI = np.array([[self.voxel_dims[0], 0.0, 0.0, min(x_arr)],
                                        [0.0, self.voxel_dims[1], 0, min(y_arr)],
                                        [0.0, 0.0, self.voxel_dims[2], min(z_arr)],
                                        [0.0, 0.0, 0.0, 1.0]])

        if self.N_voxels[0] * self.N_voxels[1] * self.N_voxels[2] * 6 != voxel_arr_c11.shape[0] + voxel_arr_c21.shape[0] + voxel_arr_c31.shape[0] + \
                voxel_arr_c22.shape[0] + voxel_arr_c32.shape[0] + voxel_arr_c33.shape[0]:
            logging.critical("Error while processing DTI data, maybe not all slices were extracted")
            raise SystemExit

        # we actually store it in the class as well, so this is redundant
        np.save(os.environ['PATIENTDIR'] + '/MRI_DTI_derived_data/affine_DTI', self.affine_DTI)

        Tensor_array = np.zeros((voxel_arr_c11.shape[0], 6), float)
        Tensor_array[:, :] = np.vstack(
            (voxel_arr_c11[:], voxel_arr_c21[:], voxel_arr_c31[:], voxel_arr_c22[:], voxel_arr_c32[:], voxel_arr_c33[:])).T

        np.save(os.environ['PATIENTDIR'] + '/MRI_DTI_derived_data/Tensor_array_DTI', Tensor_array)

        del Tensor_array, voxel_arr_c11, voxel_arr_c21, voxel_arr_c31, voxel_arr_c22, voxel_arr_c32, voxel_arr_c33

        logging.critical("----- Processing of DTI data took {} seconds -----\n".format(time.time() - start_DTI_processing))


    def obtain_DTI_class(self, inp_dict, Segm_param):
        if inp_dict["DTI_processed"] == 0:  # 1 if DTI data were already processed by the platform and corresp. meta data were created

            logging.critical("--- processing provided DTI data")
            print("\n Processing the tensor data...")
            self.process_DTI(inp_dict)
            # check relative shift of DTI from segmented MRI
            rel_DTI_coords = [self.first_vox_coords[i] - Segm_param.first_vox_coords[i] for i in range(3)]

            eps = 0.1 * min(Segm_param.voxel_dims)
            if any(x < -1 * eps for x in rel_DTI_coords):
                logging.critical("DTI data is outside of the segmented MRI data, please truncate the former.")
                raise Exception('exit')

            '''Save meta data for the future simulations with the current MRI/DTI data set'''
            with open(os.environ['PATIENTDIR'] + '/MRI_DTI_derived_data/DTI_meta.file', "wb") as f:
                pickle.dump(self.__dict__, f, pickle.HIGHEST_PROTOCOL)
            DTI_param = self
        else:
            with open(os.environ['PATIENTDIR'] + '/MRI_DTI_derived_data/DTI_meta.file', "rb") as f:
                tmp_dict = pickle.load(f)
            self.__dict__.update(tmp_dict)
            DTI_param = self
            logging.critical("--- DTI meta data were loaded\n")

        return DTI_param
