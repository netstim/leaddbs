# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 18:31:06 2021
@author: Konstantin Butenko
Neuron_array class with various methods for neuron allocation and filtering
Supports VTA, custom allocation and meta-data preparation for PAM
The script is not written in the most elegant manner, but to ensure a good understanding
IMPORTANT: '_O' in the name refers to the arrays centered around O(0,0,0)
           '_MRI' to coordinates in the MRI space
           '_PO' to the positive octant coordinates, i.e. the MRI space shifted to start in O(0,0,0). OSS-DBS processes data in PO to avoid confusion in coordinate-to-voxel indexing
"""

import h5py
import numpy as np

import logging

from dolfin import Mesh, SubMesh, Point, MeshFunction
from pandas import read_csv
from Axon_files.axon import Axon
import time
import os
import pickle
from scipy.spatial import distance


class Neuron_array(object):
    def __init__(self, input_dict,
                 MRI_param):  # for better readability, resaving the input dictionary to more explicit variables

        self.MRI_first_voxel = MRI_param.first_vox_coords  # coordinates of the first voxel (center), will be used to shift to the positive octant (PO)
        self.impl_coordinates = np.array(
            [input_dict['Implantation_coordinate_X'], input_dict['Implantation_coordinate_Y'],
             input_dict['Implantation_coordinate_Z']])  # in the MRI space, mm

        self.neuron_model = input_dict[
            'Axon_Model_Type']  # for now only 2 axons, but can be modified for other neuron models (see Butenko's doctoral thesis)
        self.brain_model = input_dict['Brain_shape_name']

        # if a neuron array is imported, just use self.pattern to check consistency
        self.pattern = {'fiber_diameter': input_dict['diam_fib'],
                        # list if multiple datasets imported from .h5, otherwise a float number
                        'num_Ranvier': input_dict['n_Ranvier'],
                        # list if multiple datasets imported from .h5, otherwise a integer number
                        'name': input_dict['pattern_model_name'],
                        # '0' or 0 if not provided
                        }
        if input_dict['Name_prepared_neuron_array'] == 0:  # changed from '0' to 0 in Dict_corrector.py
            self.imp_name = None
        else:
            self.imp_name = input_dict['Name_prepared_neuron_array']

        if input_dict['Neuron_model_array_prepared'] == 1:
            self.Type = 'Imported'
        elif input_dict['Global_rot'] == 1:  # grid allocation of neurons for VTA according
            self.Type = 'VTA'
            self.VTA_structure = {
                'Seed_coordinates': [input_dict['x_seed'], input_dict['y_seed'], input_dict['z_seed']],
                # usually at the electrode's tip or the last contact
                'N_seeding_steps': [input_dict['x_steps'], input_dict['y_steps'], input_dict['z_steps']],
                # how many additional axons will be placed along the axes (not counting the seeded one)
                'Dist_seeding_steps': [input_dict['x_step'], input_dict['y_step'], input_dict['z_step']],
                # distance between axons (center nodes) along the axes
                'X_angles_glob': input_dict['alpha_array_glob'],
                # list containing rotational angles around X axis, already converted to radians in Dict_corrector.py
                'Y_angles_glob': input_dict['beta_array_glob'],
                'Z_angles_glob': input_dict['gamma_array_glob'], }
        else:
            self.Type = 'Custom'  # custom allocation of neurons according to entries 'X_coord_old', ... and 'YZ_angles', ... in GUI_inp_dict.py
            self.custom_structure = {
                'Center_coordinates': [input_dict['X_coord_old'], input_dict['Y_coord_old'], input_dict['Z_coord_old']],
                # list of lists [[x1,x2,x2],[y1,y2,y3],...], usually at the electrode's tip or the last contact
                'X_angles_loc': input_dict['YZ_angles'],
                # list containing rotational angles around X axis. IMPORTANT: the neurons are first centered at O(0,0,0), rotated, and then translated to their center coordinates
                'Y_angles_loc': input_dict['ZX_angles'],  # already converted to radians in Dict_corrector.py
                'Z_angles_loc': input_dict['XY_angles'],
            }

    def rotate_globally(self, point_coords, inx_angle):

        """ rotates a 3D point defined by point_coords (x,y,z) around the global Cartesian axes (angles in radians)
            point_coords is usually a center node of a seeded axon
            Input  : class,  array-like shape (1,3) - Cartesian point coordinates
            Returns: a rotated 3D point (1D numpy array) """

        # passed this way for better readability
        alpha, beta, gamma = [self.VTA_structure['X_angles_glob'][inx_angle],
                              self.VTA_structure['Y_angles_glob'][inx_angle],
                              self.VTA_structure['Z_angles_glob'][inx_angle]]

        alpha_matrix = np.array([[1, 0, 0],
                                 [0, np.cos(alpha), -1 * np.sin(alpha)],
                                 [0, np.sin(alpha), np.cos(alpha)]])

        beta_matrix = np.array([[np.cos(beta), 0, np.sin(beta)],
                                [0, 1, 0],
                                [-1 * np.sin(beta), 0, np.cos(beta)]])

        gamma_matrix = np.array([[np.cos(gamma), -1 * np.sin(gamma), 0],
                                 [np.sin(gamma), np.cos(gamma), 0],
                                 [0, 0, 1]])

        xyz_alpha = np.array(point_coords)
        xyz_beta = alpha_matrix.dot(xyz_alpha)
        xyz_gamma = beta_matrix.dot(xyz_beta)

        point_coords_rotated = gamma_matrix.dot(xyz_gamma)

        return point_coords_rotated

    def place_neuron(self, alpha, beta, gamma, x_loc, y_loc,
                     z_loc):  # x_loc, y_loc, z_loc are coordinates of seeds in PO space

        """ takes a pattern model (self.pattern['Array_coord_pattern_O'], 2D numpy array (x,y,z)),
            rotates it around the global Cartesian axes (angles in radians) and shifts its center to (x_loc,y_loc,z_loc)
            Input  : class,  3 angles (radiants), 3 coordinates
            Returns: a rotated and shifted 3D pattern of points (2D numpy array) """

        # passed this way for better readability
        Arr_coord = self.pattern['Array_coord_pattern_O']

        Neuron_coord = np.zeros((Arr_coord.shape[0], 3), float)

        alpha_matrix = np.array([[1, 0, 0],
                                 [0, np.cos(alpha), -1 * np.sin(alpha)],
                                 [0, np.sin(alpha), np.cos(alpha)]])

        beta_matrix = np.array([[np.cos(beta), 0, np.sin(beta)],
                                [0, 1, 0],
                                [-1 * np.sin(beta), 0, np.cos(beta)]])

        gamma_matrix = np.array([[np.cos(gamma), -1 * np.sin(gamma), 0],
                                 [np.sin(gamma), np.cos(gamma), 0],
                                 [0, 0, 1]])

        for inx in range(Neuron_coord.shape[0]):
            xyz_alpha = np.array([Arr_coord[inx, 0], Arr_coord[inx, 1], Arr_coord[inx, 2]])
            [Neuron_coord[inx, 0], Neuron_coord[inx, 1], Neuron_coord[inx, 2]] = alpha_matrix.dot(xyz_alpha)

            xyz_beta = np.array([Neuron_coord[inx, 0], Neuron_coord[inx, 1], Neuron_coord[inx, 2]])
            [Neuron_coord[inx, 0], Neuron_coord[inx, 1], Neuron_coord[inx, 2]] = beta_matrix.dot(xyz_beta)

            xyz_gamma = np.array([Neuron_coord[inx, 0], Neuron_coord[inx, 1], Neuron_coord[inx, 2]])
            [Neuron_coord[inx, 0], Neuron_coord[inx, 1], Neuron_coord[inx, 2]] = gamma_matrix.dot(xyz_gamma)

        Neuron_coord[:, 0] = np.round(Neuron_coord[:, 0] + x_loc, 6)
        Neuron_coord[:, 1] = np.round(Neuron_coord[:, 1] + y_loc, 6)
        Neuron_coord[:, 2] = np.round(Neuron_coord[:, 2] + z_loc, 6)

        return Neuron_coord

    def mark_to_remove(self, array, index, nsegm=0):

        'simply markes (with -100000000.0) the part of the array (contains all compartments of all neurons) that belongs to the neuron that has to be removed'
        'inx is the index in array that violated the seeding rules, nsegm is the number of compartments per neuron'

        if nsegm == 0:  # if not passed from a list
            nsegm = self.pattern['num_segments']

        inx_neuron = int(index / nsegm)
        inx_start = inx_neuron * nsegm
        array[inx_start: inx_start + nsegm, :] = -100000000.0
        index = inx_start + nsegm

        return array, index, inx_neuron

    # creates a VTA grid based on settings in GUI_inp_dict.py when 'Global_rot' = 1
    # returns a 2D numpy array of middle (seeeding) nodes of axons centered around (0,0,0)
    def create_structured_array(self):

        """ creates a VTA grid based on settings in GUI_inp_dict.py when 'Global_rot' = 1
            returns a 2D numpy array of middle (seeeding) nodes of axons centered around (0,0,0)
            Input  : class
            Adds to self: seeding nodes for VTA axons (a 2D numpy array) """

        # contain coordinates of the VTA axons centered at (0,0,0)
        self.N_models_in_plane = int((self.VTA_structure['N_seeding_steps'][0] + 1) * (
                self.VTA_structure['N_seeding_steps'][1] + 1) * (self.VTA_structure['N_seeding_steps'][2] + 1))
        self.N_total_of_axons = int((self.VTA_structure['N_seeding_steps'][0] + 1) * (
                self.VTA_structure['N_seeding_steps'][1] + 1) * (
                                            self.VTA_structure['N_seeding_steps'][2] + 1) * len(
            self.VTA_structure["X_angles_glob"]))  # +1 because we have the initial axon
        self.VTA_seeds_O = np.zeros((self.N_total_of_axons, 3), float)

        # computed as x_start_point=0-(d["x_step"]*(d["x_steps"])/2)
        start_point = [-1 * self.VTA_structure['Dist_seeding_steps'][0] * self.VTA_structure['N_seeding_steps'][0] / 2,
                       -1 * self.VTA_structure['Dist_seeding_steps'][1] * self.VTA_structure['N_seeding_steps'][1] / 2,
                       -1 * self.VTA_structure['Dist_seeding_steps'][2] * self.VTA_structure['N_seeding_steps'][2] / 2]

        # seed axons for one direction
        x_one_dir = []
        y_one_dir = []
        z_one_dir = []

        # iterate over axes
        for x_ind in range(self.VTA_structure['N_seeding_steps'][
                               0] + 1):  # number of steps corresponds to additional axons along this axis besides the initial one centered at (0,0,0) along Y-axis
            for y_ind in range(self.VTA_structure['N_seeding_steps'][1] + 1):
                for z_ind in range(self.VTA_structure['N_seeding_steps'][2] + 1):
                    x_one_dir.append(start_point[0] + x_ind * self.VTA_structure['Dist_seeding_steps'][0])
                    y_one_dir.append(start_point[1] + y_ind * self.VTA_structure['Dist_seeding_steps'][1])
                    z_one_dir.append(start_point[2] + z_ind * self.VTA_structure['Dist_seeding_steps'][2])

                    # rotate (copy) axons around (0,0,0) to create VTAs of different directionality
        gl_ind = 0
        for inx_angle in range(len(self.VTA_structure[
                                       "X_angles_glob"])):  # d["alpha_array_glob"] - rotation around X, d["beta_array_glob"] - around Y. These angle lists should have the same length.
            for inx in range(len(x_one_dir)):
                Rotated_point = self.rotate_globally([x_one_dir[inx], y_one_dir[inx], z_one_dir[inx]], inx_angle)

                self.VTA_seeds_O[gl_ind, :] = Rotated_point

                gl_ind += 1

    def get_neuron_morhology(self, fib_diam,
                             N_Ranvier):  # pass fib_diam and N_Ranvier explicitly, because they might be from a list

        """
            Input  : class
            Returns: dictionary with neuron morphology, number of segments (compartments) per neuron, number of compartments per section """

        if self.neuron_model == 'McIntyre2002':

            # load the morphology
            param_ax = {
                'centered': True,
                'diameter': fib_diam  # float, not a list here! Will throw an error if uneligible fiber diameter
            }
            a = Axon(param_ax)
            nr = Axon.get_axonparams(
                a)  # we need here distances between compartments: nr["ranvier_length"], nr["para1_length"], nr["para2_length"], nr["inter_length"], nr["deltax"]

            n_comp = int(((nr["ranvier_nodes"] - 1) + nr["inter_nodes"] + nr["para1_nodes"] + nr["para2_nodes"]) / (
                    nr["ranvier_nodes"] - 1))
            n_segm = int((N_Ranvier - 1) * n_comp + 1)  # overall number of compartments on one axon incl. internodal
        elif self.neuron_model == 'McIntyre2002_ds':

            # load the morphology
            param_ax = {
                'centered': True,
                'diameter': fib_diam  # float, not a list here! Will throw an error if uneligible fiber diameter
            }
            a = Axon(param_ax)
            nr = Axon.get_axonparams(
                a)  # we need here distances between compartments: nr["ranvier_length"], nr["para1_length"], nr["para2_length"], nr["inter_length"], nr["deltax"]

            if fib_diam >= 5.7:
                n_comp = 3   # node -- -- internodal -- -- -- -- internodal -- -- node
            else:
                n_comp = 2   # mode -- -- -- internodal -- -- -- node

            n_segm = int(
                (N_Ranvier - 1) * n_comp + 1)  # overall number of compartments on one axon incl. internodal
        elif self.neuron_model == 'Reilly2016':
            n_comp = 2
            n_segm = int((
                                 N_Ranvier - 1) * n_comp + 1)  # Reilly's (Carnevale's implementation) model has 1 internodal compartment per section
            nr = {}
            nr['deltax'] = fib_diam * 200.0  # from 1 to 2 micrometers
            if fib_diam > 10.0 or fib_diam < 5.0:
                logging.critical("Wrong fiber diameter for Reilly2016, exiting")
                raise SystemExit
        else:
            logging.critical("The neuron model {} is not implemented, exiting".format(self.neuron_model))
            raise SystemExit

        return nr, n_segm, n_comp

    def generate_pattern(self):

        """ generates a neuron pattern centered around (0,0,0)
            used for VTA and sometimes custom neuron arrays
            Input  : class
            Adds to self: self.pattern['num_segments'], self.pattern['Array_coord_pattern_O'], stores it in a .csv file
            self.pattern['Array_coord_pattern_O']: coordinates of neuron compartments (segments) (2D numpy array), center at O(0,0,0) """

        nr, self.pattern['num_segments'], n_comp = self.get_neuron_morhology(self.pattern['fiber_diameter'],
                                                                             self.pattern['num_Ranvier'])

        self.pattern['Array_coord_pattern_O'] = np.zeros((self.pattern['num_segments'], 3),
                                                         float)  # _O refers to that fact that the pattern model is centered at O(0,0,0)

        # first compartment (not the middle!)
        self.pattern['Array_coord_pattern_O'][0, 0] = 0.0
        self.pattern['Array_coord_pattern_O'][0, 1] = -0.001 * nr["deltax"] * (self.pattern[
                                                                                   'num_Ranvier'] - 1) / 2.0  # deltax is in Âµm, therefore scaling (OSS-DBS operates with mm)
        self.pattern['Array_coord_pattern_O'][0, 2] = 0.0

        loc_inx = 1  # because first node (with index 0) was already seeded
        if self.neuron_model == 'McIntyre2002':

            for inx in range(1, self.pattern['num_segments']):
                if self.pattern[
                    'fiber_diameter'] >= 5.7:  # if larger than 3.0, 6 central compartments are required, otherwise 3
                    if loc_inx == 0:
                        l_step = (nr["para1_length"] + nr["ranvier_length"]) / 2000  # switch to mm from Âµm
                    if loc_inx == 1 or loc_inx == 11:
                        l_step = (nr["ranvier_length"] + nr["para1_length"]) / 2000  # switch to mm from Âµm
                    if loc_inx == 2 or loc_inx == 10:
                        l_step = (nr["para1_length"] + nr["para2_length"]) / 2000  # switch to mm from Âµm
                    if loc_inx == 3 or loc_inx == 9:
                        l_step = (nr["para2_length"] + nr["inter_length"]) / 2000  # switch to mm from Âµm
                    if loc_inx == 4 or loc_inx == 5 or loc_inx == 6 or loc_inx == 7 or loc_inx == 8:
                        l_step = nr["inter_length"] / 1000  # switch to mm from Âµm
                else:
                    if loc_inx == 0:
                        l_step = (nr["para1_length"] + nr["ranvier_length"]) / 2000  # switch to mm from Âµm
                    if loc_inx == 1 or loc_inx == 8:
                        l_step = (nr["ranvier_length"] + nr["para1_length"]) / 2000  # switch to mm from Âµm
                    if loc_inx == 2 or loc_inx == 7:
                        l_step = (nr["para1_length"] + nr["para2_nodes"]) / 2000  # switch to mm from Âµm
                    if loc_inx == 3 or loc_inx == 6:
                        l_step = (nr["para2_nodes"] + nr["inter_length"]) / 2000  # switch to mm from Âµm
                    if loc_inx == 4 or loc_inx == 5:
                        l_step = nr["inter_length"] / 1000  # switch to mm from Âµm

                loc_inx += 1
                self.pattern['Array_coord_pattern_O'][inx, :] = [0.0, self.pattern['Array_coord_pattern_O'][
                    inx - 1, 1] + l_step, 0.0]  # pattern along Y-axis

                if inx % n_comp == 0:
                    loc_inx = 0

        elif self.neuron_model == 'McIntyre2002_ds':
            for inx in range(1, self.pattern['num_segments']):

                if self.pattern['fiber_diameter'] >= 5.7:  # if larger than 3.0, 6 central compartments are required, otherwise 3
                    if loc_inx == 0:
                        l_step = nr["para1_length"] / 1000 + nr["para2_length"] / 1000 + (nr["ranvier_length"] + nr["inter_length"]) / 2000
                    if loc_inx == 1 or loc_inx == 3:
                        l_step = (nr["ranvier_length"] + nr["inter_length"]) / 2000 + nr["para1_length"] / 1000 + nr["para2_length"] / 1000  # switch to mm from micro m
                    if loc_inx == 2:
                        l_step = 5 * nr["inter_length"] / 2000
                else:
                    l_step = (nr["para1_length"] + nr["para2_length"] + nr["inter_length"])/ 1000 + (nr["ranvier_length"] + nr["inter_length"]) / 2000

                loc_inx += 1
                self.pattern['Array_coord_pattern_O'][inx, :] = [0.0, self.pattern['Array_coord_pattern_O'][
                    inx - 1, 1] + l_step, 0.0]  # pattern along Y-axis

                if inx % n_comp == 0:
                    loc_inx = 0

            # the downsamling is complete, reassign here
            if self.neuron_model == 'McIntyre2002_ds':
                self.neuron_model = 'McIntyre2002'

        elif self.neuron_model == 'Reilly2016':
            l_step = nr["deltax"] / 2000.0  # linear change of the internodal distance from 1 to 2 mm
            for inx in range(1, self.pattern['num_segments']):
                self.pattern['Array_coord_pattern_O'][inx, :] = [0.0, self.pattern['Array_coord_pattern_O'][
                    inx - 1, 1] + l_step, 0.0]  # pattern along Y-axis
        else:
            logging.critical("The neuron model is not implemented, exiting")
            raise SystemExit

        self.pattern['Array_coord_pattern_O'] = np.round(self.pattern['Array_coord_pattern_O'], 8)
        np.savetxt(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/' + str(self.pattern['name']),
                   self.pattern['Array_coord_pattern_O'],
                   delimiter=" ")

    def allocate_neurons_PO(self):

        """ allocates neurons for VTA or according to the customized input. If imported, just stores in Neuron_model_arrays/All_neuron_models.csv
            also computes the extent of the region of interest ROI: the distance from the implantation point to the furthest compartment
            Input  : class
            Adds to self: self.ROI_radius, self.VTA_coord_PO or self.custom_coord_PO for VTA and Custom, respectively
            self.VTA_coord_PO: coordinates of neuron compartments (segments) that form VTA (2D numpy array), in PO, stores in Neuron_model_arrays/All_neuron_models.csv
            self.custom_coord_PO: coordinates of custom allocated neuron compartments (segments) (2D numpy array), in PO, stores in Neuron_model_arrays/All_neuron_models.csv """

        if self.Type == 'Custom':  # if manual placement
            total_number_of_compartments = self.pattern['num_segments'] * len(
                self.custom_structure['Center_coordinates'][0]) * len(self.custom_structure["X_angles_loc"]) * len(
                self.custom_structure["Y_angles_loc"]) * len(self.custom_structure["Z_angles_loc"])
            self.custom_coord_PO = np.zeros((total_number_of_compartments, 3), float)
            loc_ind = 0

            # coordinates of the manual allocation shifted to PO
            X_locations_PO = np.array(self.custom_structure['Center_coordinates'][0]) - self.MRI_first_voxel[0]
            Y_locations_PO = np.array(self.custom_structure['Center_coordinates'][1]) - self.MRI_first_voxel[1]
            Z_locations_PO = np.array(self.custom_structure['Center_coordinates'][2]) - self.MRI_first_voxel[2]

            for inx in range(len(X_locations_PO)):  # goes through all seeding (central) nodes of neurons (axons)
                for X_angle in self.custom_structure["X_angles_loc"]:
                    for Y_angle in self.custom_structure["Y_angles_loc"]:
                        for Z_angle in self.custom_structure["Z_angles_loc"]:
                            self.custom_coord_PO[loc_ind:loc_ind + self.pattern['num_segments'], :] = self.place_neuron(
                                X_angle, Y_angle, Z_angle, X_locations_PO[inx], Y_locations_PO[inx],
                                Z_locations_PO[inx])
                            loc_ind = loc_ind + self.pattern['num_segments']

            np.savetxt(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/All_neuron_models.csv', self.custom_coord_PO,
                       delimiter=" ")
            self.ROI_radius = max(distance.cdist(self.custom_coord_PO, self.impl_coords_PO, 'euclidean'))[0]

        elif self.Type == 'VTA':  # if placement in the ordered array (as for VTA)

            total_number_of_compartments = self.pattern['num_segments'] * self.VTA_seeds_PO.shape[0]
            self.VTA_coord_PO = np.zeros((total_number_of_compartments, 3), float)

            loc_ind = 0

            for inx in range(
                    self.VTA_seeds_PO.shape[0]):  # goes through all seeding (central) nodes of neurons (axons)
                angle_X = self.VTA_structure["X_angles_glob"][int(inx / self.N_models_in_plane)]
                angle_Y = self.VTA_structure["Y_angles_glob"][int(inx / self.N_models_in_plane)]
                angle_Z = self.VTA_structure["Z_angles_glob"][int(inx / self.N_models_in_plane)]
                self.VTA_coord_PO[loc_ind:loc_ind + self.pattern['num_segments'], :] = self.place_neuron(angle_X,
                                                                                                         angle_Y,
                                                                                                         angle_Z,
                                                                                                         self.VTA_seeds_PO[
                                                                                                             inx, 0],
                                                                                                         self.VTA_seeds_PO[
                                                                                                             inx, 1],
                                                                                                         self.VTA_seeds_PO[
                                                                                                             inx, 2])

                loc_ind = loc_ind + self.pattern['num_segments']

            np.savetxt(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/All_neuron_models.csv', self.VTA_coord_PO,
                       delimiter=" ")
            self.ROI_radius = max(distance.cdist(self.VTA_coord_PO, self.impl_coords_PO, 'euclidean'))[0]

        elif self.Type == 'Imported':
            # we already have everything in self.imp_proc_coord_PO
            # self.imp_proc_coord_PO        # here already merged from all projections !!! (extract each and

            self.ROI_radius = max(distance.cdist(self.imp_proc_coord_PO, self.impl_coords_PO, 'euclidean'))[0]

            np.savetxt(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/All_neuron_models.csv', self.imp_proc_coord_PO,
                       delimiter=" ")

        logging.critical(
            "Initial neuron models can be visualized from Neuron_model_arrays/All_neuron_models.csv in Paraview")
        logging.critical("Max distance from a compartment in the neuron array to the implantation site: {}".format(
            float(self.ROI_radius)))

    def build_neuron_models(self):

        """ a manager function that calls other methods and changes coordinate spaces, can be deprecated in the future """

        if self.Type == 'VTA':  # for VTA

            self.create_structured_array()  # creates a self.VTA_seeds_O

            # the center node of the center neuron is located at the seed (the seeding coordinate is in the positive octant space)
            self.VTA_seeds_PO = np.zeros((self.VTA_seeds_O.shape[0], 3), float)
            self.VTA_seeds_PO[:, 0] = self.VTA_seeds_O[:, 0] + self.VTA_structure['Seed_coordinates'][0] - \
                                      self.MRI_first_voxel[0]
            self.VTA_seeds_PO[:, 1] = self.VTA_seeds_O[:, 1] + self.VTA_structure['Seed_coordinates'][1] - \
                                      self.MRI_first_voxel[1]
            self.VTA_seeds_PO[:, 2] = self.VTA_seeds_O[:, 2] + self.VTA_structure['Seed_coordinates'][2] - \
                                      self.MRI_first_voxel[2]

        if self.pattern['name'] == '0' or self.pattern['name'] == 0:  # we can build axon model pattern
            self.pattern['name'] = 'default_pattern.csv'  # default name
            self.generate_pattern()  # creates self.pattern['Array_coord_pattern_O']
        else:
            # load from the provided file
            Array_coord_load_get = read_csv(os.environ['PATIENTDIR'] + '/' + self.pattern['name'], delimiter=' ',
                                            header=None)
            self.pattern['Array_coord_pattern_O'] = Array_coord_load_get.values
            self.pattern['num_segments'] = self.pattern['Array_coord_pattern_O'].shape[
                0]  # all segments should be in the pattern model

        # implantation_coordinates in PO
        Imp_X_PO = self.impl_coordinates[0] - self.MRI_first_voxel[0]
        Imp_Y_PO = self.impl_coordinates[1] - self.MRI_first_voxel[1]
        Imp_Z_PO = self.impl_coordinates[2] - self.MRI_first_voxel[2]
        self.impl_coords_PO = np.array([[Imp_X_PO, Imp_Y_PO, Imp_Z_PO]])  # 2D, otherwise distance cdist won't work

        # also computes ROI_radius
        # definition of ROI is required for adaptive mesh refinement. ROI is a sphere encompassing all neuron models
        self.allocate_neurons_PO()

        logging.critical("Initial neuron models and corresponding meta data were created")

    def process_external_array(self):

        """ processes imported neuron array, filters out neurons
            that are outside of the computational domain (if brain geometry was not imported)
            and calls allocate_neurons_PO()
            Input  : class
            Adds to self: self.imp_proc_coord_PO, self.pattern['num_segments']            'Neuron_model_arrays/Processed_neuron_array_MRI_space.h5' , 'Neuron_model_arrays/All_neuron_models_by_populations.h5'
            self.imp_proc_coord_PO: coordinates of neurons (compartments/segments) that are inside the computational domain (2D numpy array), in PO,
            stored in 'Neuron_model_arrays/All_neuron_models_by_populations.h5' with different projections in separate datasets,
            'Neuron_model_arrays/Processed_neuron_array_MRI_space.h5' same in the NRI space
            self.pattern['num_segments']: a list(!) of comparments per neuron for each projection"""

        # shift to positive octant space, x_shift is -1*[0,3] entry in the affine matrix of the MRI data
        # implantation_coordinates in PO
        Imp_X_PO = self.impl_coordinates[0] - self.MRI_first_voxel[0]
        Imp_Y_PO = self.impl_coordinates[1] - self.MRI_first_voxel[1]
        Imp_Z_PO = self.impl_coordinates[2] - self.MRI_first_voxel[2]
        self.impl_coords_PO = np.array([[Imp_X_PO, Imp_Y_PO, Imp_Z_PO]])  # 2D, otherwise distance,cdist won't work

        if self.brain_model == 'Brain_substitute.brep':  # cut models if use brain approximation
            mesh_brain_sub = Mesh(os.environ['PATIENTDIR'] + "/Meshes/Mesh_brain_substitute_max_ROI.xml")

        if not isinstance(self.pattern['fiber_diameter'],
                          list):  # for simplicity, we will treat it as a list with a single entry
            self.pattern['fiber_diameter'] = [self.pattern['fiber_diameter']]

        if not isinstance(self.pattern['num_Ranvier'], list):
            self.pattern['num_Ranvier'] = [self.pattern['num_Ranvier']]

        if self.imp_name[-4:] == '.csv':

            # load the external array containing all neurons (coordinates of compartments) in the MRI space
            Array_coord_get = read_csv(os.environ['PATIENTDIR'] + '/' + self.imp_name, delimiter=' ', header=None)  #
            self.imp_coord_MRI = Array_coord_get.values
            Array_coord_in_MRI = self.imp_coord_MRI

            lst = [self.imp_name[:-4]]

        # in .h5 datasets with different morphologies can be stored
        elif self.imp_name[-3:] == '.h5':

            hf = h5py.File(os.environ['PATIENTDIR'] + '/' + self.imp_name, 'r')
            lst = list(hf.keys())  # names of datasets

        self.pattern['num_segments'] = np.zeros(len(lst), int)
        result_total = []  # appending might be too slow
        population_index = 0

        for i in lst:
            if self.imp_name[-3:] == '.h5':
                Array_coord_in_MRI = hf.get(i)  # one projection
                Array_coord_in_MRI = np.array(Array_coord_in_MRI)

            __, self.pattern['num_segments'][population_index], __ = self.get_neuron_morhology(
                self.pattern['fiber_diameter'][population_index],
                self.pattern['num_Ranvier'][population_index])  # add a consistency check

            if self.brain_model == 'Brain_substitute.brep':  # cut models if use a brain approximation

                n_models_before = int(Array_coord_in_MRI.shape[0] / self.pattern['num_segments'][population_index])
                logging.critical("Initial number of neuron models in {}: {}".format(str(i), n_models_before))
                points_outside = 0

                for inx in range(Array_coord_in_MRI.shape[0]):
                    pnt = Point(Array_coord_in_MRI[inx, 0], Array_coord_in_MRI[inx, 1], Array_coord_in_MRI[inx, 2])
                    if not (mesh_brain_sub.bounding_box_tree().compute_first_entity_collision(
                            pnt) < mesh_brain_sub.num_cells() * 100):  # if one point of the neural model is absent, the whole model is disabled
                        points_outside += 1
                        Array_coord_in_MRI, inx, __ = self.mark_to_remove(Array_coord_in_MRI, inx,
                                                                          nsegm=self.pattern['num_segments'][
                                                                              population_index])  # will also shift inx to the end of the neuron

                Array_coord_in_MRI = Array_coord_in_MRI[
                    ~np.all(Array_coord_in_MRI == -100000000.0, axis=1)]  # deletes all -10000000 enteries
                n_models_after = int(Array_coord_in_MRI.shape[0] / self.pattern['num_segments'][population_index])
                if int(n_models_before - n_models_after) != 0:
                    logging.critical(
                        "In {} {} models were outside of the approximating geometrical domain".format(str(i),
                                                                                                      n_models_before - n_models_after))

            Array_coord_in_MRI = np.round(Array_coord_in_MRI, 8)
            result_total.append(Array_coord_in_MRI)

            hf2 = h5py.File(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/Processed_neuron_array_MRI_space.h5', 'a')
            hf2.create_dataset(i, data=Array_coord_in_MRI)
            hf2.close()

            hf3 = h5py.File(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/All_neuron_models_by_populations.h5',
                            'a')  # resave to .h5
            dataset_imp_proc_coord_PO = Array_coord_in_MRI - self.MRI_first_voxel  # shift every dataset to PO
            hf3.create_dataset(i, data=dataset_imp_proc_coord_PO)
            hf3.close()

            population_index = population_index + 1

        # the downsamling is complete, reassign here
        if self.neuron_model == 'McIntyre2002_ds':
            self.neuron_model = 'McIntyre2002'

        if self.imp_name[-3:] == '.h5':
            hf.close()

        self.imp_proc_coord_MRI = np.concatenate(result_total)
        del result_total, Array_coord_in_MRI

        # shift to PO for allocaiton in the computational domain
        self.imp_proc_coord_PO = np.zeros((self.imp_proc_coord_MRI.shape[0], 3), float)
        self.imp_proc_coord_PO[:, 0] = self.imp_proc_coord_MRI[:, 0] - self.MRI_first_voxel[0]
        self.imp_proc_coord_PO[:, 1] = self.imp_proc_coord_MRI[:, 1] - self.MRI_first_voxel[1]
        self.imp_proc_coord_PO[:, 2] = self.imp_proc_coord_MRI[:, 2] - self.MRI_first_voxel[2]

        del self.imp_proc_coord_MRI  # it was saved (in h5)

        # for imported, will just save in a .csv file
        self.allocate_neurons_PO()

    def adjust_neuron_models(self, Domains, MRI_param):
        """ Adjust neuron arrays removing those that intersect with the electrode, encapsulation layer or CSF
            Stores the remaining neuron compartments (2D numpy array) in PO space
            in Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv and returns number of remaining models
            (int if VTA or custom, list if imported (per projection))
            If neuron array was imported stores projections separately in 'Neuron_model_arrays/Vert_of_Neural_model_NEURON_by_populations.h5'
            (if one projection or from .csv, stores just one dataset in .h5) """

        import os
        start_neuron_models = time.time()

        if self.Type != 'Imported':  # if the neuron array was created internally

            # or we can take it from the class: self.custom_coord_PO or self.VTA_coord_PO
            Array_coord_get = read_csv(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/All_neuron_models.csv',
                                       delimiter=' ', header=None)  #
            Array_coord = Array_coord_get.values

        elif self.Type == 'Imported':  # if the neuron array was provided

            hf = h5py.File(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/All_neuron_models_by_populations.h5', 'r')
            lst = list(hf.keys())
            List_of_arrays = []
            List_of_empty = []
            List_of_placed = []

            for i in lst:
                Array_coord = hf.get(i)
                Array_coord = np.array(Array_coord)

                if Array_coord.size == 0:
                    List_of_empty.append(i)
                    List_of_arrays.append(0)
                else:
                    List_of_arrays.append(Array_coord)
                    List_of_placed.append(Array_coord)

                del Array_coord
            hf.close()

        from CSF_refinement_new import get_CSF_voxels  # in the vicinity of the neuron array
        if self.Type != 'Imported':
            voxel_array_CSF_shifted = get_CSF_voxels(MRI_param, Array_coord, self.Type)
        else:
            voxel_array_CSF_shifted = get_CSF_voxels(MRI_param, List_of_placed, self.Type)

        points_csf, points_encap, points_outside = (0, 0, 0)

        # now we will folter out unphysiological neurins
        mesh = Mesh(os.environ['PATIENTDIR'] + "/Meshes/Mesh_unref.xml")
        subdomains_assigned = MeshFunction('size_t', mesh,
                                           os.environ['PATIENTDIR'] + "/Meshes/Mesh_unref_physical_region.xml")

        # also neurons should not pass through the encapsulation layer and floating contacts
        subdomains_enc = MeshFunction("size_t", mesh, mesh.topology().dim())
        subdomains_enc.set_all(0)

        #for i in range(len(Domains.Encup_index)):
        #    subdomains_enc.array()[
        #        subdomains_assigned.array() == Domains.Encup_index[i]] = 1  # neuron models cannot be located in CSF

        # encap at the electrode array only (defined by ROI)
        subdomains_enc.array()[subdomains_assigned.array() == Domains.Encup_index[1]] = 1  # neuron models cannot be located in encap

        # encap outside of ROI
        subdomains_enc_out = MeshFunction("size_t", mesh, mesh.topology().dim())
        subdomains_enc_out.set_all(0)
        subdomains_enc_out.array()[
            subdomains_assigned.array() == Domains.Encup_index[0]] = 1  # neuron models cannot be located in encap

        if isinstance(Domains.Float_contacts, list):
            for i in range(len(Domains.Float_contacts)):
                subdomains_enc.array()[subdomains_assigned.array() == Domains.Float_contacts[
                    i]] = 1  # neuron models cannot be located in floating conductors
        else:
            subdomains_enc.array()[
                subdomains_assigned.array() == Domains.Float_contacts] = 1  # neuron models cannot be located in floating conductors

        submesh_encup = SubMesh(mesh, subdomains_enc, 1)
        submesh_encup_out = SubMesh(mesh, subdomains_enc_out, 1)

        # if multiple pathways, those will be lists of lists
        self.neurons_idx_encap = []  # IMPORTANT: if the neuron does not intersect with encap. but with the electrode (sparse sampling), it will be treated as outside of the domain
        self.neurons_idx_csf = []

        inx = 0
        if self.Type != 'Imported':

            while inx < Array_coord.shape[0]:
                # print(inx)
                pnt = Point(Array_coord[inx, 0], Array_coord[inx, 1], Array_coord[inx, 2])

                if (submesh_encup.bounding_box_tree().compute_first_entity_collision(
                        pnt) < submesh_encup.num_cells()):  # this is a condition to check whether the point is inside encap. layer or floating conductor
                    points_encap = points_encap + 1
                    Array_coord, inx, i_neuron_encap = self.mark_to_remove(Array_coord,
                                                                           inx)  # will also shift inx to the end of the neuron
                    self.neurons_idx_encap.append(i_neuron_encap)
                else:
                    if (submesh_encup_out.bounding_box_tree().compute_first_entity_collision(
                            pnt) < submesh_encup_out.num_cells()):  # this is a condition to check whether the point is inside encap. layer or floating conductor
                        points_outside = points_outside + 1
                        Array_coord, inx, __ = self.mark_to_remove(Array_coord,inx)  # will also shift inx to the end of the neuron

                    elif not (mesh.bounding_box_tree().compute_first_entity_collision(
                            pnt) < mesh.num_cells() * 100):  # if one point of the neural model is absent, the whole model is disabled
                        points_outside = points_outside + 1
                        Array_coord, inx, __ = self.mark_to_remove(Array_coord,
                                                                   inx)  # will also shift inx to the end of the neuron
                    else:
                        # finally checks whether the neuron compartment is inside the CSF voxel
                        check1_1 = (voxel_array_CSF_shifted[:, 0] - Array_coord[inx, 0] <= MRI_param.voxel_dims[0])
                        check1_2 = (voxel_array_CSF_shifted[:, 1] - Array_coord[inx, 1] <= MRI_param.voxel_dims[1])
                        check1_3 = (voxel_array_CSF_shifted[:, 2] - Array_coord[inx, 2] <= MRI_param.voxel_dims[2])
                        check2_1 = (voxel_array_CSF_shifted[:, 0] >= Array_coord[
                            inx, 0])  # we could just check the sign
                        check2_2 = (voxel_array_CSF_shifted[:, 1] >= Array_coord[inx, 1])
                        check2_3 = (voxel_array_CSF_shifted[:, 2] >= Array_coord[inx, 2])

                        check3 = np.logical_and(np.logical_and(check1_1, check2_1),
                                                np.logical_and(np.logical_and(check1_2, check2_2),
                                                               np.logical_and(check1_3, check2_3)))
                        a = np.where((check3 == (True)))
                        if str(a) != '(array([], dtype=int64),)':
                            points_csf = points_csf + 1

                            Array_coord, inx, i_neuron_csf = self.mark_to_remove(Array_coord,
                                                                                 inx)  # will also shift inx to the end of the neuron
                            self.neurons_idx_csf.append(i_neuron_csf)
                        else:
                            inx += 1

            logging.critical(
                "Neurons in CSF, encapsulation layer (and floating conductors) and outside (and intersecting with the electrode): {0}, {1}, {2}".format(
                    points_csf, points_encap, points_outside))
            inx = 0

            del voxel_array_CSF_shifted

            Array_coord = Array_coord[~np.all(Array_coord == -100000000.0,
                                              axis=1)]  # deletes all marked (in self.mark_to_remove with -100000000.0) enteries
            Array_coord = np.round(Array_coord, 8)

            np.savetxt(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', Array_coord,
                       delimiter=" ")
            logging.critical(
                "Adjusted neuron models can be visualized from Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv in Paraview")

            N_models = int(Array_coord.shape[0] / self.pattern['num_segments'])

            if len(self.neurons_idx_csf) > 0.25 * N_models:
                logging.critical("!========================================================!")
                logging.critical("WARNING: too many neurons are in CSF")
                logging.critical("Check segmask and image normalizations")
                logging.critical("!========================================================!")

            np.savetxt(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/Adjusted_neuron_array_info.csv',
                       np.array([N_models, points_csf, points_encap, points_outside]), delimiter=" ")
            logging.critical("Number of placed neuron models: {}".format(N_models))

        if self.Type == 'Imported':

            list_of_connections = []

            # clean-up
            if os.path.isfile(
                    os.environ['PATIENTDIR'] + '/Neuron_model_arrays/Vert_of_Neural_model_NEURON_by_populations.h5'):
                os.remove(
                    os.environ['PATIENTDIR'] + '/Neuron_model_arrays/Vert_of_Neural_model_NEURON_by_populations.h5')

            number_of_points_filtered = 0
            N_models = np.zeros(len(List_of_arrays), int)  # each entry will contain number of models per projection

            for i in range(len(List_of_arrays)):

                sublist_idx_encap = []
                sublist_idx_csf = []

                if not (type(List_of_arrays[i]) is np.ndarray):  # check for empty
                    hf3 = h5py.File(
                        os.environ['PATIENTDIR'] + '/Neuron_model_arrays/Vert_of_Neural_model_NEURON_by_populations.h5',
                        'a')
                    hf3.create_dataset(lst[i], data=0)
                    hf3.close()

                    self.neurons_idx_csf.append([-1])
                    self.neurons_idx_encap.append([-1])

                else:
                    Array_coord = List_of_arrays[i]
                    while inx < Array_coord.shape[0]:

                        pnt = Point(Array_coord[inx, 0], Array_coord[inx, 1], Array_coord[inx, 2])

                        if (submesh_encup.bounding_box_tree().compute_first_entity_collision(
                                pnt) < mesh.num_cells()):  # this is a condition to check whether the point is inside encap. layer or floating conductor
                            points_encap += 1
                            Array_coord, inx, i_neuron_encap = self.mark_to_remove(Array_coord, inx,
                                                                                   nsegm=self.pattern['num_segments'][
                                                                                       i])  # will also shift inx to the end of the neuron
                            sublist_idx_encap.append(i_neuron_encap)
                        else:
                            if not (mesh.bounding_box_tree().compute_first_entity_collision(
                                    pnt) < mesh.num_cells() * 100):  # if one point of the neural model is absent, the whole model is disabled
                                points_outside += 1
                                Array_coord, inx, __ = self.mark_to_remove(Array_coord, inx,
                                                                           nsegm=self.pattern['num_segments'][i])
                            else:  # finally checks whether the neuron compartment is inside the CSF voxel
                                check1_1 = (voxel_array_CSF_shifted[:, 0] - Array_coord[inx, 0] <= MRI_param.voxel_dims[0])
                                check1_2 = (voxel_array_CSF_shifted[:, 1] - Array_coord[inx, 1] <= MRI_param.voxel_dims[1])
                                check1_3 = (voxel_array_CSF_shifted[:, 2] - Array_coord[inx, 2] <= MRI_param.voxel_dims[2])
                                check2_1 = (voxel_array_CSF_shifted[:, 0] >= Array_coord[inx, 0])
                                check2_2 = (voxel_array_CSF_shifted[:, 1] >= Array_coord[inx, 1])
                                check2_3 = (voxel_array_CSF_shifted[:, 2] >= Array_coord[inx, 2])
                                # print(check1_1)

                                check3 = np.logical_and(np.logical_and(check1_1, check2_1),
                                                        np.logical_and(np.logical_and(check1_2, check2_2),
                                                                       np.logical_and(check1_3, check2_3)))
                                a = np.where((check3 == (True)))
                                if str(a) != '(array([], dtype=int64),)':
                                    points_csf += 1
                                    Array_coord, inx, i_neuron_csf = self.mark_to_remove(Array_coord, inx,
                                                                                         nsegm=
                                                                                         self.pattern['num_segments'][
                                                                                             i])
                                    sublist_idx_csf.append(i_neuron_csf)
                                else:
                                    inx += 1

                    Array_coord = Array_coord[~np.all(Array_coord == -100000000.0,
                                                      axis=1)]  # deletes all marked (in self.mark_to_remove with -100000000.0) enteries

                    inx = 0

                    Array_coord = np.round(Array_coord, 8)
                    N_models[i] = int(Array_coord.shape[0] / self.pattern['num_segments'][i])

                    hf3 = h5py.File(
                        os.environ['PATIENTDIR'] + '/Neuron_model_arrays/Vert_of_Neural_model_NEURON_by_populations.h5',
                        'a')
                    hf3.create_dataset(lst[i], data=Array_coord)
                    hf3.close()

                    np.savetxt(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/' + lst[i] + '.csv', Array_coord,
                               delimiter=" ")
                    list_of_connections.append(lst[i])

                    number_of_points_filtered = number_of_points_filtered + Array_coord.shape[0]

                    self.neurons_idx_csf.append(sublist_idx_csf)
                    self.neurons_idx_encap.append(sublist_idx_encap)

                    if len(sublist_idx_csf) > 0.25 * N_models[i]:
                        logging.critical("!========================================================!")
                        logging.critical("WARNING: too many neurons of {} are in CSF".format(lst[i]))
                        logging.critical("Check segmask and image normalizations")
                        logging.critical("!========================================================!")

            hf = h5py.File(
                os.environ['PATIENTDIR'] + '/Neuron_model_arrays/Vert_of_Neural_model_NEURON_by_populations.h5', 'r')
            lst = list(hf.keys())
            result_total = []
            for i in lst:
                if not (i in List_of_empty):
                    # print(i)
                    a = hf.get(i)
                    a = np.array(a)
                    result_total.append(a)

            Array_coord_total = np.concatenate(result_total)
            hf.close()

            np.savetxt(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/Adjusted_neuron_array_info.csv', N_models,
                       delimiter=" ")
            np.savetxt(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv',
                       Array_coord_total, delimiter=" ")

            logging.critical(
                "Neurons in CSF, encapsulation layer (and floating conductors) and outside (and intersecting with the electrode): {0}, {1}, {2}".format(
                    points_csf, points_encap, points_outside))
            logging.critical(
                "Adjusted neuron models can be visualized from Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv in Paraview")

            # might not work like this
            # from Visualization_files.Paraview_connections_processed import show_connections
            # show_connections(list_of_connections)

            logging.critical("Number of placed neuron models per population: {}".format(N_models))
            del Array_coord_total, voxel_array_CSF_shifted

        if self.Type == 'Imported':
            del List_of_arrays
        del Array_coord

        minutes = int((time.time() - start_neuron_models) / 60)
        secnds = int(time.time() - start_neuron_models) - minutes * 60
        logging.critical("----- Adjustment of the neuron models took {} min {} sec -----\n".format(minutes, secnds))

        self.N_models = N_models
