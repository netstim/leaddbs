"""
@author: Konstantin Butenko

This script uses fiber trajectories to allocate axon models

"""

import numpy as np
import os
import h5py
import sys
import json
import scipy
import math
from scipy.io import savemat

# adopted from nibabel
from nibabel_SequenceArray import ArraySequence

class AxonModels:

    """ Model to represent axons for simulation in OSS-DBS

    """

    def __init__(self, stim_dir, hemis_idx, description_file):

        """

        Parameters
        ----------

        stim_directory: str, full path to the stimulation or result folder where allocated axons are stored
        hemis_idx: int, hemisphere ID (0 - right, 1 - left)
        description_file: str, full path to oss-dbs_parameters.mat or a .json file that contains the following parameters:

             pathway_mat_file: list of full paths to pathways files in lead-dbs format (could be just one)
             axon_diams_all: list of diameters in micrometers for all provided fibers, one per pathway
             axon_lengths_all: list of axon lengths in mm, one per pathway
             centering_coordinates: list of lists, 3-D coordinates used to center axons on fibers (e.g. active contacts)
             axon_model: str, NEURON model ('MRG2002', 'MRG2002_DS' (downsampled), 'McNeal1976' (classic McNeal's))
             combined_h5_file: str, full path to the file where axons are stored
             projection_names: list of str, optional
             connectome_name: str, optional

        Notes
        -----
        oss-dbs_parameters.mat is created via Lead-DBS. For .json parameters, see _import_custom_neurons()

        """

        self.description_file = description_file

        # make sure we are in ext_libs/OSS-DBS/Axon_Processing
        abspath = os.path.abspath(__file__)
        dname = os.path.dirname(abspath)
        os.chdir(dname)

        # not necessary to use env. variable, but just for consistency for now
        os.environ['STIMDIR'] = stim_dir

        # Lead-DBS input
        if self.description_file[-22:] == 'oss-dbs_parameters.mat':
            self._import_leaddbs_neurons(int(hemis_idx))
        elif self.description_file[-4:] == 'json':
            self.projection_names = None
            self.connectome_name = 'MyTracts'
            self._import_custom_neurons()
        else:
            print("Unsupported input format, see AxonModels in Axon_allocation.py")

    def _import_leaddbs_neurons(self, hemis_idx):

        """ Import Lead-DBS description for axon models from oss-dbs_parameters.mat

        Parameters
        ----------
        hemis_idx: int, hemisphere ID (0 - right, 1 - left)

        """

        # load .mat of different versions (WON'T WORK THIS WAY ATM!)
        try:
            file_inp = h5py.File(self.description_file, mode='r')
        except:
            print("Please, save oss-dbs_parameters using 'save(oss-dbs_parameters_path, 'settings', '-v7.3')' ")
            raise SystemExit
            file_inp = scipy.io.loadmat(self.description_file)


        # try to read from .mat
        if 'neuronModel' in file_inp['settings']:
            array_ascii = file_inp['settings']['neuronModel'][:]
            list_ascii = []
            for i in range(array_ascii.shape[0]):
                list_ascii.append(array_ascii[i][0])
            # list_ascii = map(lambda s: s.strip(), list_ascii)
            self.axon_model = ''.join(chr(i) for i in list_ascii)
            if self.axon_model not in ['MRG2002', 'MRG2002_DS', 'McNeal1976']:
                print('The selected NEURON models is not recognized, check oss-dbs_parameters.mat')
                raise SystemExit
        else:
            # Assume Reilly's (McNeals1976) by default
            self.axon_model = 'McNeal1976'

        #HEMIS_OUTPUT_PATHS = ["Results_rh", "Results_lh"]
        #output_folder = os.path.join(os.environ['STIMDIR'], HEMIS_OUTPUT_PATHS[hemis_idx])


        # connectome name within Lead-DBS (e.g. 'Multi-Tract: PetersenLUIC')
        array_ascii = file_inp['settings']['connectome'][:]
        list_ascii = []
        for i in range(array_ascii.shape[0]):
            list_ascii.append(array_ascii[i][0])
        # list_ascii = map(lambda s: s.strip(), list_ascii)
        self.connectome_name = ''.join(chr(i) for i in list_ascii)


        # 'Multi-tract' connectomes contain multiple pathways (projections) in separate .mat files
        if 'Multi-Tract' in self.connectome_name:
            # this file is pre-filtered connectome assembled in one file in Lead-DBS
            self.pathway_mat_file = [
                os.environ['STIMDIR'] + '/' + self.connectome_name.rsplit(' ', 1)[1] + '/data' + str(hemis_idx + 1) + '.mat']

            self.projection_names = []
            for i in range(len(file_inp['settings']['connectomeTractNames'][0])):
                ext_string = file_inp[file_inp['settings']['connectomeTractNames'][0][i]]
                list_ascii = []
                for j in range(ext_string.shape[0]):
                    list_ascii.append(ext_string[j][0])
                projection_name = ''.join(chr(i) for i in list_ascii)
                self.projection_names.append(projection_name)
        else:
            self.projection_names = ['default']
            # this file is pre-filtered connectome in Lead-DBS
            self.pathway_mat_file = [os.environ['STIMDIR'] + '/' + self.connectome_name + '/data' + str(hemis_idx + 1) + '.mat']


        # check which contacts are active to seed axons close to them
        # for StimSets check across all of them
        stimSets = bool(file_inp['settings']['stimSetMode'][0][0])  # if StimSets are used, create a dummy ampl_vector
        if stimSets:
            stim_protocols = np.genfromtxt(os.environ['STIMDIR'] + '/Current_protocols_' + str(hemis_idx) + '.csv',
                                           dtype=float, delimiter=',', names=True)

            total_contacts = len(list(stim_protocols[0]))
            total_protocols = stim_protocols.shape[0]

            protocols_array = np.zeros((total_protocols, total_contacts), float)
            ampl_vector = list(stim_protocols[0])  # just initialize

            for j in range(total_protocols):
                protocols_array[j, :] = list(stim_protocols[j])
                for i in range(total_contacts):
                    if not math.isnan(protocols_array[j, i]):
                        ampl_vector[i] = 1.0  # you do not need a value, just substitute NaN
        else:
            ampl_vector = list(file_inp['settings']['Phi_vector'][:, hemis_idx])

        self.centering_coordinates = []
        for i in range(len(ampl_vector)):
            if not (math.isnan(ampl_vector[i])):
                a_ref = file_inp['settings']['contactLocation'][hemis_idx][0]
                b = file_inp[a_ref]
                self.centering_coordinates.append(b[:, i])

        # hardcoded name for axons pre-filtered by Lead-DBS
        self.combined_h5_file = os.environ['STIMDIR'] + '/Allocated_axons'
        self.output_directory = self.combined_h5_file.rsplit('/', 1)[0]
        print("output_dictionary: ", self.output_directory)

        # morphology set in Lead-DBS
        self.axon_lengths_all = list(file_inp['settings']['axonLength'][:][0][:])
        self.axon_diams_all = list(file_inp['settings']['fiberDiameter'][:][0][:])

    def _import_custom_neurons(self):

        """ Import custom description for axon models from a .json dictionary

        """

        ## Example json input

        # custom_dict = {
        #     'pathway_mat_file': ['/home/konstantin/Documents/GitHub/leaddbs/connectomes/dMRI_MultiTract/Petersen/SMA_hdp_left.mat',
        #                          '/home/konstantin/Documents/GitHub/leaddbs/connectomes/dMRI_MultiTract/Petersen/SMA_hdp_right.mat',
        #                          '/home/konstantin/Documents/GitHub/leaddbs/connectomes/dMRI_MultiTract/Petersen/gpe2stn_sm_left.mat'],
        #
        #     # axon diameter and length is the same for all axons within the pathway
        #     'axon_diams_all': [5.7,5.7,3.0],
        #     'axon_lengths_all':[20.0,20.0,10.0],
        #
        #     # in this case, we just have some STN coordinates for left and right in MNI
        #     'centering_coordinates': [[7.5838, -18.3984, 1.8932],[-7.5838, -18.3984, 1.8932]],
        #     'axon_model': 'McNeal1976',
        #     'combined_h5_file': '/home/konstantin/Documents/dataset/all_tracts'
        # }

        with open(self.description_file, 'r') as fp:
            custom_dict = json.load(fp)

        self.pathway_mat_file = custom_dict['pathway_mat_file']
        self.axon_diams_all = custom_dict['axon_diams_all']
        self.axon_lengths_all = custom_dict['axon_lengths_all']
        self.centering_coordinates = custom_dict['centering_coordinates']
        self.axon_model = custom_dict['axon_model']
        self.combined_h5_file = custom_dict['combined_h5_file']
        # strip extention if provided
        if self.combined_h5_file[-3:] == '.h5':
            self.combined_h5_file = self.combined_h5_file[:-3]
        self.output_directory = self.combined_h5_file.rsplit('/', 1)[0]

        if "projection_names" in custom_dict:
            self.projection_names = custom_dict['projection_names']
        else:
            self.projection_names = [pathway_file.rsplit('/', 1)[1][0:-4] for pathway_file in self.pathway_mat_file]

        if "connectome_name" in custom_dict:
            self.connectome_name = custom_dict['connectome_name']


    def convert_fibers_to_axons(self):

        """ Seed axons iterating over all pathways

        """

        # within a projection (pathway), number of nodes of Ranvier per axon is fixed
        n_Ranviers_per_projection_all = np.zeros(len(self.axon_lengths_all), int)
        n_Neurons_all = np.zeros(len(self.axon_lengths_all), int)
        orig_n_Neurons_all = np.zeros(len(self.axon_lengths_all), int)

        # iterate over projections (fibers) and seed axons
        for i in range(len(self.axon_diams_all)):

            # various geometric parameters for a single axon
            axon_morphology = get_axon_morphology(self.axon_model, self.axon_diams_all[i], self.axon_lengths_all[i])

            # multiple .mat files (manual input)
            if len(self.pathway_mat_file) > 1:
                n_Ranviers_per_projection_all[i], n_Neurons_all[i], orig_n_Neurons_all[i] = self.deploy_axons_fibers(self.pathway_mat_file[i], self.projection_names[i],axon_morphology,False)

            # multiple pathways in one .mat file (Lead-DBS dMRI_MultiTract connectome)
            elif 'Multi-Tract' in self.connectome_name:
                n_Ranviers_per_projection_all[i], n_Neurons_all[i], orig_n_Neurons_all[i] = self.deploy_axons_fibers(self.pathway_mat_file[0], self.projection_names[i],axon_morphology,True)

            # one .mat file without pathway differentiation (Lead-DBS dMRI connectome)
            else:
                n_Ranviers_per_projection_all[i], n_Neurons_all[i], orig_n_Neurons_all[i] = self.deploy_axons_fibers(self.pathway_mat_file[0], self.projection_names[i],
                                                               axon_morphology,False)

            print(n_Neurons_all[i], " axons seeded for ", self.projection_names[i], " with ", n_Ranviers_per_projection_all[i], " nodes of Ranvier\n")

        # only add axon diameters for seeded axons
        self.axon_diams = []
        n_Ranviers_per_projection = []
        n_Neurons = []
        orig_n_Neurons = []
        for i in range(len(self.axon_diams_all)):
            if n_Ranviers_per_projection_all[i] != 0:
                self.axon_diams.append(float(self.axon_diams_all[i]))
                n_Ranviers_per_projection.append(int(n_Ranviers_per_projection_all[i]))
                n_Neurons.append(int(n_Neurons_all[i]))
                orig_n_Neurons.append(int(orig_n_Neurons_all[i]))

        self._save_axon_parameters_in_json(n_Ranviers_per_projection, n_Neurons, orig_n_Neurons)

    def _save_axon_parameters_in_json(self, n_Ranviers_per_projection, n_Neurons, orig_n_Neurons):

        """ Save minimally required axon description in a .json file

        Parameters
        ----------
        n_Ranviers_per_projection: list, number of nodes of Ranvier for axons of each pathway (one entry per pathway)
        n_Neurons: list, number of neurons seeded per pathway
        orig_n_Neurons: list, number of neurons per pathway as defined in the connectome (before Kuncel pre-filtering)
        """

        # dictionary to store axon parameters
        axon_dict = {
            'n_Ranvier': n_Ranviers_per_projection,
            'axon_diams': self.axon_diams,
            'Axon_Model_Type': self.axon_model,
            'Name_prepared_neuron_array': self.combined_h5_file + '.h5',
            'Neuron_model_array_prepared': True,
            'N_seeded_neurons': n_Neurons,
            'N_orig_neurons': orig_n_Neurons,
            'connectome_name': self.connectome_name
        }

        with open(self.output_directory + '/Allocated_axons_parameters.json', 'w') as save_as_dict:
            json.dump(axon_dict, save_as_dict)

        # np.savetxt(name_of_directory + '/' + self.combined_h5_file + '_N_nodes.csv', n_Ranviers_per_projection,
        #                delimiter=" ")

    def deploy_axons_fibers(self, pathway_file, projection_name, axon_morphology, multiple_projections_per_file=False):

        """ Convert streamlines (fibers) to axons and store in OSS-DBS supported format

        Parameters
        ----------
        pathway_file,: str, full path to .mat file containing fiber descriptions (Lead-DBS format)
        projection_name: str, pathway name
        axon_morphology: dict, geometric description of a single axon, see get_axon_morphology
        multiple_projections_per_file: bool, optional, flag if pathway_file contains multiple pathways

        Returns
        -------
        int, number of nodes of Ranvier for axons in this pathway. Returns 0 if failed to see (fiber is too short)
        int, number of axons seeded for the pathway

        Notes
        -----
        Pathways are stored as separate groups in the specified .h5 file. Axons are stored in separate 2-D datasets.
        For Paraview visualization, use axon_array_2D_<projection_name>
        For Lead-DBS visualization, use <projection_name>_axons.mat

        """

        # fallback for non hdf5, TBD!
        try:
            file = h5py.File(pathway_file, mode='r')
        except:
            file = scipy.io.loadmat(pathway_file)

        if multiple_projections_per_file == False:
            # fiber_array has 4 columns (x,y,z,fiber_index), raws - all points
            fiber_array = file['fibers'][:]
        else:
            fiber_array = file[projection_name]['fibers'][:]

        if fiber_array.ndim == 1:
            print(projection_name, 'projection is empty, check settings for fib. diameter and axon length')
            return 0,0,0  # no nodes were seeded
        else:
            
            # flip check
            if fiber_array.shape[1] == 4 and fiber_array.shape[0] != 4:
                fiber_array = fiber_array.T
                idx_shape_inx = 0
            else:
                idx_shape_inx = 1
            
            if multiple_projections_per_file == False:
                if 'origNum' in file:
                    orig_N_fibers = int(file['origNum'][0][0])
                else:
                    orig_N_fibers = int(file['idx'][:].shape[idx_shape_inx])
            else:
                if 'origNum' in file[projection_name]:
                    orig_N_fibers = int(file[projection_name]['origNum'][0][0])
                else:
                    orig_N_fibers = int(file[projection_name]['idx'][:].shape[idx_shape_inx])

        # covert fiber table to nibabel streamlines
        streamlines = convert_fibers_to_streamlines(fiber_array)

        # resample streamlines to nodes of Ranvier
        streamlines_resampled, excluded_streamlines_idx = resample_fibers_to_Ranviers(streamlines, axon_morphology)

        # truncate streamlines to match selected axon length
        # axons are seeded on the segment closest to active contacts or other ROI, see self.centering_coordinates
        streamlines_axons = place_axons_on_streamlines(streamlines_resampled, axon_morphology, self.centering_coordinates)

        # streamlines_axons already contain the position of Ranvier nodes. Now we get internodal compartments
        # and store all coordinates in a 3D array: compartment index, spatial axis, axon index
        axon_array = np.zeros((axon_morphology['n_segments'], 3, len(streamlines_axons)), float)

        # 2-D version for Paraview visualization
        axon_array_2D = np.zeros((axon_morphology['n_segments'] * len(streamlines_axons), 4), float)

        # save axons as separate datasets within groups that correspond to pathways
        hf = h5py.File(self.combined_h5_file + '.h5', 'a')
        g = hf.create_group(projection_name)

        # get local coordinates for internodal compartments
        local_comp_coords = get_local_compartment_coords(axon_morphology)
        glob_ind = 0
        for inx_axn in range(len(streamlines_axons)):
            inx_comp = 0
            for inx in range(axon_morphology['n_Ranviers'] - 1):

                # compartments are seeded along the internodal vector
                internodal_vector_normalized = normalized(streamlines_axons[inx_axn][inx+1] - streamlines_axons[inx_axn][inx])

                # positions for nodes of Ranvier are known
                axon_array[inx_comp, :, inx_axn] = streamlines_axons[inx_axn][inx]

                # now place the compartments until the next node
                for loc_comp_inx in range(1,axon_morphology['n_comp']):
                    axon_array[inx_comp + loc_comp_inx, :, inx_axn] = axon_array[inx_comp,:, inx_axn] + local_comp_coords[loc_comp_inx-1] * \
                                                                  internodal_vector_normalized[0][:]

                inx_comp = inx_comp + axon_morphology['n_comp']

            # last node of Ranview
            axon_array[-1, :, inx_axn] = streamlines_axons[inx_axn][-1]

            axon_array_2D[glob_ind:glob_ind + axon_morphology['n_segments'], :3] = axon_array[:, :, inx_axn]
            axon_array_2D[glob_ind:glob_ind + axon_morphology['n_segments'], 3] = inx_axn + 1  # because in Matlab they start from 1

            g.create_dataset('axon' + str(inx_axn), data=axon_array[:, :, inx_axn])

            glob_ind = glob_ind + axon_morphology['n_segments']

        np.savetxt(self.output_directory + '/' + 'axon_array_2D_' + projection_name + '.csv', axon_array_2D,
                   delimiter=" ")

        hf.close()

        mdic = {"fibers": axon_array_2D, "ea_fibformat": "1.0"}
        savemat(self.combined_h5_file + '_' + projection_name + "_axons.mat", mdic)

        return axon_morphology['n_Ranviers'], len(streamlines_axons), orig_N_fibers

def convert_fibers_to_streamlines(fibers):
    """ Convert Lead-DBS fibers to Nibabel streamlines

    Parameters
    ----------
    fibers,: fiber descriptions (Lead-DBS format)

    Returns
    -------
    list, streamlines stored as ArraySequence(), i.e. list that describes each fiber in a sublist

    """

    from nibabel_SequenceArray import ArraySequence
    streamlines = ArraySequence()

    # yes, indexing starts with one in those .mat files
    N_streamlines = int(fibers[3, :].max())

    k = 0
    i_previous = 0
    for i in range(N_streamlines):
        loc_counter = 0
        while ((i + 1) == fibers[3, k]):  # this is not optimal, you need to extract a pack by np.count?
            k += 1
            loc_counter += 1
            if (k == fibers[3, :].shape[0]):
                break

        stream_line = fibers[:3, i_previous:i_previous + loc_counter].T
        i_previous = k
        streamlines.append(stream_line)

    return streamlines

def get_axon_morphology(axon_model, axon_diam, axon_length=None, n_Ranviers=None):

    """ Get geometric description of a single axon

    Parameters
    ----------
     axon_model: str, NEURON model ('MRG2002', 'MRG2002_DS' (downsampled), 'McNeal1976' (classic McNeal's))
     axon_diam: float, diameter in micrometers for all fibers in the pathway
     axon_length: float, optional, axon lengths in mm for all fibers in the pathway. If not specified, provide n_Ranviers
     n_Ranviers: int, optional, number of nodes of Ranvier per axon. If not specified, provide axon_length.

    Returns
    -------
    dict

    """

    axon_morphology = {
        'axon_model': axon_model,
        'axon_diam': axon_diam,
    }

    if 'MRG2002' in axon_model:

        from Axon_files.axon import Axon
        param_ax = {
            'centered': True,
            'diameter': axon_diam
        }
        a = Axon(param_ax)
        nr = Axon.get_axonparams(a)

        axon_morphology['ranvier_length'], axon_morphology['para1_length'], axon_morphology['para2_length'], axon_morphology['node_step'] = (
            nr["ranvier_length"] * 1e-3, nr["para1_length"] * 1e-3, nr["para2_length"] * 1e-3, nr["deltax"] * 1e-3)

        if axon_model == 'MRG2002_DS':
            # downsampled version
            if axon_diam >= 5.7:
                # node -- -- internodal -- -- -- -- internodal -- -- node
                axon_morphology['n_comp'] = 3
                axon_morphology['inter_length'] = (axon_morphology['node_step']  - axon_morphology['para1_length'] * 2 -
                                                   axon_morphology['para2_length'] * 2) / 6
            else:
                # node -- -- -- internodal -- -- -- node
                axon_morphology['n_comp'] = 2
                axon_morphology['inter_length'] = (axon_morphology['node_step']  - axon_morphology['para1_length'] * 2 -
                                                   axon_morphology['para2_length'] * 2) / 3
        else:
            axon_morphology['n_comp'] = int(((nr["ranvier_nodes"] - 1) + nr["inter_nodes"] + nr["para1_nodes"] + nr["para2_nodes"]) / (
                    nr["ranvier_nodes"] - 1))

            if axon_diam >= 5.7:
                axon_morphology['inter_length'] = (axon_morphology['node_step'] - axon_morphology['para1_length'] * 2 -
                                                   axon_morphology['para2_length'] * 2) / 6
            else:
                axon_morphology['inter_length'] = (axon_morphology['node_step'] - axon_morphology['para1_length'] * 2 -
                                                   axon_morphology['para2_length'] * 2) / 3

        # check what was provided, axon_length takes precedence
        if axon_length is not None:
            axon_morphology['axon_length'] = axon_length
            axon_morphology['n_Ranviers'] = int(axon_length / axon_morphology['node_step'])
        else:
            axon_morphology['n_Ranviers'] = n_Ranviers
            axon_morphology['axon_length'] = n_Ranviers * axon_morphology['node_step']

        # always odd number of nodes of Ranvier!
        if axon_morphology['n_Ranviers'] % 2 == 0:
            axon_morphology['n_Ranviers'] -= 1
            axon_morphology['axon_length'] = axon_morphology['n_Ranviers'] * axon_morphology['node_step']

        axon_morphology['n_para1']=  nr["para1_nodes"]*(axon_morphology['n_Ranviers'] - 1) / (21 - 1)
        axon_morphology['n_para2'] = nr["para2_nodes"] * (axon_morphology['n_Ranviers'] - 1) / (21 - 1)

        axon_morphology['n_segments'] = int((axon_morphology['n_Ranviers'] - 1) * axon_morphology['n_comp'] + 1)  # overall number of points on Axon incl. internodal
        #axon_morphology['n_total'] = (axon_morphology['n_Ranviers'] - 1) * axon_morphology['n_comp'] + 1  # total incl. Ranvier

        # additional params for NEURON model, see axon.py
        axon_morphology['axon_d'], axon_morphology['node_d'], axon_morphology['para1_d'], axon_morphology['para2_d'], axon_morphology['lamellas'] = (nr["axon_diameter"],nr["node_diameter"],nr["para1_diameter"],nr["para2_diameter"],nr["lamellas"])

    elif axon_model == 'McNeal1976':
        # node -- -- -- internodal -- -- -- node
        axon_morphology['n_comp'] = 2  # only nodes and one internodal per segment
        axon_morphology['node_step'] = axon_diam * 0.2  # from 1 to 2 mm
        # check what was provided, axon_length takes precedence
        if axon_length is not None:
            axon_morphology['axon_length'] = axon_length
            axon_morphology['n_Ranviers'] = int(axon_length / axon_morphology['node_step'])
        else:
            axon_morphology['n_Ranviers'] = n_Ranviers
            axon_morphology['axon_length'] = n_Ranviers * axon_morphology['node_step']

        # always odd number of nodes of Ranvier!
        if axon_morphology['n_Ranviers'] % 2 == 0:
            axon_morphology['n_Ranviers'] -= 1
            axon_morphology['axon_length'] = axon_morphology['n_Ranviers'] * axon_morphology['node_step']

        axon_morphology['n_segments'] = int((axon_morphology['n_Ranviers'] - 1) * axon_morphology['n_comp'] + 1)
        #axon_morphology['n_total'] = (axon_morphology['n_Ranviers'] - 1) * 2 + 1  # one internodal per segment + the last Ranvier
    else:
        print("The neuron model is not implemented")
        raise SystemExit

    return axon_morphology


def normalized(vector, axis=-1, order=2):
    """ Get L2 norm of a vector

    """

    l2 = np.atleast_1d(np.linalg.norm(vector, order, axis))
    l2[l2 == 0] = 1
    return vector / np.expand_dims(l2, axis)

def place_axons_on_streamlines(streamlines_resampled, axon_morphology, centering_coordinates):

    """ Allocate axons on the streamlines with seeding points defined by centering_coordinates

    Parameters
    ----------
     streamlines_resampled: list, streamlines sampled by nodes of Ranvier, stored as ArraySequence()
     axon_morphology: dict, geometric description of a single axon, see get_axon_morphology
     centering_coordinates: list of lists, 3-D coordinates used to center axons on fibers (e.g. active contacts)

    Returns
    -------
    list, axons (truncated streamlines), stored as ArraySequence()

    """

    from scipy import spatial
    axons_ROI_centered = ArraySequence()

    for inx_axn in range(len(streamlines_resampled)):

        single_streamline_ROI_centered = np.zeros((axon_morphology['n_Ranviers'], 3), float)

        A = streamlines_resampled[inx_axn]
        distance_list = []
        index_list = []
        for j in range(len(centering_coordinates)):
            distance, index = spatial.KDTree(A).query(
                centering_coordinates[j])  # distance is a local index of closest node of Ranvier on the axon
            distance_list.append(distance)
            index_list.append(index)

        index = index_list[distance_list.index(min(distance_list))]  # index of the closest point as assigned as index

        loc_index = 0
        # choose where to start seeding the axon
        if index < int(axon_morphology['n_Ranviers'] / 2):
            # axon---fiber---fiber---fiber---fiber---#
            for i in range(0, int(axon_morphology['n_Ranviers'])):
                single_streamline_ROI_centered[loc_index, :] = A[i]
                loc_index += 1
        elif index + int(axon_morphology['n_Ranviers'] / 2) + 1 > A.shape[0]:
            # fiber---fiber---fiber---fiber---axon---#
            for i in range(A.shape[0] - axon_morphology['n_Ranviers'], A.shape[0]):
                single_streamline_ROI_centered[loc_index, :] = A[i]
                loc_index += 1
        else:
            # ---fiber---fiber---axon---fiber---fiber---#
            if axon_morphology['n_Ranviers'] % 2 == 0:
                for i in range(index - int(axon_morphology['n_Ranviers'] / 2),
                               index + int(axon_morphology['n_Ranviers'] / 2)):
                    single_streamline_ROI_centered[loc_index, :] = A[i]
                    loc_index += 1
            else:
                for i in range(index - int(axon_morphology['n_Ranviers'] / 2),
                               index + int(axon_morphology['n_Ranviers'] / 2) + 1):
                    single_streamline_ROI_centered[loc_index, :] = A[i]
                    loc_index += 1

        axons_ROI_centered.append(single_streamline_ROI_centered)

    if len(axons_ROI_centered) != len(streamlines_resampled):
        print("Failed to sample some axons!")
        raise SystemExit

    return axons_ROI_centered

def resample_fibers_to_Ranviers(streamlines, axon_morphology):

    """ Get streamlines resampled by nodes of Ranvier for a specific axonal morphology

    Parameters
    ----------
     streamlines: list, arbitrary sampled streamlines, stored as ArraySequence()
     axon_morphology: dict, geometric description of a single axon, see get_axon_morphology

    Returns
    -------
    list, resampled streamlines, stored as ArraySequence()

    """

    # resampling to nodes of Ranvier for arbitrary fiber length
    from Arbitrary_streamline_to_Ranviers import length_fiber
    lengths_streamlines_filtered = list(length_fiber(streamlines))
    streamlines_resampled = ArraySequence()

    from Arbitrary_streamline_to_Ranviers import resample_streamline_for_Ranvier
    excluded_streamlines = []
    #total_points = 0
    for streamline_index in range(len(lengths_streamlines_filtered)):
        n_Ranvier_this_axon = int(lengths_streamlines_filtered[streamline_index] / axon_morphology['node_step'])
        streamline_resampled = resample_streamline_for_Ranvier(streamlines[streamline_index],
                                                               n_Ranvier_this_axon * axon_morphology['node_step'],
                                                               n_Ranvier_this_axon)
        if len(streamline_resampled) < axon_morphology['n_Ranviers']:
            print("streamline ", streamline_index, " is too short")
            excluded_streamlines.append(streamline_index)
        else:
            streamlines_resampled.append(streamline_resampled)
            #total_points = total_points + len(streamline_resampled)

    return streamlines_resampled, excluded_streamlines

def get_local_compartment_coords(axon_morphology):

    """ Get 1-D coordinates of internodal compartments relative to the node at 0.0

    Parameters
    ----------
     axon_morphology: dict, geometric description of a single axon, see get_axon_morphology

    Returns
    -------
    Nx1 numpy.ndarray

    """

    loc_coords = np.zeros(axon_morphology['n_comp'] - 1, float)
    loc_pos = 0.0  # just for clarity

    if axon_morphology['axon_model'] == 'MRG2002' and axon_morphology['axon_diam'] >= 5.7:
        # only internodal compartments. The distances will be computed from the node of Ranvier using loc_pos

        for inx_loc in np.arange(1, axon_morphology['n_comp']):
            inx_loc = int(inx_loc)
            if inx_loc == 1:
                loc_pos = (axon_morphology['ranvier_length'] + axon_morphology['para1_length']) / 2

            if inx_loc == 2 or inx_loc == 10:
                loc_pos = loc_pos + (
                        axon_morphology['para1_length'] + axon_morphology['para2_length']) / 2
            if inx_loc == 3 or inx_loc == 9:
                loc_pos = loc_pos + (
                        axon_morphology['para2_length'] + axon_morphology['inter_length']) / 2
            if inx_loc == 4 or inx_loc == 5 or inx_loc == 6 or inx_loc == 7 or inx_loc == 8:
                loc_pos = loc_pos + (axon_morphology['inter_length']) / 1

        loc_coords[inx_loc-1] = loc_pos

    elif axon_morphology['axon_model'] == 'MRG2002' and axon_morphology['axon_diam'] < 5.7:
        for inx_loc in np.arange(1, axon_morphology['n_comp']):
            if inx_loc == 1:
                loc_pos = (axon_morphology['ranvier_length'] + axon_morphology['para1_length']) / 2
            if inx_loc == 2 or inx_loc == 7:
                loc_pos = loc_pos + (
                        axon_morphology['para1_length'] + axon_morphology['para2_length']) / 2
            if inx_loc == 3 or inx_loc == 6:
                loc_pos = loc_pos + (
                        axon_morphology['para2_length'] + axon_morphology['inter_length']) / 2
            if inx_loc == 4 or inx_loc == 5:
                loc_pos = loc_pos + (axon_morphology['inter_length']) / 1  # switch to mm from µm

        loc_coords[inx_loc-1] = loc_pos

    elif axon_morphology['axon_model'] == 'MRG2002_DS' and axon_morphology['axon_diam'] >= 5.7:
        # node -- -- internodal -- -- -- -- internodal -- -- node
        for inx_loc in np.arange(1,axon_morphology['n_comp']):  # only internodal compartments. The distances will be computed from the node of Ranvier using loc_pos
            if inx_loc == 1:
                loc_pos = (axon_morphology['ranvier_length'] + axon_morphology['inter_length']) / 2 + axon_morphology['para1_length'] + axon_morphology['para2_length']
            elif inx_loc == 2:
                loc_pos = loc_pos + 5 * axon_morphology['inter_length']
            else:
                print('wrong number of compartments')

        loc_coords[inx_loc - 1] = loc_pos

    elif axon_morphology['axon_model'] == 'MRG2002_DS' and axon_morphology['axon_diam'] < 5.7:
        # mode -- -- -- internodal -- -- -- node
        loc_coords[0] = 0.5 * axon_morphology['ranvier_length'] + 1.5 * axon_morphology['inter_length'] + axon_morphology['para1_length'] + axon_morphology['para2_length']

    elif axon_morphology['axon_model'] == 'McNeal1976':
        # mode -- -- -- internodal -- -- -- node
        loc_coords[0] = axon_morphology['node_step'] * 0.5

    return loc_coords

if __name__ == '__main__':

    """ Call to allocate axon model for OSS-DBS simulations

    Parameters
    ----------
    stim_dir: str, path to stimulation folder where oss-dbs_parameter.mat is stored
    hemis_idx: int, hemisphere ID (0 - right, 1 - left)
    NEURON model: str, optional, possible options :'MRG2002', 'MRG2002_DS' (downsampled), 'McNeal1976'
     (classic McNeal's)
        
    """

    # get stim directory
    # we can also retrieve it from the path for oss-dbs_parameters.mat
    stim_dir = sys.argv[1:][0]

    hemis_idx = sys.argv[1:][1]
    oss_dbs_parameters_path = sys.argv[1:][2]

    # process Lead-DBS input
    axons_for_PAM = AxonModels(stim_dir, hemis_idx, oss_dbs_parameters_path)
    axons_for_PAM.convert_fibers_to_axons()
