import h5py
import os
import neuron as n
import numpy as np
from pandas import read_csv
import time as tm
import multiprocessing as mp
import logging
from scipy.io import savemat
import json


class NeuronStimulation:

    """Pathway activation modeling based on extracellular voltage

    Parameters
    ----------
    pathway_dict: dict, parameters of the simulated pathway
    signal_dict: dict, parameters of the simulated signal
    folder_to_save: str, folder for outputs
    connectome_name: str, optional, name of the connectome for Lead-DBS reference
    scaling_vector: Nx1 numpy array, optional, used for N-contacts (fields) superposition
    scaling_index: int, optional, superposition index for post-processing

    """

    def __init__(self, pathway_dict, signal_dict, folder_to_save, scaling_vector=None, scaling_index=None):

        self.folder_to_save = folder_to_save

        self.connectome_name = pathway_dict['connectome_name']
        self.pathway_name = pathway_dict['pathway_name']  # str
        self.axonModel = pathway_dict['Axon_Model_Type']  # str
        self.axonDiam = pathway_dict['axon_diam']    # int, in um
        self.nRanvier = pathway_dict['n_Ranvier']    # int, number of nodes of Ranvier
        self.N_neurons = pathway_dict['N_seeded_neurons']  # int, number of neurons in the pathway
        self.orig_N_neurouns = pathway_dict['N_orig_neurons']  # int, number of neurons per pathway as defined in the connectome (before Kuncel pre-filtering)

        self.timeStep = signal_dict['time_step']   # float, in ms
        self.nSteps = signal_dict['N_time_steps']  # int, total number of steps
        self.stepsPerMs = int(1.0 / self.timeStep)
        # we always simulate only one pulse from DBS. If more needed, we should just copy the array and make sure the time vectors are in order
        self.n_pulse = 1
        self.scaling_vector = scaling_vector   # superposition-based solution scaling (contact-wise vector!)
        self.scaling = signal_dict['scaling']  # float, simple scaling of the whole solution (not contact-wise!)
        self.scaling_index = scaling_index

        from axon_allocation import get_axon_morphology

        if self.axonModel == 'McIntyre2002_ds':
            axonNEURONModel = 'McIntyre2002'
        else:
            axonNEURONModel = self.axonModel

        self.axon_morphology = get_axon_morphology(axonNEURONModel, self.axonDiam, None, self.nRanvier)

        # check actual number of n_segments in case downsampeld
        ax_mh = get_axon_morphology(self.axonModel, self.axonDiam, None, self.nRanvier)
        self.n_segments_actual = ax_mh['n_segments']


    def run_NEURON(self, v_ext):

        """Load a NEURON model and run simulation

        Parameters
        ----------
        v_ext: NxM numpy.ndarray, extracellular el. potential distribution (mV) in space (on the neuron) and time (DBS signal)

        Returns
        -------
        spike: bool, flag for induced action potential
        """

        # neuron simulation
        if self.axonModel == 'McNeal1976':
            v_init = -70.0
            n.h('{load_file("init_B5_extracellular.hoc")}')
        elif self.axonModel == "McIntyre2002" or self.axonModel == "McIntyre2002_ds":
            v_init = -80.0
            n.h('{load_file("axon4pyfull.hoc")}')
            n.h.deletenodes()
            n.h.createnodes()
            n.h.dependent_var()
            n.h.initialize()
        else:
            logging.critical("The NEURON model is not supported")
            raise SystemExit

        n.h.setupAPWatcher_0()  # 'left' end of axon
        n.h.setupAPWatcher_1()  # 'right' end of axon

        n.h.dt = self.timeStep
        n.h.tstop = self.nSteps * self.timeStep
        n.h.n_pulse = self.n_pulse
        n.h.v_init = v_init

        for i in range(v_ext.shape[0]):
            n.h.wf[i] = n.h.Vector(v_ext[i, :])  # feed the potential in time for compartment i to the NEURON model

        n.h.stimul()
        n.h.run()
        spike = bool(n.h.stoprun)  # 1 - activated, otherwise 0

        return spike

    def get_v_ext(self,v_time_sol):

        """ Convert el. potential computed by OSS-DBS to extracellular el. potential taking into account scaling parameters

        Parameters
        ----------
        v_time_sol: NxM numpy.ndarray, (abstract) el. potential distribution (in V) in space (on the neuron) and time (DBS signal)

        Returns
        -------
        v_ext: NxM numpy.ndarray, extracellular el. potential distribution (mV) in space (on the neuron) and time (DBS signal)

        """

        v_ext = np.zeros_like(v_time_sol, float)
        if self.scaling_vector is None:
            v_ext = v_time_sol * 1000.0 * self.scaling  # convert to mV
        else:
            for contact_i in range(len(self.scaling_vector)):
                v_ext = v_ext + v_time_sol * 1000.0 * self.scaling * self.scaling_vector[contact_i]  # convert to mV

        return v_ext


    def modify_hoc_file(self):
        """ Directly insert neuron parameters to .hoc files.
            This function should be substituted by passing parameters to .hoc directly.

        """

        if self.axonModel == "McIntyre2002" or self.axonModel == "McIntyre2002_ds":
        
            if self.axonDiam >= 5.7:
                axoninter = (self.nRanvier - 1) * 6
            else:
                axoninter = (self.nRanvier - 1) * 3

            v_init = -80.0
            #from Axon_files.Parameter_insertion_python3 import paste_to_hoc_python3
            from Axon_files.Parameter_insertion_python3 import paste_to_hoc_python3
            paste_to_hoc_python3(self.nRanvier, self.axon_morphology['n_para1'], self.axon_morphology['n_para2'],
                                 axoninter, self.axon_morphology['n_segments'], v_init, self.axonDiam, self.axon_morphology['para1_length']*1e3,
                                 self.axon_morphology['para2_length']*1e3, self.axon_morphology['ranvier_length']*1e3, self.axon_morphology['node_d'],
                                 self.axon_morphology['axon_d'], self.axon_morphology['para1_d'], self.axon_morphology['para2_d'],
                                 self.axon_morphology['node_step']*1e3, self.axon_morphology["lamellas"], self.stepsPerMs)

        elif self.axonModel == "McNeal1976":
            n_internodal = self.axon_morphology['n_segments'] - self.nRanvier 
            v_init = -70.0
            from Axon_files.McNeal1976.Parameter_insertion_python3_McNeal import paste_to_hoc_python3
            #from Axon_files.Reilly2016.Parameter_insertion_python3_Reilly import paste_to_hoc_python3
            paste_to_hoc_python3(self.nRanvier, n_internodal, self.axon_morphology['n_segments'],
                                 v_init, self.stepsPerMs)


    def create_leaddbs_outputs(self, Axon_Lead_DBS):

        """ Export axons with activation state in Lead-DBS supported format

        Parameters
        ----------
        Axon_Lead_DBS: NxM numpy.ndarray, geometry, index and activation status of neurons (equivalent of connectome.fibers format in Lead-DBS)

        """

        mdic = {"fibers": Axon_Lead_DBS, "ea_fibformat": "1.0",
                "connectome_name": self.connectome_name}  # For Lead-DBS .mat files

        if self.scaling_index is None:
            if self.pathway_name is None:
                savemat(self.folder_to_save + "/Axon_state.mat", mdic)
            else:
                savemat(self.folder_to_save + "/Axon_state_" + self.pathway_name + ".mat", mdic)
        else:
            if self.pathway_name is None:
                savemat(self.folder_to_save + "/Axon_state_" + str(self.scaling_index) + ".mat", mdic)
            else:
                savemat(self.folder_to_save  + "/Axon_state_" + self.pathway_name + "_" + str(
                    self.scaling_index) + ".mat", mdic)


    def create_paraview_outputs(self, Axon_Lead_DBS):

        """ Export axons with activation state in Paraview supported format

        Parameters
        ----------
        Axon_Lead_DBS: NxM numpy.ndarray, geometry, index and activation status of neurons (equivalent of connectome.fibers format in Lead-DBS)

        """

        if self.scaling_index is None:
            if self.pathway_name is None:
                np.savetxt(self.folder_to_save + "/Axon_state.csv", Axon_Lead_DBS, delimiter=",", header='x-pt,y-pt,z-pt,idx,status')
            else:
                np.savetxt(self.folder_to_save + "/Axon_state_" + self.pathway_name + ".csv", Axon_Lead_DBS, delimiter=",", header='x-pt,y-pt,z-pt,idx,status')
        else:
            if self.pathway_name is None:
                np.savetxt(self.folder_to_save + "/Axon_state_" + str(self.scaling_index) + ".csv", Axon_Lead_DBS, delimiter=",", header='x-pt,y-pt,z-pt,idx,status')
            else:
                np.savetxt(self.folder_to_save  + "/Axon_state_" + self.pathway_name + "_" + str(
                    self.scaling_index) + ".csv", Axon_Lead_DBS, delimiter=",", header='x-pt,y-pt,z-pt,idx,status')


    def store_axon_statuses(self, percent_activated, percent_damaged, percent_csf):

        """ store PAM results

        Parameters
        ----------
        percent_activated: float, percent of the original(!) number of neurons that are activated for the particular stimulation
        percent_damaged: float, percent of the original(!) number of neurons that are 'damaged' for the particular electrode placement
        percent_csf: float, percent of the original(!) number of neurons that intersect with CSF for the particular brain segmentation

        Note
        ----------
        For activation state of particular neuron see 'fiberActivation*' files as those restore original(!) indices as in the connectome.
        """

        summary_dict = {
            'percent_activated': percent_activated,
            'percent_damaged': percent_damaged,
            'percent_csf': percent_csf,
        }

        if self.scaling_index is None:
            if self.pathway_name is None:
                with open(self.folder_to_save + '/Pathway_status.json', 'w') as save_as_dict:
                    json.dump(summary_dict, save_as_dict)
            else:
                summary_dict['pathway_name'] = self.pathway_name
                with open(self.folder_to_save + '/Pathway_status_' + self.pathway_name + '.json', 'w') as save_as_dict:
                    json.dump(summary_dict, save_as_dict)
        else:
            summary_dict['scaling_index'] = str(self.scaling_index)
            if self.pathway_name is None:
                with open(self.folder_to_save + '/Pathway_status_' + str(self.scaling_index) + '.json', 'w') as save_as_dict:
                    json.dump(summary_dict, save_as_dict)
            else:
                summary_dict['pathway_name'] = self.pathway_name
                with open(self.folder_to_save + '/Pathway_status_' + self.pathway_name + "_" + str(self.scaling_index) + '.json', 'w') as save_as_dict:
                    json.dump(summary_dict, save_as_dict)


    def check_pathway_activation(self, pathway_dataset):

        """ Parallelized probing of action potentials at all stimulated neurons (axons) described by supplied pathway datagroups

        Parameters
        ----------
        pathway_dataset: h5 group, contains datasets that describe geometries for all neuron (axons)
        pre_status: Nx1 numpy.ndarray, vector of initial neuron (axon) statuses (0 - available for probing, -1 - intersected with electrode, -2 - inside CSF)

        Note
        ----------
        Creates 'Axon_states_*' .mat and .csv files for visualization in Lead-DBS and Paraview, respectively.
        Also stores summary statistics in 'Pathway_status_*.json'
        """

        # use half of CPUs
        N_proc = mp.cpu_count() / 2
        self.modify_hoc_file()

        # check pre-status
        pre_status = pathway_dataset['Status']

        # initialize outputs
        Axon_Lead_DBS = np.zeros((self.N_neurons * self.n_segments_actual, 5), float)
        List_of_activated = []
        List_of_not_activated = []
        Activated_models = 0

        neuron_index = 0
        while neuron_index < self.N_neurons:

            # run parallel processing
            proc = []
            j_proc = 0  # counter for processes
            output = mp.Queue()
            while j_proc < N_proc and neuron_index < self.N_neurons:

                # get neuron geometry and field solution
                neuron = pathway_dataset['axon' + str(neuron_index)]
                Axon_Lead_DBS[neuron_index * self.n_segments_actual:(neuron_index + 1) * self.n_segments_actual,
                :3] = np.array(neuron['Points[mm]'])

                # add index
                Axon_Lead_DBS[neuron_index * self.n_segments_actual:(neuron_index + 1) * self.n_segments_actual,
                3] = neuron_index + 1  # because Lead-DBS numbering starts from 1

                # temp binary solution. True - potential probed
                # if not(pre_status[neuron_index]):
                #     Axon_Lead_DBS[neuron_index * self.n_segments_actual:(neuron_index + 1) * self.n_segments_actual, 4] = -1
                #     neuron_index += 1
                #     continue

                # check which neurons were flagged with CSF and electrode intersection, skip probing of those
                if pre_status[neuron_index] != 0:
                    Axon_Lead_DBS[neuron_index * self.n_segments_actual:(neuron_index + 1) * self.n_segments_actual, 4] = pre_status[neuron_index]
                    neuron_index += 1
                    continue

                neuron_time_sol = np.array(neuron['Potential[V]'])

                processes = mp.Process(target=self.get_axon_status, args=(neuron_index, neuron_time_sol, output))
                proc.append(processes)

                j_proc += 1
                neuron_index += 1

            for p in proc:
                p.start()
            for p in proc:
                p.join()

            # check the status of batch processed neurons
            neurons_idxs_stat = [output.get() for p in proc]
            for n_idx_stat in neurons_idxs_stat:   #  n_idx_stat is a list[neuron index, status (1 or 0)]
                if n_idx_stat[1] == 1:
                    Activated_models += 1
                    List_of_activated.append(n_idx_stat[0])
                else:
                    List_of_not_activated.append(n_idx_stat[0])

        # iterate over all neurons initially placed by OSS-DBS and assign their statuses
        for neuron_index in range(self.N_neurons):
            if neuron_index in List_of_activated:
                Axon_Lead_DBS[neuron_index * self.n_segments_actual:(neuron_index + 1) * self.n_segments_actual, 4] = 1
            elif neuron_index in List_of_not_activated:
                Axon_Lead_DBS[neuron_index * self.n_segments_actual:(neuron_index + 1) * self.n_segments_actual, 4] = 0
            else:
                # the status was already assigned
                continue

        self.create_leaddbs_outputs(Axon_Lead_DBS)
        self.create_paraview_outputs(Axon_Lead_DBS)

        percent_activated = np.round(Activated_models/float(self.orig_N_neurouns)*100,2)
        percent_damaged = np.round(np.sum(pre_status == -1)/float(self.orig_N_neurouns)*100,2)
        percent_csf = np.round(np.sum(pre_status == -2)/float(self.orig_N_neurouns)*100,2)

        #percent_activated = np.round(Activated_models/float(self.orig_N_neurouns)*100,2)
        #percent_damaged = np.round((len(pre_status)-np.count_nonzero(pre_status))/float(self.orig_N_neurouns)*100,2)
        #percent_csf = np.round(np.sum(pre_status == -2)/float(self.orig_N_neurouns)*100,2)

        print("\n\nPathway ",self.pathway_name, ": ")
        print("Activated neurons: ", percent_activated, " %")
        print("Neurons with status -1: ", percent_damaged, " %")
        print("Neurons with status -2: ", percent_csf, " %")

        self.store_axon_statuses(percent_activated, percent_damaged, percent_csf)

    def get_axon_status(self, neuron_index, v_time_sol, output):

        """ Probe action potential at a neuron (axon)

        Parameters
        ----------
        neuron_index: int, index of neuron in the pathway starting from 0
        v_time_sol: NxM numpy.ndarray, (abstract) el. potential distribution (in V) in space (on the neuron) and time (DBS signal)

        Returns
        -------
        list, neuron index in the pathway and its activation status (1 or 0)

        """

        if self.axonModel == 'McIntyre2002_ds':
            v_time_sol = self.upsample_voltage(v_time_sol)

        v_ext = self.get_v_ext(v_time_sol)
        spike = self.run_NEURON(v_ext)

        return output.put([neuron_index, spike])


    def upsample_voltage(self, v_time_sol):

        # ToDo: estimate ratios directly from the morphology

        """ Upsample el. potential distribution for downsampled neurons by interpolating

        Parameters
        ----------
        v_time_sol: LxM numpy.ndarray, (abstract) el. potential distribution (in V) in space (on the neuron) and time (DBS signal)

        Returns
        -------
        v_time_sol_full: NxM numpy.ndarray

        """

        # let's interpolate voltage between node - center_l - center_r - node
        # assume 11 segments
        # n_segments_ds = ((n_segments_full - 1) / 11) * 3 +1

        v_time_sol_full = np.zeros((self.axon_morphology['n_segments'],v_time_sol.shape[1]), float)

        if self.axonDiam >= 5.7:
            # fill out nodes first
            for k in np.arange(0, self.axon_morphology['n_segments'], 11):
                z = int(k / 11) * 3
                v_time_sol_full[k, :] = v_time_sol[z,:]

            # now two segments in between
            for k in np.arange(3, self.axon_morphology['n_segments'], 11):
                z = int(k / 11) * 3 + 1
                v_time_sol_full[k, :] = v_time_sol[z,:]

            for k in np.arange(8, self.axon_morphology['n_segments'], 11):
                z = int(k / 11) * 3 + 2
                v_time_sol_full[k, :] = v_time_sol[z,:]

            # node -- -- intern -- -- -- -- intern -- -- node  ->  node-para1-para2-intern-intern-intern-intern-intern-intern-para2-para1-node
            internodal_length = (self.axon_morphology['node_step'] - self.axon_morphology['ranvier_length'] - (self.axon_morphology['para2_length'] + self.axon_morphology['para1_length'])*2) / 6.0
            dist_node_internode = (self.axon_morphology['ranvier_length'] + internodal_length) * 0.5 + self.axon_morphology['para1_length'] + self.axon_morphology['para2_length']
            ratio_1 = (self.axon_morphology['ranvier_length'] + self.axon_morphology['para1_length']) * 0.5 / dist_node_internode
            ratio_2 = ratio_1 + (self.axon_morphology['para1_length'] + self.axon_morphology['para2_length']) * 0.5 / dist_node_internode
            dist_intern_intern = internodal_length * 5
            ratio_3 = internodal_length / dist_node_internode
            ratio_4 = ratio_3 + internodal_length / dist_node_internode

            # now interpolate to the rest
            list_interp = [[1, 2], [4, 5, 6, 7], [9, 10]]  # local indices of interpolated segments
            for interv in range(len(list_interp)):
                for j in np.arange(0, self.axon_morphology['n_segments'] - 1, 11):
                    if interv == 0:
                        v_time_sol_full[j + 1, :] = (1-ratio_1)  * v_time_sol_full[j, :] + (ratio_1) * v_time_sol_full[j + 3, :]
                        v_time_sol_full[j + 2, :] = (1-ratio_2)  * v_time_sol_full[j, :] + (ratio_2) * v_time_sol_full[j + 3, :]
                    elif interv == 1:
                        v_time_sol_full[j + 4, :] = (1-ratio_3) * v_time_sol_full[j + 3, :] + ratio_3 * v_time_sol_full[j + 8, :]
                        v_time_sol_full[j + 5, :] = (1-ratio_4) * v_time_sol_full[j + 3, :] + ratio_4 * v_time_sol_full[j + 8, :]
                        v_time_sol_full[j + 6, :] = ratio_4 * v_time_sol_full[j + 3, :] + (1-ratio_4) * v_time_sol_full[j + 8, :]
                        v_time_sol_full[j + 7, :] = ratio_3 * v_time_sol_full[j + 3, :] + (1-ratio_3) * v_time_sol_full[j + 8, :]
                    else:
                        v_time_sol_full[j + 9, :] = ratio_2 * v_time_sol_full[j + 8, :] + (1-ratio_2) * v_time_sol_full[j + 11, :]
                        v_time_sol_full[j + 10, :] = ratio_1 * v_time_sol_full[j + 8, :] + (1-ratio_1) * v_time_sol_full[j + 11, :]
        else:
            # let's interpolate voltage between node - center - node
            # assume 8 segments
            # fill out nodes first
            for k in np.arange(0, self.axon_morphology['n_segments'], 8):
                z = int(k / 8) * 2
                v_time_sol_full[k, :] = v_time_sol[z,:]

            # now the center between nodes
            for k in np.arange(4, self.axon_morphology['n_segments'], 8):
                z = int(k / 8) * 2 + 1
                v_time_sol_full[k, :] = v_time_sol[z,:]

            # node -- -- -- internodal -- -- -- node  ->  node-para1-para2-intern-intern-intern-para2-para1-node
            dist_node_internode = self.axon_morphology['node_step'] / 2.0
            internodal_length = (self.axon_morphology['node_step'] - self.axon_morphology['ranvier_length'] - (self.axon_morphology['para2_length'] + self.axon_morphology['para1_length'])*2) / 3.0
            ratio_1 = (self.axon_morphology['ranvier_length'] + self.axon_morphology['para1_length']) * 0.5/ dist_node_internode
            ratio_2 = ratio_1 + (self.axon_morphology['para1_length'] + self.axon_morphology['para2_length']) * 0.5/ dist_node_internode
            ratio_3 = ratio_2 + (
                        self.axon_morphology['para2_length'] + internodal_length) * 0.5 / dist_node_internode

            # now interpolate to the rest
            list_interp = [[1, 2, 3], [5, 6, 7]]  # local indices of interpolated segments
            for interv in range(len(list_interp)):
                for j in np.arange(0, self.axon_morphology['n_segments'] - 1, 8):
                    if interv == 0:  # ratios based on intercompartment distances
                        v_time_sol_full[j + 1, :] = (1-ratio_1) * v_time_sol_full[j, :] + (ratio_1) * v_time_sol_full[j + 4, :]
                        v_time_sol_full[j + 2, :] = (1-ratio_2) * v_time_sol_full[j, :] + (ratio_2) * v_time_sol_full[j + 4, :]
                        v_time_sol_full[j + 3, :] = (1-ratio_3) * v_time_sol_full[j, :] + (ratio_3) * v_time_sol_full[j + 4, :]
                    else:
                        v_time_sol_full[j + 5, :] = (ratio_3) * v_time_sol_full[j + 4, :] + (1-ratio_3) * v_time_sol_full[j + 8, :]
                        v_time_sol_full[j + 6, :] = (ratio_2) * v_time_sol_full[j + 4, :] + (1-ratio_2) * v_time_sol_full[j + 8, :]
                        v_time_sol_full[j + 7, :] = (ratio_1) * v_time_sol_full[j + 4, :] + (1-ratio_1) * v_time_sol_full[j + 8, :]

        return v_time_sol_full

    # def restore_downsampled(self, n_segments):
    #
    #     if self.axonDiam >= 5.7:
    #         n_segments_full = ((n_segments - 1) / 3) * 11 + 1
    #     else:
    #         n_segments_full = ((n_segments - 1) / 2) * 8 + 1
    #
    #     if n_segments_full.is_integer():
    #         n_segments_full = int(n_segments_full)
    #     else:
    #         logging.critical('Mismatch between the downsampled and the full model')
    #         raise SystemExit
    #
    #     return n_segments_full
