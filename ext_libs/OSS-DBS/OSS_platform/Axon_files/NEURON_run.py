# -*- coding: utf-8 -*-
"""

@author: Konstantin Butenko

run_simulation_with_NEURON() prepares a parallelized solution of cable equations in NEURON
solve_parallel_NEURON() calls external NEURON scripts

At the moment, two cable models are supported
    the commonly used double cable equation by McIntyre, Grill and Richardson
    and Carnevale's (Reilly2016) implementation of classic McNeal's model

Returns the states of neurons in response to the extracellular stimulation
This is the concluding script of the OSS-DBS pipeline, various outputs are generated.
"""

import h5py
import os
import neuron as n
import numpy as np
from pandas import read_csv
import matplotlib.pyplot as plt
import time as tm
import multiprocessing as mp
import logging
from scipy.io import savemat


def solve_parallel_NEURON(d, last_point, N_index_glob, N_index, n_segments, fib_diam, S_vector, output):

    """ For the given neuron (defined by N_index and last_point) loads the distribution of the electric potential
        in space and time and imports it to the NEURON model. If the model responses with an action potential the
        spike variable is set to 0.5 """

    # t_step in seconds here
    tstop = (d['t_step'] * 1000.0) * d["t_steps_trunc"]  # convert to ms
    dt = d['t_step'] * 1000.0
    n_pulse = 1  # we always simulate only one pulse from DBS. If need more, we should just copy the array and make sure the time vectors are in order

    # upsample McIntyre's model
    if d['downsampled_axon'] == True:
        if d["Axon_Model_Type"] == "McIntyre2002" and fib_diam >= 5.7:
            n_segments_true = ((n_segments - 1) / 3) * 11 + 1
        elif d["Axon_Model_Type"] == "McIntyre2002" and fib_diam < 5.7:
            n_segments_true = ((n_segments - 1) / 2) * 8 + 1

        if n_segments_true.is_integer():
            n_segments_true = int(n_segments_true)
            V_art = np.zeros((n_segments_true, d["t_steps_trunc"]), float)
        else:
            logging.critical('Mismatch between the downsampled and the full model')
            raise SystemExit
    else:
        V_art = np.zeros((n_segments, d["t_steps_trunc"]), float)


    if type(S_vector) is not np.ndarray and not isinstance(S_vector, list):
        axon_in_time = np.load(os.environ['LGFDIR'] + '/Axons_in_time/Signal_t_conv' + str(
            n_segments - 1 + N_index * n_segments + last_point) + '.npy')  #to distinguish axons in different populations, we indexed them with the global index of the last compartment
        if d['downsampled_axon'] == False:
            for i in range(n_segments):
                V_art[i, :] = axon_in_time[i, :d["t_steps_trunc"]] * (1000) * d["Ampl_scale"]  # convert to mV
        else:
            # let's interpolate voltage between node - center_l - center_r - node
            # assume 11 segments
            # n_segments_new = ((n_segments - 1) / 11) * 3 +1

            if fib_diam >= 5.7:
                # fill out nodes first
                for k in np.arange(0, n_segments_true, 11):
                    z = int(k / 11) * 3
                    V_art[k, :] = axon_in_time[z, :d["t_steps_trunc"]] * (1000) * d["Ampl_scale"]  # convert to mV

                # now two segments in between
                for k in np.arange(3, n_segments_true, 11):
                    z = int(k / 11) * 3 + 1
                    V_art[k, :] = axon_in_time[z, :d["t_steps_trunc"]] * (1000) * d["Ampl_scale"]  # convert to mV

                for k in np.arange(8, n_segments_true, 11):
                    z = int(k / 11) * 3 + 2
                    V_art[k, :] = axon_in_time[z, :d["t_steps_trunc"]] * (1000) * d["Ampl_scale"]  # convert to mV

                # now interpolate to the rest
                list_interp = [[1, 2], [4, 5, 6, 7], [9, 10]]  # local indices of interpolated segments
                for interv in range(len(list_interp)):
                    for j in np.arange(0, n_segments_true - 1, 11):
                        if interv == 0:
                            V_art[j + 1, :] = 0.962 * V_art[j, :] + 0.038 * V_art[j + 3, :]
                            V_art[j + 2, :] = 0.77 * V_art[j, :] + 0.23 * V_art[j + 3, :]
                        elif interv == 1:
                            V_art[j + 4, :] = 0.80 * V_art[j + 3, :] + 0.20 * V_art[j + 8, :]
                            V_art[j + 5, :] = 0.60 * V_art[j + 3, :] + 0.40 * V_art[j + 8, :]
                            V_art[j + 6, :] = 0.40 * V_art[j + 3, :] + 0.60 * V_art[j + 8, :]
                            V_art[j + 7, :] = 0.20 * V_art[j + 3, :] + 0.80 * V_art[j + 8, :]
                        else:
                            V_art[j + 9, :] = 0.23 * V_art[j + 8, :] + 0.77 * V_art[j + 11, :]
                            V_art[j + 10, :] = 0.038 * V_art[j + 8, :] + 0.962 * V_art[j + 11, :]

            else:
                # let's interpolate voltage between node - center - node
                # assume 8 segments
                # fill out nodes first
                for k in np.arange(0, n_segments_true, 8):
                    z = int(k / 8) * 2
                    V_art[k, :] = axon_in_time[z, :d["t_steps_trunc"]] * (1000) * d["Ampl_scale"]  # convert to mV

                # now the center between nodes
                for k in np.arange(4, n_segments_true, 8):
                    z = int(k / 8) * 2 + 1
                    V_art[k, :] = axon_in_time[z, :d["t_steps_trunc"]] * (1000) * d["Ampl_scale"]  # convert to mV

                # now interpolate to the rest
                list_interp = [[1, 2, 3], [5, 6, 7]]  # local indices of interpolated segments
                for interv in range(len(list_interp)):
                    for j in np.arange(0, n_segments_true - 1, 8):
                        if interv == 0: # ratios based on intercompartment distances
                            V_art[j + 1, :] = 0.985 * V_art[j, :] + 0.015 * V_art[j + 4, :]
                            V_art[j + 2, :] = 0.90 * V_art[j, :] + 0.10 * V_art[j + 4, :]
                            V_art[j + 3, :] = 0.60 * V_art[j, :] + 0.40 * V_art[j + 4, :]
                        else:
                            V_art[j + 5, :] = 0.40 * V_art[j + 4, :] + 0.60 * V_art[j + 8, :]
                            V_art[j + 6, :] = 0.10 * V_art[j + 4, :] + 0.90 * V_art[j + 8, :]
                            V_art[j + 7, :] = 0.015 * V_art[j + 4, :] + 0.985 * V_art[j + 8, :]

    else:
        for j in range(len(S_vector)):
            if S_vector[j] == None:
                S_vector[j] = 0.0  # Important: Here we put it to 0 to drop it's contribution. It does not affect the field solution, because it was not treated as a ground in FFEM

        for contact_i in range(len(S_vector)):
            if S_vector[contact_i] == 0.0:
                continue
            else:
                axon_in_time = np.load(os.environ['LGFDIR'] + '/Axons_in_time/Signal_t_conv' + str(
                    n_segments - 1 + N_index * n_segments + last_point) + "_" + str(contact_i) + '.npy')

                if d['downsampled_axon'] == False:
                    for i in range(n_segments):
                        # print(S_vector[contact_i])
                        V_art[i, :] = V_art[i, :] + axon_in_time[i, :d["t_steps_trunc"]] * (1000) * d["Ampl_scale"] * S_vector[contact_i]  # convert to mV
                else:
                    # let's interpolate voltage between node - center_l - center_r - node
                    # assume 11 segments
                    # n_segments_new = ((n_segments - 1) / 11) * 3 +1

                    V_art_temp = np.zeros((n_segments_true, d["t_steps_trunc"]), float)

                    if fib_diam >= 5.7:
                        # fill out nodes first
                        for k in np.arange(0, n_segments_true, 11):
                            z = int(k / 11) * 3
                            V_art[k, :] = V_art[k, :] + axon_in_time[z, :d["t_steps_trunc"]] * (1000) * d[
                                "Ampl_scale"] * S_vector[contact_i] # convert to mV
                                
                            V_art_temp[k,:] = axon_in_time[z, :d["t_steps_trunc"]] * (1000) * d["Ampl_scale"] * S_vector[contact_i] # convert to mV

                        # now two segments in between
                        for k in np.arange(3, n_segments_true, 11):
                            z = int(k / 11) * 3 + 1
                            V_art[k, :] = V_art[k, :] + axon_in_time[z, :d["t_steps_trunc"]] * (1000) * d[
                                "Ampl_scale"] * S_vector[contact_i] # convert to mV
                                
                            V_art_temp[k,:] = axon_in_time[z, :d["t_steps_trunc"]] * (1000) * d["Ampl_scale"] * S_vector[contact_i] # convert to mV

                        for k in np.arange(8, n_segments_true, 11):
                            z = int(k / 11) * 3 + 2
                            V_art[k, :] = V_art[k, :] + axon_in_time[z, :d["t_steps_trunc"]] * (1000) * d[
                                "Ampl_scale"] * S_vector[contact_i]  # convert to mV
                                
                            V_art_temp[k,:] = axon_in_time[z, :d["t_steps_trunc"]] * (1000) * d["Ampl_scale"] * S_vector[contact_i] # convert to mV


                        # now interpolate to the rest
                        list_interp = [[1, 2], [4, 5, 6, 7], [9, 10]]  # local indices of interpolated segments
                        for interv in range(len(list_interp)):
                            for j in np.arange(0, n_segments_true - 1, 11):
                                if interv == 0:
                                    V_art[j + 1, :] = V_art[j + 1, :] + 0.962 * V_art_temp[j, :] + 0.038 * V_art_temp[j + 3, :]
                                    V_art[j + 2, :] = V_art[j + 2, :] + 0.77 * V_art_temp[j, :] + 0.23 * V_art_temp[j + 3, :]
                                elif interv == 1:
                                    V_art[j + 4, :] = V_art[j + 4, :] + 0.80 * V_art_temp[j + 3, :] + 0.20 * V_art_temp[j + 8, :]
                                    V_art[j + 5, :] = V_art[j + 5, :] + 0.60 * V_art_temp[j + 3, :] + 0.40 * V_art_temp[j + 8, :]
                                    V_art[j + 6, :] = V_art[j + 6, :] + 0.40 * V_art_temp[j + 3, :] + 0.60 * V_art_temp[j + 8, :]
                                    V_art[j + 7, :] = V_art[j + 7, :] + 0.20 * V_art_temp[j + 3, :] + 0.80 * V_art_temp[j + 8, :]
                                else:
                                    V_art[j + 9, :] = V_art[j + 9, :] + 0.23 * V_art_temp[j + 8, :] + 0.77 * V_art_temp[j + 11, :]
                                    V_art[j + 10, :] = V_art[j + 10, :] + 0.038 * V_art_temp[j + 8, :] + 0.962 * V_art_temp[j + 11, :]

                    else:
                        # let's interpolate voltage between node - center - node
                        # assume 8 segments
                        # fill out nodes first
                        for k in np.arange(0, n_segments_true, 8):
                            z = int(k / 8) * 2
                            V_art[k, :] = V_art[k, :] + axon_in_time[z, :d["t_steps_trunc"]] * (1000) * d[
                                "Ampl_scale"] * S_vector[contact_i]  # convert to mV
                                
                            V_art_temp[k,:] = axon_in_time[z, :d["t_steps_trunc"]] * (1000) * d["Ampl_scale"] * S_vector[contact_i] # convert to mV

                        # now the center between nodes
                        for k in np.arange(4, n_segments_true, 8):
                            z = int(k / 8) * 2 + 1
                            V_art[k, :] = V_art[k, :] + axon_in_time[z, :d["t_steps_trunc"]] * (1000) * d[
                                "Ampl_scale"] * S_vector[contact_i]  # convert to mV
                                
                            V_art_temp[k,:] = axon_in_time[z, :d["t_steps_trunc"]] * (1000) * d["Ampl_scale"] * S_vector[contact_i] # convert to mV

                        # now interpolate to the rest
                        list_interp = [[1, 2, 3], [5, 6, 7]]  # local indices of interpolated segments
                        for interv in range(len(list_interp)):
                            for j in np.arange(0, n_segments_true - 1, 8):
                                if interv == 0:  # ratios based on intercompartment distances
                                    V_art[j + 1, :] = V_art[j + 1, :] + 0.985 * V_art_temp[j, :] + 0.015 * V_art_temp[j + 4, :]
                                    V_art[j + 2, :] = V_art[j + 2, :] + 0.90 * V_art_temp[j, :] + 0.10 * V_art_temp[j + 4, :]
                                    V_art[j + 3, :] = V_art[j + 3, :] + 0.60 * V_art_temp[j, :] + 0.40 * V_art_temp[j + 4, :]
                                else:
                                    V_art[j + 5, :] = V_art[j + 5, :] + 0.40 * V_art_temp[j + 4, :] + 0.60 * V_art_temp[j + 8, :]
                                    V_art[j + 6, :] = V_art[j + 6, :] + 0.10 * V_art_temp[j + 4, :] + 0.90 * V_art_temp[j + 8, :]
                                    V_art[j + 7, :] = V_art[j + 7, :] + 0.015 * V_art_temp[j + 4, :] + 0.985 * V_art_temp[j + 8, :]


    ##  only if we want to save potential in time on axons
    #np.save(os.environ['PATIENTDIR'] + '/Field_on_axons_in_time/'+str(population_name)+'axon_'+str(N_index_glob), V_art)

    if d['Axon_Model_Type'] == "Reilly2016":
        n.h('{load_file("init_B5_extracellular.hoc")}')
    elif d['Axon_Model_Type'] == "McIntyre2002":
        n.h('{load_file("axon4pyfull.hoc")}')
        n.h.deletenodes()
        n.h.createnodes()
        n.h.dependent_var()
        n.h.initialize()
    else:
        logging.critical("The NEURON model is not supported")
        raise SystemExit

    n.h.setupAPWatcher_0() # 'left' end of axon
    n.h.setupAPWatcher_1() # 'right' end of axon

    n.h.dt = dt
    n.h.tstop = tstop
    n.h.n_pulse = n_pulse
    n.h.v_init = d['v_init']

    for i in range(0,V_art.shape[0]):
        n.h.wf[i] = n.h.Vector(V_art[i,:])        # feed the potential in time for compartment i to the NEURON model

    n.h.stimul()
    n.h.run()
    spike = n.h.stoprun-0.5

    if spike == 0.5:
        return output.put([N_index_glob, N_index])
    else:
        return output.put([N_index_glob, -1])

def run_simulation_with_NEURON(d, Neuron_models, shift_to_MRI_space, population_index=-1, last_point=0, S_vector=None, scaling_index=None):

    if population_index == -1:      # in this case -1 will refer to the only entry in the lists
        Neuron_models.neurons_idx_encap = [Neuron_models.neurons_idx_encap]
        Neuron_models.neurons_idx_csf = [Neuron_models.neurons_idx_csf]

    # note: these should have been lists before
    if not isinstance(Neuron_models.pattern['num_Ranvier'], list) and not isinstance(Neuron_models.pattern['num_Ranvier'], (np.ndarray, np.generic)): # in this case -1 will refer to the only entry in the lists
        Neuron_models.pattern['num_Ranvier'] = [Neuron_models.pattern['num_Ranvier']]

    if not isinstance(Neuron_models.N_models, list) and not isinstance(Neuron_models.N_models, (np.ndarray, np.generic)):               # in this case -1 will refer to the only entry in the lists
        Neuron_models.N_models = [Neuron_models.N_models]

    # choose the output folder depending on the simulated hemisphere
    if d['Stim_side'] == 0:
        stim_folder = 'Results_rh/'
    else:
        stim_folder = 'Results_lh/'

    #  load compartments of the relevant neuron models (for multiple protocols, might make sense to cache it)
    if population_index == -1:            # only one population is simulated
        population_name = ''
        Vert_full_get = read_csv(os.environ['PATIENTDIR']+'/Neuron_model_arrays/All_neuron_models.csv', delimiter=' ', header=None)     # get all neuron models
        Vert_full = Vert_full_get.values
        Vert_full = np.round(Vert_full,8)

        Vert_get = read_csv(os.environ['PATIENTDIR']+'/Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', delimiter=' ', header=None)    # get only physiologically correct neuron models
        Vert = Vert_get.values
        Vert = np.round(Vert,8)
    else:
        hf = h5py.File(os.environ['PATIENTDIR']+'/Neuron_model_arrays/All_neuron_models_by_populations.h5', 'r')
        lst = list(hf.keys())
        population_name = str(lst[population_index])+'/'
        Vert_full = hf.get(lst[population_index])
        Vert_full = np.array(Vert_full)
        hf.close()
        Vert_full = np.round(Vert_full,8)

        hf2 = h5py.File(os.environ['PATIENTDIR']+'/Neuron_model_arrays/Vert_of_Neural_model_NEURON_by_populations.h5', 'r')
        lst = list(hf2.keys())
        Vert = hf2.get(lst[population_index])
        Vert = np.array(Vert)
        hf2.close()
        Vert = np.round(Vert,8)


    #  Here we assume that all axons have the same number of nodes of Ranvier (and hence the length) and the same morphology
    if d["Axon_Model_Type"] == "McIntyre2002":
        from Axon_files.axon import Axon
        param_ax={
        'centered': True,
        'diameter': Neuron_models.pattern['fiber_diameter'][population_index]
        }
        ax=Axon(param_ax)
        axon_dict=Axon.get_axonparams(ax)
        paranodes1 = axon_dict["para1_nodes"]*(Neuron_models.pattern['num_Ranvier'][population_index] - 1) / (21 - 1)
        paranodes2 = axon_dict["para2_nodes"]*(Neuron_models.pattern['num_Ranvier'][population_index] - 1) / (21 - 1)
        if axon_dict["fiberD"] > 3.0:
            axoninter = (Neuron_models.pattern['num_Ranvier'][population_index] - 1) * 6
        else:
            axoninter = (Neuron_models.pattern['num_Ranvier'][population_index] - 1) * 3
        n_segments = int(Neuron_models.pattern['num_Ranvier'][population_index] + paranodes1 + paranodes2 + axoninter)

        #passing through n.h. does not work sometimes, so we do insert the parameters straight to the file
        from Axon_files.Parameter_insertion_python3 import paste_to_hoc_python3, paste_paraview_vis_python3
        paste_to_hoc_python3(Neuron_models.pattern['num_Ranvier'][population_index],paranodes1,paranodes2,axoninter,n_segments,d['v_init'],axon_dict["fiberD"],axon_dict["para1_length"],axon_dict["para2_length"],axon_dict["ranvier_length"],axon_dict["node_diameter"],axon_dict["axon_diameter"],axon_dict["para1_diameter"],axon_dict["para2_diameter"],axon_dict["deltax"],axon_dict["lamellas"],int(1.0/(d['t_step'] * 1000.0)))
    else:
        axoninter = Neuron_models.pattern['num_Ranvier'][population_index] - 1
        n_segments = int((Neuron_models.pattern['num_Ranvier'][population_index] - 1) * 2) + 1       # lumped internodal segment

        #passing through n.h. does not work sometimes, so we do insert the parameters straight to the file
        from Axon_files.Reilly2016.Parameter_insertion_python3_Reilly import paste_to_hoc_python3
        paste_to_hoc_python3(Neuron_models.pattern['num_Ranvier'][population_index],axoninter,n_segments,d['v_init'],int(1.0/(d['t_step'] * 1000.0)))

    if d['downsampled_axon'] == True:
        # reassign here to downsampled
        n_segments = Neuron_models.pattern['num_segments'][population_index]

    number_neurons_initially = int(Vert_full.shape[0] / n_segments)

    #  stores indices of neurons accounting for those removed by OSS-DBS (but not by Kuncel-VTA!)
    neuron_global_index_array = np.zeros(Neuron_models.N_models[population_index],int)

    # Various output structures

    # Status of neurons NOT removed during Kuncel-VTA (-1 - was not placed, 0 - was not activated, 1 - was activated)
    Vert_full_status = np.zeros(number_neurons_initially, int)

    # Contains info whether the placed(!) neuron was activated (highlights compartments, mostly for visualization purposes)
       
    Nodes_status = np.zeros((Neuron_models.N_models[population_index] * n_segments, 4), float)
    Nodes_status[:,:3] = Vert[:,:]
    Axon_Lead_DBS = np.zeros((number_neurons_initially*n_segments,5), float)  # same but for Lead-DBS (also needs the axon index)

    # x1,y1,z1,x2,y2,z2,status. Holds info only about placed(!) axons (start and end). Important: coordinates will be in the initial MRI space!
    connection_status_MRI = np.zeros((Neuron_models.N_models[population_index], 7), float)

    List_of_activated = []
    List_of_not_activated = []   # not activated but tested in NEURON (i.e. placed)
    Activated_models, int_counter, neuron_index = (0, 0, 0)  # various counters

    # to check progress
    neurons_quart = [int(Neuron_models.N_models[population_index]/4.0), int(2*Neuron_models.N_models[population_index]/4.0), int(3*Neuron_models.N_models[population_index]/4.0)]

    # run NEURON simulations in parallel
    while neuron_index < Neuron_models.N_models[population_index]:
        proc = []
        j_proc = 0 #counter for processes
        output = mp.Queue()
        while j_proc < d["number_of_processors"] and neuron_index < Neuron_models.N_models[population_index]:

            # recover the original index of the neuron (alternatively, we could use a look-up list pointing to original indices)
            first_axon_point=np.array([Vert[neuron_index*n_segments,0],Vert[neuron_index*n_segments,1],Vert[neuron_index*n_segments,2]])
            second_axon_point=np.array([Vert[neuron_index*n_segments+1,0],Vert[neuron_index*n_segments+1,1],Vert[neuron_index*n_segments+1,2]])
            last_axon_point=np.array([Vert[neuron_index*n_segments+n_segments-1,0],Vert[neuron_index*n_segments+n_segments-1,1],Vert[neuron_index*n_segments+n_segments-1,2]])

            inx_first=np.flatnonzero((Vert_full == first_axon_point).all(1))
            inx_second=np.flatnonzero((Vert_full == second_axon_point).all(1))
            inx_last=np.flatnonzero((Vert_full == last_axon_point).all(1))

            # assuming we do not have axons that start (first two points) and end in the same points
            for j in inx_first:
                for j_second in inx_second:
                    if j_second-j == 1:
                        for j_last in inx_last:
                            if j_last-j == n_segments-1:
                                inx_first_true = j
                                break

            neuron_global_index_array[neuron_index] = int(inx_first_true/n_segments)

            processes = mp.Process(target = solve_parallel_NEURON, args=(d,last_point,neuron_global_index_array[neuron_index],neuron_index,n_segments,Neuron_models.pattern['fiber_diameter'][population_index],S_vector,output))
            proc.append(processes)

            j_proc += 1
            neuron_index += 1
            
            if Neuron_models.N_models[population_index] > 500 and neuron_index in neurons_quart: # progress only for large pathways
                logging.critical("{}% of neuron models were processed".format(int(neuron_index * 100 / Neuron_models.N_models[population_index]) + 1))
            
        for p in proc:
            p.start()
        for p in proc:
            p.join()

        # returns list, where activated models have their corresponding numbers, others just -1
        Activated_numbers = [output.get() for p in proc]

        # n_mdls[0] - global index of the neuron model (before adjusting for encap, CSF, etc). n_mdls[1] - local
        for n_mdls in Activated_numbers:            #  n_mdls is list[N_glob,N_loc]!
            if n_mdls[1] != -1:                     #  -1 points to the NOT activated models
                Nodes_status[n_segments*n_mdls[1]:(n_segments*n_mdls[1]+n_segments),3] = 1.0
                Activated_models += 1
                List_of_activated.append(n_mdls[0])
            else:
                List_of_not_activated.append(int(n_mdls[0]))  #

            int_counter += 1
    logging.critical("{} models were activated".format(Activated_models))

    num_removed = 0
    # iterate over all neurons initially placed by OSS-DBS and assign their states
    for axon_i in range(number_neurons_initially):
        Axon_Lead_DBS[axon_i*n_segments:n_segments*(axon_i+1),:3] = Vert_full[axon_i*n_segments:n_segments*(axon_i+1), :3] - shift_to_MRI_space
        Axon_Lead_DBS[axon_i*n_segments:n_segments*(axon_i+1),3] = axon_i+1   # because Lead-DBS number them from 1
        if axon_i in List_of_activated:
            Vert_full_status[axon_i] = 1
            Axon_Lead_DBS[axon_i*n_segments:n_segments*(axon_i+1),4] = 1
        elif axon_i in List_of_not_activated:
            Vert_full_status[axon_i] = 0
            Axon_Lead_DBS[axon_i*n_segments:n_segments*(axon_i+1),4] = 0
        else:
            # check if axon_i is in the lists of removed
            if axon_i in Neuron_models.neurons_idx_encap[population_index]:
                Vert_full_status[axon_i] = -1     # intersected with the electrode / encap
                Axon_Lead_DBS[axon_i*n_segments:n_segments*(axon_i+1), 4] = -1
            elif axon_i in Neuron_models.neurons_idx_csf[population_index]:
                Vert_full_status[axon_i] = -2     # traversed CSF
                Axon_Lead_DBS[axon_i*n_segments:n_segments*(axon_i+1), 4] = -2
            else:
                Vert_full_status[axon_i] = -3     # outside of the domain/ lost?
                Axon_Lead_DBS[axon_i*n_segments:n_segments*(axon_i+1), 4] = -3


    loc_ind_start = 0
    for i in range(Neuron_models.N_models[population_index]):
        connection_status_MRI[i,:3] = Nodes_status[loc_ind_start,:3] - shift_to_MRI_space
        connection_status_MRI[i,3:6] = Nodes_status[loc_ind_start+n_segments-1,:3] - shift_to_MRI_space
        connection_status_MRI[i,6] = Nodes_status[loc_ind_start,3]
        loc_ind_start = loc_ind_start+n_segments

    # #for Paraview visualization purposes in MRI space (only activated models)
    # Nodes_status_MRI_space=np.zeros((Neuron_models.N_models*n_segments,4),float)
    # Nodes_status_MRI_space[:,:3]=Nodes_status[:,:3]-shift_to_MRI_space
    # Nodes_status_MRI_space[:,3]=Nodes_status[:,3]
    # Nodes_status_MRI_space_only_activated = np.delete(Nodes_status_MRI_space, np.where(Nodes_status_MRI_space[:,3] == 0.0)[0], axis=0)


    # Simple way to check activations (only those that were not exluded by Kuncel-VTA)
    Summary_status = np.zeros(5,float)    # Activated, non-activated , 'damaged'
    Summary_status[0] = len(List_of_activated)
    Summary_status[1] = len(List_of_not_activated)
    Summary_status[2] = len(Neuron_models.neurons_idx_encap[population_index])
    Summary_status[3] = len(Neuron_models.neurons_idx_csf[population_index])
    Summary_status[4] = number_neurons_initially-np.sum(Summary_status[:4])

    List_of_activated = np.asarray(List_of_activated)

#    # Outdated: we need to add a specific output for VATs (something with contour extraction)
#    if d["Name_prepared_neuron_array"]==None:
#        hf = h5py.File('Field_solutions/Activation/VAT_Neuron_array.h5', 'a')
#        hf.create_dataset('VAT_Neuron_array_'+str(Activated_models), data=Nodes_status_MRI_space_only_activated)
#        hf.close()
#    else:
#        hf = h5py.File('Field_solutions/Activation/Activation_in_'+d["Name_prepared_neuron_array"][:-3]+'.h5', 'a')
#        hf.create_dataset(str(lst[population_index])+'_'+str(Activated_models), data=Nodes_status_MRI_space_only_activated)
#        hf.close()

    mdic = {"fibers": Axon_Lead_DBS, "ea_fibformat": "1.0", "connectome_name": d["Name_prepared_neuron_array"][:-3]}  # For Lead-DBS .mat files
    if population_index == -1:
        logging.critical("{}% activation (including damaged neurons)\n".format(np.round(Activated_models/float(number_neurons_initially)*100,2)))

        if scaling_index == None:
            np.savetxt(os.environ['PATIENTDIR']+'/'+stim_folder+'Last_run.csv', List_of_activated, delimiter=" ")
            np.save(os.environ['PATIENTDIR']+'/'+stim_folder+'Connection_status_MRI', connection_status_MRI)
            savemat(os.environ['PATIENTDIR'] + "/" + stim_folder + "Axon_state.mat", mdic)
            np.savetxt(os.environ['PATIENTDIR']+'/'+stim_folder+'Neurons_status.csv', Vert_full_status, delimiter=" ")  #Ningfei prefers .csv
            np.savetxt(os.environ['PATIENTDIR'] + '/' + stim_folder + 'Summary_status.csv', Summary_status, delimiter=" ")
            #np.savetxt(os.environ['PATIENTDIR']+'/'+stim_folder+'Activation_Neuron_Array_'+str(Activated_models)+'_MRI.csv', Nodes_status_MRI_space_only_activated, delimiter=" ")
        else:
            np.savetxt(os.environ['PATIENTDIR'] + '/' + stim_folder + 'Last_run_' + str(scaling_index) + '.csv', List_of_activated, delimiter=" ")
            np.save(os.environ['PATIENTDIR'] + '/' + stim_folder + 'Connection_status_MRI_' + str(scaling_index), connection_status_MRI)
            savemat(os.environ['PATIENTDIR'] + "/" + stim_folder + "Axon_state_" + str(scaling_index) + ".mat", mdic)
            np.savetxt(os.environ['PATIENTDIR'] + '/' + stim_folder + 'Neurons_status_' + str(scaling_index) + '.csv', Vert_full_status, delimiter=" ")  # Ningfei prefers .csv
            np.savetxt(os.environ['PATIENTDIR'] + '/' + stim_folder + 'Summary_status_' + str(scaling_index) + '.csv',
                       Summary_status, delimiter=" ")
            #np.savetxt(os.environ['PATIENTDIR'] + '/' + stim_folder + 'Activation_Neuron_Array_' + str(Activated_models) + '_' + str(scaling_index) + '.csv', Nodes_status_MRI_space_only_activated, delimiter=" ")            np.savetxt(os.environ['PATIENTDIR'] + '/' + stim_folder + 'Neurons_status.csv', Vert_full_status, delimiter=" ")
    else:
        logging.critical("{}% activation in {} (including damaged neurons)\n".format(np.round(Activated_models/float(number_neurons_initially)*100,2),lst[population_index]))

        if scaling_index == None:
            savemat(os.environ['PATIENTDIR'] + "/" + stim_folder + "Axon_state_" + str(lst[population_index]) + ".mat", mdic)
            np.savetxt(os.environ['PATIENTDIR']+'/'+stim_folder+'Last_run_in_'+str(lst[population_index])+'.csv', List_of_activated, delimiter=" ")
            np.save(os.environ['PATIENTDIR']+'/'+stim_folder+'Connection_status_MRI_'+str(lst[population_index]), connection_status_MRI)
            #np.savetxt(os.environ['PATIENTDIR']+'/'+stim_folder+'Activation_'+d["Name_prepared_neuron_array"][:-3]+'___'+str(lst[population_index])+'_'+str(Activated_models)+'_MRI.csv', Nodes_status_MRI_space_only_activated, delimiter=" ")

            hf = h5py.File(os.environ['PATIENTDIR'] + '/' + stim_folder + 'Summary_status.h5', 'a')
            hf2 = h5py.File(os.environ['PATIENTDIR'] + '/' + stim_folder + 'Neurons_status.h5', 'a')
        else:
            savemat(os.environ['PATIENTDIR'] + "/" + stim_folder + "Axon_state_" + str(lst[population_index]) + '_' + str(scaling_index) + ".mat", mdic)
            np.savetxt(os.environ['PATIENTDIR']+'/'+stim_folder+'Last_run_in_'+str(lst[population_index])+ '_' + str(scaling_index) + '.csv', List_of_activated, delimiter=" ")
            np.save(os.environ['PATIENTDIR']+'/'+stim_folder+'Connection_status_MRI_'+str(lst[population_index])+ '_' + str(scaling_index) , connection_status_MRI)
            #np.savetxt(os.environ['PATIENTDIR'] + '/' + stim_folder + 'Activation_' + d["Name_prepared_neuron_array"][:-3] + '___' + str(lst[population_index]) + '_' + str(Activated_models) + '_' + str(scaling_index) + '.csv', Nodes_status_MRI_space_only_activated, delimiter=" ")

            hf = h5py.File(os.environ['PATIENTDIR'] + '/' + stim_folder + 'Summary_status_' + str(scaling_index) + '.h5', 'a')
            hf2 = h5py.File(os.environ['PATIENTDIR'] + '/' + stim_folder + 'Neurons_status_' + str(scaling_index) + '.h5', 'a')

        hf.create_dataset(str(lst[population_index]), data=Summary_status)
        hf.close()

        hf2.create_dataset(str(lst[population_index]), data=Vert_full_status)
        hf2.close()

    return Activated_models

