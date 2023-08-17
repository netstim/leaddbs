'''
    By K. Butenko
    This script solves an optimization problem by adjusting the scaling (S_vector) for the unit current solutions
        1) V(P_neuron)=S1*V1(P_neuron)+S2*V2(P_neuron)+...   , where V1 is the solution for Contact1-Ground, the rest are floating
        2) Computes activation profiles for the scaling in NEURON (see test_scaling(...))
        3) Computes distances from these profiles to the optimal profiles for each symptom (see compute_global_score(...))
        4) Computes the global scores as the sum of weighted distances across all symptoms (see compute_global_score(...))
    The main output are currents, activation profiles and distances over opt. iterations stored in .csv files
'''

import numpy as np
import time
import pickle
import subprocess
import importlib
import os
import json
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=FutureWarning)
    import h5py

import logging
import time
import shutil

def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


def test_scaling(S_vector,d,MRI_param,Xs_signal_norm,Neuron_models,FR_vector_signal,t_vector,name_sorted_solution,inx_start_octv,scaling_index,VTA_edge):

    ''' Compute pathay activation for current protocol defined in S_vector'''

    logging.critical("----- Estimating neuron activity -----")
    start_neuron=time.time()

    oss_plat_cont = os.getcwd()

    #==================== Prepare NEURON model =====================#
    # we need to copy Axon_files to the patient folder to ensure stability
    src = oss_plat_cont + "/Axon_files"
    dst = os.environ['PATIENTDIR'] + "/Axon_files"

    if os.path.isdir(dst):
        os.chdir(os.environ['PATIENTDIR'])
        os.system('rm -fr Axon_files')
    os.makedirs(dst)

    copytree(src, dst, symlinks=False, ignore=None)
    os.chdir(dst)  # we now operate in Axon_files/ in the stim folder

    if d["Axon_Model_Type"] == 'McIntyre2002':
        with open(os.devnull, 'w') as FNULL: subprocess.call('nocmodl axnode.mod', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        with open(os.devnull, 'w') as FNULL: subprocess.call('nrnivmodl axnode', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    elif d["Axon_Model_Type"] == 'Reilly2016':
        os.chdir("Reilly2016/")
        if sys.platform == 'win32':
            logging.critical(
                "If using Windows, please make sure your precompile Reilly2016 and comment out the next line")
            raise SystemExit
        with open(os.devnull, 'w') as FNULL: subprocess.call('nrnivmodl', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    # ==================== Run NEURON model =====================#
    from Axon_files.NEURON_run import run_simulation_with_NEURON
    if isinstance(d["n_Ranvier"],list) and len(d["n_Ranvier"])>1:
        Number_of_activated = 0
        last_point=0
        for i in range(len(d["n_Ranvier"])):
            Number_of_activated_population = run_simulation_with_NEURON(d, Neuron_models, np.array(MRI_param.MRI_shift), population_index=i, last_point=last_point, S_vector=S_vector, scaling_index=scaling_index)
            Number_of_activated = Number_of_activated+Number_of_activated_population
            
            last_point=Neuron_models.pattern['num_segments'][i]*Neuron_models.N_models[i]+last_point
    else:
        Number_of_activated = run_simulation_with_NEURON(d, Neuron_models, np.array(MRI_param.MRI_shift), S_vector=S_vector, scaling_index=scaling_index)

    os.chdir(oss_plat_cont)

    minutes=int((time.time() - start_neuron)/60)
    secnds=int(time.time() - start_neuron)-minutes*60
    logging.critical("----- NEURON calculations took {} min {} sec -----\n".format(minutes, secnds))

    # ==================== Post-process the output =====================#
    if d['Stim_side'] == 0:
        stim_folder = 'Results_rh/'
    else:
        stim_folder = 'Results_lh/'

    # retrieve activations across populations
    if isinstance(d["n_Ranvier"],list):
        activation_in_populations=np.zeros(len(d["n_Ranvier"]),int)
        hf = h5py.File(os.environ['PATIENTDIR']+'/'+stim_folder+'Neurons_status_'+str(scaling_index)+'.h5', 'r')
        lst = list(hf.keys())
        for i in range(len(lst)):
            Axon_status = hf.get(lst[i])
            num_activ_in_population = 0
            Axon_status = np.array(Axon_status)
            for j in range(Axon_status.shape[0]):
                if Axon_status[j] == 1:
                    num_activ_in_population += 1
            activation_in_populations[i] = num_activ_in_population
        hf.close()
    else:
        activation_in_populations = Number_of_activated

    if os.path.isdir(os.environ['PATIENTDIR']+'/Field_solutions/Activation'):     # we always re-run NEURON simulation
        os.system('rm -fr '+os.environ['PATIENTDIR']+'/Field_solutions/Activation')
        os.makedirs(os.environ['PATIENTDIR']+'/Field_solutions/Activation')

    # if os.path.isdir(os.environ['PATIENTDIR']+'/Axons_in_time'):     # we always re-run NEURON simulation
    #     os.system('rm -fr '+os.environ['PATIENTDIR']+'/Axons_in_time')
    #     os.makedirs(os.environ['PATIENTDIR']+'/Axons_in_time')

    # store activations over StimSets
    if d['Current_sets'] == True:
        iter_results = []
        iter_results.append(scaling_index)
        if type(activation_in_populations) is np.ndarray:
            for result in activation_in_populations:
                iter_results.append(result)
        else:
            iter_results.append(activation_in_populations)
        iter_results = np.vstack((iter_results)).T

        if scaling_index == 0:
            if d['Stim_side'] == 0:
                with open(os.environ['PATIENTDIR'] + '/Activations_over_StimSets_rh.csv', 'w') as f_handle:
                    np.savetxt(f_handle, iter_results)
            else:
                with open(os.environ['PATIENTDIR'] + '/Activations_over_StimSets_lh.csv', 'w') as f_handle:
                    np.savetxt(f_handle, iter_results)
        else:
            if d['Stim_side'] == 0:
                with open(os.environ['PATIENTDIR'] + '/Activations_over_StimSets_rh.csv', 'a') as f_handle:
                    np.savetxt(f_handle, iter_results)
            else:
                with open(os.environ['PATIENTDIR'] + '/Activations_over_StimSets_lh.csv', 'a') as f_handle:
                    np.savetxt(f_handle, iter_results)

    return activation_in_populations


def compute_global_score(S_vector, *args):
    d, profile_dict, Soft_SE_dict, SE_dict, MRI_param, Xs_signal_norm, N_array, FR_vector_signal, t_vector, name_sorted_solution, inx_start_octv = args

    # read this way for consistency
    # the function will work only for a proper Lead-DBS import (connectome folder, oss-dbs_parameters.mat)
    # get all pathways that survived Kuncel(!) pre-filtering and original(!) number of fibers
    from Pathways_Stats import get_simulated_pathways
    Pathways, axons_in_path = get_simulated_pathways(d['Stim_side'])

    # IMPORTANT: if fibers start and end at the same spot, they will be considered as the same (because of local to global mapping)
    activation_profile_abs = test_scaling(S_vector,d,MRI_param,Xs_signal_norm,N_array,FR_vector_signal,t_vector,name_sorted_solution,inx_start_octv,0,0)
    activation_profile = activation_profile_abs / axons_in_path

    # scalar value if only one .mat for the whole connectome

    # first check strict side-effects, penalize if the threshold is exceeded
    for key in SE_dict:
        if d['Stim_side'] == 0 and not("_rh" in key):
            continue
        elif d['Stim_side'] == 1 and not("_lh" in key):
            continue

        activ_threshold_profile = list(SE_dict[key].keys())
        for i in range(len(activ_threshold_profile)):

            if not (activ_threshold_profile[i] in Pathways):
                print(activ_threshold_profile[i], " was not in the training set. Perhaps, it is too far from the electrode")
            else:
                inx = Pathways.index(activ_threshold_profile[i])
                if np.any(activation_profile[inx] > SE_dict[key][activ_threshold_profile[i]][0]):  # [0] - activation rate, [1] - weight of the pathway
                    return 100000000.0

    # load fixed symptom weights
    fp = open(d['symptom_weights_file'])
    symptom_weights = json.load(fp)
    symptom_weights = symptom_weights['fixed_symptom_weights']
    remaining_weights = 1.0
    N_fixed = 0
    # is it iterating only across the correct side?
    for key in symptom_weights:
        remaining_weights = remaining_weights - symptom_weights[key]
        N_fixed += 1

    # get symptom-wise difference between activation and target profiles
    from Optim_strategies import get_symptom_distances
    [__, symptom_diff, symptoms_list] = get_symptom_distances(activation_profile, profile_dict, Soft_SE_dict,
                                               symptom_weights, Pathways, d['Stim_side'], score_symptom_metric=d['similarity_metric'])
    # symptom_diff is in the symptom space, not pathway! So it will have a different dimensionality

    # compute weight for non-fixed as the equal distribution of what remained
    if len(symptoms_list) != N_fixed:
        rem_weight = remaining_weights / (len(symptoms_list) - N_fixed)
    else:
        rem_weight = 0.0

    # for now, the global score is a simple summation (exactly like in Optim_strategies.py)
    loc_ind = 0
    for key in symptoms_list:
        if d['Stim_side'] == 0 and not ("_rh" in key):
            continue
        elif d['Stim_side'] == 1 and not ("_lh" in key):
            continue

        if key in symptom_weights:
            global_score = symptom_weights[key] * symptom_diff[loc_ind]
        else:
            global_score = rem_weight * symptom_diff[loc_ind]

        loc_ind += 1

    # do clean-up in Results_
    if os.path.isdir(os.environ['PATIENTDIR']+'/Results_rh'):
        os.system('rm -fr '+os.environ['PATIENTDIR']+'/Results_rh')
        os.makedirs(os.environ['PATIENTDIR']+'/Results_rh')

    if os.path.isdir(os.environ['PATIENTDIR']+'/Results_lh'):
        os.system('rm -fr '+os.environ['PATIENTDIR']+'/Results_lh')
        os.makedirs(os.environ['PATIENTDIR']+'/Results_lh')

    # store info about current protocols
    current_data = []
    for i in range(len(S_vector)):
        current_data.append(S_vector[i])
    current_data_array = np.vstack((current_data)).T

    # store activation profiles (only simulated)
    activation_data = []
    for act in activation_profile:
        activation_data.append(act)
    activation_data_array = np.vstack((activation_data)).T

    # store non-weighted symptom distances and global score
    score_data = symptom_diff.tolist()
    score_data.append(global_score)
    score_data_array = np.vstack((score_data)).T

    # you need a pre clean-up
    if d['Stim_side'] == 0:
        with open(os.environ['PATIENTDIR']+'/NB_0/current_over_iterations.csv', 'a') as f_handle:
            np.savetxt(f_handle, current_data_array)
        with open(os.environ['PATIENTDIR']+'/NB_0/activation_over_iterations.csv', 'a') as f_handle:
            np.savetxt(f_handle, activation_data_array)
        with open(os.environ['PATIENTDIR']+'/NB_0/scores_over_iterations.csv', 'a') as f_handle:
            np.savetxt(f_handle, score_data_array)
    else:
        with open(os.environ['PATIENTDIR']+'/NB_1/current_over_iterations.csv', 'a') as f_handle:
            np.savetxt(f_handle, current_data_array)
        with open(os.environ['PATIENTDIR']+'/NB_1/activation_over_iterations.csv', 'a') as f_handle:
            np.savetxt(f_handle, activation_data_array)
        with open(os.environ['PATIENTDIR']+'/NB_1/scores_over_iterations.csv', 'a') as f_handle:
            np.savetxt(f_handle, score_data_array)

    return global_score
