'''
    By K. Butenko
    This script
    2) solves optimization problem by adjusting the scalings for the solutions
        2.1) V(P_neuron)=S1*V1(P_neuron)+S2*V2(P_neuron)+...              , where V1 is the solution for Contact1-Ground, the rest are floating
        2.2) manages IFFT
        2.3) run NEURON,get the activation rate, compare with the desired rate
    3) if the the best solution over all steps is found, the system solved with adjusted potentials to compute the currents:
        BC1=S1*1+S2*V2(C1)+S3*V3(C1)+... , where V2(C1) is the floating potential on Contact1 when Contact2-Ground is solved

'''

import numpy as np
import time
import pickle
import subprocess
import importlib
import os
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


def conduct_unit_IFFT(d, Xs_signal_norm, N_models, N_segm, FR_vector_signal, t_vector, name_sorted_solution,
                      inx_start_octv):
    # stores sunit solution (el. potential on axons in space and time) over contacts

    from IFFT_contact_ground import convolute_signal_with_field_and_compute_ifft
    if d["spectrum_trunc_method"] == 'Octave Band Method':
        if isinstance(d["n_Ranvier"], list):  # if different populations
            last_point = 0
            hf = h5py.File(os.environ['PATIENTDIR'] + '/' + d["Name_prepared_neuron_array"], 'r')
            lst_population_names = list(hf.keys())
            hf.close()
            for i in range(len(d["n_Ranvier"])):
                # print("in ",lst_population_names[i]," population")
                last_point = convolute_signal_with_field_and_compute_ifft(d, Xs_signal_norm, N_models[i], N_segm[i],
                                                                          FR_vector_signal, t_vector,
                                                                          name_sorted_solution,
                                                                          inx_st_oct = inx_start_octv, dif_axons=True,
                                                                          last_point = last_point)
        else:
            convolute_signal_with_field_and_compute_ifft(d, Xs_signal_norm, N_models, N_segm, FR_vector_signal,
                                                         t_vector, name_sorted_solution, inx_st_oct=inx_start_octv,
                                                         dif_axons=False, last_point=0)
    else:
        logging.scale("Spectrum truncation with Octave Band Method has to be enabled")
        raise SystemExit

    return True


def test_scaling(S_vector,d,Xs_signal_norm,N_models,N_segm,FR_vector_signal,t_vector,A,name_sorted_solution,inx_start_octv,VTA_IFFT,scaling_index,VTA_parameters):

    if VTA_IFFT == 1:

        VTA_edge,VTA_full_name,VTA_resolution = VTA_parameters

        from Parallel_IFFT_on_VTA_array import scale_and_get_IFFT_on_VTA_array
        Max_signal_for_point = scale_and_get_IFFT_on_VTA_array(S_vector,d["number_of_processors"],name_sorted_solution,d,FR_vector_signal,Xs_signal_norm,t_vector,N_segm,inx_start_octv)
        from VTA_from_array import get_VTA_scaled
        get_VTA_scaled(d,VTA_full_name,Max_signal_for_point,N_segm,VTA_edge,VTA_resolution,scaling_index)
        return True

    logging.critical("----- Estimating neuron activity -----")
    start_neuron=time.time()

    oss_plat_cont = os.getcwd()

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
        from Axon_files.NEURON_direct_run_scaled import run_simulation_with_NEURON
    elif d["Axon_Model_Type"] == 'Reilly2016':
        os.chdir("Reilly2016/")
        logging.critical("Please, precompile Reilly2016 and comment out the next line")
        #raise SystemExit
        with open(os.devnull, 'w') as FNULL: subprocess.call('nrnivmodl', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        from Axon_files.Reilly2016.NEURON_Reilly2016_scaled import run_simulation_with_NEURON

    if isinstance(d["n_Ranvier"],list) and len(d["n_Ranvier"])>1:
        Number_of_activated = 0
        last_point=0
        for i in range(len(d["n_Ranvier"])):
            Number_of_activated_population = run_simulation_with_NEURON(d, S_vector,last_point,i,d["diam_fib"][i],d["n_Ranvier"][i],N_models[i],d["Ampl_scale"],d["number_of_processors"],scaling_index,d["Name_prepared_neuron_array"])
            Number_of_activated = Number_of_activated+Number_of_activated_population

            #if d["Axon_Model_Type"] == 'Reilly2016':
            #    os.chdir("Reilly2016/")
            last_point=N_segm[i]*N_models[i]+last_point

        #if d["Axon_Model_Type"] == 'Reilly2016':
        #    os.chdir("..")
    else:
        if isinstance(d["diam_fib"],list):
            d["diam_fib"]=d["diam_fib"][0]
            d["n_Ranvier"]=d["n_Ranvier"][0]
        Number_of_activated = run_simulation_with_NEURON(d,S_vector,0,-1,d["diam_fib"],d["n_Ranvier"],N_models[0],d["Ampl_scale"],d["number_of_processors"],scaling_index)

    os.chdir(oss_plat_cont)

    minutes=int((time.time() - start_neuron)/60)
    secnds=int(time.time() - start_neuron)-minutes*60
    logging.critical("----- NEURON calculations took {} min {} sec -----\n".format(minutes, secnds))

    if d['Stim_side']==0:
        stim_folder='Results_rh/'
    else:
        stim_folder='Results_lh/'


    if isinstance(d["n_Ranvier"],list):
        activation_in_populations=np.zeros(len(d["n_Ranvier"]),int)
        hf = h5py.File(os.environ['PATIENTDIR']+'/'+stim_folder+'Network_status_'+str(scaling_index)+'.h5', 'r')
        lst=list(hf.keys())
        for i in range(len(lst)):
            Axon_status=hf.get(lst[i])
            num_activ_in_population=0
            Axon_status=np.array(Axon_status)
            for j in range(Axon_status.shape[0]):
                if Axon_status[j]==1:
                    num_activ_in_population+=1
            activation_in_populations[i]=num_activ_in_population
        hf.close()
    else:
        activation_in_populations=Number_of_activated

    if os.path.isdir(os.environ['PATIENTDIR']+'/Field_solutions/Activation'):     # we always re-run NEURON simulation
        os.system('rm -fr '+os.environ['PATIENTDIR']+'/Field_solutions/Activation')
        os.makedirs(os.environ['PATIENTDIR']+'/Field_solutions/Activation')

    # if os.path.isdir(os.environ['PATIENTDIR']+'/Axons_in_time'):     # we always re-run NEURON simulation
    #     os.system('rm -fr '+os.environ['PATIENTDIR']+'/Axons_in_time')
    #     os.makedirs(os.environ['PATIENTDIR']+'/Axons_in_time')


    return activation_in_populations


def compute_similarity(S_vector, *args):
    d, Xs_signal_norm, N_models, N_segm, FR_vector_signal, t_vector, A, name_sorted_solution, inx_start_octv = args

    activation_profile = test_scaling(S_vector,d,Xs_signal_norm,N_models,N_segm,FR_vector_signal,t_vector,A,name_sorted_solution,inx_start_octv,d["Full_Field_IFFT"],0,0)
    # scalar value if only one .mat for the whole connectome

    # do clean-up in Results_
    if os.path.isdir(os.environ['PATIENTDIR']+'/Results_rh'):
        os.system('rm -fr '+os.environ['PATIENTDIR']+'/Results_rh')
        os.makedirs(os.environ['PATIENTDIR']+'/Results_rh')

    if os.path.isdir(os.environ['PATIENTDIR']+'/Results_lh'):
        os.system('rm -fr '+os.environ['PATIENTDIR']+'/Results_lh')
        os.makedirs(os.environ['PATIENTDIR']+'/Results_lh')

    if len(d['optimal_profile']) == 1:
        logging.critical("Computed activation rate - Optimal rate for {}: {}".format(S_vector, activation_profile - d['optimal_profile'][0]))
        distance = abs(activation_profile - d['optimal_profile'][0])
    else:
        from scipy.spatial.distance import canberra, cityblock, euclidean, braycurtis, cosine
        if d['similarity_metric'] == 'Canberra':
            distance = canberra(activation_profile, d['optimal_profile'], w = d['profile_weighting'])
        elif d['similarity_metric'] == 'Manhattan':
            distance = cityblock(activation_profile, d['optimal_profile'], w = d['profile_weighting'])
        elif d['similarity_metric'] == 'Euclidean':
            distance = euclidean(activation_profile, d['optimal_profile'], w = d['profile_weighting'])
        elif d['similarity_metric'] == 'Cosine':
            distance = cosine(activation_profile, d['optimal_profile'], w = d['profile_weighting'])
        elif d['similarity_metric'] == 'Bray-Curtis':
            distance = braycurtis(activation_profile, d['optimal_profile'], w = d['profile_weighting'])
        else:
            logging.critical("Metric {} is not supported".format(d['similarity_metric']))
            raise SystemExit
        logging.critical("{} distance for {}: {}".format(d['similarity_metric'], S_vector, distance))

    # store info of iterations
    optim_data = []
    for i in range(len(S_vector)):
        optim_data.append(S_vector[i])
    optim_data.append(distance)
    optim_data_array = np.vstack((optim_data)).T
    with open(os.environ['PATIENTDIR']+'/current_optim_iterations.csv', 'a') as f_handle:
        np.savetxt(f_handle, optim_data_array)

    return distance


def find_activation(current_comb,d,Xs_signal_norm,N_models,N_segm,FR_vector_signal,t_vector,A,name_sorted_solution,inx_start_octv,scaling_index,VTA_param=0):

    import time
    start_current_run=time.time()

    activation=test_scaling(current_comb,d,Xs_signal_norm,N_models,N_segm,FR_vector_signal,t_vector,A,name_sorted_solution,inx_start_octv,d["Full_Field_IFFT"],scaling_index,VTA_param)

    minutes=int((time.time() - start_current_run)/60)
    secnds=int(time.time() - start_current_run)-minutes*60
    logging.critical("----- Solved for current protocol in {} min {} sec -----\n".format(minutes, secnds))

    return activation


