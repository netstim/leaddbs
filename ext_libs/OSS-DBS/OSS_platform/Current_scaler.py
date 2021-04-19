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

import time

import shutil

def test_scaling(S_vector,d,Xs_signal_norm,N_models,N_segm,FR_vector_signal,t_vector,A,name_sorted_solution,inx_start_octv,VTA_IFFT,scaling_index,VTA_parameters):

    #d,Xs_signal_norm,N_models,N_segm,FR_vector_signal,t_vector,A,name_sorted_solution,inx_start_octv=args

    if VTA_IFFT==1:

        VTA_edge,VTA_full_name,VTA_resolution=VTA_parameters

        from Parallel_IFFT_on_VTA_array import scale_and_get_IFFT_on_VTA_array
        Max_signal_for_point=scale_and_get_IFFT_on_VTA_array(S_vector,d["number_of_processors"],name_sorted_solution,d,FR_vector_signal,Xs_signal_norm,t_vector,N_segm,inx_start_octv)
        from VTA_from_array import get_VTA_scaled
        get_VTA_scaled(d,VTA_full_name,Max_signal_for_point,N_segm,VTA_edge,VTA_resolution,scaling_index)
        return True

    from IFFT_with_scaling_factor import convolute_signal_with_field_and_compute_ifft
    if d["spectrum_trunc_method"]=='Octave Band Method':
        if isinstance(d["n_Ranvier"],list):             #if different populations
            last_point=0
            hf = h5py.File(os.environ['PATIENTDIR']+'/'+d["Name_prepared_neuron_array"], 'r')
            lst_population_names=list(hf.keys())
            hf.close()
            for i in range(len(d["n_Ranvier"])):
                #print("in ",lst_population_names[i]," population")
                last_point=convolute_signal_with_field_and_compute_ifft(S_vector,d,Xs_signal_norm,N_models[i],N_segm[i],FR_vector_signal,t_vector,A,name_sorted_solution,inx_st_oct=inx_start_octv,dif_axons=True,last_point=last_point)
        else:
            convolute_signal_with_field_and_compute_ifft(S_vector,d,Xs_signal_norm,N_models,N_segm,FR_vector_signal,t_vector,A,name_sorted_solution,inx_st_oct=inx_start_octv,dif_axons=False,last_point=0)
    else:
        print("Not implemented")
        raise SystemExit


    print("----- Estimating neuron activity -----")
    start_neuron=time.time()

    if d["Axon_Model_Type"] == 'McIntyre2002':
        os.chdir("Axon_files/")
        with open(os.devnull, 'w') as FNULL: subprocess.call('nocmodl axnode.mod', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        with open(os.devnull, 'w') as FNULL: subprocess.call('nrnivmodl axnode', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        from Axon_files.NEURON_direct_run_scaled import run_simulation_with_NEURON
    elif d["Axon_Model_Type"] == 'Reilly2016':
        os.chdir("Axon_files/Reilly2016/")
        print("Please, precompile Reilly2016 and comment out the next line")
        raise SystemExit
        #with open(os.devnull, 'w') as FNULL: subprocess.call('nrnivmodl', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        from Axon_files.Reilly2016.NEURON_Reilly2016_scaled import run_simulation_with_NEURON

    if isinstance(d["n_Ranvier"],list) and len(d["n_Ranvier"])>1:
        Number_of_activated=0
        last_point=0
        for i in range(len(d["n_Ranvier"])):
            Number_of_activated_population=run_simulation_with_NEURON(last_point,i,d["diam_fib"][i],1000*d["t_step"],1000.0/d["freq"],d["n_Ranvier"][i],N_models[i],d["v_init"],t_vector.shape[0],d["Ampl_scale"],d["number_of_processors"],d['Stim_side'],scaling_index,d["Name_prepared_neuron_array"])
            Number_of_activated=Number_of_activated+Number_of_activated_population

            os.chdir("Axon_files/")
            if d["Axon_Model_Type"] == 'Reilly2016':
                os.chdir("Reilly2016/")
            last_point=N_segm[i]*N_models[i]+last_point
        os.chdir("..")
        if d["Axon_Model_Type"] == 'Reilly2016':
            os.chdir("..")
    else:
        if isinstance(d["diam_fib"],list):
            d["diam_fib"]=d["diam_fib"][0]
            d["n_Ranvier"]=d["n_Ranvier"][0] 
        Number_of_activated=run_simulation_with_NEURON(0,-1,d["diam_fib"],1000*d["t_step"],1000.0/d["freq"],d["n_Ranvier"],N_models,d["v_init"],t_vector.shape[0],d["Ampl_scale"],d["number_of_processors"],d['Stim_side'],scaling_index)

    if isinstance(d["n_Ranvier"],list) and len(d["n_Ranvier"])>1:
        with open(os.devnull, 'w') as FNULL: subprocess.call('python Visualization_files/Paraview_connections_activation.py', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)


    minutes=int((time.time() - start_neuron)/60)
    secnds=int(time.time() - start_neuron)-minutes*60
    print("----- NEURON calculations took ",minutes," min ",secnds," s -----\n")

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


    if os.path.isdir(os.environ['PATIENTDIR']+'/Axons_in_time'):     # we always re-run NEURON simulation
        os.system('rm -fr '+os.environ['PATIENTDIR']+'/Axons_in_time')
        os.makedirs(os.environ['PATIENTDIR']+'/Axons_in_time')


    return activation_in_populations

#def adapt_with_annealing(lb,ub,max_iterations,args_all):
#
#    from scipy.optimize import dual_annealing
#
#    if os.path.isfile('Best_scaling_yet.csv'):
#        #initial_scaling=np.genfromtxt('Best_scaling_yet.csv', delimiter=' ')
#        #initial_scaling=[-7.818776816129684448e-01,-1.069913357496261597e+00,2.497921139001846313e-01,-1.049459293484687805e+00]
#        #initial_scaling=[3.374655917286872864e-01,-1.232426971197128296e+00,-3.297385610640048981e-01,-1.211972907185554504e+00]
#        #initial_scaling=[1.126232415437698364e+00,1.097056448459625244e+00,3.340458869934082031e-01,1.117510512471199036e+00] #50,293 dif
#        initial_scaling=[1.316384516656398773e+00,2.556199058890342712e-01,8.183348253369331360e-01,-5.758354514837265015e-01] #48,437 dif
#        res = dual_annealing(test_scaling, bounds=list(zip(lb, ub)),args=args_all,maxfun=max_iterations, seed=19234394,x0=initial_scaling, no_local_search=True) #number of iterations by counter
#    else:
#        res = dual_annealing(test_scaling, bounds=list(zip(lb, ub)),args=args_all,maxfun=max_iterations, initial_temp=25000.0, visit=2.85, seed=int(np.random.random()*1000), no_local_search=True) #number of iterations by counter
#
#    print("Optimization results: ",res.x,res.fun)
#    return(res.x,res.fun)

def find_activation(current_comb,d,Xs_signal_norm,N_models,N_segm,FR_vector_signal,t_vector,A,name_sorted_solution,inx_start_octv,scaling_index,VTA_param=0):

    import time
    start_current_run=time.time()

    activation=test_scaling(current_comb,d,Xs_signal_norm,N_models,N_segm,FR_vector_signal,t_vector,A,name_sorted_solution,inx_start_octv,d["Full_Field_IFFT"],scaling_index,VTA_param)

    minutes=int((time.time() - start_current_run)/60)
    secnds=int(time.time() - start_current_run)-minutes*60
    print("----- Solved for current protocol in: ",minutes," min ",secnds," s -----\n")

    return activation


