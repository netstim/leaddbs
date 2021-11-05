#    Copyright (C) 2019 Konstantin Butenko, konstantin.butenko@uni-rostock.de
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.



'''
#@author: K.Butenko (konstantin.butenko@uni-rostock.de), routines for Salome prepared in cooperation with P.Trieu. Script "tissue dielectrics.py" and axon files are adapted from FanPy (https://arxiv.org/abs/1705.10478)
GUI tree is developed in cooperation with S.A.Udongwo
The platform is in beta test and should be used with caution and a proper amount of healthy suspicion.
'''

'''This is the main script that governs the simulation flow

It contains functions:

run_master_study(), which is useful when the platform is run interatively to solve UQ or optimization problems.
It will evaluate which truncation method and with how many frequencies we need to use for a particular setup.
It can also estimate the difference in potential for different profile salome scripts to estimate the initial meshing quality (profile scripts should be created manually).
Check out the function before you use it. By default, it tries to match absolutely the result obtained with the full spectrum, which is normally unnecessary if you have thousands of neurons.

run_full_model(master_dict) governs the simulation flow taking master_dict as the input that will modify the input dictionary (defined in GUI_inp_dict.py)

'''
print("\nOSS-DBS by K.Butenko --- version 0.4")
print("Butenko K, Bahls C, Schroeder M, Koehling R, van Rienen U (2020) 'OSS-DBS: Open-source simulation platform for deep brain stimulation with a comprehensive automated modeling.' PLoS Comput Biol 16(7): e1008023. https://doi.org/10.1371/journal.pcbi.1008023")
print("____________________________________\n")

import numpy as np
import time
import pickle
import subprocess
import importlib
import os
import warnings
import json
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py

import shutil
import logging

def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)

def run_full_model(master_dict):

    start_simulation_run = time.time()

    import os
    os.environ['PATIENTDIR'] = '/opt/Patient' # Use fixed mount path for docker

    #===================Load and update input dictionary======================#
    #should be loaded this way for iterative studies (some simulation state variables change during a run)

    # should be changed: load directly from the patient folder a json

    import sys
    sys.path.insert(1, os.environ['PATIENTDIR'])
    import GUI_inp_dict
    importlib.reload(GUI_inp_dict)
    from GUI_inp_dict import d

    # all implicit parameters (not controlled in OSS-DBS GUI) should be passed like this
    with open(os.environ['PATIENTDIR'] + '/Lead_DBS_input.json', 'r') as fp:
        lead_dict = json.load(fp)
    fp.close()
    # IMPORTANT: do not just update the dictionary, some parameters are explicitly changed upstream
    d['stretch'] = lead_dict['stretch']
    logging.critical('Electrode array stretch: {}'.format(d['stretch']))
    d['StimSets'] = lead_dict['StimSets']

    log_file = open(os.environ['PATIENTDIR']+'/log_file_hemi_' + str(d["Stim_side"]) + '.log', 'w')

    from Dict_corrector import rearrange_Inp_dict
    d=rearrange_Inp_dict(d)             #misc. transformation of parameters to the platform's format
    d.update(master_dict)               #modifies the user provided input dictionary (e.g. for UQ study), check run_master_study() function . Warning: this does not work update the encap. layer properties and the solver during adaptive mesh refiment, because these data are realoaded from the original dictionary

    # this is a block dedicated to a new module 'Optimizer'
    # this is the internal loop of optimization (current protocol for a given position)
    d['Optimizer'] = 0
    if d['Optimizer'] == 1:  # run simulated annealing (no local search) from scipy on the unit current solutions, i.e., just recompute the NEURON model
                             # but beware that the mesh is not specifically refined for all protocols

        # just an example how it will look like (will be imported and transformed from Lead-DBS)
        d['num_iterations'] = 5
        d['min_bound_per_contact'] = [-0.003, -0.003, -0.003, -0.003]   # in A! same length as the electrode
        d['max_bound_per_contact'] = [0.003, 0.003, 0.003, 0.003]  # same length as the electrode
        d['optimal_profile'] = [0.5, 0.4, 0.3, 0.4, 0.5]    # activation rates (1.0 = 100%) as many as .mat files you have for the connectome
                                                            # if one value, provide it still in list

        # this snippet is only relevant if you have a profile, not a single .mat
        d['profile_weighting'] = [1.0, 1.0, 1.0, 1.0, 1.0]  # weighting for pathways in the profile, not for symptoms!
        d['similarity_metric'] = 'Canberra'  # Criterion to evaluate how similar two activation profiles are. Supported: Bray-Curtis, Canberra, Manhattan, Euclidean, Cosine

        if len(d['min_bound_per_contact']) != len(d['Phi_vector']):
            logging.critical('Error: the length of the bound vector does not match the number of contacts')
            raise SystemExit

        if d["Full_Field_IFFT"] == 1:
            logging.critical("Optimization is yet not supported for VTA from E-field/Rattay's function")
            raise SystemExit

        d["current_control"], cc_multicontact, d["Skip_mesh_refinement"], d["external_grounding"], d["EQS_core"] = (1, 1, 1, True, "QS")
        logging.critical("At the moment, optimization is limited to QS formulation of current-controlled mode without mesh refinement")
        logging.critical("When running optimization, grounding is fixed to the casing (but optimization might create 'pseudo-grounding' contacts)")

        d['Phi_vector'] = [1.0] * len(d['Phi_vector'])  # unit vector

        # do we need an initial guess?

    # This is a block dedicated to a new functionality: charge-balancing
    d['Charge_balancing'] = False  # only required for symmetric balancing, but can be also used for low amplitude
    if d['Charge_balancing'] == True:
        d['Balancing_type'] = 'Symmetric' # Symmetric - Counter pulse is identical to the DBS pulse and follows right after
                                          # Low_amplitude - always a rectangular counter pulse that has a pulse width = 1/(DBS_Freq) - DBS_pulse_width, and the amplitude is adjusted accordingly to have the same charge

    if d['StimSets'] == 1 and os.path.isfile(os.environ['PATIENTDIR']+'/Current_protocols_'+str(d['Stim_side'])+'.csv'):
        stim_protocols = np.genfromtxt(os.environ['PATIENTDIR']+'/Current_protocols_'+str(d['Stim_side'])+'.csv', dtype=float, delimiter=',', names=True)
        if stim_protocols.size:
            d['Current_sets']=True
            d["Skip_mesh_refinement"]=1
            logging.critical("When testing different current set, adaptive refinement is unavailable, make sure the mesh is prerefined")
            d["EQS_core"]="QS"
            logging.critical("When testing different current set, only QS formulation is currently available")
            cc_multicontact=True
            d["spectrum_trunc_method"]="Octave Band Method"
            logging.critical("When testing different current set, only Octave Band Method is currently available")
            d['Phi_vector']=[1.0] * len(d['Phi_vector'])          # unit vector
            logging.critical("When testing different current sets, grounding is fixed to the casing (but you can imitate grounding by assigning a value of -1.0*sum(Icontacts)) to one of the contacts")
            d["external_grounding"]=True
            d["current_control"]=1

            if d["Full_Field_IFFT"] == 1:
                logging.critical("Field superposition is yet not supported for VTA from E-field/Rattay's function")
                raise SystemExit

            import math

            Currents_to_check=[]
            for i in range(stim_protocols.shape[0]):
                stim_prot=list(stim_protocols[i])
                for j in range(len(stim_prot)):
                    if math.isnan(stim_prot[j]):
                        stim_prot[j]=None
                    elif stim_prot[j]==0.0:
                        logging.critical('0.0 always refers to grounding in OSS-DBS. Please, type "passive" or "float" for contacts that do not deliver currents.')
                        raise SystemExit
                    else:
                        stim_prot[j]=stim_prot[j]*0.001            # Lead-DBS stores in mA
                if len(d['Phi_vector']) != len(stim_prot):
                    logging.critical("Current protocols do not match the number of contacts on the electrode, exiting")
                    raise SystemExit
                Currents_to_check.append(stim_prot)
        else:
            d['Current_sets']=False
    else:
        d['Current_sets']=False

    #=========Check the simulation setup and state, load the corresponding data=========#

    if d["Stim_side"] == 0:
        logging.critical("Processing right hemisphere\n")
    else:
        logging.critical("Processing left hemisphere\n")

    if d["current_control"]==1 and d["CPE_activ"]==1:
        d["CPE_activ"]=0
        logging.critical("Disabling CPE for current-controlled simulation")

    Phi_vector_active_non_zero = [x for x in d["Phi_vector"] if (x is not None) and (x != 0.0)]
    cc_multicontact = False
    if d["current_control"] == 1 and len(Phi_vector_active_non_zero) > 1:       #multicontact current-controlled case
        cc_multicontact = True

    from Sim_state import check_state
    check_state(d)      #will switch the state of simulation depending on the corresp. data in the input_dictionary, will also manage project folders

    if d["voxel_arr_MRI"] == 0 and d["voxel_arr_DTI"] == 1:
        logging.critical("MRI data is new, the DTI data will be reprocessed")
        d["voxel_arr_DTI"] = 0

    #loading of meta data depending on the simulatation setup and state
    import os

    if d["Init_mesh_ready"]==1:
        with open(os.environ['PATIENTDIR']+'/Meshes/Mesh_ind.file', "rb") as f:
            Domains = pickle.load(f)

    #========================Formating MRI and DTI data=======================#
    from MRI_DTI_prep_new import obtain_MRI_class
    MRI_param=obtain_MRI_class(d)       #also creates 'MRI_DTI_derived_data/Tissue_array_MRI.csv' and meta data

    anisotrop,DTI_param=(0,0)          #initizialization, 0 by default for isotropic

    if d["DTI_data_name"]!=0:	#0 if we do not provide DTI data
        anisotrop=1         #the study will use tensor conductivity
        from MRI_DTI_prep_new import obtain_DTI_class
        DTI_param=obtain_DTI_class(d,MRI_param) #also creates 'MRI_DTI_derived_data/Tensor_array_DTI.csv' and meta data

    #========================Geometry generation==============================#

    if d["Full_Field_IFFT"] == 1:       # for this case, we use a neuron array that matches the VTA array in dimensions
        d['Axon_Model_Type']='Reilly2016'
        d['x_seed'],d['y_seed'],d['z_seed']=(d['Implantation_coordinate_X'],d['Implantation_coordinate_Y']+3.0,d['Implantation_coordinate_Z']+5.0)  # it makes sense to shift it a bit from the tip
        d['diam_fib']=5.0
        d['n_Ranvier']=22
        d['x_step'],d['y_step'],d['z_step']=(1.0,1.0,1.0)
        d['x_steps'],d['y_steps'],d['z_steps']=(20.0,0.0,20.0)  #we assume that Z-axis is ventra-dorsal in the MRI
        d['Global_rot']=1
        d['alpha_array_glob']=[0]
        d['beta_array_glob']=[0]
        d['gamma_array_glob']=[0]
        if d['Electrode_type']=="AA_rodent_monopolar" or d['Electrode_type']=="SR_rodent":    #rodent VTA
            d['Axon_Model_Type']='Reilly2016'
            d['x_seed'],d['y_seed'],d['z_seed']=(d['Implantation_coordinate_X'],d['Implantation_coordinate_Y'],d['Implantation_coordinate_Z'])  # it makes sense to shift it a bit from the tip
            d['diam_fib']=5.0
            d['n_Ranvier']=3
            d['x_step'],d['y_step'],d['z_step']=(0.1,0.1,0.1)
            d['x_steps'],d['y_steps'],d['z_steps']=(20.0,0.0,20.0)  #we assume that Z-axis is ventra-dorsal in the MRI
            d['Global_rot']=1
            d['alpha_array_glob']=[0]
            d['beta_array_glob']=[0]
            d['gamma_array_glob']=[0]


    if d["Init_mesh_ready"] == 0:

        if d["Brain_shape_name"] == 0:   #Creates a brain approximation (ellisploid)
            from CAD_Salome import build_brain_approx
            x_length,y_length,z_length = build_brain_approx(d, MRI_param)      #also creates 'brain_subsitute.brep'
            Brain_shape_name = 'Brain_substitute.brep'

        if d["Brain_shape_name"] != 0:
            Brain_shape_name = d["Brain_shape_name"]

        #==================Initial neuron array generation====================#
        from Neural_array_processing import Neuron_array
        N_array = Neuron_array(d, MRI_param)

        if d["Init_neuron_model_ready"]==0 and d["Neuron_model_array_prepared"]==0:
            logging.critical("----- Creating initial neuron array -----")

            N_array.build_neuron_models() #builds a pattern model, if not provided, then builds a neuron array and stores in 'Neuron_model_arrays/All_neuron_models.csv'
            with open(os.environ['PATIENTDIR']+'/Neuron_model_arrays/Neuron_array_class.file', "wb") as f:
                pickle.dump(N_array, f, pickle.HIGHEST_PROTOCOL)
        elif d["Neuron_model_array_prepared"]==1 and d["Init_neuron_model_ready"]==0:
            logging.critical("----- Creating initial neuron array from a provided neuron array -----")
            N_array.process_external_array()   # adjusts the prepared neuron array to the computational domain (only for brain substitutes!), then stores the nueron array in 'Neuron_model_arrays/All_neuron_models.csv'
            with open(os.environ['PATIENTDIR']+'/Neuron_model_arrays/Neuron_array_class.file', "wb") as f:
                pickle.dump(N_array, f, pickle.HIGHEST_PROTOCOL)
        else:
            logging.critical("--- Initial neuron array was loaded\n")
            with open(os.environ['PATIENTDIR']+'/Neuron_model_arrays/Neuron_array_class.file', "rb") as f:
                N_array = pickle.load(f)

        #if brain substitute is used, it will be enlarged to encompass previosly defined neuron array (if necessary)
        if Brain_shape_name=='Brain_substitute.brep':
            needs_a_rebuid=0
            if N_array.ROI_radius > (x_length/2):
                d['Approximating_Dimensions'][0] = N_array.ROI_radius*2+0.1
                logging.critical("increasing length along x to encompass the neuron array\n")
                needs_a_rebuid=1
            if N_array.ROI_radius > (y_length/2):
                d['Approximating_Dimensions'][1] = N_array.ROI_radius*2+0.1
                logging.critical("increasing length along y to encompass the neuron array\n")
                needs_a_rebuid=1
            if N_array.ROI_radius > (z_length/2):
                d['Approximating_Dimensions'][2] = N_array.ROI_radius*2+0.1
                logging.critical("increasing length along z to encompass the neuron array\n")
                needs_a_rebuid=1
            if needs_a_rebuid == 1:
                x_length,y_length,z_length=build_brain_approx(d,MRI_param)      #also creates 'brain_subsitute.brep'
                logging.critical("\n")
            if N_array.ROI_radius > min((x_length/2),(y_length/2),(z_length/2)):
                logging.critical("ROI_radius: ", N_array.ROI_radius)
                logging.critical("ROI is still bigger than the computational domain.")
                raise SystemExit

        #===================Final geometry generation=========================#
        from CAD_Salome import build_final_geometry
        Domains = build_final_geometry(d,MRI_param,Brain_shape_name,N_array.ROI_radius,cc_multicontact)       #creates and discretizes the geometry with the implanted electrode, encapsulation layer and ROI, converts to the approp. format. The files are stored in Meshes/

    #===============Adjusting neuron array====================================#

    if d["Adjusted_neuron_model_ready"] == 0:

        if d["Init_mesh_ready"]==1:
            with open(os.environ['PATIENTDIR']+'/Neuron_model_arrays/Neuron_array_class.file', "rb") as f:
                N_array = pickle.load(f)

        # subtracts neurons from the previously define All_neuron_models.csv, if the are non-physical (inside encap. layer or CSF, outside of the domain, intersect with the electrode geometry.)
        N_array.adjust_neuron_models(Domains, MRI_param)
        with open(os.environ['PATIENTDIR']+'/Neuron_model_arrays/Neuron_array_class.file', "wb") as f:
            pickle.dump(N_array, f, pickle.HIGHEST_PROTOCOL)
    else:
        logging.critical("--- Adjusted neuron array was loaded\n")
        with open(os.environ['PATIENTDIR']+'/Neuron_model_arrays/Neuron_array_class.file', "rb") as f:
            N_array = pickle.load(f)

    number_of_points = int(np.sum(N_array.pattern['num_segments'] * N_array.N_models))


    if d["Full_Field_IFFT"] == 1:       # rename later

        from VTA_from_array import create_VTA_array,resave_as_verts,get_VTA
        VTA_edge,VTA_full_name,VTA_resolution= create_VTA_array(d['x_seed'],d['y_seed'],d['z_seed'],d['Electrode_type'])
        arrays_shape = resave_as_verts(VTA_full_name)
        number_of_points=sum(arrays_shape)
        VTA_parameters=[VTA_edge,VTA_full_name,VTA_resolution]

#=============================Signal creation=================================#
    #in case there was an update of the signal
    Phi_vector=[x for x in d["Phi_vector"] if x is not None] #now we don't need None values (floating potentials), Contacts contain only the active ones
    Domains.fi=Phi_vector
    with open(os.environ['PATIENTDIR']+'/Meshes/Mesh_ind.file', "wb") as f:
        pickle.dump(Domains, f, pickle.HIGHEST_PROTOCOL)


    #logging.critical('Domains fi: ',Domains.fi)
    if d["current_control"]==1 and cc_multicontact==False:     #if multicontact, then it will be scaled in parallel comp.
        A=max(Domains.fi[:], key=abs)
        Phi_max=1      #not needed
    else:
        A=1.0 #for vc we will apply Phi_vector as boundary conditions with appropriate amplitude
        Phi_max=max(Domains.fi[:], key=abs)

    if d["signal_generation_ready"]==0:
        logging.critical("----- Generating DBS signal -----")

        from Signal_generator import generate_signal
        [t_vector,signal_out,Xs_signal_norm,FR_vector_signal]=generate_signal(d,A,Phi_max,cc_multicontact)

        Xs_storage=np.zeros((np.real(Xs_signal_norm).shape[0],2),float)

        Xs_storage[:,0]=np.real(Xs_signal_norm)
        Xs_storage[:,1]=np.imag(Xs_signal_norm)
        np.savetxt(os.environ['PATIENTDIR']+'/Stim_Signal/Xs_storage_full.csv', Xs_storage, delimiter=" ")
        np.savetxt(os.environ['PATIENTDIR']+'/Stim_Signal/t_vector.csv', t_vector, delimiter=" ")
        np.savetxt(os.environ['PATIENTDIR']+'/Stim_Signal/FR_vector_signal.csv', FR_vector_signal, delimiter=" ")
    else:
        logging.critical("--- DBS signal is taken from the previous simulation\n")
        Xs_recovered = np.genfromtxt(os.environ['PATIENTDIR']+'/Stim_Signal/Xs_storage_full.csv', delimiter=' ')
        Xs_signal_norm=np.vectorize(complex)(Xs_recovered[:,0],Xs_recovered[:,1])
        FR_vector_signal = np.genfromtxt(os.environ['PATIENTDIR']+'/Stim_Signal/FR_vector_signal.csv', delimiter=' ')
        t_vector = np.genfromtxt(os.environ['PATIENTDIR']+'/Stim_Signal/t_vector.csv', delimiter=' ')


#=============================CSF_refinement==================================#

    if d["Parallel_comp_ready"]==0:
        if d["external_grounding"]==1:
            d["external_grounding"]==True
            from Ext_ground_preref import refine_external_ground
            refine_external_ground(Domains)

    '''if we want to skip CSF and adaptive mesh refinement'''
    if d["Skip_mesh_refinement"]==1:
        from CSF_refinement_new import Dummy_CSF
        Dummy_CSF()         #will resave 'Meshes/Mesh_unref.xml' as 'Results_adaptive/mesh_adapt.xml.gz' and the same for subdomains and boundaries
        d["CSF_mesh_ready"]=1
        d["Adapted_mesh_ready"]=1
        logging.critical("CSF and adaptive mesh refinement was skipped\n")
    else:
        #from Signal_generator import pick_refinement_freqs
        if d["refinement_frequency"][0]==-1:   # frequencies basing on the power spectrum distribution
            from Signal_generator import pick_refinement_freqs
            ref_freqs,Xs_on_ref_freqs=pick_refinement_freqs(FR_vector_signal,Xs_signal_norm,d["num_ref_freqs"],A)  #A is needed to switch back to unit pulse. We use unit pulse in adaptive to avoid confusion (math modules always scale the signal to the required value)
        else:
            ref_freqs=d["refinement_frequency"]

    if d["CSF_mesh_ready"]==0:
        #from CSF_refinement_MPI import launch_CSF_refinement
        from CSF_refinement_new import launch_CSF_refinement
        Scaling_CSF=launch_CSF_refinement(d,MRI_param,DTI_param,Domains,anisotrop,cc_multicontact,ref_freqs)

        from Visualization_files.paraview_find_arrayname import get_Para_Array_name
        from Parameter_insertion import insert_f_name
    elif d["Skip_mesh_refinement"]==0:      #should be changed
        Scaling_CSF=np.genfromtxt(os.environ['PATIENTDIR']+'/CSF_ref/Scaling_CSF.csv', delimiter=' ')


#====================Mesh adaptive algorithm==================================#

    if d["Adapted_mesh_ready"]==0:
        from Mesh_adaption_hybrid import mesh_adapter
        Ampl_on_vert=mesh_adapter(MRI_param,DTI_param,Scaling_CSF,Domains,d,anisotrop,cc_multicontact,ref_freqs)     #also saves adapted mesh
        np.savetxt(os.environ['PATIENTDIR']+'/Results_adaptive/Ampl_on_vert.csv', Ampl_on_vert, delimiter=" ")    #just to check

#===================Truncate the frequency spectrum===========================#
    if d["spectrum_trunc_method"]!='No Truncation':
        from Truncation_of_spectrum import get_freqs_for_calc
        FR_vector_signal_new,add_trunc_data=get_freqs_for_calc(d,FR_vector_signal,Xs_signal_norm,t_vector)
        if d["spectrum_trunc_method"]=='Octave Band Method':
            inx_start_octv=add_trunc_data
        else:
            Xs_signal_norm_new=add_trunc_data

#==========Calculate freq in parallel and rearrange field array===============#

    if d['Current_sets'] == True or d['Optimizer'] == 1:

        if d["Parallel_comp_ready"]==0:

            if ["Parallel_comp_interrupted"]==1:
                import os
                if not (os.path.isfile(os.environ['PATIENTDIR']+'/Field_solutions/Phi_real_scaled_'+str(d["freq"])+'Hz.pvd') or os.path.isfile(os.environ['PATIENTDIR']+'Field_solutions/Phi_real_unscaled_'+str(d["freq"])+'Hz.pvd')):     #to make sure that there were interrupted computations
                    logging.critical("There were no previous computations, 'Parallel_comp_interrupted' is put to 0")
                    ["Parallel_comp_interrupted"] == 0

            from Parallel_unit_current_calc import calculate_in_parallel
            '''calculate_in_parallel will save a sorted_solution array if IFFT is pointwise or will save the whole field for each frequency in Field_solutions_functions/ if full_IFFT is requested'''

            if d["spectrum_trunc_method"]=='No Truncation':
                logging.critical("----- Calculating electric field in the frequency spectrum -----")
                calculate_in_parallel(d,FR_vector_signal,Domains,MRI_param,DTI_param,anisotrop,number_of_points,cc_multicontact)

            if d["spectrum_trunc_method"]=='High Amplitude Method' or d["spectrum_trunc_method"]=='Cutoff Method' or d["spectrum_trunc_method"]=='Octave Band Method':
                logging.critical("----- Calculating electric field in the truncated frequency spectrum -----")
                calculate_in_parallel(d,FR_vector_signal_new,Domains,MRI_param,DTI_param,anisotrop,number_of_points,cc_multicontact)
        else:
            logging.critical("--- Results of calculations in the frequency spectrum were loaded\n")

        if d["spectrum_trunc_method"]=='No Truncation':
            name_sorted_solution=os.environ['PATIENTDIR']+'/Field_solutions/sorted_solution_per_contact.csv'
        else:
            name_sorted_solution=os.environ['PATIENTDIR']+'/Field_solutions/sorted_solution_per_contact_'+str(d["spectrum_trunc_method"])+'_'+str(d["trunc_param"])+'.csv'

        d['number_of_processors'] = d['number_of_processors'] * 2  # use multithreading
        logging.critical("Switching to multithreading")

        if d["IFFT_ready"] == 0 and d["Full_Field_IFFT"] != 1:

            if os.path.isdir(os.environ['PATIENTDIR']+'/Axons_in_time'):     # we always re-run NEURON simulation
                os.system('rm -fr '+os.environ['PATIENTDIR']+'/Axons_in_time')
                os.makedirs(os.environ['PATIENTDIR']+'/Axons_in_time')

            from Current_scaler import conduct_unit_IFFT
            conduct_unit_IFFT(d, Xs_signal_norm, N_array.N_models,N_array.pattern['num_segments'], FR_vector_signal, t_vector, name_sorted_solution,
                              inx_start_octv)

        if d["Full_Field_IFFT"] == 1:       # rename later
            N_array.pattern['num_segments'] = arrays_shape
            logging.critical('VTA approximation for current sets/optimization has to be enabled, contact the developers')
        else:
            VTA_parameters=0


        if d['Current_sets'] == True:
            it_num = 0
            for current_comb in Currents_to_check:
                from Current_scaler import find_activation
                Activation=find_activation(current_comb,d,Xs_signal_norm,N_array.N_models,N_array.pattern['num_segments'],FR_vector_signal,t_vector,A,name_sorted_solution,inx_start_octv,it_num,VTA_param=VTA_parameters)
                it_num += 1
        elif d['Optimizer'] == 1:
            logging.critical('Running optimization of current protocol')
            from scipy.optimize import dual_annealing
            from Current_scaler import compute_similarity

            args_all = [d, Xs_signal_norm, N_array.N_models,
                        N_array.pattern['num_segments'], FR_vector_signal, t_vector, A, name_sorted_solution,
                        inx_start_octv]
            if os.path.isfile(os.environ['PATIENTDIR']+'/Best_scaling_yet.csv'):
                initial_scaling = np.genfromtxt(os.environ['PATIENTDIR']+'/Best_scaling_yet.csv', delimiter = ' ')
                res = dual_annealing(compute_similarity,bounds=list(zip(d['min_bound_per_contact'],d['max_bound_per_contact'])), args = args_all, x0 = initial_scaling, maxfun = d['num_iterations'], seed = 42, visit = 2.62, no_local_search = True)
            else:
                res = dual_annealing(compute_similarity,bounds=list(zip(d['min_bound_per_contact'],d['max_bound_per_contact'])), args = args_all, maxfun = d['num_iterations'], seed = 42, visit = 2.62, no_local_search = True)

            # if you have outer loop for electrode placement optimization, you would want to compare res.fun against the previous best solution here
            np.savetxt(os.environ['PATIENTDIR']+'Best_scaling_yet.csv', res.x, delimiter = ' ')
            opt_res = np.append(res.x,res.fun)
            logging.critical("Current optimization results: {}".format(' '.join(map(str, list(opt_res)))))

        if d["Stim_side"] == 0:
            subprocess.call(['touch', os.environ['PATIENTDIR']+'/success_rh.txt'])
        else:
            subprocess.call(['touch', os.environ['PATIENTDIR']+'/success_lh.txt'])

        return True

    if d["Parallel_comp_ready"] == 0:
        if ["Parallel_comp_interrupted"] == 1:
            if not (os.path.isfile(os.environ['PATIENTDIR']+'/Field_solutions/Phi_real_scaled_'+str(d["freq"])+'Hz.pvd') or os.path.isfile(os.environ['PATIENTDIR']+'/Field_solutions/Phi_real_unscaled_'+str(d["freq"])+'Hz.pvd')):     #to make sure that there were interrupted computations
                logging.critical("There were no previous computations, 'Parallel_comp_interrupted' is put to 0")
                ["Parallel_comp_interrupted"] == 0

        from Parallel_field_calc import calculate_in_parallel
        '''calculate_in_parallel will save a sorted_solution array if IFFT is pointwise or will save the whole field for each frequency in Field_solutions_functions/ if full_IFFT is requested'''

        if d["spectrum_trunc_method"]=='No Truncation':
            logging.critical("----- Calculating electric field in the frequency spectrum -----")
            calculate_in_parallel(d,FR_vector_signal,Domains,MRI_param,DTI_param,anisotrop,number_of_points,cc_multicontact)

        if d["spectrum_trunc_method"]=='High Amplitude Method' or d["spectrum_trunc_method"]=='Cutoff Method' or d["spectrum_trunc_method"]=='Octave Band Method':
            logging.critical("----- Calculating electric field in the truncated frequency spectrum -----")
            calculate_in_parallel(d,FR_vector_signal_new,Domains,MRI_param,DTI_param,anisotrop,number_of_points,cc_multicontact)
    else:
        logging.critical("--- Results of calculations in the frequency spectrum were loaded\n")

    if d["spectrum_trunc_method"]=='No Truncation':
        name_sorted_solution=os.environ['PATIENTDIR']+'/Field_solutions/sorted_solution.csv'
    else:
        name_sorted_solution=os.environ['PATIENTDIR']+'/Field_solutions/sorted_solution_'+str(d["spectrum_trunc_method"])+'_'+str(d["trunc_param"])+'.csv'


        '''sorted_solution structure is'''
        '''''''''x y z Phi_r Phi_im Freq'''
        '''pnt1'''                '''130'''
        '''...'''                 '''260'''
        '''pnt1'''                '''1300000'''
        '''pnt2'''                '''130'''
        '''...'''
        '''Points are defined by the neuron array and might be duplicated. The order of points depends on the neural models'''

        '''If we enable Full_Field_IFFT, sorted_solution will be not be created'''

#=============================Full Field IFFT=================================#
    if d["IFFT_ready"] == 1:
        logging.critical("--- Results of IFFT (FFEM) were loaded\n")

    #if IFFT_on_VTA_array == 1:
    if d["Full_Field_IFFT"]  == 1:
        #Max_signal_for_point=ifft_on_VTA_array(name_sorted_solution,d,FR_vector_signal,Xs_signal_norm,t_vector,d["T"],inx_start_octv,arrays_shape)

        from Parallel_IFFT_on_VTA_array import get_IFFT_on_VTA_array
        Max_signal_for_point=get_IFFT_on_VTA_array(d['number_of_processors'],name_sorted_solution,d,FR_vector_signal,Xs_signal_norm,t_vector,d["T"],inx_start_octv,arrays_shape)

        get_VTA(d,VTA_full_name,Max_signal_for_point,arrays_shape,VTA_edge,VTA_resolution)
        #home_dir=os.path.expanduser("~")
        if d["Stim_side"]==0:
            subprocess.call(['touch', os.environ['PATIENTDIR']+'/success_rh.txt'])
        else:
            subprocess.call(['touch', os.environ['PATIENTDIR']+'/success_lh.txt'])

        return True

        d["IFFT_ready"] = 1               #modification of dictionary

    '''This will be enabled after the next update'''
    # if d["IFFT_ready"] == 0 and d["Full_Field_IFFT"] == 1:
    #     from Full_IFFT_field_function import get_field_in_time

    #     VTA_size = get_field_in_time(d,FR_vector_signal,Xs_signal_norm,t_vector)        # also uses data from Field_solutions_functions/ and if spectrum truncation is applied, than data from Stim_Signal/
    #
    #     d["IFFT_ready"] = 1               #modification of dictionary

    #     if d["VTA_from_divE"]==True or d["VTA_from_E"]==True:           #else it will just jump to NEURON_direct_run.py
    #         minutes=int((time.time() - start_simulation_run)/60)
    #         secnds=int(time.time() - start_simulation_run)-minutes*60
    #         total_seconds=time.time() - start_simulation_run
    #         logging.critical("---Simulation run took ",minutes," min ",secnds," s ")

    #         return None,VTA_size

#==============Truncation of the already computed spectrum====================#

    if d["Truncate_the_obtained_full_solution"] == 1 and d["IFFT_ready"] == 0:
        if d["spectrum_trunc_method"]=='High Amplitude Method' or d["spectrum_trunc_method"] == 'Cutoff Method':
            logging.critical("----- Conducting IFFT truncating already obtained full solution -----")
            from Field_IFFT_on_different_axons import convolute_signal_with_field_and_compute_ifft
            if isinstance(d["n_Ranvier"],list):             #if different populations
                last_point=0
                hf = h5py.File(d["Name_prepared_neuron_array"], 'r')
                lst_population_names=list(hf.keys())
                hf.close()
                for i in range(len(d["n_Ranvier"])):
                    logging.critical("in {} population".format(lst_population_names[i]))
                    last_point=convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm,N_array.N_models[i],N_array.pattern['num_segments'][i],FR_vector_signal,t_vector,A,os.environ['PATIENTDIR']+'/Field_solutions/sorted_solution.csv',dif_axons=True,last_point=last_point)
            else:
                convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm,N_array.N_models,N_array.pattern['num_segments'],FR_vector_signal,t_vector,A,os.environ['PATIENTDIR']+'/Field_solutions/sorted_solution.csv')
        else:
            logging.critical("Truncation of the obtained full solution is only for high. ampl and cutoff methods")

#==============================Parall. IFFT===================================#

    if d["IFFT_ready"]==0 and d["Truncate_the_obtained_full_solution"]!=1:
        logging.critical("----- Conducting signal scaling and IFFT -----")
        from Field_IFFT_on_different_axons import convolute_signal_with_field_and_compute_ifft

        if d["spectrum_trunc_method"]=='No Truncation':
            if isinstance(d["n_Ranvier"],list):             #if different populations
                last_point=0
                hf = h5py.File(d["Name_prepared_neuron_array"], 'r')
                lst_population_names=list(hf.keys())
                hf.close()
                for i in range(len(d["n_Ranvier"])):
                    logging.critical("in {} population".format(lst_population_names[i]))
                    last_point=convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm,N_array.N_models[i],N_array.pattern['num_segments'][i],FR_vector_signal,t_vector,A,name_sorted_solution,dif_axons=True,last_point=last_point)
            else:
                convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm,N_array.N_models,N_array.pattern['num_segments'],FR_vector_signal,t_vector,A,name_sorted_solution)

        if (d["spectrum_trunc_method"]=='High Amplitude Method' or d["spectrum_trunc_method"]=='Cutoff Method') and d["Truncate_the_obtained_full_solution"]==0:
            if isinstance(d["n_Ranvier"],list):             #if different populations
                last_point=0
                hf = h5py.File(d["Name_prepared_neuron_array"], 'r')
                lst_population_names=list(hf.keys())
                hf.close()
                for i in range(len(d["n_Ranvier"])):
                    logging.critical("in {} population".format(lst_population_names[i]))
                    last_point=convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm_new,N_array.N_models[i],N_array.pattern['num_segments'][i],FR_vector_signal_new,t_vector,A,name_sorted_solution,dif_axons=True,last_point=last_point)
            else:
                convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm_new,N_array.N_models,N_array.pattern['num_segments'],FR_vector_signal_new,t_vector,A,name_sorted_solution)

        if d["spectrum_trunc_method"]=='Octave Band Method':
            if isinstance(d["n_Ranvier"],list):             #if different populations
                last_point=0
                hf = h5py.File(os.environ['PATIENTDIR']+'/'+d["Name_prepared_neuron_array"], 'r')
                lst_population_names=list(hf.keys())
                hf.close()
                for i in range(len(d["n_Ranvier"])):
                    logging.critical("in {} population".format(lst_population_names[i]))
                    last_point=convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm,N_array.N_models[i],N_array.pattern['num_segments'][i],FR_vector_signal,t_vector,A,name_sorted_solution,inx_st_oct=inx_start_octv,dif_axons=True,last_point=last_point)
            else:
                convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm,N_array.N_models,N_array.pattern['num_segments'],FR_vector_signal,t_vector,A,name_sorted_solution,inx_st_oct=inx_start_octv,dif_axons=False,last_point=0)



    '''Just to compute impedance in time, on if CPE is added or current-control mode'''

    if (d["CPE_activ"]==1 or d["current_control"]==1) and cc_multicontact==False:# and d["IFFT_ready"]==0:        #modify later
        from Field_IFFT_on_different_axons import compute_Z_ifft
        logging.critical("----- Calculating impedance -----\n")
        if d["spectrum_trunc_method"]=='No Truncation' or d["Truncate_the_obtained_full_solution"]==1:
            Imp_in_time=compute_Z_ifft(d,Xs_signal_norm,t_vector,A)
        elif d["spectrum_trunc_method"]=='Octave Band Method':
            Imp_in_time=compute_Z_ifft(d,Xs_signal_norm,t_vector,A,i_start_octv=inx_start_octv)
        else:
            Imp_in_time=compute_Z_ifft(d,Xs_signal_norm_new,t_vector,A)
#===========================NEURON model simulation===========================#

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

    logging.critical("----- Estimating neuron activity -----")
    start_neuron = time.time()

    if d["Axon_Model_Type"] == 'McIntyre2002':
        with open(os.devnull, 'w') as FNULL: subprocess.call('nocmodl axnode.mod', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        with open(os.devnull, 'w') as FNULL: subprocess.call('nrnivmodl', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        from Axon_files.NEURON_direct_run import run_simulation_with_NEURON
    elif d["Axon_Model_Type"] == 'Reilly2016':
        os.chdir("Reilly2016/")
        with open(os.devnull, 'w') as FNULL: subprocess.call('nrnivmodl', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        from Axon_files.Reilly2016.NEURON_Reilly2016 import run_simulation_with_NEURON

    if isinstance(d["n_Ranvier"],list):
        Number_of_activated=0
        last_point=0
        for i in range(len(d["n_Ranvier"])):
            Number_of_activated_population=run_simulation_with_NEURON(d, last_point,i,d["diam_fib"][i],d["n_Ranvier"][i],N_array.N_models[i],d["Ampl_scale"],d["number_of_processors"],d["Name_prepared_neuron_array"])
            Number_of_activated=Number_of_activated+Number_of_activated_population

            last_point=N_array.pattern['num_segments'][i]*N_array.N_models[i]+last_point
    else:
        Number_of_activated = run_simulation_with_NEURON(d, 0,-1,d["diam_fib"],d["n_Ranvier"],N_array.N_models[0],d["Ampl_scale"],d["number_of_processors"])

    os.chdir(oss_plat_cont)

    #if isinstance(d["n_Ranvier"],list) and len(d["n_Ranvier"])>1:
        #with open(os.devnull, 'w') as FNULL: subprocess.call('python Visualization_files/Paraview_connections_activation.py', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        #if d['Show_paraview_screenshots']==1:
        #    subprocess.call('xdg-open "Images/Axon_activation.png"',shell=True)
    #else:
        #subprocess.call('python Visualization_files/Paraview_csv_activation.py', shell=True)
        #if d['Show_paraview_screenshots']==1:
         #   subprocess.call('xdg-open os.environ['PATIENTDIR']+"/Images/Activated_neurons.png"',shell=True)

    minutes=int((time.time() - start_neuron)/60)
    secnds=int(time.time() - start_neuron)-minutes*60
    logging.critical("----- NEURON calculations took {} min {} sec ----- \n".format(minutes, secnds))
    #logging.critical("----- NEURON calculations took ",minutes," min ",secnds," s -----\n")

    minutes=int((time.time() - start_simulation_run)/60)
    secnds=int(time.time() - start_simulation_run)-minutes*60
    total_seconds=time.time() - start_simulation_run
    logging.critical("----- Simulation run took {} min {} sec -----".format(minutes, secnds))
    #logging.critical("---Simulation run took ",minutes," min ",secnds," s ")

    #home_dir=os.path.expanduser("~")
    if d["Stim_side"]==0:
        subprocess.call(['touch', os.environ['PATIENTDIR']+'/success_rh.txt'])
    else:
        subprocess.call(['touch', os.environ['PATIENTDIR']+'/success_lh.txt'])

    return True


#===========================Master Study======================================#


'''Master study is useful for optimization and UQ studies'''
'''It will evaluate which truncation method and with how many frequencies we need to use for a particular setup'''
'''If adaptive mesh refinement is disabled, it will also estimate the difference in potential for 2 profile salome scripts (profile scripts should be created manually)'''

def run_master_study():

    '''Computations on initial mesh generated in SALOME'''

    Salome_profile_list=['SNEX100','SNEX100_UQ2','SNEX100_UQ3']     #create 'SNEX100_UQ2_profile','SNEX100_UQ3_profile' and 'SNEX100_UQ4_profile' from 'SNEX100_profile' (in Electrode_files/) by varying mesh parameters
    Salome_best_profile='SNEX100_UQ4'   #name of the script with a highest mesh requirements (which might be too strict for iterative studies). Will be used as a benchmark.

    master_dict={'Electrode_type':Salome_best_profile}
    total_seconds_bench,Number_of_activated_bench=run_full_model(master_dict)

    Result_diff_profiles=np.zeros((len(Salome_profile_list),2),float)

    for i in range(len(Salome_profile_list)):
        master_dict={'Electrode_type':Salome_profile_list[i], 'voxel_arr_MRI':1, 'voxel_arr_DTI':1, 'Init_neuron_model_ready':1}        #some steps were already done
        Result_diff_profiles[i,:]=run_full_model(master_dict)


    logging.critical("Number of activated neurons, the benchmark and the profiles : ", Number_of_activated_bench,Result_diff_profiles[:,1])
    logging.critical("Time in seconds, the benchmark and the profiles : ", total_seconds_bench,Result_diff_profiles[:,0])
    logging.critical("Choose the best option for 'Electrode_type'. Comment out '''Computations on initial mesh generated in SALOME''' section")

    raise SystemExit



    '''Spectrum truncation'''

    Number_of_activated_high_ampl=0               #initizialization
    Number_of_activated_octaves=0               #initizialization

    #first, run for the full spectrum
    master_dict={"spectrum_trunc_method":'No Truncation', "Truncate_the_obtained_full_solution":0}
    Run_time,Number_of_activated_benchmark=run_full_model(master_dict)

    #secondly, extract results from the full solution on freqs. with high ampl. of FFT
    trunc_param_N_freqs=25      #we will start with truncation to 25 frequencies
    while Number_of_activated_high_ampl!=Number_of_activated_benchmark:  #here in the next version we could give a percentage of deviation
        master_dict={"spectrum_trunc_method":'High Amplitude Method', "Truncate_the_obtained_full_solution":1, "trunc_param":trunc_param_N_freqs}
        Run_time_high_ampl,Number_of_activated_high_ampl=run_full_model(master_dict)
        trunc_param_N_freqs=trunc_param_N_freqs+25

    #now run the same study, but without extracting from the full solution
    master_dict={"spectrum_trunc_method":'High Amplitude Method', "Truncate_the_obtained_full_solution":0, "trunc_param":trunc_param_N_freqs-25}
    Run_time_high_ampl,Number_of_activated_high_ampl=run_full_model(master_dict)

    def func_cutoff_oct(cutoff_oct,*data_pass):     #defines the frequency after which octaves will be placed depending on the total number of freqs (the final total number can deviate by 1 frequency)
        from Inp_dict import d
        Sim_time=1.0/d["freq"]      #one pulse per simulation
        freq_max=d["freq"]*Sim_time/d["t_step"]

        return int(np.log2((freq_max-cutoff_oct)*np.sqrt(2.0)/d["freq"])+cutoff_oct/d["freq"])+1-N_freqs_oct

    from scipy.optimize import fsolve

    N_freqs_oct=10  #we will start with truncation to 10 frequencies
    while Number_of_activated_octaves!=Number_of_activated_benchmark:  #here in the next version we could give a percentage of deviation
        data_pass=(N_freqs_oct)
        cutoff_octaves=fsolve(func_cutoff_oct,10,args=data_pass)

        master_dict={"spectrum_trunc_method":'Octave Band Method', "Truncate_the_obtained_full_solution":0, "trunc_param":cutoff_octaves}
        Run_time_oct,Number_of_activated_octaves=run_full_model(master_dict)
        N_freqs_oct=N_freqs_oct+5
        if Run_time_oct>Run_time_high_ampl:
            logging.critical("Octave method is already slower than high amplitude method., no need to continue.")
            break


    if Run_time_oct<Run_time_high_ampl and Run_time_oct<Run_time:
        logging.critical("Octave truncation is the fastest")
    elif Run_time_oct>Run_time_high_ampl and Run_time_high_ampl<Run_time:
        logging.critical("High amplitude truncation is the fastest")
    else:
        logging.critical("Truncation did not accelerate the computations")

    logging.critical("Time for Full_run,High_ampl,Octaves: ", Run_time,Run_time_high_ampl,Run_time_oct)
    logging.critical("Activation for Full_run,High_ampl,Octaves: ", Number_of_activated_benchmark,Number_of_activated_high_ampl,Number_of_activated_octaves)


#master_dict = {}
#run_full_model(master_dict)

from datetime import datetime

# datetime object containing current date and time
now = datetime.now()
date_and_time = now.strftime("%d-%m-%Y___%H-%M-%S")  # EUROPEAN format

import logging
logging.basicConfig(filename='/opt/Patient' + '/complete_log_' + date_and_time + '.log', format='[%(asctime)s]:%(message)s', level=logging.ERROR)


logf = open("/opt/Patient/last_error.log", "w")

try:
    master_dict = {}
    run_full_model(master_dict)
except Exception:
    logging.exception("Fatal error in main loop")