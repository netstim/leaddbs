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
print("\n\n\nOSS-DBS by K.Butenko --- version 0.3")
print("Manuscript is submitted")
print("____________________________________\n")

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



def run_full_model(master_dict):
    
    start_simulation_run=time.time() 
    
    #===================Load and update input dictionary======================#
    #should be loaded this way for iterative studies (some simulation state variables change during a run)
    import GUI_inp_dict
    importlib.reload(GUI_inp_dict)
    from GUI_inp_dict import d
    
    from Dict_corrector import rearrange_Inp_dict
    d=rearrange_Inp_dict(d)             #misc. transformation of parameters to the platform's format
    d.update(master_dict)               #modifies the user provided input dictionary (e.g. for UQ study), check run_master_study() function . Warning: this does not work update the encap. layer properties and the solver during adaptive mesh refiment, because these data are realoaded from the original dictionary
   
    #=========Check the simulation setup and state, load the corresponding data=========#
    if d["current_control"]==1 and d["CPE_activ"]==1:
        d["CPE_activ"]=0
        print("Disabling CPE for current-controlled simulation")

    Phi_vector_active_non_zero=[x for x in d["Phi_vector"] if (x is not None) and (x!=0.0)]
    cc_multicontact=False
    if d["current_control"]==1 and len(Phi_vector_active_non_zero)>1:       #multicontact current-controlled case
        cc_multicontact=True

    from Sim_state import check_state
    check_state(d)      #will switch the state of simulation depending on the corresp. data in the input_dictionary, will also manage project folders 
       
    if d["voxel_arr_MRI"]==0 and d["voxel_arr_DTI"]==1:     
        print("MRI data is new, the DTI data will be reprocessed")
        d["voxel_arr_DTI"]==0
    
    #loading of meta data depending on the simulatation setup and state
    if d["Init_neuron_model_ready"]==1:     
        [ranvier_nodes, para1_nodes, para2_nodes, inter_nodes, ranvier_length, para1_length, para2_length, inter_length, deltax, diam_fib,n_Ranvier,ROI_radius,N_segm]=np.genfromtxt('Neuron_model_arrays/Neuron_model_misc.csv', delimiter=' ')
        param_axon=[ranvier_nodes, para1_nodes, para2_nodes, inter_nodes, ranvier_length, para1_length, para2_length, inter_length, deltax, diam_fib]            

        if d["Neuron_model_array_prepared"]==0:    
            with open('Neuron_model_arrays/Neuron_param_class.file', "rb") as f:
                Neuron_param = pickle.load(f)
        if d["Neuron_model_array_prepared"]==1:
            Neuron_param=0          #not needed for a pre-defined neuron array
            if d["Name_prepared_neuron_array"][-3:]=='.h5':     
                n_segments_fib_diam_array=np.load('Neuron_model_arrays/Neuron_populations_misc.npy')
                N_segm=n_segments_fib_diam_array[:,0]
                N_segm=N_segm.astype(int)
    if d["Init_mesh_ready"]==1:
        with open('Meshes/Mesh_ind.file', "rb") as f:
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
    if d["Init_mesh_ready"]==0:
                        
        if d["Brain_shape_name"]==0:   #Creates a brain approximation (ellisploid)
            from CAD_Salome import build_brain_approx
            x_length,y_length,z_length=build_brain_approx(d,MRI_param)      #also creates 'brain_subsitute.brep'
            Brain_shape_name='Brain_substitute.brep'
            
        if d["Brain_shape_name"]!=0:
            Brain_shape_name=d["Brain_shape_name"]
    
        #==================Initial neuron array generation====================#
        if d["Init_neuron_model_ready"]==0 and d["Neuron_model_array_prepared"]==0:
            print("----- Creating initial neuron array -----")
            
            from Neuron_models_arangement_new import build_neuron_models
            Neuron_param,param_axon,ROI_radius,N_segm=build_neuron_models(d,MRI_param)     #builds a pattern model, if not provided, then builds a neuron array and stores in 'Neuron_model_arrays/All_neuron_models.csv'. Also, creates corresp. meta data.
            
        if d["Neuron_model_array_prepared"]==1 and d["Init_neuron_model_ready"]==0:
            print("----- Creating initial neuron array from a provided neuron array -----")
            from Neuron_models_arangement_new import create_meta_data_for_predefined_models, cut_models_by_domain
            
            cut_models_by_domain(d,Brain_shape_name,d["Name_prepared_neuron_array"])      #to adjust the prepared neuron array to the computational domain (only for brain substitutes!)
            ROI_radius,param_axon,N_segm=create_meta_data_for_predefined_models(d,MRI_param)  #shifts coordinates of the provided neuron array to the positive octant coord. and stores in 'Neuron_model_arrays/All_neuron_models.csv'. Also creates corresp. meta data.
            Neuron_param=0  #This class is irrelevant if we have a full neuron array predefined
            
        if d["Neuron_model_array_prepared"]==1 and d["Init_neuron_model_ready"]==1:
            Neuron_param=0  #This class is irrelevant if we have a full neuron array predefined
            print("--- Initial neuron array was loaded\n")

        #if brain substitute is used, it will be enlarged to encompass previosly defined neuron array (if necessary)
        if Brain_shape_name=='Brain_substitute.brep':
            needs_a_rebuid=0
            if ROI_radius > (x_length/2):
                d['Approximating_Dimensions'][0]=ROI_radius*2+0.1
                print("increasing length along x to encompass the neuron array\n")
                needs_a_rebuid=1
            if ROI_radius > (y_length/2):
                d['Approximating_Dimensions'][1]=ROI_radius*2+0.1
                print("increasing length along y to encompass the neuron array\n")
                needs_a_rebuid=1
            if ROI_radius > (z_length/2):
                d['Approximating_Dimensions'][2]=ROI_radius*2+0.1
                print("increasing length along z to encompass the neuron array\n")
                needs_a_rebuid=1
            if needs_a_rebuid==1:
                x_length,y_length,z_length=build_brain_approx(d,MRI_param)      #also creates 'brain_subsitute.brep'
                print("\n")
            if ROI_radius > min((x_length/2),(y_length/2),(z_length/2)):               
                print("ROI_radius: ",ROI_radius)
                print("ROI is still bigger than the computational domain.")
                raise SystemExit
 
        #===================Final geometry generation=========================#
        from CAD_Salome import build_final_geometry
        Domains=build_final_geometry(d,MRI_param,Brain_shape_name,ROI_radius,cc_multicontact)       #creates and discretizes the geometry with the implanted electrode, encapsulation layer and ROI, converts to the approp. format. The files are stored in Meshes/
    
    #===============Adjusting neuron array====================================#
    
    if d["Adjusted_neuron_model_ready"]==0:
        from Neuron_models_arangement_new import adjust_neuron_models 
        N_models=adjust_neuron_models(d,MRI_param,Domains,Neuron_param,param_axon)       #subtracts neurons from the previously define All_neuron_models.csv, if the are non-physical (inside encap. layer or CSF, outside of the domain, intersect with the electrode geometry.)      
    else:
        if d["Neuron_model_array_prepared"]==1:
            if d["Name_prepared_neuron_array"][-3:]=='.h5':     #if imported with h5
                N_models = np.genfromtxt('Neuron_model_arrays/Adjusted_neuron_array_info.csv', delimiter=' ')
                N_models=N_models.astype(int)
            else:
                N_models,points_csf,points_encap_and_float_contacts,points_outside=np.genfromtxt('Neuron_model_arrays/Adjusted_neuron_array_info.csv', delimiter=' ')
                N_models=int(N_models)                
        else:
            N_models,points_csf,points_encap_and_float_contacts,points_outside=np.genfromtxt('Neuron_model_arrays/Adjusted_neuron_array_info.csv', delimiter=' ')
            N_models=int(N_models)
        print("--- Neuron array meta data were loaded\n")
        
    number_of_points=int(np.sum(N_segm*N_models))

    #subprocess.call('python Paraview_InitMesh_and_Neurons.py', shell=True)
    if d['Show_paraview_screenshots']==1:
        subprocess.call('xdg-open "Images/InitMesh_and_Neurons.png"',shell=True)


#=============================Signal creation=================================#
    #in case there was an update of the signal
    Phi_vector=[x for x in d["Phi_vector"] if x is not None] #now we don't need None values (floating potentials), Contacts contain only the active ones
    Domains.fi=Phi_vector
    with open('Meshes/Mesh_ind.file', "wb") as f:
        pickle.dump(Domains, f, pickle.HIGHEST_PROTOCOL)

    if d["current_control"]==1 and cc_multicontact==False:     #if multicontact, then it will be scaled in parallel comp.
        A=max(Domains.fi[:], key=abs)
        Phi_max=1      #not needed
    else:
        A=1.0 #for vc we will apply Phi_vector as boundary conditions with appropriate amplitude
        Phi_max=max(Domains.fi[:], key=abs)

    if d["signal_generation_ready"]==0:
        print("----- Generating DBS signal -----")
        
        from Signal_generator import generate_signal
        [t_vector,signal_out,Xs_signal_norm,FR_vector_signal]=generate_signal(d,A,Phi_max,cc_multicontact)
                        
        Xs_storage=np.zeros((np.real(Xs_signal_norm).shape[0],2),float)
        
        Xs_storage[:,0]=np.real(Xs_signal_norm)
        Xs_storage[:,1]=np.imag(Xs_signal_norm)
        np.savetxt('Stim_Signal/Xs_storage_full.csv', Xs_storage, delimiter=" ")
        np.savetxt('Stim_Signal/t_vector.csv', t_vector, delimiter=" ")
        np.savetxt('Stim_Signal/FR_vector_signal.csv', FR_vector_signal, delimiter=" ")
    else:
        print("--- DBS signal is taken from the previous simulation\n")
        Xs_recovered = np.genfromtxt('Stim_Signal/Xs_storage_full.csv', delimiter=' ')
        Xs_signal_norm=np.vectorize(complex)(Xs_recovered[:,0],Xs_recovered[:,1])
        FR_vector_signal = np.genfromtxt('Stim_Signal/FR_vector_signal.csv', delimiter=' ')
        t_vector = np.genfromtxt('Stim_Signal/t_vector.csv', delimiter=' ')       


#=============================CSF_refinement==================================#

    '''if we want to skip CSF and adaptive mesh refinement'''
    if d["Skip_mesh_refinement"]==1:
        from CSF_refinement_new import Dummy_CSF
        Dummy_CSF()         #will resave 'Meshes/Mesh_unref.xml' as 'Results_adaptive/mesh_adapt.xml.gz' and the same for subdomains and boundaries
        d["CSF_mesh_ready"]=1
        d["Adapted_mesh_ready"]=1
        print("CSF and adaptive mesh refinement was skipped\n")
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
        #vtu_file = "CSF_ref/Field_r_1.0000000.vtu"
        #f_name="Visualization_files/Paraview_field_CSF_ref.py"
        #f_nam = get_Para_Array_name(vtu_file)
        #insert_f_name(f_nam,f_name)
        
    elif d["Skip_mesh_refinement"]==0:      #should be changed
        Scaling_CSF=np.genfromtxt('CSF_ref/Scaling_CSF.csv', delimiter=' ')

        
    import os
    #with open(os.devnull, 'w') as FNULL: subprocess.call('python Visualization_files/Paraview_field_CSF_ref.py', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    with open(os.devnull, 'w') as FNULL: subprocess.call('python Visualization_files/Paraview_CSFref.py', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    
    if d['Show_paraview_screenshots']==1 and d["Skip_mesh_refinement"]==0:
        subprocess.call('xdg-open "Images/CSF_ref.png"',shell=True)
        #subprocess.call('xdg-open "Images/Field_Adapted_CSF.png"',shell=True)
    
    
#====================Mesh adaptive algorithm==================================#

    if d["Adapted_mesh_ready"]==0:
        from Mesh_adaption_hybrid import mesh_adapter        
        Ampl_on_vert=mesh_adapter(MRI_param,DTI_param,Scaling_CSF,Domains,d,anisotrop,cc_multicontact,ref_freqs)     #also saves adapted mesh
        np.savetxt('Results_adaptive/Ampl_on_vert.csv', Ampl_on_vert, delimiter=" ")    #just to check
     
    with open(os.devnull, 'w') as FNULL: subprocess.call('python Visualization_files/Paraview_adapted.py', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    if d['Show_paraview_screenshots']==1 and d["Skip_mesh_refinement"]==0:  
        subprocess.call('xdg-open "Images/Adapted_mesh.png"',shell=True)
                      
#===================Truncate the frequency spectrum===========================#
    if d["spectrum_trunc_method"]!='No Truncation':

        from Truncation_of_spectrum import get_freqs_for_calc
        FR_vector_signal_new,add_trunc_data=get_freqs_for_calc(d,FR_vector_signal,Xs_signal_norm,t_vector)
        if d["spectrum_trunc_method"]=='Octave Band Method':
            inx_start_octv=add_trunc_data
        else:
            Xs_signal_norm_new=add_trunc_data  

#==========Calculate freq in parallel and rearrange field array===============#
    if d["Parallel_comp_ready"]==0:
        
        if ["Parallel_comp_interrupted"]==1:
            import os
            if not (os.path.isfile('Field_solutions/Phi_real_scaled_'+str(d["freq"])+'Hz.pvd') or os.path.isfile('Field_solutions/Phi_real_unscaled_'+str(d["freq"])+'Hz.pvd')):     #to make sure that there were interrupted computations 
                print("There were no previous computations, 'Parallel_comp_interrupted' is put to 0")
                ["Parallel_comp_interrupted"]==0
        
        from Parallel_field_calc import calculate_in_parallel
        '''calculate_in_parallel will save a sorted_solution array if IFFT is pointwise or will save the whole field for each frequency in Field_solutions_functions/ if full_IFFT is requested'''
        
        if d["spectrum_trunc_method"]=='No Truncation':   
            print("----- Calculating electric field in the frequency spectrum -----")
            calculate_in_parallel(d,FR_vector_signal,Domains,MRI_param,DTI_param,anisotrop,number_of_points,cc_multicontact)
            
        if d["spectrum_trunc_method"]=='High Amplitude Method' or d["spectrum_trunc_method"]=='Cutoff Method' or d["spectrum_trunc_method"]=='Octave Band Method':
            print("----- Calculating electric field in the truncated frequency spectrum -----")
            calculate_in_parallel(d,FR_vector_signal_new,Domains,MRI_param,DTI_param,anisotrop,number_of_points,cc_multicontact)
    else:
        print("--- Results of calculations in the frequency spectrum were loaded\n")

    if d["spectrum_trunc_method"]=='No Truncation':
        name_sorted_solution='Field_solutions/sorted_solution.csv'  
    else:
        name_sorted_solution='Field_solutions/sorted_solution_'+str(d["spectrum_trunc_method"])+'_'+str(d["trunc_param"])+'.csv'
  

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
        print("--- Results of IFFT (FFEM) were loaded\n")
    
    if d["IFFT_ready"] == 0 and d["Full_Field_IFFT"] == 1:
        from Full_IFFT_field_function import get_field_in_time
        
        VTA_size = get_field_in_time(d,FR_vector_signal,Xs_signal_norm,t_vector)        # also uses data from Field_solutions_functions/ and if spectrum truncation is applied, than data from Stim_Signal/
                                                                    # if VTA_from_NEURON is enabled, will save pointwise solutions in Points_in_time/ 
        d["IFFT_ready"] = 1               #modification of dictionary
        
        if d["VTA_from_divE"]==True or d["VTA_from_E"]==True:           #else it will just jump to NEURON_direct_run.py
            minutes=int((time.time() - start_simulation_run)/60)
            secnds=int(time.time() - start_simulation_run)-minutes*60
            total_seconds=time.time() - start_simulation_run
            print("---Simulation run took ",minutes," min ",secnds," s ") 
            
            return None,VTA_size

#==============Truncation of the already computed spectrum====================#
        
    if d["Truncate_the_obtained_full_solution"] == 1 and d["IFFT_ready"] == 0:
        if d["spectrum_trunc_method"]=='High Amplitude Method' or d["spectrum_trunc_method"] == 'Cutoff Method':
            print("----- Conducting IFFT truncating already obtained full solution -----")
            from Field_IFFT_on_different_axons import convolute_signal_with_field_and_compute_ifft
            if isinstance(d["n_Ranvier"],list):             #if different populations
                last_point=0
                hf = h5py.File(d["Name_prepared_neuron_array"], 'r')
                lst_population_names=list(hf.keys())
                hf.close()
                for i in range(len(d["n_Ranvier"])):
                    print("in ",lst_population_names[i]," population")
                    last_point=convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm,N_models[i],N_segm[i],FR_vector_signal,t_vector,A,'Field_solutions/sorted_solution.csv',dif_axons=True,last_point=last_point)
            else:
                convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm,N_models,N_segm,FR_vector_signal,t_vector,A,'Field_solutions/sorted_solution.csv')       
        else:        
            print("Truncation of the obtained full solution is only for high. ampl and cutoff methods")
    
#==============================Parall. IFFT===================================#      
    
    if d["IFFT_ready"]==0 and d["Truncate_the_obtained_full_solution"]!=1:
        print("----- Conducting signal scaling and IFFT -----")
        from Field_IFFT_on_different_axons import convolute_signal_with_field_and_compute_ifft
 
        if d["spectrum_trunc_method"]=='No Truncation':
            if isinstance(d["n_Ranvier"],list):             #if different populations
                last_point=0
                hf = h5py.File(d["Name_prepared_neuron_array"], 'r')
                lst_population_names=list(hf.keys())
                hf.close()
                for i in range(len(d["n_Ranvier"])):
                    print("in ",lst_population_names[i]," population")
                    last_point=convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm,N_models[i],N_segm[i],FR_vector_signal,t_vector,A,name_sorted_solution,dif_axons=True,last_point=last_point)
            else:
                convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm,N_models,N_segm,FR_vector_signal,t_vector,A,name_sorted_solution)  
                       
        if (d["spectrum_trunc_method"]=='High Amplitude Method' or d["spectrum_trunc_method"]=='Cutoff Method') and d["Truncate_the_obtained_full_solution"]==0:
            if isinstance(d["n_Ranvier"],list):             #if different populations
                last_point=0
                hf = h5py.File(d["Name_prepared_neuron_array"], 'r')
                lst_population_names=list(hf.keys())
                hf.close()
                for i in range(len(d["n_Ranvier"])):
                    print("in ",lst_population_names[i]," population")
                    last_point=convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm_new,N_models[i],N_segm[i],FR_vector_signal_new,t_vector,A,name_sorted_solution,dif_axons=True,last_point=last_point)
            else:
                convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm_new,N_models,N_segm,FR_vector_signal_new,t_vector,A,name_sorted_solution)  

        if d["spectrum_trunc_method"]=='Octave Band Method':
            if isinstance(d["n_Ranvier"],list):             #if different populations
                last_point=0
                hf = h5py.File(d["Name_prepared_neuron_array"], 'r')
                lst_population_names=list(hf.keys())
                hf.close()
                for i in range(len(d["n_Ranvier"])):
                    print("in ",lst_population_names[i]," population")
                    last_point=convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm,N_models[i],N_segm[i],FR_vector_signal,t_vector,A,name_sorted_solution,inx_st_oct=inx_start_octv,dif_axons=True,last_point=last_point)
            else:
                convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm,N_models,N_segm,FR_vector_signal,t_vector,A,name_sorted_solution,inx_st_oct=inx_start_octv,dif_axons=False,last_point=0)  


        
    '''Just to compute impedance in time, on if CPE is added or current-control mode'''

    if (d["CPE_activ"]==1 or d["current_control"]==1) and cc_multicontact==False and d["spectrum_trunc_method"]!='Octave Band Method' and d["IFFT_ready"]==0:        #modify later
        from Field_IFFT_on_different_axons import compute_Z_ifft
        print("-----Calculating impedance over time-----\n")
        if d["spectrum_trunc_method"]=='No Truncation' or d["Truncate_the_obtained_full_solution"]==1:
            Imp_in_time=compute_Z_ifft(d,Xs_signal_norm,t_vector,A)
        else:
            Imp_in_time=compute_Z_ifft(d,Xs_signal_norm_new,t_vector,A)
        
    
#===========================NEURON model simulation===========================#

    print("----- Estimating neuron activity -----")
    start_neuron=time.time()

    if d["Axon_Model_Type"] == 'McIntyre2002': 
        os.chdir("Axon_files/")
        with open(os.devnull, 'w') as FNULL: subprocess.call('nocmodl axnode.mod', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        with open(os.devnull, 'w') as FNULL: subprocess.call('nrnivmodl', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        from Axon_files.NEURON_direct_run import run_simulation_with_NEURON
    elif d["Axon_Model_Type"] == 'Reilly2016':
        os.chdir("Axon_files/Reilly2016/")
        with open(os.devnull, 'w') as FNULL: subprocess.call('nrnivmodl', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        from Axon_files.Reilly2016.NEURON_Reilly2016 import run_simulation_with_NEURON
        
    if isinstance(d["n_Ranvier"],list) and len(d["n_Ranvier"])>1:        
        Number_of_activated=0
        last_point=0
        for i in range(len(d["n_Ranvier"])):
            Number_of_activated_population=run_simulation_with_NEURON(last_point,i,d["diam_fib"][i],1000*d["t_step"],1000.0/d["freq"],d["n_Ranvier"][i],N_models[i],d["v_init"],t_vector.shape[0],d["Ampl_scale"],d["number_of_processors"])
            Number_of_activated=Number_of_activated+Number_of_activated_population
            os.chdir("Axon_files/")
            if d["Axon_Model_Type"] == 'Reilly2016':
                os.chdir("Reilly2016/")
            last_point=N_segm[i]*N_models[i]+last_point
        os.chdir("..")     
        if d["Axon_Model_Type"] == 'Reilly2016':
            os.chdir("..")
    else:
        Number_of_activated=run_simulation_with_NEURON(0,-1,d["diam_fib"],1000*d["t_step"],1000.0/d["freq"],d["n_Ranvier"],N_models,d["v_init"],t_vector.shape[0],d["Ampl_scale"],d["number_of_processors"])
    
    if isinstance(d["n_Ranvier"],list) and len(d["n_Ranvier"])>1:
        with open(os.devnull, 'w') as FNULL: subprocess.call('python Visualization_files/Paraview_connections_activation.py', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        if d['Show_paraview_screenshots']==1:
            subprocess.call('xdg-open "Images/Axon_activation.png"',shell=True)
    else:
        subprocess.call('python Visualization_files/Paraview_csv_activation.py', shell=True)
        if d['Show_paraview_screenshots']==1:
            subprocess.call('xdg-open "Images/Activated_neurons.png"',shell=True)

    minutes=int((time.time() - start_neuron)/60)
    secnds=int(time.time() - start_neuron)-minutes*60
    print("----- NEURON calculations took ",minutes," min ",secnds," s -----\n") 

    minutes=int((time.time() - start_simulation_run)/60)
    secnds=int(time.time() - start_simulation_run)-minutes*60
    total_seconds=time.time() - start_simulation_run
    print("---Simulation run took ",minutes," min ",secnds," s ")  
    
    return total_seconds,Number_of_activated


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
        
        
    print("Number of activated neurons, the benchmark and the profiles : ", Number_of_activated_bench,Result_diff_profiles[:,1])
    print("Time in seconds, the benchmark and the profiles : ", total_seconds_bench,Result_diff_profiles[:,0])
    print("Choose the best option for 'Electrode_type'. Comment out '''Computations on initial mesh generated in SALOME''' section")
    
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
            print("Octave method is already slower than high amplitude method., no need to continue.")
            break
        
        
    if Run_time_oct<Run_time_high_ampl and Run_time_oct<Run_time:
        print("Octave truncation is the fastest")       
    elif Run_time_oct>Run_time_high_ampl and Run_time_high_ampl<Run_time:
        print("High amplitude truncation is the fastest")       
    else:
        print("Truncation did not accelerate the computations")
                
    print("Time for Full_run,High_ampl,Octaves: ", Run_time,Run_time_high_ampl,Run_time_oct)
    print("Activation for Full_run,High_ampl,Octaves: ", Number_of_activated_benchmark,Number_of_activated_high_ampl,Number_of_activated_octaves)
        

#run_master_study()    #in case we want to find optimal spectrum truncation method and initial mesh settings.

master_dict={}          #you can implement UQ or optimization by adding a function that will create master_dict with entries to be optimized. Name of entries should be the same as in GUI_inp_dict.py
Run_time_high_ampl,Number_of_activated_high_ampl=run_full_model(master_dict)
