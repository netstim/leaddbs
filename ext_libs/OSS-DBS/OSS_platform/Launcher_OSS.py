'''
#@author: Konstantin Butenko (konstantin.butenko@uni-rostock.de), routines for Salome prepared in cooperation with P.Trieu. Script "tissue dielectrics.py" and axon files are adapted from FanPy (https://arxiv.org/abs/1705.10478)

The platform is in alpha test and should be used with caution and a proper amount of healthy suspicion.
'''

        
'''This is the main script that governs the simulation flow

It contains functions: 

run_master_study(), which is useful for optimization and UQ studies
It will evaluate which truncation method and with how many frequencies we need to use for a particular setup
If the adataptive mesh alg. are not usable in the study, it will also estimate the difference in potential for 2 profile salome scripts (profile scripts should be created manually)

run_full_model(master_dict) governs the simulation flow taking master_dict as the input file

In the section Uncertainty Quantification you will find an example how to conduct the analysis using the platform

Optimization algorithms are to be implemented'''
print("")
print("OSS platform for DBS by K.Butenko --- version 0.2")
print("Manuscript is in preparation")
print("____________________________________")

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



#def run_full_model(encap_thickness_UQ,encap_cond_UQ,y_t_UQ):     #example of the header for UQ a analysis
def run_full_model(master_dict):
    
 
    start_simulation_run=time.time() 
    
    #should be loaded in this way for UQ/optimization studies
    import Inp_dict 
    importlib.reload(Inp_dict)
    from Inp_dict import d as d_init    
    from GUI_inp_dict import d
    
    '''For the master study'''
   
    
    from Dict_corrector import rearrange_Inp_dict
    d=rearrange_Inp_dict(d)

    d.update(master_dict)               #WARNING: modifies the user provided input, check run_master_study() function on the bottom
    
       
    
    '''=============Should be modified by user depending on problem========='''
    
    #Example for a UQ study
    
    #'''Uncomment, if necessary'''
    #encap_thickness_UQ,encap_cond_UQ,y_t_UQ=(0.4,0.05,0.0)              #if testing with no uncertainty is required
    
    #'''Below assign the optimized parameters to the input dictionary parameters, e.g.'''
    #d["encap_thickness"]=encap_thickness_UQ    
    #d["Yt"]=d["Yt"]+y_t_UQ   #here y_t_UQ is a parameter of missplacement
    
    #UQ_sim_parameter=np.vstack((encap_thickness_UQ,encap_cond_UQ,y_t_UQ)).T       #to store data on the study
    #print("UQ param (thickness, cond, y_t): ",encap_thickness_UQ,encap_cond_UQ,y_t_UQ)
    
    '''===============Store data on the study==============================='''
#    f=open('UQ_stuff/UQ_sim_parameters.csv','ab')
#    np.savetxt(f, UQ_sim_parameter, delimiter=" ")
#    f.close()
#
#    srvs_array = np.genfromtxt('UQ_stuff/UQ_sim_parameters.csv', delimiter=' ')
#    if srvs_array.ndim==1:
#        Study_step=0
#    else:
#        Study_step=srvs_array.shape[0]-1
#    print("Study_step: ",Study_step)
    
    
    '''===========Check the simulation state, load the corresponding data============'''
    print(" ")
    if d["current_control"]==1 and d["CPE_activ"]==1:
        d["CPE_activ"]=0
        print("Disabling CPE for current-controlled simulation")

    Phi_vector_active_non_zero=[x for x in d["Phi_vector"] if (x is not None) and (x!=0.0)]
    cc_multicontact=False
    if d["current_control"]==1 and len(Phi_vector_active_non_zero)>1:       #multicontact current-controlled case
        cc_multicontact=True

    from Sim_state import check_state
    check_state(d)      #will switch the state of simulation depending on the corresp. data in the input_dictionary 
    
    if d["voxel_arr_MRI"]==0 and d["voxel_arr_DTI"]==1:
        print("MRI data is new, the DTI data will be reprocessed")
        d["voxel_arr_DTI"]==0
    
    if d["Init_neuron_model_ready"]==1:
        [ranvier_nodes, para1_nodes, para2_nodes, inter_nodes, ranvier_length, para1_length, para2_length, inter_length, deltax, diam_fib,n_Ranvier,ROI_radius,N_segm]=np.genfromtxt('Neuron_model_arrays/Neuron_model_misc.csv', delimiter=' ')
        param_axon=[ranvier_nodes, para1_nodes, para2_nodes, inter_nodes, ranvier_length, para1_length, para2_length, inter_length, deltax, diam_fib]            


        if d["Neuron_model_array_prepared"]==0:    #otherwise we do not need Neuron_param class
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
    MRI_param=obtain_MRI_class(d)       #also creates 'MRI_DTI_files/voxel_array_MRI.csv' and meta data
       
    anisotrop=0             #initizialization, 0 by default for isotropic
    DTI_param=0             #initizialization, 0 by default for isotropic
    
    if d["DTI_data_name"]!=0:	#0 if we do not provide DTI data. DTI nd MRI data must be in the same coordinate space and DTI box should not exceed MRI box!!!
        anisotrop=1         #the study will use tensor conductivity
        from MRI_DTI_prep_new import obtain_DTI_class
        DTI_param=obtain_DTI_class(d,MRI_param) #also creates 'MRI_DTI_files/voxel_array_DTI.csv' and meta data
        
    
    #========================Geometry generation==============================#
    if d["Init_mesh_ready"]==0:
                        
        if d["Brain_shape_name"]==0:   #Creates a brain approximation, ellipsoid or box (add box!)
            from CAD_Salome import build_brain_approx
            x_length,y_length,z_length=build_brain_approx(d,MRI_param)      #also creates 'brain_subsitute.brep'
            Brain_shape_name='Brain_substitute.brep'
            
        if d["Brain_shape_name"]!=0:
            Brain_shape_name=d["Brain_shape_name"]
    
        #=====================Initial Neuron models===========================#
        if d["Init_neuron_model_ready"]==0 and d["Neuron_model_array_prepared"]==0:
            print("----- Creating initial neuron array -----")
            
            from Neuron_models_arangement_new import build_neuron_models
            Neuron_param,param_axon,ROI_radius,N_segm=build_neuron_models(d,MRI_param)     #builds a pattern model, if not provided, then builds a neuron array and stores in 'Neuron_model_arrays/All_neuron_models.csv'. Also, creates corresp. meta data.
            
        if d["Neuron_model_array_prepared"]==1 and d["Init_neuron_model_ready"]==0:
            print("----- Creating initial neuron array from a provided neuron array-----")
            from Neuron_models_arangement_new import create_meta_data_for_predefined_models, cut_models_by_domain
            
            cut_models_by_domain(d,Brain_shape_name,d["Name_prepared_neuron_array"])      #to adjust the prepared neuron array to the comp. domain (only for brain substitutes!)
            ROI_radius,param_axon,N_segm=create_meta_data_for_predefined_models(d,MRI_param)  #shifts coordinates op the provided neuron array to the positive octant coord. and stores in 'Neuron_model_arrays/All_neuron_models.csv'. Also creates corresp. meta data.
            '''This class is irrelevant if we have a full neuron array predefined'''
            Neuron_param=0
            
        if d["Neuron_model_array_prepared"]==1 and d["Init_neuron_model_ready"]==1:
            Neuron_param=0
            print("--- Initial neuron array was loaded")
            print(" ")

        if Brain_shape_name=='Brain_substitute.brep':
            needs_a_rebuid=0
            if ROI_radius > (x_length/2):
                d["x_length"]=ROI_radius*2+0.1
                print("changing d['x_length'] to encompass the neuron array")
                needs_a_rebuid=1
            if ROI_radius > (y_length/2):
                d["y_length"]=ROI_radius*2+0.1
                print("changing d['y_length'] to encompass the neuron array")
                needs_a_rebuid=1
            if ROI_radius > (z_length/2):
                d["z_length"]=ROI_radius*2+0.1
                print("changing d[z_length'] to encompass the neuron array")
                needs_a_rebuid=1
            if needs_a_rebuid==1:
                print(" ")
                x_length,y_length,z_length=build_brain_approx(d,MRI_param)      #also creates 'brain_subsitute.brep'
                print(" ")
            if ROI_radius > min((x_length/2),(y_length/2),(z_length/2)):               
                print("ROI_radius: ",ROI_radius)
                print("ROI is still bigger than the computational domain.")
                raise SystemExit
 
        #========================Final geometry generation==================================#
        from CAD_Salome import build_final_geometry
        Domains=build_final_geometry(d,MRI_param,Brain_shape_name,ROI_radius,cc_multicontact)       #creates and discretizes the geometry with the implanted electrode, encapsulation layer and ROI, converts to the approp. format. The files are stored in Meshes/
    
    #===============Updated neuron model==========================================#
    
    if d["Adjusted_neuron_model_ready"]==0:
        from Neuron_models_arangement_new import adjust_neuron_models 
        N_models=adjust_neuron_models(d,MRI_param,Domains,Neuron_param,param_axon)       #subtracts neurons from the previously define All_neuron_models.csv, if the are non-physical (inside encap. layer or CSF, outside of the domain, intersect with the electrode geometry.) 
        #The final neuron array is stored in 'Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv'
        import os
        direct = os.getcwd()
        with open(os.devnull, 'w') as FNULL: subprocess.call('salome -t python '+ 'Visualization_files/Paraview_InitMesh_and_Neurons.py' +' --ns-port-log='+direct+'/salomePort.txt', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        #################### Using this script to close Salome python script
        port_file = open(direct+'/salomePort.txt','r')
        killPort = int(port_file.readline())
        port_file.close()
        #Kill the session with the specified port:
        subprocess.call('/opt/SALOME-8.3.0-UB16.04/BINARIES-UB16.04/KERNEL/bin/salome/killSalomeWithPort.py %s' % killPort,shell=True)   
     
     

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
        print("--- Neuron array meta data were loaded")
        print(" ")
        

   
    #subprocess.call('python Paraview_InitMesh_and_Neurons.py', shell=True)
    if d['Show_paraview_screenshots']==1:
        subprocess.call('xdg-open "Images/InitMesh_and_Neurons.png"',shell=True)
     
    number_of_points=np.sum(N_segm*N_models)
    number_of_points=int(number_of_points)

    #dielectric properties of encap. are loaded directly during mesh adaption. To be changed in the later version.  
    #Cond_perm_encap=np.array([d["encap_tissue_type"],d["encap_scaling_cond"],d["encap_scaling_perm"]])        
    #np.savetxt('Cond_perm_encap.csv', Cond_perm_encap, delimiter=" ")   
    
    #=========================CSF_refinement=====================================#

    '''if we want to skip CSF and adaptive mesh refinement'''
    if d["Skip_mesh_refinement"]==1:
        from CSF_refinement_new import Dummy_CSF
        Dummy_CSF()         #will resave 'Meshes/Mesh_unref.xml' as 'Results_adaptive/mesh_adapt.xml.gz' and the same for subdomains and boundaries
        d["CSF_mesh_ready"]=1
        d["Adapted_mesh_ready"]=1
        print("CSF and adaptive mesh refinement was skipped")
        print(" ")
   

    if d["CSF_mesh_ready"]==0:   
        from CSF_refinement_new import launch_CSF_refinement
      
        Scaling_CSF=launch_CSF_refinement(d,MRI_param,DTI_param,Domains,anisotrop,cc_multicontact)  
        
        from Visualization_files.paraview_find_arrayname import get_Para_Array_name
        from Parameter_insertion import insert_f_name
        vtu_file = "CSF_ref/Field_r_1.0000000.vtu"
        f_name="Visualization_files/Paraview_field_CSF_ref.py"
        f_nam = get_Para_Array_name(vtu_file)
        insert_f_name(f_nam,f_name)
        
    elif d["Skip_mesh_refinement"]==0:
        Scaling_CSF=np.genfromtxt('CSF_ref/Scaling_CSF.csv', delimiter=' ')

        
    import os
    with open(os.devnull, 'w') as FNULL: subprocess.call('python Visualization_files/Paraview_field_CSF_ref.py', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    with open(os.devnull, 'w') as FNULL: subprocess.call('python Visualization_files/Paraview_CSFref.py', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    
    if d['Show_paraview_screenshots']==1 and d["Skip_mesh_refinement"]==0:
        subprocess.call('xdg-open "Images/CSF_ref.png"',shell=True)
        subprocess.call('xdg-open "Images/Field_Adapted_CSF.png"',shell=True)
    
    
    #====================Mesh adaptive algorithm==================================#

    if d["Adapted_mesh_ready"]==0:
        
        from Mesh_adaption import adapt_mesh
        
        Ampl_on_vert=adapt_mesh(MRI_param,DTI_param,Scaling_CSF,Domains,d,anisotrop,cc_multicontact)     #also saves adapted mesh
        np.savetxt('Results_adaptive/Ampl_on_vert.csv', Ampl_on_vert, delimiter=" ")    #just to check
     
    #subprocess.call('python Visualization_files/Paraview_adapted.py', shell=True)
    with open(os.devnull, 'w') as FNULL: subprocess.call('python Visualization_files/Paraview_adapted.py', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    if d['Show_paraview_screenshots']==1 and d["Skip_mesh_refinement"]==0:  
        subprocess.call('xdg-open "Images/Adapted_mesh.png"',shell=True)
            
    
    #=========================Signal creation=====================================#
    
    Phi_vector=[x for x in d["Phi_vector"] if x is not None] #now we don't need None values (floating potentials), Contacts contain only the active ones
    Domains.fi=Phi_vector

    if d["current_control"]==1 and not (cc_multicontact==True):     #if multicontact, then it is already scaled in parallel comp.
        A=max(Domains.fi[:], key=abs)
        Phi_max=-1      #not needed
    else:
        A=1.0 #for vc we will apply Phi_vector as boundary conditions with appropriate amplitude
        Phi_max=max(Domains.fi[:], key=abs)

             
    if d["signal_generation_ready"]==0:
        print("----- Generating DBS signal -----")
        print(" ")
        
        from Signal_generator import generate_signal
        [t_vector,signal_out,Xs_signal_norm,FR_vector_signal]=generate_signal(d,A,Phi_max)
                        
        Xs_storage=np.zeros((np.real(Xs_signal_norm).shape[0],2),float)
        
        Xs_storage[:,0]=np.real(Xs_signal_norm)
        Xs_storage[:,1]=np.imag(Xs_signal_norm)
        np.savetxt('Stim_Signal/Xs_storage_full.csv', Xs_storage, delimiter=" ")

        np.savetxt('Stim_Signal/t_vector.csv', t_vector, delimiter=" ")
        np.savetxt('Stim_Signal/FR_vector_signal.csv', FR_vector_signal, delimiter=" ")
    else:
        Xs_recovered = np.genfromtxt('Stim_Signal/Xs_storage_full.csv', delimiter=' ')
        Xs_signal_norm=np.vectorize(complex)(Xs_recovered[:,0],Xs_recovered[:,1])
        FR_vector_signal = np.genfromtxt('Stim_Signal/FR_vector_signal.csv', delimiter=' ')
        t_vector = np.genfromtxt('Stim_Signal/t_vector.csv', delimiter=' ')  
        
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
            if not os.path.isfile('Field_solutions/Phi_real_'+str(d["freq"])+'Hz.pvd'):     #to make sure that there were interrupted computations 
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
    if d["IFFT_ready"] == 0 and d["Full_Field_IFFT"] == 1:
        from Full_IFFT_field_function import get_field_in_time
        
        #here is hardwired for my UQ study!!!
        VTA_size = get_field_in_time(d,FR_vector_signal,Xs_signal_norm,t_vector)        # also uses data from Field_solutions_functions/ and if spectrum truncation is applied, than data from Stim_Signal/
                                                                    # if VTA_from_NEURON is enabled, will save pointwise solutions in Points_in_time/ 
        d["IFFT_ready"] = 1               #modification of d
        
        if d["VTA_from_divE"]==True or d["VTA_from_E"]==True:
            #np.savetxt('UQ_stuff/VTA_'+str(Study_step)+'.csv', np.array([VTA_size,VTA_SN,VTA_STN,VTA_EPN]), delimiter=" ")
            #print("Simulation run is completed")
            minutes=int((time.time() - start_simulation_run)/60)
            secnds=int(time.time() - start_simulation_run)-minutes*60
            total_seconds=time.time() - start_simulation_run
            print("---Simulation run took ",minutes," min ",secnds," s ") 
            
            return None,VTA_size

        #else it will just jump to  NEURON_sep_points.py                
    
          
#==============Truncation of the already computed spectrum====================#
        
    if d["Truncate_the_obtained_full_solution"] == 1 and d["IFFT_ready"] == 0:
        if d["spectrum_trunc_method"]=='High Amplitude Method' or d["spectrum_trunc_method"] == 'Cutoff Method':
            print("----- Conducting IFFT truncating already obtained full solution -----")
            
            
            if isinstance(d["n_Ranvier"],list):             #if different populations
                from Field_IFFT_on_different_axons import convolute_signal_with_field_and_compute_ifft
                last_point=0
                for i in range(len(d["n_Ranvier"])):
                    last_point=convolute_signal_with_field_and_compute_ifft(N_models[i],last_point,d,Xs_signal_norm,N_segm[i],FR_vector_signal,t_vector,A,'Field_solutions/sorted_solution.csv')
            else:
                from Field_IFFT_on_axons import convolute_signal_with_field_and_compute_ifft
                convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm,N_models,N_segm,FR_vector_signal,t_vector,A,'Field_solutions/sorted_solution.csv')
        
        else:        
            print("Truncation of the obtained full solution is only for high. ampl and cutoff methods")
        #manipulate with solution_sort, TBD
    
    #=============================Parall. IFFT====================================#      
    
    if d["IFFT_ready"]==0 and d["Truncate_the_obtained_full_solution"]!=1:
        print("----- Conducting signal scaling and IFFT -----")
        
        #scale_signal will save pointwise solutions in Points_in_time/ 

        if d["spectrum_trunc_method"]=='No Truncation':
            if isinstance(d["n_Ranvier"],list):             #if different populations
                from Field_IFFT_on_different_axons import convolute_signal_with_field_and_compute_ifft
                last_point=0
                for i in range(len(d["n_Ranvier"])):
                    last_point=convolute_signal_with_field_and_compute_ifft(N_models[i],last_point,d,Xs_signal_norm,N_segm[i],FR_vector_signal,t_vector,A,name_sorted_solution)
            else:
                from Field_IFFT_on_axons import convolute_signal_with_field_and_compute_ifft
                convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm,N_models,N_segm,FR_vector_signal,t_vector,A,name_sorted_solution)  
                       
        if (d["spectrum_trunc_method"]=='High Amplitude Method' or d["spectrum_trunc_method"]=='Cutoff Method') and d["Truncate_the_obtained_full_solution"]==0:

            if isinstance(d["n_Ranvier"],list):             #if different populations
                from Field_IFFT_on_different_axons import convolute_signal_with_field_and_compute_ifft
                last_point=0
                for i in range(len(d["n_Ranvier"])):
                    last_point=convolute_signal_with_field_and_compute_ifft(N_models[i],last_point,d,Xs_signal_norm_new,N_segm[i],FR_vector_signal_new,t_vector,A,name_sorted_solution)
            else:
                from Field_IFFT_on_axons import convolute_signal_with_field_and_compute_ifft
                convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm_new,N_models,N_segm,FR_vector_signal_new,t_vector,A,name_sorted_solution)  

        if d["spectrum_trunc_method"]=='Octave Band Method':

            if isinstance(d["n_Ranvier"],list):             #if different populations
                from Field_IFFT_on_different_axons import convolute_signal_with_field_and_compute_ifft
                last_point=0
                for i in range(len(d["n_Ranvier"])):
                    last_point=convolute_signal_with_field_and_compute_ifft(N_models[i],last_point,d,Xs_signal_norm,N_segm[i],FR_vector_signal,t_vector,A,name_sorted_solution,inx_start_octv)
            else:
                from Field_IFFT_on_axons import convolute_signal_with_field_and_compute_ifft
                convolute_signal_with_field_and_compute_ifft(d,Xs_signal_norm,N_models,N_segm,FR_vector_signal,t_vector,A,name_sorted_solution,inx_start_octv)  


        
    '''Just to compute impedance in time, on if CPE is added or current-control mode'''
    #print(d["CPE_activ"],d["current_control"],d["spectrum_trunc_method"],d["IFFT_ready"])
    if (d["CPE_activ"]==1 or d["current_control"]==1) and cc_multicontact==False and d["spectrum_trunc_method"]!='Octave Band Method' and d["IFFT_ready"]==0:        #modify later
        from Field_IFFT_on_axons import compute_Z_ifft
        print("-----Calculating impedance over time-----")
        print(" ")
        if d["spectrum_trunc_method"]=='No Truncation' or d["Truncate_the_obtained_full_solution"]==1:
            Imp_in_time=compute_Z_ifft(d,Xs_signal_norm,t_vector,A)
        else:
            Imp_in_time=compute_Z_ifft(d,Xs_signal_norm_new,t_vector,A)
        
    
    #=============================Neural model====================================#

    print("-----Estimating neuron activity-----")

   
   
    os.chdir("Axon_files/")
    with open(os.devnull, 'w') as FNULL: subprocess.call('nocmodl axnode.mod', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    with open(os.devnull, 'w') as FNULL: subprocess.call('nrnivmodl', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    #'''We have to call Neuron_sep_points like this because the NEURON code is for python 2.7'''
    

    #print(os.getcwd())
    if isinstance(d["n_Ranvier"],list):
        start_neuron=time.time()
        from Axon_files.NEURON_direct_run import run_simulation_with_NEURON
        #from Parameter_insertion import paste_dif_neuron_model_param
        last_point=0
        for i in range(len(d["n_Ranvier"])):
            run_simulation_with_NEURON(last_point,i,d["diam_fib"][i],1000*d["t_step"],1000.0/d["freq"],d["n_Ranvier"][i],N_models[i],d["v_init"],t_vector.shape[0],d["Ampl_scale"],d["number_of_processors"])
            os.chdir("Axon_files/")
            last_point=N_segm[i]*N_models[i]+last_point
         
        minutes=int((time.time() - start_neuron)/60)
        secnds=int(time.time() - start_neuron)-minutes*60
        print("----- NEURON calculations took ",minutes," min ",secnds," s -----")   


    else:
        
        from Axon_files.NEURON_direct_run import run_simulation_with_NEURON 
        run_simulation_with_NEURON(0,-1,d["diam_fib"],1000*d["t_step"],1000.0/d["freq"],d["n_Ranvier"],N_models,d["v_init"],t_vector.shape[0],d["Ampl_scale"],d["number_of_processors"])
        os.chdir("Axon_files/")
    
    os.chdir("..")

#=====================Example for UQ==========================================#
    if os.path.isfile('Field_solutions/Activation/Last_run.csv'):
        Activated_neurons = np.genfromtxt('Field_solutions/Activation/Last_run.csv', delimiter=' ')
        #Number_of_activated=Activated_neurons.shape[0]
        if Activated_neurons.size == 1:
            Number_of_activated = 1
        else:
            Number_of_activated=Activated_neurons.shape[0]
    else:
        Number_of_activated = 0
        Activated_neurons=0

#=============================================================================#
       
#    print("Number_of_activated: ", Number_of_activated)
#        
#        
#    '''For the particular study'''        
#    N_models_total=900     #number of initially located models 
#    
#    model_activ_array=np.zeros(int(N_models_total),float)
#    for i in range(int(N_models_total)):
#        if np.any(Activated_neurons==i):
#            model_activ_array[i]=1.0
            
    #with open('UQ_stuff/Active_models'+str(Study_step)+'.csv','w') as f_handle:
    #    np.savetxt(f_handle,model_activ_array)
#            
#    models_vector=np.arange(int(N_models))
   
    
#    STN_axons=300.0
#    EPN_axons=300.0
#    SN_axons=300.0
#    
#    
#    STN_threshold=STN_axons
#    EPN_threshold=STN_axons+EPN_axons
#    SN_threshold=STN_axons+EPN_axons+SN_axons  
#    
#    Overall_activ=0.0
#    STN_activ=0.0
#    EPN_activ=0.0
#    SN_activ=0.0
#    #model_activ_array=np.zeros(int(N_models),float)
#    for i in xrange(int(N_models_total)):     #here we should have 
#        if np.any(Activated_neurons==i):
#            Overall_activ=Overall_activ+1.0
#            if i<STN_threshold:
#                STN_activ=STN_activ+1.0
#            if i>=STN_threshold and i<EPN_threshold:
#                STN_activ=STN_activ+1.0
#            if i>=EPN_threshold and i<SN_threshold:
#                STN_activ=STN_activ+1.0
    
    #return models_vector,model_activ_array
    
    #activ_vector=(Overall_activ, STN_activ, EPN_activ, SN_activ)
    #step_UQ
    
    #number_of_loc=np.arange(4.0)
    minutes=int((time.time() - start_simulation_run)/60)
    secnds=int(time.time() - start_simulation_run)-minutes*60
    total_seconds=time.time() - start_simulation_run
    print("---Simulation run took ",minutes," min ",secnds," s ")  
    
    #return None,Number_of_activated
    return total_seconds,None


#==========================Uncertainty Quantification====================================))======#
#import uncertainpy as un
#import chaospy as cp                 # To create distributions
##from LauncherUQ import run_full_model
##subprocess.call('python LauncherUQ.py', shell=True)         #do not forget to use appropriate Input_dic
##from Inp_dictUQ import d		#input dictionary
##locals().update(d)

#import os
#if not os.path.isdir('UQ_stuff'):
    #os.makedirs('UQ_stuff')


#
##run_full_model(0.2)
#
##model = un.Model(run=run_full_model,labels=["Model number", "Status"])
#model = un.Model(run=run_full_model,labels=["Location", "Number_of_activated_axons"])
#
##encap_thickness_dist = cp.Truncnorm(lo=0.1, up=0.4, mu=0.2, sigma=0.1)  #better have beta?
#encap_thickness_dist = cp.Beta(11, 61.5, lower=0.0, upper=1.0)   #better have beta? CHECK THIS!
#
##encap_thickness_dist = cp.Beta(3.0, 4.2, lo=0.05, up=0.3)   #better have beta? CHECK THIS!
#
##encap_cond_dist = cp.Beta(1.1, 5, lower=0.05, upper=0.2)   #EMBC first edition
#encap_cond_dist = cp.Beta(18.0, 40, lower=0.00, upper=0.2)   #better have beta? CHECK THIS!
#
#y_t_dist=cp.Truncnorm(lower=-0.45, upper=0.45, mu=0, sigma=0.2)
#
#parameters = {"encap_thickness_UQ": encap_thickness_dist,"encap_cond_UQ": encap_cond_dist,"y_t_UQ": y_t_dist}
#
## Set up the uncertainty quantification
#UQ = un.UncertaintyQuantification(model=model,
#                                  parameters=parameters,CPUs=None)
#
##'''Visualization of distribution'''
##X = cp.Truncnorm(lower=-0.45, upper=0.45, mu=0, sigma=0.2)    #better have beta? CHECK THIS!
##Y = cp.Iid(X, 1000)
##import matplotlib.pyplot as plt
##plt.plot(Y.sample(),np.arange(1000), 'ro')     #1000 probes
##plt.show()
#
#
## Perform the uncertainty quantification using
## polynomial chaos with point collocation (by default)
#data = UQ.quantify()
#
#
#variance=data['run_full_model'].variance
#print("Variance: ",variance)
#mean=data['run_full_model'].mean
#print("Mean: ",mean)
#
#Sobol=data['run_full_model'].sobol_first
#print(Sobol)
#
#P5=data['run_full_model'].percentile_5
#print(P5)
#
#P95=data['run_full_model'].percentile_95
#print(P95)

#===========================Master Study=======================================================#



#Number_of_activated=run_full_model()

'''Master study is useful for optimization and UQ studies'''
'''It will evaluate which truncation method and with how many frequencies we need to use for a particular setup'''
'''If the adataptive mesh alg. are not usable in the study, it will also estimate the difference in potential for 2 profile salome scripts (profile scripts should be created manually)'''

def run_master_study():

    '''Computations on initial mesh generated in SALOME'''
    
    Salome_profile_list=['SNEX100','SNEX100_UQ2','SNEX100_UQ3']
    Salome_best_profile='SNEX100_UQ4'   #name of the script with a highest mesh requirements (which might be too strict for iterative studies). Will be used as a benchmark.
    
    master_dict={'Electrode_type':Salome_best_profile}
    #total_seconds_bench,Number_of_activated_bench=run_full_model(master_dict)
    #total_seconds_bench=4416.0
    
    
    Result_diff_profiles=np.zeros((len(Salome_profile_list),2),float)
    
    for i in range(len(Salome_profile_list)):
        master_dict={'Electrode_type':Salome_profile_list[i], 'voxel_arr_MRI':1, 'voxel_arr_DTI':1, 'Init_neuron_model_ready':1}        #some steps were already done
        Result_diff_profiles[i,:]=run_full_model(master_dict)
        
        
    print("Number of activated neurons, the benchmark and the profiles : ", Number_of_activated_bench,Result_diff_profiles[:,1])
    print("Time is seconds, the benchmark and the profiles : ", total_seconds_bench,Result_diff_profiles[:,0])
    print("Choose the best option and put it in input dictionary. Comment out '''Computations on initial mesh generated in SALOME''' section")
    
    raise SystemExit
    
    
    
    '''Spectrum truncation'''
    
    Number_of_activated_high_ampl=0               #initizialization
    Number_of_activated_octaves=0               #initizialization
    
    #run_full_model(spectrum_truncation_method,truncate_obtained_full,trunc_param=0)
    #first, run for the full spectrum
    master_dict={"spectrum_trunc_method":'No Truncation', "Truncate_the_obtained_full_solution":0}
    Run_time,Number_of_activated_benchmark=run_full_model(master_dict)
   
    #"trunc_param"
    #secondly, extract results from the full solution on freqs. with high ampl. of FFT
    trunc_param_N_freqs=40      #will start with truncation to 40 frequencies
    while Number_of_activated_high_ampl!=Number_of_activated_benchmark:  #here maybe we could give a percentage of deviation
        master_dict={"spectrum_trunc_method":'High Amplitude Method', "Truncate_the_obtained_full_solution":1, "trunc_param":trunc_param_N_freqs}
        Run_time_high_ampl,Number_of_activated_high_ampl=run_full_model(master_dict)
        trunc_param_N_freqs=trunc_param_N_freqs+25
        
    #now run the same study, but without extracting from the full solution 
    master_dict={"spectrum_trunc_method":'High Amplitude Method', "Truncate_the_obtained_full_solution":0, "trunc_param":trunc_param_N_freqs-25}
    Run_time_high_ampl,Number_of_activated_high_ampl=run_full_model(master_dict)
#        
#    #trunc_param_octave_freq=10000           #compute manually, give an example  
#    #if repit. rate is 130 Hz, max freq is 500000 and cutoff is 10000
#    #(130*octave_scale)/np.sqrt(2.0)   octave_scale=2**x,     x-number of freqs after cut off
#    #130*2**x/np.sqrt(2.0)>490000       x>log2(490000*np.sqrt(2.0)/130)
  
        
    def func_cutoff_oct(cutoff_oct,*data_pass):     #defines the frequency after which octaves will be placed depending on the total number of freqs (the final total number can deviate by 1 frequency)
        from Inp_dict import d
        Sim_time=1.0/d["freq"]      #one pulse per simulation
        freq_max=d["freq"]*Sim_time/d["t_step"]
        
        return int(np.log2((freq_max-cutoff_oct)*np.sqrt(2.0)/d["freq"])+cutoff_oct/d["freq"])+1-N_freqs_oct

    #N_freqs_oct=int(np.log2((freq_max-cutoff_oct)*np.sqrt(2.0)/d["freq"])+cutoff_oct/d["freq"])+1
    from scipy.optimize import fsolve

    N_freqs_oct=40
    while Number_of_activated_octaves!=Number_of_activated_benchmark:  #here maybe we could give a percentage of deviation
        data_pass=(N_freqs_oct)
        cutoff_octaves=fsolve(func_cutoff_oct,10,args=data_pass)
        
        master_dict={"spectrum_trunc_method":'Octave Band Method', "Truncate_the_obtained_full_solution":0, "trunc_param":cutoff_octaves}
        Run_time_oct,Number_of_activated_octaves=run_full_model(master_dict)
        N_freqs_oct=N_freqs_oct+5
        if Run_time_oct>Run_time_high_ampl:
            print("Octave method is already slower than high ampl., no need to continue.")
            break
        
        
    if Run_time_oct<Run_time_high_ampl and Run_time_oct<Run_time:
        print("Octave truncation is the fastest")       
    elif Run_time_oct>Run_time_high_ampl and Run_time_high_ampl<Run_time:
        print("High amplitude truncation is the fastest")       
    else:
        print("Truncation did not accelerate the computations")
                
    print("Time for Full_run,High_ampl,Octaves: ", Run_time,Run_time_high_ampl,Run_time_oct)
    print("Activation for Full_run,High_ampl,Octaves: ", Number_of_activated_benchmark,Number_of_activated_high_ampl,Number_of_activated_octaves)
        
#run_master_study() 
#Run_time_high_ampl,Number_of_activated_high_ampl=run_full_model(1,0,50)
#master_dict={"spectrum_trunc_method":3, "Truncate_the_obtained_full_solution":0, "trunc_param":5000.0}
master_dict={}
Run_time_high_ampl,Number_of_activated_high_ampl=run_full_model(master_dict)
