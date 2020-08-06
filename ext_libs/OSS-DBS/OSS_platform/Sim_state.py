#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Butenko K.

Function check_state(inp_dict) will shift the simulation flow to the last stage marked by the user in the input dictionary
It will also manage the output folders accordingly using manage_folders(inp_dict)


"""



#import numpy as np
import os
import shutil

def copy_rename(old_file_name, new_file_name):
        src_dir= os.curdir
        dst_dir= os.path.join(os.curdir , "subfolder")
        src_file = os.path.join(src_dir, old_file_name)
        shutil.copy(src_file,dst_dir)
        
        dst_file = os.path.join(dst_dir, old_file_name)
        new_dst_file_name = os.path.join(dst_dir, new_file_name)
        os.rename(dst_file, new_dst_file_name)

def manage_folders(d):
    ##print('Some results from previous simulations will be deleted')
    ##warning=str(raw_input('Enter STOP to exit\n'))
    #warning='go'
    #if warning=='STOP' or warning=='Stop' or warning=='stop':
    #    print("exiting")
    #    raise Exception('exit')

    if os.path.isdir('/opt/Patient/Tensors') and d["Parallel_comp_ready"]!=1:
        shutil.rmtree('/opt/Patient/Tensors')

    if not os.path.isdir('/opt/Patient/Tensors'):
        os.makedirs('/opt/Patient/Tensors')
        
    if not os.path.isdir('/opt/Patient/Images'):
        os.makedirs('/opt/Patient/Images')
    elif d["Init_neuron_model_ready"]==0:   # a totally new simulation, old images can be deleted
        shutil.rmtree('/opt/Patient/Images')
        os.makedirs('/opt/Patient/Images')
        
    if d["voxel_arr_MRI"]==0 and d["voxel_arr_DTI"]==0:
        if os.path.isdir('/opt/Patient/MRI_DTI_derived_data'):
            shutil.rmtree('/opt/Patient/MRI_DTI_derived_data')
        os.makedirs('/opt/Patient/MRI_DTI_derived_data')
    if d["Init_mesh_ready"]!=1:
        if os.path.isdir('/opt/Patient/Meshes'):
            shutil.rmtree('/opt/Patient/Meshes')
        os.makedirs('/opt/Patient/Meshes')
    if d["CSF_mesh_ready"]!=1:
        if os.path.isdir('/opt/Patient/CSF_ref'):
            shutil.rmtree('/opt/Patient/CSF_ref')
        os.makedirs('/opt/Patient/CSF_ref')
    if d["Adapted_mesh_ready"]!=1:
        if os.path.isdir('/opt/Patient/Results_adaptive'):
            shutil.rmtree('/opt/Patient/Results_adaptive')
        os.makedirs('/opt/Patient/Results_adaptive')
    if d["signal_generation_ready"]!=1:
        if os.path.isdir('/opt/Patient/Stim_Signal'):
            shutil.rmtree('/opt/Patient/Stim_Signal')
        os.makedirs('/opt/Patient/Stim_Signal')
    if d["Parallel_comp_ready"]!=1 and d["Parallel_comp_interrupted"]!=1:
        if os.path.isdir('/opt/Patient/Field_solutions'):
            shutil.rmtree('/opt/Patient/Field_solutions')
        os.makedirs('/opt/Patient/Field_solutions')
        os.makedirs('/opt/Patient/Field_solutions/Activation')
        os.makedirs('/opt/Patient/Field_solutions/Animation_files')
        if os.path.isdir('/opt/Patient/Field_solutions_functions'):
            shutil.rmtree('/opt/Patient/Field_solutions_functions')
        os.makedirs('/opt/Patient/Field_solutions_functions')
    if d["IFFT_ready"]!=1:
        if os.path.isdir('/opt/Patient/Points_in_time'):
            shutil.rmtree('/opt/Patient/Points_in_time')
        os.makedirs('/opt/Patient/Points_in_time')
        if os.path.isdir('/opt/Patient/Animation_Field_in_time'):
            shutil.rmtree('/opt/Patient/Animation_Field_in_time')
        os.makedirs('/opt/Patient/Animation_Field_in_time')
    if (d["Init_neuron_model_ready"]==0 and d["Adjusted_neuron_model_ready"]==0):
        if os.path.isdir('/opt/Patient/Neuron_model_arrays'):
            shutil.rmtree('/opt/Patient/Neuron_model_arrays')
        os.makedirs('/opt/Patient/Neuron_model_arrays')

    if os.path.isdir('/opt/Patient/Field_solutions/Activation'):     # we always re-run NEURON simulation
        shutil.rmtree('/opt/Patient/Field_solutions/Activation')
        os.makedirs('/opt/Patient/Field_solutions/Activation')
    
    return True


def check_state(d):
      
    if d["IFFT_ready"]==1:
        d["voxel_arr_MRI"]=1     
        d["Init_mesh_ready"]=1        
        d["Init_neuron_model_ready"]=1
        d["Adjusted_neuron_model_ready"]=1
        d["CSF_mesh_ready"]=1
        d["Adapted_mesh_ready"]=1
        d["Parallel_comp_ready"]=1
        d["signal_generation_ready"]=1
    if d["Parallel_comp_ready"]==1:
        d["voxel_arr_MRI"]=1     
        d["Init_mesh_ready"]=1        
        d["Init_neuron_model_ready"]=1
        d["Adjusted_neuron_model_ready"]=1
        d["CSF_mesh_ready"]=1
        d["Adapted_mesh_ready"]=1
        d["signal_generation_ready"]=1
#    if d["signal_generation_ready"]==1:        
#        d["Adapted_mesh_ready"]=1
#        d["voxel_arr_MRI"]=1     
#        d["Init_mesh_ready"]=1        
#        d["Init_neuron_model_ready"]=1
#        d["Adjusted_neuron_model_ready"]=1
#        d["CSF_mesh_ready"]=1
    if d["Adapted_mesh_ready"]==1:
        d["voxel_arr_MRI"]=1     
        d["Init_mesh_ready"]=1        
        d["Init_neuron_model_ready"]=1
        d["Adjusted_neuron_model_ready"]=1
        d["CSF_mesh_ready"]=1
    if d["CSF_mesh_ready"]==1:
        d["voxel_arr_MRI"]=1     
        d["Init_mesh_ready"]=1        
        d["Init_neuron_model_ready"]=1
        d["Adjusted_neuron_model_ready"]=1
    if d["Adjusted_neuron_model_ready"]==1:
        d["voxel_arr_MRI"]=1
        d["Init_neuron_model_ready"]=1        
        d["Init_mesh_ready"]=1
    if d["Init_mesh_ready"]==1:
        d["Init_neuron_model_ready"]=1        
        d["voxel_arr_MRI"]=1        
    if d["Init_neuron_model_ready"]==1:
        d["voxel_arr_MRI"]=1     
        
    manage_folders(d)    
    
    print("Folders were adjusted\n")

    return True