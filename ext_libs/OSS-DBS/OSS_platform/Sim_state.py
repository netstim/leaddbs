#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Butenko K.

Function check_state(inp_dict) will shift the simulation flow to the last stage marked by the user in the input dictionary
It will also manage the output folders accordingly using manage_folders(inp_dict)

"""

# import numpy as np
import os
import shutil
import logging


def copy_rename(old_file_name, new_file_name):
    src_dir = os.curdir
    dst_dir = os.path.join(os.curdir, "subfolder")
    src_file = os.path.join(src_dir, old_file_name)
    shutil.copy(src_file, dst_dir)

    dst_file = os.path.join(dst_dir, old_file_name)
    new_dst_file_name = os.path.join(dst_dir, new_file_name)
    os.rename(dst_file, new_dst_file_name)


def manage_folders(d):
    ##print('Some results from previous simulations will be deleted')
    ##warning=str(raw_input('Enter STOP to exit\n'))
    # warning='go'
    # if warning=='STOP' or warning=='Stop' or warning=='stop':
    #    print("exiting")
    #    raise Exception('exit')

    if os.path.isdir(os.environ['PATIENTDIR'] + '/Tensors') and d["Parallel_comp_ready"] != 1:
        shutil.rmtree(os.environ['PATIENTDIR'] + '/Tensors')

    if not os.path.isdir(os.environ['PATIENTDIR'] + '/Tensors'):
        os.makedirs(os.environ['PATIENTDIR'] + '/Tensors')

    if not os.path.isdir(os.environ['PATIENTDIR'] + '/Images'):
        os.makedirs(os.environ['PATIENTDIR'] + '/Images')
    elif d["Init_neuron_model_ready"] == 0:  # a totally new simulation, old images can be deleted
        shutil.rmtree(os.environ['PATIENTDIR'] + '/Images')
        os.makedirs(os.environ['PATIENTDIR'] + '/Images')

    if d["Segm_MRI_processed"] == 0 and d["DTI_processed"] == 0:
        if os.path.isdir(os.environ['PATIENTDIR'] + '/MRI_DTI_derived_data'):
            shutil.rmtree(os.environ['PATIENTDIR'] + '/MRI_DTI_derived_data')
        os.makedirs(os.environ['PATIENTDIR'] + '/MRI_DTI_derived_data')
    if d["Init_mesh_ready"] != 1:
        if os.path.isdir(os.environ['PATIENTDIR'] + '/Meshes'):
            shutil.rmtree(os.environ['PATIENTDIR'] + '/Meshes')
        os.makedirs(os.environ['PATIENTDIR'] + '/Meshes')
    if d["CSF_mesh_ready"] != 1:
        if os.path.isdir(os.environ['PATIENTDIR'] + '/CSF_ref'):
            shutil.rmtree(os.environ['PATIENTDIR'] + '/CSF_ref')
        os.makedirs(os.environ['PATIENTDIR'] + '/CSF_ref')
    if d["Adapted_mesh_ready"] != 1:
        if os.path.isdir(os.environ['PATIENTDIR'] + '/Results_adaptive'):
            shutil.rmtree(os.environ['PATIENTDIR'] + '/Results_adaptive')
        os.makedirs(os.environ['PATIENTDIR'] + '/Results_adaptive')
    if d["Signal_generated"] != 1:
        if os.path.isdir(os.environ['PATIENTDIR'] + '/Stim_Signal'):
            shutil.rmtree(os.environ['PATIENTDIR'] + '/Stim_Signal')
        os.makedirs(os.environ['PATIENTDIR'] + '/Stim_Signal')
    if d["Parallel_comp_ready"] != 1 and d["Parallel_comp_interrupted"] != 1:
        if os.path.isdir(os.environ['PATIENTDIR'] + '/Field_solutions'):
            shutil.rmtree(os.environ['PATIENTDIR'] + '/Field_solutions')
        os.makedirs(os.environ['PATIENTDIR'] + '/Field_solutions')
        os.makedirs(os.environ['PATIENTDIR'] + '/Field_solutions/Activation')
        os.makedirs(os.environ['PATIENTDIR'] + '/Field_solutions/Animation_files')
        if os.path.isdir(os.environ['PATIENTDIR'] + '/Field_solutions_functions'):
            shutil.rmtree(os.environ['PATIENTDIR'] + '/Field_solutions_functions')
        os.makedirs(os.environ['PATIENTDIR'] + '/Field_solutions_functions')
    if d["IFFT_ready"] != 1:
        if os.path.isdir(os.environ['LGFDIR'] + '/Axons_in_time'):
            os.system('rm -fr ' + os.environ['LGFDIR'] + '/Axons_in_time')
        os.makedirs(os.environ['LGFDIR'] + '/Axons_in_time')
        if os.path.isdir(os.environ['PATIENTDIR'] + '/Animation_Field_in_time'):
            shutil.rmtree(os.environ['PATIENTDIR'] + '/Animation_Field_in_time')
        os.makedirs(os.environ['PATIENTDIR'] + '/Animation_Field_in_time')
    if d["Init_neuron_model_ready"] == 0 and d["Adjusted_neuron_model_ready"] == 0:
        if os.path.isdir(os.environ['PATIENTDIR'] + '/Neuron_model_arrays'):
            shutil.rmtree(os.environ['PATIENTDIR'] + '/Neuron_model_arrays')
        os.makedirs(os.environ['PATIENTDIR'] + '/Neuron_model_arrays')
    if os.path.isdir(os.environ['PATIENTDIR'] + '/Field_solutions/Activation'):  # we always re-run NEURON simulation
        shutil.rmtree(os.environ['PATIENTDIR'] + '/Field_solutions/Activation')
        os.makedirs(os.environ['PATIENTDIR'] + '/Field_solutions/Activation')
    if d['Stim_side'] == 0:
        if os.path.isdir(os.environ['PATIENTDIR'] + '/Results_rh'):
            shutil.rmtree(os.environ['PATIENTDIR'] + '/Results_rh')
        os.makedirs(os.environ['PATIENTDIR'] + '/Results_rh')
    if d['Stim_side'] == 1:
        if os.path.isdir(os.environ['PATIENTDIR'] + '/Results_lh'):
            shutil.rmtree(os.environ['PATIENTDIR'] + '/Results_lh')
        os.makedirs(os.environ['PATIENTDIR'] + '/Results_lh')

    return True


def check_state(d):
    if d['number_of_processors'] == 0:
        physical_cores = os.popen("""lscpu -b -p=Core,Socket | grep -v '^#' | sort -u | wc -l""").read()[:-1]
        d['number_of_processors'] = int(
            physical_cores)  # this option is only active if Docker App is used (on macOS and Windows)
        logging.critical("Number of cores available for Docker: {}".format(d['number_of_processors']))

    logging.critical("Number of processors used: {}".format(d['number_of_processors']))

    if d["IFFT_ready"] == 1:
        d["Segm_MRI_processed"] = 1
        d["Init_mesh_ready"] = 1
        d["Init_neuron_model_ready"] = 1
        d["Adjusted_neuron_model_ready"] = 1
        d["CSF_mesh_ready"] = 1
        d["Adapted_mesh_ready"] = 1
        d["Parallel_comp_ready"] = 1
        d["Signal_generated"] = 1
    if d["Parallel_comp_ready"] == 1:
        d["Segm_MRI_processed"] = 1
        d["Init_mesh_ready"] = 1
        d["Init_neuron_model_ready"] = 1
        d["Adjusted_neuron_model_ready"] = 1
        d["CSF_mesh_ready"] = 1
        d["Adapted_mesh_ready"] = 1
        d["Signal_generated"] = 1
    #    if d["Signal_generated"]==1:
    #        d["Adapted_mesh_ready"]=1
    #        d["Segm_MRI_processed"]=1
    #        d["Init_mesh_ready"]=1
    #        d["Init_neuron_model_ready"]=1
    #        d["Adjusted_neuron_model_ready"]=1
    #        d["CSF_mesh_ready"]=1
    if d["Adapted_mesh_ready"] == 1:
        d["Segm_MRI_processed"] = 1
        d["Init_mesh_ready"] = 1
        d["Init_neuron_model_ready"] = 1
        d["Adjusted_neuron_model_ready"] = 1
        d["CSF_mesh_ready"] = 1
    if d["CSF_mesh_ready"] == 1:
        d["Segm_MRI_processed"] = 1
        d["Init_mesh_ready"] = 1
        d["Init_neuron_model_ready"] = 1
        d["Adjusted_neuron_model_ready"] = 1
    if d["Adjusted_neuron_model_ready"] == 1:
        d["Segm_MRI_processed"] = 1
        d["Init_neuron_model_ready"] = 1
        d["Init_mesh_ready"] = 1
    if d["Init_mesh_ready"] == 1:
        d["Init_neuron_model_ready"] = 1
        d["Segm_MRI_processed"] = 1
    if d["Init_neuron_model_ready"] == 1:
        d["Segm_MRI_processed"] = 1

    manage_folders(d)

    logging.critical("Folders were adjusted\n")

    return True
