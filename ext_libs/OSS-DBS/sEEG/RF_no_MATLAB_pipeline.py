#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 14:43:20 2026

@author: forel
"""

import os
import shutil
import subprocess
import numpy as np
import pandas as pd
#from scipy.io import loadmat
import json

AMPLITUDE = 1.0 # in mA
HEMI_SUFFIX = ['_rh','_lh']

from ossdbs.electrodes.defaults import default_electrode_parameters

def check_electrode_availability(reco_electrode):
    
    """
    Check if OSS supports the electrode model

    Args:
        reco_electrode: str, electrode name as defined in ea_resolve_elspec.m
    Returns:
        str, OSS label for the electrode
        
    """

        
    electrode_names = {

        "ELAINE Rat Electrode": "MicroProbesRodentElectrode",
        "JR MicroElectrode":"MicroElectrode",  # not sure about the name in Lead-DBS
        "LS MicroElectrode":"MicroElectrode",  # I just changed the parameters
    }

    for lead in electrode_names.keys():
        if lead == reco_electrode:
            electrode_name = electrode_names[lead]
            return electrode_name
    
    return False


if __name__ == '__main__':
    
    electrode_coords = pd.read_csv('/home/forel/Documents/data/LS/contact_coords.csv',delimiter=',')
    
    tail_coords = np.zeros((electrode_coords.shape[0],3),float)
    
    for el_i in range(38,tail_coords.shape[0]):
        micro_tip = np.array([electrode_coords['x_micro'][el_i],electrode_coords['y_micro'][el_i],electrode_coords['z_micro'][el_i]])
        contact = np.array([electrode_coords['x_macro'][el_i],electrode_coords['y_macro'][el_i],electrode_coords['z_macro'][el_i]])
    
        directions = contact - micro_tip
        unit_directions = directions / np.linalg.norm(directions)
        
        # 1 mm shifted upwards
        tail_coords[el_i,:] = 1.0*unit_directions + contact
        
        InputDict = {
            'stimFolder': '/home/forel/Documents/data/LS/stim_' + str(el_i),
            'path2segmask': '/home/forel/Documents/GitHub/leaddbs/templates/space/MNI152NLin2009bAsym/segmask.nii',
            'StimMode': 'CC',
            'hemi_idx': 0, # 0 for RH, 1 for LH
            'el_type': 'LS MicroElectrode',
            'remove_el': True,
            # coordinates
            'contact_coords_x': [contact[0]],
            'contact_coords_y': [contact[1]],
            'contact_coords_z': [contact[2]],
            # provide tail coord for monopolar contacts
            'tail_coords': list(tail_coords[el_i,:]),
            'first_contact_offset': 0.5,  # distance from the very tip of the electrode to the center of the first contact
            # output parameters
            'nii_resolution': 0.5,
            'nii_voxel_number': 71, # per edge
            
        }

        if os.path.exists(InputDict['stimFolder']) and os.path.isdir(InputDict['stimFolder']):
            shutil.rmtree(InputDict['stimFolder'])
        os.mkdir(InputDict['stimFolder'])
    
        if InputDict['StimMode'] == 'CC':
            current_controlled = 1
        else:
            current_controlled = 0
    
        # ToDo: check if the electrode is directional
        oss_electrode = check_electrode_availability(InputDict['el_type'])  
        elec_params = default_electrode_parameters[oss_electrode]
    
        # pass one electrode at a time
        contact_coords = np.array([InputDict['contact_coords_x'],InputDict['contact_coords_y'],InputDict['contact_coords_z']])
    
        # for unilateral implantations to make sure the pipeline does not break, simply copy the values from one side to another
        # they will not be actually used!
        # if contact_coords_rh == None:
        #     contact_coords_rh, marker_head_coords_rh, marker_y_coords_rh, el_type_rh = (contact_coords_lh, marker_head_coords_lh, marker_y_coords_lh, el_type_lh)
        #     el_type = el_type_lh
        # elif contact_coords_lh == None:
        #     contact_coords_lh, marker_head_coords_lh, marker_y_coords_lh, el_type_lh = contact_coords_rh, marker_head_coords_rh, marker_y_coords_rh, el_type_rh
        #     el_type = el_type_rh
        
        
        ### Update the default OSS-DBS dictionary and run simulations
        from custom_dict import oss_params
        
        # iterate over contacts changing Phi and StimCenter
        N_contacts = contact_coords.shape[1]
        
        hemi_idx = InputDict['hemi_idx']
        first_contact = contact_coords[:,0]
        if N_contacts != 1:
            last_contact = contact_coords[:,-1]
            actual_span = np.linalg.norm(last_contact - first_contact)
        else:
            actual_span = 0.0
        
        for cnt_i in range(N_contacts):
    
            hemi_sim_folder = os.path.join(InputDict['stimFolder'],'OSS_sim_files' + HEMI_SUFFIX[hemi_idx])
            if os.path.exists(hemi_sim_folder) and os.path.isdir(hemi_sim_folder):
                shutil.rmtree(hemi_sim_folder)
            os.mkdir(hemi_sim_folder)
            
            stim_protocol = np.zeros(N_contacts,float)
            stim_protocol[cnt_i] = AMPLITUDE
                    
            # grounding parameters
            if current_controlled:
                # switch from mA to A
                grounded_current = -0.001 * np.round(np.sum(stim_protocol), 9)  # could be 0
                case_grounding = True # if no current is actually grounded, this is used to reference voltages 
            else:
                grounded_current = 0.0 # not relevant
                print("Voltage-controlled mode is used, case grounding is ON")
                case_grounding = True
            
            # get tip position from the head and tail markers
    
            if N_contacts == 1:
                directions = np.array(InputDict['tail_coords']) - first_contact
            else:
                directions = last_contact - first_contact
                
            unit_directions = directions / np.linalg.norm(directions)
            
            offset = InputDict['first_contact_offset']
            tip_pos = first_contact - offset * unit_directions
            
            oss_params["BrainRegion"]["Center"]["x[mm]"] = first_contact[0]
            oss_params["BrainRegion"]["Center"]["y[mm]"] = first_contact[1]
            oss_params["BrainRegion"]["Center"]["z[mm]"] = first_contact[2]
            
            oss_params["BrainRegion"]["Dimension"]["x[mm]"] = 50.0 + np.abs(unit_directions[0]) * actual_span * 2.0
            oss_params["BrainRegion"]["Dimension"]["y[mm]"] = 50.0 + np.abs(unit_directions[1]) * actual_span * 2.0
            oss_params["BrainRegion"]["Dimension"]["z[mm]"] = 50.0 + np.abs(unit_directions[2]) * actual_span * 2.0
            
            oss_params["Electrodes"][0]["Name"] = oss_electrode
            
            oss_params["Electrodes"][0]["Rotation[Degrees]"] = 0.0 # NO ROTATION ASSUMED!
            oss_params["Electrodes"][0]["Direction"]["x[mm]"] = unit_directions[0]
            oss_params["Electrodes"][0]["Direction"]["y[mm]"] = unit_directions[1]
            oss_params["Electrodes"][0]["Direction"]["z[mm]"] = unit_directions[2]
            
            oss_params["Electrodes"][0]["TipPosition"]["x[mm]"] = tip_pos[ 0]
            oss_params["Electrodes"][0]["TipPosition"]["y[mm]"] = tip_pos[1]
            oss_params["Electrodes"][0]["TipPosition"]["z[mm]"] = tip_pos[2]
            
            oss_params["PointModel"]["Lattice"]["Active"] = True
            oss_params["PointModel"]["Lattice"]["Center"]["x[mm]"] = first_contact[0]
            oss_params["PointModel"]["Lattice"]["Center"]["y[mm]"] = first_contact[1]
            oss_params["PointModel"]["Lattice"]["Center"]["z[mm]"] = first_contact[2]
            oss_params["PointModel"]["Lattice"]["Shape"] = {"x": InputDict['nii_voxel_number'], "y": InputDict['nii_voxel_number'], "z": InputDict['nii_voxel_number']}
            oss_params["PointModel"]["Lattice"]["PointDistance[mm]"] = InputDict['nii_resolution']
            oss_params["PointModel"]["Lattice"]["ExportField"] = True
            oss_params["PointModel"]["Lattice"]["CollapseVTA"] = InputDict['remove_el']
            
            oss_params["MaterialDistribution"]["MRIPath"] = InputDict['path2segmask']
            oss_params["StimulationSignal"]["CurrentControlled"] = current_controlled
            oss_params["Surfaces"][0]["Active"] = case_grounding
            oss_params["Surfaces"][0]["Current[A]"] = grounded_current
    
            oss_params["OutputPath"] = os.path.join(InputDict['stimFolder'],'Results_protocol_' + str(cnt_i) + HEMI_SUFFIX[hemi_idx])
    
            # contact dictionary is specified separately
    
            cnt_range = range(N_contacts)
                
            # cntct_dicts is a list of the contacts that will go into the dict
            # for this electrode
            cntct_dicts = np.empty(len(cnt_range), dtype=object)
            cntcts_made = 0
            
            for i in cnt_range:
                # all (truly) non-active contacts are floating with 0A
                if stim_protocol[i] == 0.0:
                    floating = True
                    cntct_dicts[cntcts_made] = {
                        # Assuming one-indexed contact ids
                        "Contact_ID": i + 1,
                        "Active": False,
                        "Current[A]": 0.0,
                        "Voltage[V]": 0.0,
                        "Floating": True,
                    }
                else:
                    # for current-controlled, we have a pseudo non-active contact
                    if current_controlled:
                        cntct_dicts[cntcts_made] = {
                            "Contact_ID": i + 1,  # OSS numbers contacts starting from 1!
                            "Active": False,
                            "Current[A]": stim_protocol[i]*0.001,
                            "Voltage[V]": 0.0,
                            "Floating": True,
                            "MaxMeshSizeEdge": 0.08545132017764238  # hardcoded for now,
                        }
                    else:
                        cntct_dicts[cntcts_made] = {
                            "Contact_ID": i + 1,  # OSS numbers contacts starting from 1!
                            "Active": True,
                            "Current[A]": 0.0,
                            "Voltage[V]": stim_protocol[i],
                            "Floating": False,
                            "MaxMeshSizeEdge": 0.08545132017764238  # hardcoded for now,
                        }
            
                cntcts_made += 1
            
            oss_params['Electrodes'][0]['Contacts'] = cntct_dicts.tolist()
    
            # save the settings
            json_settings = json.dumps(oss_params, indent=2)
            with open(os.path.join(InputDict['stimFolder'],'oss-dbs_parameters_adj.json'), "w") as outfile:
                outfile.write(json_settings)
                
            # run OSS-DBS
            with open(os.devnull, 'w') as FNULL: subprocess.call('ossdbs ' + os.path.join(InputDict['stimFolder'],'oss-dbs_parameters_adj.json'),shell=True)
    
                
        #         # clean-up
        
        #         lead_settings.add_stim_center(np.array((contact_coords_rh[:,cnt_i],contact_coords_lh[:,cnt_i])))
                
        #         # save .mat
        #         hdf5storage.savemat(mat_file_path, mat_data, format='7.3')
                
        #         with open(os.devnull, 'w') as FNULL: subprocess.call('leaddbs2ossdbs  --hemi_side ' + str(hemi_idx) + ' ' + os.path.join(stimFolder,'oss-dbs_parameters_adj.mat') + ' --output_path ' + hemi_sim_folder,shell=True)
        #         with open(os.devnull, 'w') as FNULL: subprocess.call('ossdbs ' + os.path.join(stimFolder,'oss-dbs_parameters.json'),shell=True)
        
        
                    # we should then copy it with a unique identifier and convert to the same space using MATLAB
        
        # settings = h5py.File(str(oss_settings_template), "r")
        # settings = settings['settings']
        
        # settings['Patient_folder']