#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 18:29:05 2026

@author: forel
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 19:05:08 2025

@author: forel
"""


import os
import sys
import shutil
import subprocess
import numpy as np
#from scipy.io import loadmat
import h5py
import json
from dataclasses import asdict

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
        "Abbott Directed 6172 (short)": "AbbottStJudeDirected6172",
        "Abbott Directed 6173 (long)": "AbbottStJudeDirected6173",
        "Abbott ActiveTip (6146-6149)": "AbbottStJudeActiveTip6146_6149",
        "Abbott ActiveTip (6142-6145)": "AbbottStJudeActiveTip6142_6145",
        "St. Jude ActiveTip (6146-6149)": "AbbottStJudeActiveTip6146_6149",
        "St. Jude ActiveTip (6142-6145)": "AbbottStJudeActiveTip6142_6145",
        "St. Jude Directed 6180": "AbbottStJudeDirected6172",
        "St. Jude Directed 6172 (short)": "AbbottStJudeDirected6172",
        "St. Jude Directed 6173 (long)": "AbbottStJudeDirected6173",
        "Boston Scientific Vercise": "BostonScientificVercise",
        "Boston Scientific Vercise Directed": "BostonScientificVerciseDirected",
        "Boston Scientific Vercise Cartesia HX": "BostonScientificCartesiaHX",
        "Boston Scientific Vercise Cartesia X": "BostonScientificCartesiaX",
        "ELAINE Rat Electrode": "MicroProbesRodentElectrode",
        "Medtronic 3387": "Medtronic3387",
        "Medtronic 3389": "Medtronic3389",
        "Medtronic 3391": "Medtronic3391",
        "SceneRay SR1210": "Medtronic3387",
        "SceneRay SR1200": "Medtronic3389",
        "Medtronic B33005": "MedtronicSenSightB33005",
        "Medtronic B33015": "MedtronicSenSightB33015",
        "PINS Medical L301": "PINSMedicalL301",
        "PINS Medical L302": "PINSMedicalL302",
        "PINS Medical L303": "PINSMedicalL303",
        "NeuroPace DL-344-3.5": "NeuroPaceDL344_3_5",
        "NeuroPace DL-344-10": "NeuroPaceDL344_10",
        "DIXI D08-05AM": "DixiSEEG5",
        "DIXI D08-08AM": "DixiSEEG8",
        "DIXI D08-10AM": "DixiSEEG10",
        "DIXI D08-12AM": "DixiSEEG12",
        "DIXI D08-15AM": "DixiSEEG15",
        "DIXI D08-18AM": "DixiSEEG18",
        "PMT 2102-08-091": "PMTsEEG2102_08",
        "PMT 2102-10-091": "PMTsEEG2102_10",
        "PMT 2102-12-091": "PMTsEEG2102_12",
        "PMT 2102-14-091": "PMTsEEG2102_14",
        "PMT 2102-16-091": "PMTsEEG2102_16",
        "PMT 2102-08-094": "PMTsEEG2102_08",  # same, just with permanent stylet
        "PMT 2102-10-094": "PMTsEEG2102_10",
        "PMT 2102-12-094": "PMTsEEG2102_12",
        "PMT 2102-14-094": "PMTsEEG2102_14",
        "PMT 2102-16-094": "PMTsEEG2102_16",
    }

    for lead in electrode_names.keys():
        if lead == reco_electrode:
            electrode_name = electrode_names[lead]
            return electrode_name
    
    return False

def extract_electrode_coords(filepath: str, space: str, hemi_idx: int):
    """
    Loads a MATLAB MAT-file containing electrode reconstruction data and extracts
    the 'coords_mm' field based on the specified coordinate space.

    Args:
        filepath (str): The path to the input MAT-file (e.g., 'sub-BER001_desc-reconstruction.mat').
        space (str): The desired coordinate space. Must be 'native' or 'mni'.
        hemi_idx (int): The hemisphere index, 0 - RH, 1 - LH 

    Returns:
        np.ndarray or None: The extracted electrode coordinates (coords_mm) 
                            as a NumPy array, or None if extraction fails.
    """
    if space not in ['native', 'mni']:
        print(f"Error: 'space' argument must be 'native' or 'mni', received '{space}'.")
        return None

    if not os.path.exists(filepath):
        print(f"Error: File not found at path: {filepath}")
        return None

    try:
        # Load the MATLAB file. Squeeze_me=True helps flatten single-element arrays.
        #data = loadmat(filepath, squeeze_me=True)
        data = h5py.File(str(filepath), "r")
    except Exception as e:
        print(f"Error loading MAT-file: {e}")
        return None

    # The main reconstruction structure is typically stored under a key like 'reco'.
    # Adjust this key if your MAT file uses a different root variable name.
    reco_key = 'reco' 
    if reco_key not in data:
        print(f"Error: Root variable '{reco_key}' not found in the MAT-file.")
        return None
        
    reco = data[reco_key]

    #print(reco['scrf']['coords_mm'])
    
    # general info
    if 'electrode' in reco:
        electrode_reference = reco['electrode']['dbs']
        electrode_group = data[electrode_reference[hemi_idx,0]]
        
        # just check if it is dataset (then no electrode, otherwise it is a group)
        if isinstance(electrode_group, h5py.Dataset):
            print("No electrode on this side")
            return None,None,None,None
        
        ascii_codes = electrode_group['elmodel'][:,0]
    else:
        electrode_reference = reco['props']['elmodel'][hemi_idx,0]
        electrode_group = data[electrode_reference]
        ascii_codes = electrode_group[:,0]
        
    electrode_model = "".join(np.vectorize(chr)(ascii_codes))

    try:
        if space == 'native':
            # Priority 1: reco.scrf
            # MATLAB structs loaded by scipy often require indexing like [()] or [0, 0]
            # to access the nested fields, but we try direct key access first.
            
            # Check for 'scrf' (Screen/Coregistration Refined) reconstruction
            if 'scrf' in reco:
                print("Extracting coordinates from reco.scrf (native space with brain shift correction).")
                coords_reference = reco['scrf']['coords_mm']
                coords_group = data[coords_reference[hemi_idx,0]]
                
                marker_head_reference = reco['scrf']['markers']['head']
                marker_head_group = data[marker_head_reference[hemi_idx,0]]
                
                marker_y_reference = reco['scrf']['markers']['y']
                marker_y_group = data[marker_y_reference[hemi_idx,0]]
                
                return coords_group, marker_head_group, marker_y_group, electrode_model
            
            # Priority 2: reco.native
            elif 'native' in reco:
                print("Extracting coordinates from reco.native (native space).")
                coords_reference = reco['native']['coords_mm']
                coords_group = data[coords_reference[hemi_idx,0]]
                
                marker_head_reference = reco['native']['markers']['head']
                marker_head_group = data[marker_head_reference[hemi_idx,0]]
                
                marker_y_reference = reco['native']['markers']['y']
                marker_y_group = data[marker_y_reference[hemi_idx,0]]
                
                return coords_group, marker_head_group, marker_y_group, electrode_model
            else:
                print(f"Error: 'coords_mm' not found in reco.scrf or reco.native.")
                return None
                
        elif space == 'mni':
            # Extract from reco.mni for MNI space
            if 'mni' in reco:
                print("Extracting coordinates from reco.mni (MNI space).")
                coords_reference = reco['mni']['coords_mm']
                coords_group = data[coords_reference[hemi_idx,0]]
                
                marker_head_reference = reco['mni']['markers']['head']
                marker_head_group = data[marker_head_reference[hemi_idx,0]]
                
                marker_y_reference = reco['mni']['markers']['y']
                marker_y_group = data[marker_y_reference[hemi_idx,0]]
                
                return coords_group, marker_head_group, marker_y_group, electrode_model
            else:
                print(f"Error: 'coords_mm' not found in reco.mni.")
                return None
                
    except Exception as e:
        print(f"An unexpected error occurred during data extraction: {e}")
        return None


# # Define the file path (assuming the uploaded file is accessible)
# path2subject = '/home/forel/Documents/data/CologneStimFit/derivatives/leaddbs/sub-CIR01'  # perhaps, this is not really needed and can be set to JD
# path2segmask = '/home/forel/Documents/data/CologneStimFit/derivatives/leaddbs/sub-CIR01/stimulations/native/NBstim/segmask.nii'
# stimFolder = '/home/forel/Documents/data/CologneStimFit/derivatives/leaddbs/sub-CIR01/stimulations/native/NBstim'

# space = 'mni'

if __name__ == '__main__':

    # called from MATLAB
    # sys.argv[1] - full path to the json dict
    
    # load json
    with open(sys.argv[1], 'r') as fp:
        InputDict = json.load(fp)
    fp.close()

    if InputDict['StimMode'] == 'CC':
        current_controlled = 1
    else:
        current_controlled = 0

    contact_coords_rh = np.array([InputDict['contact_coords_rh_x'],InputDict['contact_coords_rh_y'],InputDict['contact_coords_rh_z']])
    contact_coords_lh = np.array([InputDict['contact_coords_lh_x'],InputDict['contact_coords_lh_y'],InputDict['contact_coords_lh_z']])
    
    # for unilateral implantations to make sure the pipeline does not break, simply copy the values from one side to another
    # they will not be actually used!
    # if contact_coords_rh == None:
    #     contact_coords_rh, marker_head_coords_rh, marker_y_coords_rh, el_type_rh = (contact_coords_lh, marker_head_coords_lh, marker_y_coords_lh, el_type_lh)
    #     el_type = el_type_lh
    # elif contact_coords_lh == None:
    #     contact_coords_lh, marker_head_coords_lh, marker_y_coords_lh, el_type_lh = contact_coords_rh, marker_head_coords_rh, marker_y_coords_rh, el_type_rh
    #     el_type = el_type_rh
    
    # load the template settings
    if contact_coords_rh.shape[1] == 8:
        oss_settings_template = 'oss-dbs_parameters_8contacts.mat'
    elif contact_coords_rh.shape[1] == 4:
        oss_settings_template = 'oss-dbs_parameters_4contacts.mat'
    else:
        print("TBA for monocontact electrodes")
        raise SystemExit
        
        
    ### update lead_settings to utilize already implemented functionality
    from lead_settings_class import LeadSettings
    lead_settings = LeadSettings(oss_settings_template)
    
    lead_settings.add_pat_fold(InputDict['path2subject'])
    lead_settings.add_mri_name(InputDict['path2segmask'])   # run synthSeg
    lead_settings.add_dti_name('no dti')
    
    lead_settings.add_stim_set_mode(0.0)   # 
    lead_settings.add_phi_vec(np.zeros((2,contact_coords_rh.shape[1])))  # initialize
    
    if InputDict['space'] == 'native':
        lead_settings.add_y_mark_nat(np.array([InputDict['marker_y_coords_rh'],InputDict['marker_y_coords_lh']]))
        lead_settings.add_head_nat(np.array([InputDict['marker_head_coords_rh'],InputDict['marker_head_coords_lh']]))
    else:
        lead_settings.add_y_mark_mni(np.array([InputDict['marker_y_coords_rh'],InputDict['marker_y_coords_lh']]))
        lead_settings.add_head_mni(np.array([InputDict['marker_head_coords_rh'],InputDict['marker_head_coords_lh']]))
    
    lead_settings.add_imp_coord(np.array((contact_coords_rh[:,0],contact_coords_lh[:,0])))
    if contact_coords_rh.shape[1] == 1:
        # get tail
        print("TBA for monocontact electrodes")
        raise SystemExit
    elif 'DIXI' in InputDict['el_type']:
        lead_settings.add_sec_coord(np.array((contact_coords_rh[:,3],contact_coords_lh[:,3])))
    else:
        lead_settings.add_sec_coord(np.array((contact_coords_rh[:,-1],contact_coords_lh[:,-1])))
    
    # explicit update to circumvent h5py rigidity
    contactLocation_reference = lead_settings._settings["contactLocation"][0,:][0]
    lead_settings._settings[contactLocation_reference][...] = contact_coords_rh[:,:]
    
    contactLocation_reference = lead_settings._settings["contactLocation"][1,:][0]
    lead_settings._settings[contactLocation_reference][...] = contact_coords_lh[:,:]
    
    ### Update the default OSS-DBS dictionary and run simulations
    from custom_dict import oss_params
    
    # iterate over contacts changing Phi and StimCenter
    N_contacts = lead_settings.get_phi_vec().shape[1]
    for hemi_idx in range(2):
        
        contact_locations = lead_settings.get_cntct_loc(hemi_idx)
        first_contact = np.array(
            [
                contact_locations[0][0],
                contact_locations[1][0],
                contact_locations[2][0],
            ]
        )
        last_contact = np.array(
            [
                contact_locations[0][-1],
                contact_locations[1][-1],
                contact_locations[2][-1],
            ]
        )
        actual_span = np.linalg.norm(last_contact - first_contact)
        
        # add segmented rings!
        for cnt_i in range(N_contacts):
    
            hemi_sim_folder = os.path.join(InputDict['stimFolder'],'OSS_sim_files' + HEMI_SUFFIX[hemi_idx])
            if os.path.exists(hemi_sim_folder) and os.path.isdir(hemi_sim_folder):
                shutil.rmtree(hemi_sim_folder)
            os.mkdir(hemi_sim_folder)
            
            phi_vec = np.zeros((2,contact_coords_rh.shape[1]))
            phi_vec[hemi_idx,cnt_i] = AMPLITUDE
            lead_settings.add_phi_vec(phi_vec)
            stim_protocol = phi_vec[hemi_idx,:]
                    
            # grounding parameters
            if current_controlled:
                # switch from mA to A
                grounded_current = 0.001 * np.round(np.sum(stim_protocol), 9)  # could be 0
                case_grounding = True # if no current is actually grounded, this is used to reference voltages 
            else:
                grounded_current = 0.0 # not relevant
                print("Voltage-controlled mode is used, case grounding is ON")
                case_grounding = True
            
            lead_settings.add_stim_center(np.array((contact_coords_rh[:,cnt_i],contact_coords_lh[:,cnt_i])))
    
            if hemi_idx == 0:
                contact_coords = contact_coords_rh
            else:
                contact_coords = contact_coords_lh
                
                
            
            first_last_contact_coords = [list(contact_coords[:,0]),list(contact_coords[:,-1])]
            
            oss_electrode = check_electrode_availability(InputDict['el_type'])  
            elec_params = default_electrode_parameters[oss_electrode]
            
            # get tip position from the head and tail markers
            unit_directions, tip_pos, specs_array_length = lead_settings.get_tip_position(
                oss_electrode, hemi_idx
            )
            grid_center, grid_resolution = lead_settings.get_grid_parameters(
                oss_electrode, hemi_idx, unit_directions, specs_array_length
            )
            
            stretched_parameters = lead_settings.stretch_electrode(oss_electrode, hemi_idx)
            
            oss_params["BrainRegion"]["Center"]["x[mm]"] = lead_settings.get_imp_coord()[hemi_idx, 0]
            oss_params["BrainRegion"]["Center"]["y[mm]"] = lead_settings.get_imp_coord()[hemi_idx, 1]
            oss_params["BrainRegion"]["Center"]["z[mm]"] = lead_settings.get_imp_coord()[hemi_idx, 2]
            
            oss_params["BrainRegion"]["Dimension"]["x[mm]"] = 50.0 + np.abs(unit_directions[hemi_idx, 0]) * actual_span * 2.0
            oss_params["BrainRegion"]["Dimension"]["y[mm]"] = 50.0 + np.abs(unit_directions[hemi_idx, 1]) * actual_span * 2.0
            oss_params["BrainRegion"]["Dimension"]["z[mm]"] = 50.0 + np.abs(unit_directions[hemi_idx, 2]) * actual_span * 2.0
            
            oss_params["Electrodes"][0]["Name"] = oss_electrode + 'Custom'
            oss_params["Electrodes"][0]["CustomParameters"] = stretched_parameters
            
            oss_params["Electrodes"][0]["Rotation[Degrees]"] = lead_settings.get_rot_z(hemi_idx)
            oss_params["Electrodes"][0]["Direction"]["x[mm]"] = unit_directions[hemi_idx, 0]
            oss_params["Electrodes"][0]["Direction"]["y[mm]"] = unit_directions[hemi_idx, 1]
            oss_params["Electrodes"][0]["Direction"]["z[mm]"] = unit_directions[hemi_idx, 2]
            
            oss_params["Electrodes"][0]["TipPosition"]["x[mm]"] = tip_pos[hemi_idx, 0]
            oss_params["Electrodes"][0]["TipPosition"]["y[mm]"] = tip_pos[hemi_idx, 1]
            oss_params["Electrodes"][0]["TipPosition"]["z[mm]"] = tip_pos[hemi_idx, 2]
            
            oss_params["PointModel"]["Lattice"]["Active"] = not (bool(lead_settings.get_calc_axon_act()))
            oss_params["PointModel"]["Lattice"]["Center"]["x[mm]"] = grid_center[0]
            oss_params["PointModel"]["Lattice"]["Center"]["y[mm]"] = grid_center[1]
            oss_params["PointModel"]["Lattice"]["Center"]["z[mm]"] = grid_center[2]
            oss_params["PointModel"]["Lattice"]["Shape"] = {"x": 71, "y": 71, "z": 71}
            oss_params["PointModel"]["Lattice"]["PointDistance[mm]"] = grid_resolution
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