'''
    By K. Butenko
    Runs OSS-DBS to compute stimulating fields for sEEG
'''

import pandas as pd
import numpy as np
import sys
import os
import json
import subprocess
import re
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

def get_geom_definitions(contact_locations):
    
    """
    Determine the extent and the center of the comp. domain  

    Args:
        contact_locations: list of lists, 3D coordinates for the first and the last active contacts
    Returns:
        dict, dimensions of the brain approximation
        NumPy array of shape (3,1), coordinates of Lattice
        NumPy array of shape (3,1), unit vector of the lead trajectory (from distal to proximal)
        
    """
    
    # check the distance between first and last contact for brain approx. dimensions
    first_contact = np.array(
        [
            contact_locations[0][0],
            contact_locations[0][1],
            contact_locations[0][2],
        ]
    )
    last_contact = np.array(
        [
            contact_locations[-1][0],
            contact_locations[-1][1],
            contact_locations[-1][2],
        ]
    )
    actual_span = np.linalg.norm(last_contact - first_contact)
    
    if actual_span > 12.0:
        print("Large active contact span! The simulated domain will be scaled, but Lattice settings should be adjusted!.")
    
    unit_directions = (last_contact - first_contact) / actual_span
    unit_directions = unit_directions.flatten()
    
    Dimensions = {
                "x[mm]": 100.0
                + np.abs(unit_directions[0]) * actual_span * 2.0,
                "y[mm]": 100.0
                + np.abs(unit_directions[1]) * actual_span * 2.0,
                "z[mm]": 100.0
                + np.abs(unit_directions[2]) * actual_span * 2.0,
            }
    
    grid_center = (last_contact + first_contact) / 2
    
    return Dimensions, grid_center, unit_directions
    

if __name__ == '__main__':

    # called from MATLAB
    # sys.argv[1] - str, full path to the reconstruction file
    # sys.argv[2] - str, full path to the stimulation protocol
    # sys.argv[3] - 'CC' if current-controlled
    # sys.argv[4] - Electrode ID (as integer), optional     

    ''' input processing '''
    SEEG_recos = sys.argv[1]
    _,extension = os.path.splitext(SEEG_recos)
    if extension == '.tsv':
        # either we get reco from tsv (Clemens' format)
        SEEG_recos_df = pd.read_csv(SEEG_recos, sep='\t')
    elif extension == '.mat':
        print("Lead-DBS reconstruction files are currently not supported")
        raise SystemExit()
        # reconstruction (a la Garance)
        import h5py
        SEEG_recos_mat = h5py.File(str(SEEG_recos), "r")
        # create a pandas dataframe analogous to above
        
    # stimulations in any case from .csv
    # separate file for each electrode!
    SEEG_stim = sys.argv[2]
    SEEG_stim_df = pd.read_csv(SEEG_stim, sep=',')
    SEEG_stim_df = SEEG_stim_df.replace('None', np.nan)
    SEEG_stim_array = SEEG_stim_df.to_numpy()     
    
    if sys.argv[3] == 'CC':
        current_controlled = True
    else:
        current_controlled = False  

    Electrode_ID = int(sys.argv[4])
    if Electrode_ID != None:
        # if elecrode ID is specified, drop entries for the rest
        SEEG_recos_df = SEEG_recos_df[SEEG_recos_df['Electrode_ID'] == Electrode_ID]
    
    # some auto-definitions
    stim_folder = os.path.dirname(SEEG_recos)
    contacts = SEEG_stim_df.columns.tolist()

    # load default oss-dbs_parameters.json
    #with open('/home/forel/Documents/data/SEEG/code/oss-dbs_parameters.json', 'r') as fp:
    #    ossdbs_dict = json.load(fp)
        
    # we actually rebuild the electrode for each stimulation to address the bend issue
    for stim_i in range(SEEG_stim_df.shape[0]):
    
        # electrode type and ID
        electrode = SEEG_recos_df.loc[SEEG_recos_df['Electrode_ID'] == Electrode_ID, 'electrode']
        oss_electrode = check_electrode_availability(electrode.iloc[0])  
        if Electrode_ID == None:
            Electrode_ID = SEEG_recos_df['Electrode_ID'][stim_i]
        
        # ToDO: sanity check for the number of contacts in the stim. protocol file for this electrode model
        
        ''' Check which contacts are active and build the trajectory between them. If one contact is used, select the adjacent to it'''
        
        stim_protocol = np.zeros(SEEG_stim_array[stim_i,:].shape)
        for i,s in enumerate(SEEG_stim_array[stim_i,:]):
            if isinstance(s,str):
                # actual number
                stim_protocol[i] = float(s)
            else:
                # nan
                stim_protocol[i] = s
        
        cnt_active_idx = np.where(~np.isnan(stim_protocol))[0][:]
        used_contacts = [contacts[i] for i in cnt_active_idx]  
               
        cnts_first_last = []
        cnts_labels = []
        cnts_first_last.append(SEEG_recos_df[(SEEG_recos_df['name'] == used_contacts[0]) & (SEEG_recos_df['Electrode_ID'] == Electrode_ID)])        
        cnts_labels.append(used_contacts[0])
        
        flip = False
        if len(used_contacts) == 1:
            # if only one contact was used in stim, we need to find either the next or the previous contact
            # to build the trajectory
            
            # IMPORTANT: this is a hard assumption that contact labels start with 1!
            index_on_electrode = int(re.sub(r"\D", "",used_contacts[0]))
            if index_on_electrode == len(contacts): 
                # last contact on the electrode
                flip = True

                # check the previous contact coordinates from the reconstruction file                
                row_index = SEEG_recos_df[(SEEG_recos_df['name'] == used_contacts[0]) & (SEEG_recos_df['Electrode_ID'] == Electrode_ID)].index.min()
                previous_row_index = SEEG_recos_df.index[SEEG_recos_df.index.get_loc(row_index) - 1]
                if SEEG_recos_df.loc[previous_row_index,'Electrode_ID'] == Electrode_ID:
                    cnts_labels.append(SEEG_recos_df.loc[previous_row_index,'name'])
                    cnts_first_last.append(SEEG_recos_df[(SEEG_recos_df['name'] == cnts_labels[-1]) & (SEEG_recos_df['Electrode_ID'] == Electrode_ID)])
                    
            else:
                # check the next contact coordinates from the reconstruction file      
                row_index = SEEG_recos_df[(SEEG_recos_df['name'] == used_contacts[0]) & (SEEG_recos_df['Electrode_ID'] == Electrode_ID)].index.min()
                next_row_index = SEEG_recos_df.index[SEEG_recos_df.index.get_loc(row_index) + 1]
                if SEEG_recos_df.loc[next_row_index,'Electrode_ID'] == Electrode_ID:
                    cnts_labels.append(SEEG_recos_df.loc[next_row_index,'name'])
                    cnts_first_last.append(SEEG_recos_df[(SEEG_recos_df['name'] == cnts_labels[-1]) & (SEEG_recos_df['Electrode_ID'] == Electrode_ID)])
                
                #cnts_labels.append(contacts[next_cnt_inx])
                #cnts_first_last.append(SEEG_recos_df[(SEEG_recos_df['name'] == contacts[next_cnt_inx]) & (SEEG_recos_df['Electrode_ID'] == Electrode_ID)])
        else:
            cnts_first_last.append(SEEG_recos_df[(SEEG_recos_df['name'] == used_contacts[-1]) & (SEEG_recos_df['Electrode_ID'] == Electrode_ID)])
            cnts_labels.append(used_contacts[-1])
        
        # check contact coordinates
        for cnt_i in range(len(cnts_first_last)):
            if cnts_first_last[cnt_i].shape[0] > 1:
                print("Non-unique definition of the contact coordinates for ",cnts_labels[cnt_i])
                raise SystemExit()
            elif cnts_first_last[cnt_i].shape[0] == 0:
                print("Contact coordinates not found for ",cnts_labels[cnt_i])
                raise SystemExit()
        
        # these are two contact coordinates to define the trajectory
        contact_coords = [[SEEG_recos_df['x'][cnts_first_last[0].index], SEEG_recos_df['y'][cnts_first_last[0].index], SEEG_recos_df['z'][cnts_first_last[0].index]],[SEEG_recos_df['x'][cnts_first_last[1].index], SEEG_recos_df['y'][cnts_first_last[1].index], SEEG_recos_df['z'][cnts_first_last[1].index]]]
        
        '''  now get the implantation coordinates '''
        Dimensions,grid_center,unit_directions = get_geom_definitions(contact_coords)
        if flip:
            unit_directions = unit_directions * -1.0
 
        elec_params = default_electrode_parameters[oss_electrode]
        # first_contact = first_active_coords - active_index * (spacing) * unit_direction
        # this might be wrong if the first contact (active tip) has a different length
        imp_coords = np.array([contact_coords[0][0].values[0],contact_coords[0][1].values[0],contact_coords[0][2].values[0]]) - cnt_active_idx[0] * (elec_params.contact_length + elec_params.contact_spacing) * unit_directions

        # offset = from tip to the center of the first contact 
        offset = elec_params.get_center_first_contact() * 1.0
        tip_position = imp_coords - offset * unit_directions

        # grounding parameters
        if current_controlled:
            # switch from mA to A
            grounded_current = -0.001 * np.round(np.sum(stim_protocol[cnt_active_idx.astype(int)]), 9)  # could be 0
            case_grounding = True # if no current is actually grounded, this is used to reference voltages 
        else:
            grounded_current = 0.0 # not relevant
            print("Voltage-controlled mode is used, case grounding is ON")
            case_grounding = True

        custom_params = {

            "BrainRegion": {
                # center at the head marker
                "Center": {
                    "x[mm]":  grid_center[0][0],
                    "y[mm]":  grid_center[1][0],
                    "z[mm]":  grid_center[2][0],
                },
                "Dimension": {
                  "x[mm]": Dimensions['x[mm]'],
                  "y[mm]": Dimensions['y[mm]'],
                  "z[mm]": Dimensions['z[mm]']
                },
                "Shape": "Ellipsoid"
            },
            "Electrodes": [
                {
                  "Name": oss_electrode,
                  "Rotation[Degrees]": -179.18978179152586,  # not relevant
                  "Direction": {
                    "x[mm]": unit_directions[0],
                    "y[mm]": unit_directions[1],
                    "z[mm]": unit_directions[2]
                  },
                  "TipPosition": {
                    "x[mm]": tip_position[0],
                    "y[mm]": tip_position[1],
                    "z[mm]": tip_position[2]
                  },
                  "EncapsulationLayer": {
                    "Thickness[mm]": 0.0,
                    "Material": "White matter",
                    "DielectricModel": "ColeCole4",
                    "DielectricParameters": None,
                    "MaxMeshSize": 1000000.0
                  }
                }
              ],

            "Surfaces": [
                {
                    "Name": "BrainSurface",
                    "Active": case_grounding,
                    "Current[A]": grounded_current,
                    "Voltage[V]": 0.0,
                }
            ],
            "MaterialDistribution": {
                "MRIPath": os.path.join(stim_folder,'segmask.nii'),
                "DiffusionTensorActive": False,
                "DTIPath": '',
            },
            "DielectricModel": {
                "Type": "ColeCole4",
                "CustomParameters": None
            },
            "Mesh": {
              "LoadMesh": False,
              "LoadPath": "",
              "MeshingHypothesis": {
                "Type": "Default",
                "MaxMeshSize": 1000000.0,
                "MeshSizeFilename": ""
              },
              "HPRefinement": {
                "Active": False,
                "Levels": 2,
                "Factor": 0.125
              },
              "AdaptiveMeshRefinement": {
                "Active": False,
                "MaxIterations": 10,
                "ErrorTolerance": 0.1
              },
              "MaterialRefinementSteps": 1,
              "MeshSize": {
                "Edges": {},
                "Faces": {},
                "Volumes": {}
              },
              "SaveMesh": False,
              "SavePath": "mesh"
            },
            "EQSMode": False,
            "FEMOrder": 2,
            "ComputeImpedance": False,
            "StimulationSignal": {
                "Type": "Multisine",
                "ListOfFrequencies": [
                    10000
                ],
                "Frequency[Hz]": 130.0,
                "PulseWidth[us]": 60.0,   # not relevant for stim volumes
                "CounterPulseWidth[us]": 0.0,
                "InterPulseWidth[us]": 0.0,
                "SpectrumMode": "FullSpectrum",
                "CutoffFrequency": 100000000.0,
                "CurrentControlled": current_controlled
            },
            "Solver": {
              "Type": "CG",
              "Preconditioner": "local",
              "PreconditionerKwargs": {},
              "MaximumSteps": 1000,
              "Precision": 1e-10
            },
            "PointModel": {
                "Lattice": {
                    "Active": True,
                    "Center": {
                        "x[mm]": grid_center[0][0],
                        "y[mm]": grid_center[1][0],
                        "z[mm]": grid_center[2][0]
                    },
                    "Shape": {"x": 51, "y": 51, "z": 51},
                    "Direction": {
                        "x[mm]": 0,
                        "y[mm]": 0,
                        "z[mm]": 1,
                    },
                    "PointDistance[mm]": 0.4,
                    "CollapseVTA": True,  # questionable
                }
            },
            "OutputPath": os.path.join(os.path.dirname(SEEG_recos),'Results_VTA_E' + str(Electrode_ID) + '_protocol_' + str(stim_i)),
            "SaveImpedance": False,
            "ExportVTK": True,
            "TemplateSpace": False,
            "ModelSide": 0,
            "CalcAxonActivation": False,
            "ActivationThresholdVTA[V-per-m]": 200.0,
            "FailFlag": 'rh',
            "OutOfCore": False
        }
            
        # contact dictionary is specified separately
        
        # only explicity model contacts between the first and the last active
        if len(used_contacts) == 1:
            cnt_range = range(cnt_active_idx[0],cnt_active_idx[0]+1)
            print(cnt_range)
        else:
            cnt_range = range(cnt_active_idx[0],cnt_active_idx[-1]+1)
            
        # cntct_dicts is a list of the contacts that will go into the dict
        # for this electrode
        cntct_dicts = np.empty(len(cnt_range), dtype=object)
        cntcts_made = 0
        
        for i in cnt_range:
            # all (truly) non-active contacts are floating with 0A
            if np.isnan(stim_protocol[i]):
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
        
        custom_params['Electrodes'][0]['Contacts'] = cntct_dicts.tolist()

        # save the settings
        json_settings = json.dumps(custom_params, indent=2)
        with open(os.path.join(stim_folder,'oss-dbs_parameters.json'), "w") as outfile:
            outfile.write(json_settings)
            
        # run OSS-DBS
        with open(os.devnull, 'w') as FNULL: subprocess.call('ossdbs ' + os.path.join(stim_folder,'oss-dbs_parameters.json'),shell=True)
