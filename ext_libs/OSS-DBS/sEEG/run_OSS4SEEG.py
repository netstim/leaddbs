'''
    By K. Butenko
    Runs OSS-DBS to compute 'VTRs' for sEEG
'''

import pandas as pd
import numpy as np
import sys
import os
import json
import subprocess
import re

from run_OSS4SEEG_Stim_no_shift import check_electrode_availability, get_geom_definitions, extract_index
from ossdbs.electrodes.defaults import default_electrode_parameters

MAX_RF_CONTACT_INDEX = 2   # for single contact RFs, do not go above the 3rd contact for simulation to avoid extremelely large VCMs (but the lead location is adjusted to match the actual contact coordinate) 

if __name__ == '__main__':

    # called from MATLAB
    # sys.argv[1] - full path to the reconstruction file
    # sys.argv[2] - 'CC' if current-controlled
    # sys.argv[3] - Electrode ID (as integer), optional 

    ''' input processing '''
    SEEG_recos = sys.argv[1]
    _,extension = os.path.splitext(SEEG_recos)
    if extension == '.tsv' or extension == '.csv':
        # either we get reco from tsv (BIDS format)
        SEEG_recos_df = pd.read_csv(SEEG_recos)  # make sure that contact numbering in the ascending order
    elif extension == '.mat':
        print("Lead-DBS reconstruction files are currently not supported")
        raise SystemExit()
        # reconstruction (a la Garance)
        import h5py
        SEEG_recos_mat = h5py.File(str(SEEG_recos), "r")
        # create a pandas dataframe analogous to above

    # if sys.argv[2] == 'CC':
    #     current_controlled = True
    # else:
    #     current_controlled = False
    current_controlled = False   # always use VC for VTRs

    if len(sys.argv) > 3:
        Electrode_ID = int(sys.argv[3])
        SEEG_recos_df = SEEG_recos_df[SEEG_recos_df['group'] == Electrode_ID]
        Electrode_ID_given = True
    else:
        Electrode_ID_given = False

    # some auto-definitions
    stim_folder = os.path.dirname(SEEG_recos)
    Amplitude = [1.0]  # 1 V. The exact value is not important for VTRs
    contacts2simulate = SEEG_recos_df.name   # all reconstructed

    # we iterate over all contacts of all electrodes, but rebuilding the geometry everytime
    for cnt_i in range(len(contacts2simulate)):
    
        # electrode type and ID
        oss_electrode = check_electrode_availability(SEEG_recos_df['electrode'][cnt_i])  
        if not Electrode_ID_given:
            # definition from the reconstruction sheet
            Electrode_ID = SEEG_recos_df['group'][cnt_i]
    
        ''' determine two contacts (active and adjacent) to build the trajctory '''
        cnt_ID = contacts2simulate[cnt_i]  # actual label
        flip = False;
        
        # IMPORTANT: this is a hard assumption that contact labels start with 1!
        index_on_electrode = extract_index(cnt_ID) - 1  # integer index
        if index_on_electrode < 0:
            print("Numbering for contact labels is expected to start from 1!")
            raise SystemExit
            
        if cnt_i+1 == len(contacts2simulate):
            # last listed contact
            cnt_ID2 = contacts2simulate[cnt_i-1]
            flip = True
        else:
            cnt_ID2 = contacts2simulate[cnt_i+1]
            index_on_electrode2 = extract_index(cnt_ID2) - 1
            if index_on_electrode2 - index_on_electrode != 1:
                # last contact, use previous to define a trajectory
                cnt_ID2 = contacts2simulate[cnt_i-1]
                flip = True
        
        # find indices of these contacts in SEEG_recos_df
        inx = SEEG_recos_df.index[SEEG_recos_df['name'] == cnt_ID]
        # check if the second contact coordinates are present
        if (SEEG_recos_df['name'].eq(cnt_ID2)).any():
            inx2 = SEEG_recos_df.index[SEEG_recos_df['name'] == cnt_ID2]
        else:
            print("The second contact was not found, cannot reconstruct the trajectory!")         
            raise SystemExit
        
        # these are two contact coordinates to define the trajectory
        contact_coords = [[SEEG_recos_df['x'][inx].values[0], SEEG_recos_df['y'][inx].values[0], SEEG_recos_df['z'][inx].values[0]],[SEEG_recos_df['x'][inx2].values[0], SEEG_recos_df['y'][inx2].values[0], SEEG_recos_df['z'][inx2].values[0]]]
        
        '''  now get the implantation coordinates '''
        Dimensions,grid_center,unit_directions = get_geom_definitions(contact_coords)
        if flip:
            # flip if the second contact is distal
            unit_directions = unit_directions * -1.0
 
        elec_params = default_electrode_parameters[oss_electrode]
        if index_on_electrode > MAX_RF_CONTACT_INDEX:
            # "shift" the lead to avoid large VCMs
            simulated_contact_index = MAX_RF_CONTACT_INDEX
        else:
            simulated_contact_index = index_on_electrode
            
        # this might be wrong if the first contact (active tip) has a different length
        if 'BF' in oss_electrode:
            # first contact spacing is different
            imp_coords = np.array([contact_coords[0][0],contact_coords[0][1],contact_coords[0][2]]) - (simulated_contact_index * (elec_params.contact_length + elec_params.contact_spacing) + bool(simulated_contact_index) * (elec_params.first_contact_spacing-elec_params.contact_spacing)) * unit_directions  
        else:
            imp_coords = np.array([contact_coords[0][0],contact_coords[0][1],contact_coords[0][2]]) - simulated_contact_index * (elec_params.contact_length + elec_params.contact_spacing) * unit_directions         
        
        # offset = from tip to the center of the first contact 
        offset = elec_params.get_center_first_contact() * 1.0
        tip_position = imp_coords - offset * unit_directions
            
        for amp in Amplitude:
        
            custom_params = {
    
                "BrainRegion": {
                    # center at the head marker
                    "Center": {
                        "x[mm]":  grid_center[0],
                        "y[mm]":  grid_center[1],
                        "z[mm]":  grid_center[2],
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
                      "Contacts": [         # we need only one contact! Move the forth one for now!
                        {
                          "Contact_ID": simulated_contact_index + 1,
                          "Active": True,
                          "Current[A]": False,
                          "Voltage[V]": amp,
                          "Floating": False,
                          "SurfaceImpedance[Ohmm]": {
                            "real": 0.0,
                            "imag": 0.0
                          },
                          "MaxMeshSize": 1000000.0,
                          "MaxMeshSizeEdge": 0.08545132017764238  # hardcoded for now
                        }
                      ],                      
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
                        "Active": True,
                        "Current[A]": False,
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
                    "PulseTopWidth[us]": 0.0,
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
                            "x[mm]": grid_center[0],
                            "y[mm]": grid_center[1],
                            "z[mm]": grid_center[2]
                        },
                        "Shape": {"x": 71, "y": 71, "z": 71},
                        "Direction": {
                            "x[mm]": 0,
                            "y[mm]": 0,
                            "z[mm]": 1,
                        },
                        "PointDistance[mm]": 0.4,
                        "CollapseVTA": True,  # questionable
                    }
                },
                "OutputPath": os.path.join(os.path.dirname(SEEG_recos),'Results_VTR_' + str(Electrode_ID) + '_1V_' + cnt_ID),
                "SaveImpedance": False,
                "ExportVTK": True,
                "TemplateSpace": False,
                "ModelSide": 0,
                "CalcAxonActivation": False,
                "ActivationThresholdVTA[V-per-m]": 200.0,
                "FailFlag": 'rh',
                "OutOfCore": False,
                "TruncateAfterActivePartRatio": None
            }
    
            # save the settings
            json_settings = json.dumps(custom_params, indent=2)
            with open(os.path.join(stim_folder,'oss-dbs_parameters.json'), "w") as outfile:
                outfile.write(json_settings)
    
            # run OSS-DBS
            with open(os.devnull, 'w') as FNULL: subprocess.call('ossdbs ' + os.path.join(stim_folder,'oss-dbs_parameters.json'),shell=True)
            
