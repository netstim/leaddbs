"""
@author: Konstantin Butenko

# This script reads input files from Lead-DBS and updates the default dictionary of GUI.
"""

# import tables
import h5py  # works better
import numpy as np
import subprocess
# import pickle
import json
import os
import sys


# updates default dictionary
def get_input_from_LeadDBS(settings_location, index_side, cluster_run=False):  # 0 - rhs, 1 - lhs

    index_side = int(index_side)

    # these parameters will be imported from Lead-DBS
    input_dict = {
        'MRI_data_name': "name_MRI.nii.gz",  # segmented MRI data
        'DTI_data_name': "name_DTI.nii.gz",  # scaled tensor data
        'CSF_index': 0.0,  # index of the tissue in the segmented MRI data
        'WM_index': 0.0,
        'GM_index': 0.0,
        'default_material': 'GM',  # GM, WM or CSF
        'Electrode_type': 'St_Jude6148',
        'Implantation_coordinate_X': -10000000000.0,  # in mm in the MRI data space
        'Implantation_coordinate_Y': -10000000000.0,
        'Implantation_coordinate_Z': -10000000000.0,
        'Second_coordinate_X': -10000000000.0,
        'Second_coordinate_Y': -10000000000.0,
        'Second_coordinate_Z': -10000000000.0,
        'Rotation_Z': 0.0,  # rotation around the lead axis in degrees
        'current_control': 0,  # 0 - VC, 1 - CC
        'Pulse_amp': [None, None, None, None],
        # Signal vector: give an amplitude. If CC, 0.0 refers to 0 V (ground), other numbers are in A. None is for floating potentials
        'Activation_threshold_VTA': 0.0,  # threshold for Astrom VTA (V/mm), compute using the Lead-DBS function.
        'VTA_approx': 0,
        'external_grounding': False,
        'VTA_from_E': 1,
        'VTA_from_divE': 0,
        'Stim_side': 0,  # 0 - rh, 1 - lh
        'Neuron_model_array_prepared': 0,
        'stretch': 1.0,
        'number_of_processors': 0,
        'Approximating_Dimensions': [80.0, 80.0, 80.0],
        'Aprox_geometry_center': [0.0, 0.0, 0.0],
        'el_order': 2,
        'StimSets': 0,
    }

    print("\nInput from ", settings_location, "\n")

    path = os.path.normpath(settings_location)
    path.split(os.sep)

    try:
        patient_folder = path.split(os.sep)[-5]
        print('Patient: ', path.split(os.sep)[-5])
    except:
        patient_folder = "John Doe"  # stim folder was passed, patient_folder not relevant

    patient_folder = path.split(os.sep)[-5]

    file = h5py.File(str(settings_location), 'r')

    # if file.root.settings.current_control[0][0]!=file.root.settings.current_control[0][1]:
    if all(~np.isnan(file['settings']['current_control'][0])):
        if file['settings']['current_control'][0][0] != file['settings']['current_control'][0][-1]:
            print("Simultaneous use of VC and CC is not allowed for safety reasons!")
            raise SystemExit

    # IMPORTANT: it is actually not native but scrf!
    head_native = file['settings']['headNative'][:, index_side]
    y = file['settings']['yMarkerNative'][:, index_side] - head_native
    y_postop = y / np.linalg.norm(y)
    phi = np.arctan2(-y_postop[0], y_postop[1])
    input_dict["Rotation_Z"] = phi * 180.0 / np.pi

    input_dict['Stim_side'] = index_side
    #Pulse_amp = file.root.settings.Phi_vector[:,index_side]
    Pulse_amp = file['settings']['Phi_vector'][:, index_side]
    if file['settings']['current_control'][0][0] == 1 or file['settings']['current_control'][0][-1] == 1:
        input_dict['current_control'] = 1
        Pulse_amp = Pulse_amp * 0.001  # because Lead-DBS uses mA as the input

    Pulse_amp = list(Pulse_amp)
    import math
    for i in range(len(Pulse_amp)):
        if math.isnan(Pulse_amp[i]):
            Pulse_amp[i] = None
        # Fix of VC random grounding bug for Lead-DBS stim settings
        elif Pulse_amp[i] == 0.0 and input_dict['current_control'] == 0:
            Pulse_amp[i] = None  # this is allowed only for this case since we can't set contact to 0 V in Lead-DBS
            print("swapped 0V contact to None in VC")

    StimSets = int(file['settings']['stimSetMode'][0][0])
    input_dict['StimSets'] = StimSets

    if all(v is None for v in Pulse_amp) and StimSets == 0:
        print("No stimulation defined for this hemisphere")
        return -1

    input_dict['Pulse_amp'] = Pulse_amp

    # convert ascii from Matlab struct to a string
    # array_ascii=file.root.settings.MRI_data_name[:]
    array_ascii = file['settings']['MRI_data_name'][:]
    list_ascii = []
    for i in range(array_ascii.shape[0]):
        list_ascii.append(array_ascii[i][0])
    # list_ascii = map(lambda s: s.strip(), list_ascii)
    name_split = ''.join(chr(i) for i in list_ascii)
    input_dict['MRI_data_name'] = name_split.rsplit(os.sep, 1)[-1]
    stim_folder = os.path.dirname(str(settings_location))
    print("Stim folder:", stim_folder)

    # array_ascii=file.root.settings.DTI_data_name[:]
    array_ascii = file['settings']['DTI_data_name'][:]
    list_ascii = []
    if array_ascii[0] == 0:
        input_dict['DTI_data_name'] = ''
    else:
        for i in range(array_ascii.shape[0]):
            list_ascii.append(array_ascii[i][0])
        # list_ascii = map(lambda s: s.strip(), list_ascii)
        input_dict['DTI_data_name'] = ''.join(chr(i) for i in list_ascii)

    input_dict['CSF_index'] = file['settings']['CSF_index'][0][0]
    input_dict['WM_index'] = file['settings']['WM_index'][0][0]
    input_dict['GM_index'] = file['settings']['GM_index'][0][0]

    # input_dict['CSF_index']=file.root.settings.CSF_index[0][0]
    # input_dict['WM_index']=file.root.settings.WM_index[0][0]
    # input_dict['GM_index']=file.root.settings.GM_index[0][0]

    # array_ascii=file.root.settings.default_material[:]
    array_ascii = file['settings']['default_material'][:]
    list_ascii = []
    for i in range(array_ascii.shape[0]):
        list_ascii.append(array_ascii[i][0])
    # list_ascii = map(lambda s: s.strip(), list_ascii)
    default_material = ''.join(chr(i) for i in list_ascii)

    if default_material == 'GM':
        input_dict['default_material'] = 3
    elif default_material == 'WM':
        input_dict['default_material'] = 2
    elif default_material == 'CSF':
        input_dict['default_material'] = 1
    else:
        print("Unrecognized default material")

    # array_ascii=file.root.settings.Electrode_type[:]
    array_ascii = file['settings']['Electrode_type'][:]
    list_ascii = []
    for i in range(array_ascii.shape[0]):
        list_ascii.append(array_ascii[i][0])
    # list_ascii = map(lambda s: s.strip(), list_ascii)
    Electrode_type = ''.join(chr(i) for i in list_ascii)

    if Electrode_type == 'Medtronic 3389':
        input_dict['Electrode_type'] = "Medtronic3389"  # 1
        normal_array_length = 6.0  # center of the first to the center of the last
    elif Electrode_type == 'Medtronic 3387':
        normal_array_length = 9.0
        input_dict['Electrode_type'] = "Medtronic3387"  # 1
    elif Electrode_type == 'Medtronic 3391':
        input_dict['Electrode_type'] = "Medtronic3391"  # 1
        normal_array_length = 21.0
    elif Electrode_type == 'St. Jude Directed 6172 (short)' or Electrode_type == 'St. Jude Directed 6180':  # just different marker colors
        input_dict['Electrode_type'] = "St_Jude6180"  # 1
        normal_array_length = 6.0
    elif Electrode_type == 'St. Jude Directed 6173 (long)':
        normal_array_length = 9.0
        input_dict['Electrode_type'] = "St_Jude6173"  # 1
    elif Electrode_type == 'St. Jude ActiveTip (6142-6145)':  # just different tail lenghts, but it does not matter here
        input_dict['Electrode_type'] = "St_Jude6142"
        normal_array_length = 9.75
    elif Electrode_type == 'St. Jude ActiveTip (6146-6149)':  # just different tail lenghts, but it does not matter here
        input_dict['Electrode_type'] = "St_Jude6148"  # 1
        normal_array_length = 6.75
    elif Electrode_type == 'Boston Scientific Vercise':
        input_dict['Electrode_type'] = "Boston_Scientific_Vercise"  # 1
        normal_array_length = 6.0  # still the 4th level. BEWARE: Stretch is defined between the 1st and the 4th, but applied between 1st and 8th!!!
        normal_array_length = 14.0
    elif Electrode_type == 'Boston Scientific Vercise Directed':
        input_dict['Electrode_type'] = "Boston_Scientific_Vercise_Cartesia"  # 1
        normal_array_length = 6.0
    elif Electrode_type == 'PINS Medical L301':
        input_dict['Electrode_type'] = "PINS_L301"
        normal_array_length = 6.0
    elif Electrode_type == 'PINS Medical L302':
        input_dict['Electrode_type'] = "PINS_L302"
        normal_array_length = 9.0
    elif Electrode_type == 'PINS Medical L303':
        input_dict['Electrode_type'] = "PINS_L303"
        normal_array_length = 18.0
    elif Electrode_type == 'ELAINE Rat Electrode':
        input_dict['Electrode_type'] = "AA_rodent_monopolar"
        normal_array_length = 1.0  # irrelevant for monopolar electrodes
    else:
        print(
            "The electrode is not yet implemented, but we will be happy to add it. Contact us via konstantin.butenko@uni-rostock.de")

    input_dict['Implantation_coordinate_X'], input_dict['Implantation_coordinate_Y'], input_dict[
        'Implantation_coordinate_Z'] = file['settings']['Implantation_coordinate'][:, index_side]
    input_dict['Second_coordinate_X'], input_dict['Second_coordinate_Y'], input_dict['Second_coordinate_Z'] = \
        file['settings']['Second_coordinate'][:, index_side]
    input_dict['Aprox_geometry_center'] = [input_dict['Implantation_coordinate_X'],
                                           input_dict['Implantation_coordinate_Y'],
                                           input_dict['Implantation_coordinate_Z']]

    el_array_length = np.sqrt((input_dict['Implantation_coordinate_X'] - input_dict['Second_coordinate_X']) ** 2 + (
                input_dict['Implantation_coordinate_Y'] - input_dict['Second_coordinate_Y']) ** 2 + (
                                          input_dict['Implantation_coordinate_Z'] - input_dict[
                                      'Second_coordinate_Z']) ** 2)
    stretch = el_array_length / normal_array_length

    Estimate_In_Template = file['settings']['Estimate_In_Template'][0][0]

    if abs(stretch - 1.0) < 0.01 or Electrode_type == 'ELAINE Rat Electrode':  # 1% tolerance or monopolar electrodes
        input_dict["stretch"] = 1.0
    else:
        input_dict["stretch"] = stretch

    input_dict['Activation_threshold_VTA'] = file['settings']['Activation_threshold_VTA'][0][
        0]  # threshold is the same for both hemispheres
    input_dict['external_grounding'] = bool(file['settings']['Case_grounding'][:, index_side][0])
    input_dict['Neuron_model_array_prepared'] = int(
        file['settings']['calcAxonActivation'][0][0])  # external model (e.g. from fiber tractograpy)

    if input_dict['Neuron_model_array_prepared'] != 1:
        input_dict['VTA_approx'] = 1  # for now we have only these two options
    else:
        with open(stim_folder + '/Allocated_axons_parameters.json', 'r') as fp:
            neuron_dict = json.load(fp)
        input_dict.update(neuron_dict)

    input_dict['Stim_side'] = index_side
    input_dict['patient_folder'] = patient_folder

    with open(stim_folder + '/Lead_DBS_input.json', 'w') as save_as_dict:
        json.dump(input_dict, save_as_dict)

    interactive_mode = int(file['settings']['interactiveMode'][0][0])

    return stim_folder, index_side, interactive_mode, cluster_run, patient_folder, StimSets


if __name__ == '__main__':

    oss_dbs_folder = os.path.dirname(os.path.realpath(sys.argv[0]))
    os.chdir(oss_dbs_folder)

    print(sys.argv[1:])
    if len(sys.argv[1:]) > 2 and sys.argv[3] == '1':
        cluster_sim = int(sys.argv[3])
    else:
        cluster_sim = 0
    print("Cluster run:", bool(cluster_sim))

    path_to_patient, side, interactive_mode, cluster_run, patient_folder, StimSets = get_input_from_LeadDBS(
        *sys.argv[1:3], cluster_sim)

    if os.environ.get('SINGULARITY_NAME'):  # not sure about this part
        os.environ['PATIENTDIR'] = path_to_patient
        os.environ['TMPDIR'] = sys.argv[5]
    else:
        os.environ['PATIENTDIR'] = '/opt/Patient'  # Use fixed mount path for docker

    process_1 = -1
    if cluster_run == 0:
        if path_to_patient != -1:
            if sys.platform == 'linux':
                if os.environ.get('SINGULARITY_NAME'):
                    output = subprocess.run(
                        ['xterm', '-e', 'python3', 'GUI_tree_files/AppUI.py', path_to_patient, str(side),
                         str(interactive_mode), str(patient_folder), str(StimSets)])
                else:
                    output = subprocess.run(
                        ['xterm', '-e', 'python3', 'GUI_tree_files/AppUI.py', path_to_patient, str(side),
                         str(interactive_mode), str(patient_folder), str(StimSets)])
            elif sys.platform == 'darwin':
                open_terminal = 'tell application "Terminal" to do script "cd \'' + oss_dbs_folder + '\';'
                open_gui = ' '.join(
                    ['python3 GUI_tree_files/AppUI.py', path_to_patient, str(side), str(interactive_mode),
                     str(patient_folder), str(StimSets), ';exit"'])
                output = subprocess.run(['osascript', '-e', open_terminal + open_gui])
            elif sys.platform == 'win32':
                output = subprocess.run(
                    ['start', 'cmd', '/c', 'python', 'GUI_tree_files/AppUI.py', path_to_patient, str(side),
                     str(interactive_mode), str(patient_folder), str(StimSets)], shell=True)
            else:
                print("The system's OS does not support OSS-DBS")
                raise SystemExit
    else:

        dir_code = os.path.dirname(os.getcwd())  # OSS-DBS folder to be mount, NOT OSS_platform folder

        import getpass

        # /fast/users/ might need to be changed depending on the cluster 
        user_pseudo_home_mount = '/fast/users/' + getpass.getuser() + ':/root'
        user_pseudo_work = '/fast/work/users/' + getpass.getuser()

        env_vars = 'PATIENTDIR=' + os.environ['PATIENTDIR'] + ',TMPDIR=' + os.environ['TMPDIR'] + ','

        output = subprocess.run(['singularity', 'run', '--cleanenv',
                                 '-B', dir_code + ':/opt/OSS-DBS',
                                 '-B', path_to_patient + ':/opt/Patient',
                                 '-B', sys.argv[5] + ':/opt/temp',
                                 '-B', user_pseudo_work,
                                 '-B', user_pseudo_home_mount,
                                 '--home', '/root',
                                 '--env', env_vars, sys.argv[6], 'python3',
                                 '/opt/OSS-DBS/OSS_platform/Launcher_OSS_lite.py', str(cluster_run)])

        ## When launching MATLAB in Singularity?

        # if os.environ.get('SINGULARITY_NAME'):
        #    output = subprocess.run(['python3', 'Launcher_OSS_cluster.py'])
        # else:
        #    print("Only Singularity cluster solution is currently supported")
        #    raise SystemExit
