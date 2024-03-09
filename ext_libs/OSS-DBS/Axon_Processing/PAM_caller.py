import h5py
import numpy as np
import json
import os
import subprocess
import sys

def launch_PAM(neuron_folder, folder_to_save, points_h5_file, pathways_params_file, scaling, scaling_index=None):
    """
    Parameters
    ----------
    neuron_folder: str, path to folder where NEURON models stored
    folder_to_save: str, path to folder where results are stored. Lead-DBS expects <stim_folder>/Results_<hemis>
    points_h5_file: str, path to .h5 containing the time domain solution for the pathways (point model, usually oss_time_result.h5)
    pathways_params_file: str, path to .json containing parameters for the pathways (usually Allocated_axons_parameters.json)
    scaling: float, optional, scaling factor for the whole solution (different from scaling_vector)
    scaling_index: int, optional, index of the scaling factor or scaling vector

    """

    # load files
    with open(pathways_params_file, 'r') as fp:
        pathways_dict = json.load(fp)

    # get to the right NEURON folder and compile
    if pathways_dict['Axon_Model_Type'] == 'McNeal1976':

        os.chdir(neuron_folder + "/McNeal1976")
        with open(os.devnull, 'w') as FNULL:
            if sys.platform == 'win32':
                subprocess.call('mknrndll', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
            else:
                subprocess.call('nrnivmodl', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    elif pathways_dict['Axon_Model_Type'] == "MRG2002" or pathways_dict['Axon_Model_Type'] == "MRG2002_DS":

        os.chdir(neuron_folder)
        with open(os.devnull, 'w') as FNULL:
            subprocess.call('nocmodl axnode.mod', shell=True, stdout=FNULL,
                            stderr=subprocess.STDOUT)  # might not work with remote hard drives
        if sys.platform == 'win32':
            with open(os.devnull, 'w') as FNULL:
                subprocess.call('mknrndll', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        else:
            with open(os.devnull, 'w') as FNULL:
                subprocess.call('nrnivmodl', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    from Axon_files.neuron_simulation import NeuronStimulation

    # load solution
    hf = h5py.File(points_h5_file , 'r')
    pathways = list(hf.keys())
    pathways.remove('TimeSteps[s]')

    # signal parameters can be extracted from solution
    TimeSteps = np.array(hf['TimeSteps[s]'])
    signal_dict = {
     'time_step': np.round(1000.0 * (TimeSteps[1] - TimeSteps[0]), 6),   # in ms
     'scaling': scaling,       # from GUI
     'N_time_steps': TimeSteps.shape[0]  #
    }

    #scaling_vectors = [[1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0],[]]

    pathway_idx = 0
    for pathway_name in pathways:

        pathway_dataset = hf[pathway_name]

        pathway_dict = {
            'pathway_name': pathway_name,
            'Axon_Model_Type': pathways_dict['Axon_Model_Type'],
            'axon_diam': pathways_dict['axon_diams'][pathway_idx],
            'n_Ranvier': pathways_dict['n_Ranvier'][pathway_idx],
            'N_seeded_neurons': pathways_dict['N_seeded_neurons'][pathway_idx],
            'N_orig_neurons': pathways_dict['N_orig_neurons'][pathway_idx],
            'connectome_name': pathways_dict['connectome_name'],
        }

        pathwayNEURON = NeuronStimulation(pathway_dict, signal_dict, folder_to_save, None, scaling_index)
        pathwayNEURON.check_pathway_activation(pathway_dataset)

        pathway_idx += 1

    hf.close()

if __name__ == '__main__':
    """ Call to probe action potentials for a given time domain solution

    Parameters
    ----------
    neuron_folder: str, path to folder where NEURON models stored
    folder_to_save: str, path to folder where results are stored. Lead-DBS expects <stim_folder>/Results_<hemis>
    points_h5_file: str, path to .h5 containing the time domain solution for the pathways (point model)
    pathways_params_file: str, path to .json containing parameters for the pathways
    scaling: float, optional, scaling factor for the whole solution (different from scaling_vector)
    scaling_index: int, optional, index of the scaling factor or scaling vector

    """
    neuron_folder = sys.argv[1:][0]
    folder_to_save = sys.argv[1:][1]
    points_h5_file = sys.argv[1:][2]
    pathways_params_file = sys.argv[1:][3]
    if len(sys.argv[1:]) >= 5:
        scaling = float(sys.argv[1:][4])
    else:
        scaling = 1.0

    if len(sys.argv[1:]) >= 6:
        scaling_index = int(sys.argv[1:][5])
    else:
        scaling_index = None

    launch_PAM(neuron_folder, folder_to_save, points_h5_file, pathways_params_file, scaling, scaling_index)
