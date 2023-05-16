'''
    By K. Butenko
    Routines to post-process VAT-based pathway recruitment
'''



import os
import json
import numpy as np

def remove_failed_protocols(Currents_all, ActivationResults_all):

    # remove protocols for which the simulator failed, they were marked with -1.0
    rows_idx = np.where(np.any(ActivationResults_all==-1.0,axis=1))[0]
    ActivationResults = np.delete(ActivationResults_all, rows_idx, axis=0)
    # adjust the iteration vector
    ActivationResults[:,0] = np.arange(ActivationResults.shape[0])

    Currents = np.delete(Currents_all, rows_idx, axis=0)

    return Currents, ActivationResults

def get_VAT_pathways(side):

    # load simulated pathways from json
    with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/VAT_pathways.json', 'r') as fp:
        vat_paths_dict = json.load(fp)
    fp.close()
    vat_paths_dict = vat_paths_dict['vat_paths_dict']
    Pathways_fo = []
    axons_in_path_fo = []
    sort_inx = []
    for key in vat_paths_dict:
        Pathways_fo.append(key)
        axons_in_path_fo.append(vat_paths_dict[key][0])
        sort_inx.append(vat_paths_dict[key][1] - 1)  # MATLAB to Python

    # re-sort as in Activations_over_StimSets
    # this is weird
    Pathways = []
    axons_in_path = []
    for idx in range(len(sort_inx)):
        Pathways.append(Pathways_fo[sort_inx.index(idx)])
        axons_in_path.append(axons_in_path_fo[sort_inx.index(idx)])

    return Pathways, axons_in_path