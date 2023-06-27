'''
    By K. Butenko
    Functions for PathwayTune (see description in the headers)
'''
import numpy as np
import json
import os
import h5py
from scipy.stats import qmc
import csv
import pandas
import sys

def determine_el_type(el_model):

    ''' Checks electrode configuration based on the model '''

    concentric4 = ['PINS_Medical_L303','PINS_Medical_L302','PINS_Medical_L301', 'Medtronic_3389','Medtronic_3387','Medtronic_3391','St._Jude ActiveTip_(6142-6145)', 'St._Jude_ActiveTip_(6146-6149)']
    concentric8 = ['Boston_Scientific_Vercise']
    segmented8 = ['St._Jude_Directed_6172_(short)', 'St._Jude_Directed_6180','St._Jude_Directed_6173_(long)','Boston_Scientific_Vercise_Directed']
    if el_model in concentric4:
        return 'concentric4'
    elif el_model in concentric8:
        return 'concentric8'
    elif el_model in segmented8:
        return 'segmented8'
    else:
        print('The electrode model was not recognized')
        return 'Not classified electrode model'


def load_AP_from_OSSDBS(side, inters_as_stim=False):

    ''' Load activation profile and pathway labels from OSS-DBS results
        if inters_as_stim == True, fibers inside encapsulation and/or outside of the domain will be treated as activated '''

    # get all pathways that survived Kuncel(!) pre-filtering and original(!) number of fibers
    # this function will work only for a proper Lead-DBS import (connectome folder, oss-dbs_parameters.mat)
    from Pathways_Stats import get_simulated_pathways
    Pathways, axons_in_path = get_simulated_pathways(side)

    if side == 0:
        res_folder = os.environ['STIMDIR'] + '/' + 'Results_rh/'
    else:
        res_folder = os.environ['STIMDIR'] + '/' + 'Results_lh/'

    hf = h5py.File(res_folder + 'Summary_status.h5', 'r')
    lst = list(hf.keys())
    if len(lst) != len(Pathways):
        print("Number of pathways passed to OSS-DBS and returned is not matching")
        raise SystemExit

    Perc_Activation = np.zeros(len(lst), float)

    for pathway_index in range(len(lst)):
        pathway_name = str(lst[pathway_index])
        if Pathways[pathway_index] != pathway_name:
            print("Mismatch in pathway ordering, check oss-dbs_parameters.mat and OSS-DBS log")
            raise SystemExit
        else:
            Activation_in_pathway = hf.get(lst[pathway_index])
            Activation_in_pathway = np.array(Activation_in_pathway)

            if inters_as_stim == True:
                Perc_Activation[pathway_index] = (Activation_in_pathway[0] + Activation_in_pathway[2] + Activation_in_pathway[4]) / axons_in_path[pathway_index]
            else:
                Perc_Activation[pathway_index] = Activation_in_pathway[0] / axons_in_path[pathway_index]

    hf.close()

    return Perc_Activation, Pathways



def launch_weight_optimizer(activation_profile_dict, netblend_dict, fixed_symptom_weights, side, approx_pathways):

    ''' Launch ANN-based optimization that will find optimal current protocol I while adjusting weights(!):
        Global_score = W1 * DS1(I) + W2 * DS2(I) + ... , where
        W1, W2 - symptom weights, including fixed and adjusted
        DS1, DS2 - distances in pathway activation space from target profiles of symptom 1 and 2
        to the activation profile for I '''

    ## Default algorithm parameters
    #netblend_dict['similiarity_metric'] = 'Canberra'  # or Bray-Curtis, Euclidean, etc
    #netblend_dict['optim_alg'] = 'Dual Annealing'  # or PSO
    #netblend_dict['num_iterations_ANN'] = 100  # number of ANN iterations to optimize current at the given electrode position

    # better regenerate them
    from Improvement4Protocol import create_NB_dictionaries
    profile_dict, Soft_SE_dict, SE_dict = create_NB_dictionaries(side, activation_profile_dict,
                                                                 disease='spontaneous human combustion')

    # # load previously approved symptom-specific profiles
    # with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/profile_dict.json', 'r') as fp:
    #     profile_dict = json.load(fp)
    # fp.close()
    #
    # with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Soft_SE_dict.json', 'r') as fp:
    #     Soft_SE_dict = json.load(fp)
    # fp.close()
    #
    # with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/SE_dict.json', 'r') as fp:
    #     SE_dict = json.load(fp)
    # fp.close()


    if netblend_dict['optim_alg'] == 'Dual Annealing':

        # optimize in respect to all symptoms with some weights W fixed
        # solve min(W*distance(main_symptoms) + (1-W)*sum(distance(others)))
        # weights_others ~ 1/distance(others) - this will be stored in a separate file and read in

        # ANN based
        from tensorflow.keras.models import load_model
        approx_model = load_model(os.environ['STIMDIR'] + '/NB_' + str(side) + '/ANN_approved_model')

        # this part is with ANN
        from scipy.optimize import dual_annealing
        from Optim_strategies import choose_weights_minimizer
        res = dual_annealing(choose_weights_minimizer,
                             bounds=list(zip(netblend_dict['min_bound_per_contact'], netblend_dict['max_bound_per_contact'])),
                             args=[approx_model, fixed_symptom_weights, netblend_dict['similarity_metric'], profile_dict, Soft_SE_dict, SE_dict, side, approx_pathways], maxfun=netblend_dict['num_iterations_ANN'], seed=42, visit=2.62,
                             no_local_search=True)

        optimized_current = res.x  # in A!

    if netblend_dict['optim_alg'] == 'PSO':

        # ANN based
        from tensorflow.keras.models import load_model
        approx_model = load_model(os.environ['STIMDIR'] + '/NB_' + str(side) + '/ANN_approved_model')

        from pyswarms.single.global_best import GlobalBestPSO
        bounds = (netblend_dict['min_bound_per_contact'], netblend_dict['max_bound_per_contact'])
        options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9}
        optimizer = GlobalBestPSO(n_particles=20, dimensions=len(netblend_dict['min_bound_per_contact']), options=options, bounds=bounds)

        from Optim_strategies import prepare_swarm
        cost, optimized_current = optimizer.optimize(prepare_swarm, iters=netblend_dict['num_iterations_ANN'], args_to_pass=[approx_model, fixed_symptom_weights, netblend_dict['similarity_metric'], profile_dict, Soft_SE_dict, SE_dict, side, approx_pathways])

    # ============================================ Plotting ===========================================================#

    optim_stim = np.reshape(np.array(optimized_current), (-1, len(optimized_current)))
    activation_profile = approx_model.predict(optim_stim, verbose=0)
    activation_profile = activation_profile[0]  # get the actual array

    estim_weights_and_total_score = np.genfromtxt(
        os.environ['STIMDIR'] + '/NB_' + str(side) + '/Estim_weights_and_total_score.csv', dtype=float,
        delimiter=' ')

    non_weighted_symptom_dist = np.genfromtxt(
        os.environ['STIMDIR'] + '/NB_' + str(side) + '/NW_symptom_distances.csv', dtype=float,
        delimiter=' ')

    from RoutinesForResults import get_activation_prediction
    get_activation_prediction(optimized_current, activation_profile, approx_pathways, non_weighted_symptom_dist, profile_dict,
                              Soft_SE_dict, side, plot_results=True, score_symptom_metric='Canberra',
                              estim_weights_and_total_score=estim_weights_and_total_score, fixed_symptom_weights=fixed_symptom_weights)



if __name__ == '__main__':

    # called from MATLAB
    # sys.argv[1] - stim folder
    # sys.argv[2] - Activation profile dictionary
    # sys.argv[3] - side (0 - right hemisphere)

    os.environ['STIMDIR'] = sys.argv[1]
    activation_profile_dict = sys.argv[2]
    side = int(sys.argv[3])


    # load parameters from .json folder generated in previous steps
    with open(os.environ['STIMDIR'] + '/netblend_dict_file.json', 'r') as fp:
        netblend_dict = json.load(fp)
    fp.close()
    netblend_dict = netblend_dict['netblendict']

    # load Fixed_symptoms
    with open(os.environ['STIMDIR'] + '/Fixed_symptoms.json', 'r') as fp:
        fixed_symptom_weights_dict = json.load(fp)
    fp.close()
    fixed_symptom_weights_dict = fixed_symptom_weights_dict['fixed_symptom_weights']

    # load ANN approximated pathways
    with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/ANN_abs_errors.json', 'r') as fp:
        pathways_errors_dict = json.load(fp)
    fp.close()
    # IMPORTANT: the pathways' order is preserved as they were processed in ANN!
    approx_pathways = list(pathways_errors_dict.keys())

    # clean-up
    if os.path.isfile(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Estim_weights_and_total_score.csv'):
        os.remove(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Estim_weights_and_total_score.csv')
    if os.path.isfile(os.environ['STIMDIR'] + '/NB_' + str(side) + '/All_iters_estim_weights_and_total_score.csv'):
        os.remove(os.environ['STIMDIR'] + '/NB_' + str(side) + '/All_iters_estim_weights_and_total_score.csv')

    launch_weight_optimizer(activation_profile_dict, netblend_dict, fixed_symptom_weights_dict, side, approx_pathways)