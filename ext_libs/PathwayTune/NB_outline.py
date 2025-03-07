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

SIDE_SUFFIX = ['rh','lh']

def determine_el_type(el_model):

    ''' Checks electrode configuration based on the model '''

    concentric4 = ['PINS_Medical_L303','PINS_Medical_L302','PINS_Medical_L301', 'Medtronic_3389','Medtronic_3387','Medtronic_3391','St._Jude ActiveTip_(6142-6145)', 'St._Jude_ActiveTip_(6146-6149)']
    concentric8 = ['Boston_Scientific_Vercise']
    segmented8 = ['St._Jude_Directed_6172_(short)', 'St._Jude_Directed_6180','St._Jude_Directed_6173_(long)','Boston_Scientific_Vercise_Directed', 'Medtronic_B33005', 'Medtronic B33015']
    if el_model in concentric4:
        return 'concentric4'
    elif el_model in concentric8:
        return 'concentric8'
    elif el_model in segmented8:
        return 'segmented8'
    else:
        print('The electrode model was not recognized')
        return 'Not classified electrode model'

def launch_weight_optimizer(stim_folder, netblend_dict, fixed_symptom_weights, side, approx_pathways):

    ''' Launch ANN-based optimization that will find optimal current protocol I while adjusting weights(!):
        Global_score = W1 * DS1(I) + W2 * DS2(I) + ... , where
        W1, W2 - symptom weights, including fixed and adjusted
        DS1, DS2 - distances in pathway activation space from target profiles of symptom 1 and 2
        to the activation profile for I '''

    ## Default algorithm parameters
    #netblend_dict['similiarity_metric'] = 'Canberra'  # or Bray-Curtis, Euclidean, etc
    #netblend_dict['optim_alg'] = 'Dual Annealing'  # or PSO
    #netblend_dict['num_iterations_ANN'] = 100  # number of ANN iterations to optimize current at the given electrode position


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

    # # iterate over each previously approved profile
    # with open(os.path.join(stim_folder, netblend_dict['ActivProfileDict']), 'r') as fp:
    #     target_profiles = json.load(fp)
    # fp.close()

    if netblend_dict['optim_alg'] == 'Dual Annealing':

        # optimize in respect to all symptoms with some weights W fixed
        # solve min(W*distance(main_symptoms) + (1-W)*sum(distance(others)))
        # weights_others ~ 1/distance(others) - this will be stored in a separate file and read in

        # ANN based
        from tensorflow.keras.models import load_model
        approx_model = load_model(os.path.join(stim_folder,'NB_' + str(SIDE_SUFFIX[side]), 'ANN_approved_model'))

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
    stim_dir = sys.argv[1]
    side = int(sys.argv[2])

    # load parameters from .json folder generated in previous steps
    with open(os.path.join(stim_dir,'netblend_dict.json'), 'r') as fp:
        netblend_dict = json.load(fp)
    fp.close()
    netblend_dict = netblend_dict['netblendict']

    # load Fixed_symptoms
    with open(os.path.join(stim_dir,'Fixed_symptoms.json'), 'r') as fp:
        fixed_symptom_weights_dict = json.load(fp)
    fp.close()
    fixed_symptom_weights_dict = fixed_symptom_weights_dict['fixed_symptom_weights']

    # load ANN approximated pathways
    with open(os.path.join(stim_dir,'NB_' + str(SIDE_SUFFIX[side]),'ANN_abs_errors.json'), 'r') as fp:
        pathways_errors_dict = json.load(fp)
    fp.close()
    # IMPORTANT: the pathways' order is preserved as they were processed in ANN!
    approx_pathways = list(pathways_errors_dict.keys())

    # clean-up
    if os.path.isfile(os.path.join(stim_dir,'NB_' + str(SIDE_SUFFIX[side]),'Estim_weights_and_total_score.csv')):
        os.remove(os.path.join(stim_dir,'NB_' + str(SIDE_SUFFIX[side]),'Estim_weights_and_total_score.csv'))
    if os.path.isfile(os.path.join(stim_dir,'NB_' + str(SIDE_SUFFIX[side]),'All_iters_estim_weights_and_total_score.csv')):
        os.remove(os.path.join(stim_dir,'NB_' + str(SIDE_SUFFIX[side]),'All_iters_estim_weights_and_total_score.csv'))

    launch_weight_optimizer(stim_dir, netblend_dict, fixed_symptom_weights_dict, side, approx_pathways)
