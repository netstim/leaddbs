
'''
    By K. Butenko
    Functions for PathwayTune (see description in the headers)
'''

import pandas
import numpy as np
import os
from scipy.spatial.distance import canberra, cityblock, euclidean, braycurtis, cosine
import json
import copy

def get_symptom_distances(activation_profile, Target_profiles, Soft_SE_thresh, fixed_symptom_weights, approx_pathways, side, score_symptom_metric='Canberra'):

    ''' Compute distances in pathway activation space from Target_profiles of symptoms and soft-side effects
        to the activation_profile '''

    # here we can merge target profiles for symptoms and threshold profiles for soft side-effects
    Target_profiles_and_SE = copy.deepcopy(Target_profiles)
    Target_profiles_and_SE.update(Soft_SE_thresh)

    N_symptoms_side = 0
    for key in Target_profiles_and_SE:
        if side == 0 and "_rh" in key:
            N_symptoms_side += 1
        elif side == 1 and "_lh" in key:
            N_symptoms_side += 1

    symp_distances = np.zeros(N_symptoms_side, float)
    symp_inx = 0
    sum_symp_nonfixed = 0
    symptom_list = []

    for key in Target_profiles_and_SE:
        if side == 0 and not ("_rh" in key):
            continue
        elif side == 1 and not ("_lh" in key):
            continue

        target_rates = []
        predicted_rates = []
        weights_for_pathways = []

        activ_target_profile = list(Target_profiles_and_SE[key].keys())

        for i in range(len(activ_target_profile)):

            target_rates.append(Target_profiles_and_SE[key][activ_target_profile[i]][0])
            weights_for_pathways.append(Target_profiles_and_SE[key][activ_target_profile[i]][2])

            if activ_target_profile[i] in approx_pathways:

                inx = approx_pathways.index(activ_target_profile[i])

                # if the activation is below the threshold, assign the threshold (so that the distance is 0)
                if key in Soft_SE_thresh and Target_profiles_and_SE[key][activ_target_profile[i]][0] > activation_profile[inx]:
                    predicted_rates.append(Target_profiles_and_SE[key][activ_target_profile[i]][0])
                else:
                    predicted_rates.append(activation_profile[inx])

            else:  # if not a part of the approx model, assign the threshold (so that the distance is 0)
                predicted_rates.append(Target_profiles_and_SE[key][activ_target_profile[i]][0])
                #print("Percent activation was not found for pathway ", activ_target_profile[i], "assigning null distance")

        # within the symptom, weights_for_pathways should sum up to 1
        # based on the acceptance of symptom-tract val equality
        weights_for_pathways = np.array(weights_for_pathways) / sum(weights_for_pathways)

        # distance in less important pathways is less penalized WITHIN the symptom
        if score_symptom_metric == 'Canberra':
            symp_distances[symp_inx] = canberra(predicted_rates, target_rates, w=weights_for_pathways) # / len(weights_for_pathways)
        elif score_symptom_metric == 'Manhattan':
            symp_distances[symp_inx] = cityblock(predicted_rates, target_rates, w=weights_for_pathways) # / len(weights_for_pathways)
        elif score_symptom_metric == 'Euclidean':
            symp_distances[symp_inx] = euclidean(predicted_rates, target_rates, w=weights_for_pathways) # / len(weights_for_pathways)
        elif score_symptom_metric == 'Cosine':
            symp_distances[symp_inx] = cosine(predicted_rates, target_rates, w=weights_for_pathways) # / len(weights_for_pathways)
        elif score_symptom_metric == 'Bray-Curtis':
            symp_distances[symp_inx] = braycurtis(predicted_rates, target_rates, w=weights_for_pathways) # / len(weights_for_pathways)
        else:
            print("Metric is not supported")
            raise SystemExit

        symptom_list.append(key)  # symptoms / soft-side effects for the given hemisphere (defined as side)

        # also return symptom distances for non-fixed (adjusted) weights
        if key not in fixed_symptom_weights:
            sum_symp_nonfixed += symp_distances[symp_inx]

        symp_inx += 1

    return sum_symp_nonfixed, symp_distances, symptom_list

def prepare_swarm(x, args_to_pass=[]):
    """Higher-level method to do forward_prop in the
    whole swarm.

    Inputs
    ------
    x: numpy.ndarray of shape (n_particles, dimensions)
        The swarm that will perform the search

    Returns
    -------
    numpy.ndarray of shape (n_particles, )
        The computed loss for each particle
    """
    n_particles = x.shape[0]
    j = [choose_weights_minimizer(x[i], args_to_pass[:]) for i in range(n_particles)]
    return np.array(j)

def choose_weights_minimizer(stim_vector, *args):
    ''' Compute estimated weights and global score based on symptom distances for the activation profile
        for stim_vector estimated with the approximation model '''

    mode_for_SSE = 'reverse'  # 'normal', 'reverse', 'exclude'

    # quick fix
    if len(args) == 1:
        approx_model, fixed_symptom_weights, score_symptom_metric, Target_profiles, Soft_SE_thresh, SE_thresh, side, approx_pathways = args[0]
    else:
        approx_model, fixed_symptom_weights, score_symptom_metric, Target_profiles, Soft_SE_thresh, SE_thresh, side, approx_pathways = args

    # IMPORTANT: activation profile should be composed by the same order of pathways as PD_profiles
    #activation_profile = get_profiles_ANN(stim_vector)

    # estimate the activation profile
    stim_array = np.reshape(np.array(stim_vector), (-1, len(stim_vector)))
    activation_profile = approx_model.predict(stim_array, verbose=0)
    activation_profile = activation_profile[0]   # get the actual array

    # first check the strict thresholds
    for key in SE_thresh:
        if side == 0 and not ("_rh" in key):
            continue
        elif side == 1 and not ("_lh" in key):
            continue

        activ_threshold_profile = list(SE_thresh[key].keys())
        for i in range(len(activ_threshold_profile)):

            if activ_threshold_profile[i] in approx_pathways:  # else we assume that the activation is 0 (below threshold)
                inx = approx_pathways.index(activ_threshold_profile[i])
                if np.any(activation_profile[inx] > SE_thresh[key][activ_threshold_profile[i]][0]):  # [0] - activation rate, [1] - weight of the pathway
                    #print('Side-effect threshold for this pathway activation was exceeded, discarding the solution')


                    # we need to create dummy output here
                    iter_estim_weights_and_total_score = (len(Target_profiles) + len(Soft_SE_thresh) + 1) * [0.0]
                    with open(
                            os.environ['STIMDIR'] + '/NB_' + str(side) + '/All_iters_estim_weights_and_total_score.csv',
                            'a') as f_handle:
                        np.savetxt(f_handle, np.vstack((iter_estim_weights_and_total_score)).T)

                    return 1000000.0
            #else:
            #    print("Pathway ", activ_threshold_profile[i], " was not included to the approximation model")


    # now we compute symptom distances (including soft-side effects)
    sum_symp_nonfixed, symp_distances, symptoms_list = get_symptom_distances(activation_profile, Target_profiles, Soft_SE_thresh, fixed_symptom_weights, approx_pathways, side, score_symptom_metric)

    # here we can have a soft-threshold analog. to the definition above
    # but for now we just store them all in symp_distances


    null_activation_profile = np.zeros(activation_profile.shape[0])
    [__, null_symptom_diff, symptoms_list] = get_symptom_distances(null_activation_profile, Target_profiles, Soft_SE_thresh, [], approx_pathways, side, score_symptom_metric)

    # also get symptom distances for 100% activation to estimate worst case scenario for soft-side effects
    max_activation_profile = 100.0 * np.ones(activation_profile.shape[0])
    [__, max_symptom_diff, symptoms_list] = get_symptom_distances(max_activation_profile, Target_profiles, Soft_SE_thresh, [], approx_pathways, side, score_symptom_metric)

    from RoutinesForResults import get_improvement_from_distance
    Impr_pred, estim_symp_improv_dict, symptom_labels_marked = get_improvement_from_distance(Target_profiles, Soft_SE_thresh, side, symp_distances, max_symptom_diff, null_symptom_diff, estim_weights_and_total_score=0, fixed_symptom_weights=0)
    # check improvement for non-fixed symptoms and remaining symtpom weight to be optimized
    symp_inx = 0
    Impr_non_fixed = 0.0
    Rest_weight = 1.0

    # do not estimate weights for soft side-effects
    Impr_pred_no_SSE = np.where(Impr_pred < 0.0, 0.001, Impr_pred)

    # here we can merge target profiles for symptoms and threshold profiles for soft side-effects
    Target_profiles_and_SE = copy.deepcopy(Target_profiles)
    Target_profiles_and_SE.update(Soft_SE_thresh)

    for symptom in Target_profiles_and_SE:

        if side == 0 and not ("_rh" in symptom):
            continue
        elif side == 1 and not ("_lh" in symptom):
            continue

        if symptom not in fixed_symptom_weights:
            if mode_for_SSE == 'exclude':  # no weighting estimation for soft side-effects
                Impr_non_fixed = Impr_non_fixed + Impr_pred_no_SSE[symp_inx, 0]
            elif mode_for_SSE == 'normal':  # the less worsening the higher the weight
                if Impr_pred[symp_inx, 0] < 0.0:  # soft side-effect, estimate improvement as 100% worsening  - predicted worsening
                    Impr_non_fixed = Impr_non_fixed + (1.0 + Impr_pred[symp_inx, 0])
                else:
                    Impr_non_fixed = Impr_non_fixed + Impr_pred[symp_inx, 0]
            elif mode_for_SSE == 'reverse':  # larger weight for larger worsening
                if Impr_pred[symp_inx, 0] < 0.0:  # just flip the sign to assign higher weight for higher worsening
                    Impr_non_fixed = Impr_non_fixed + abs(Impr_pred[symp_inx, 0])
                else:
                    Impr_non_fixed = Impr_non_fixed + Impr_pred[symp_inx, 0]
        else:
            Rest_weight = Rest_weight - fixed_symptom_weights[symptom]

        symp_inx += 1

    # estimated symptom weights (only for adjusted) and global score
    estim_symp_weights_norm = np.zeros(len(Target_profiles_and_SE), float)
    total_estim_weighted_score = 0.0

    # weight here and not above just for clarity (but can be combined)
    estim_symp_weight_norm_dict = {}


    symp_inx = 0
    for symptom in Target_profiles_and_SE:

        if side == 0 and not ("_rh" in symptom):
            continue
        elif side == 1 and not ("_lh" in symptom):
            continue

        if symptom not in fixed_symptom_weights:
            # this is the key part
            # do not make any val assumption at this point

            if estim_symp_improv_dict[symptom] < 0.0:
                if mode_for_SSE == 'exclude':
                    estim_symp_improv_dict[symptom] = 0.0  # nullify here to avoid weight estimation for soft SE
                elif mode_for_SSE == 'normal':
                    # soft side-effect, estimate improvement as 100% worsening  - predicted worsening
                    estim_symp_improv_dict[symptom] = 1.0 + estim_symp_improv_dict[symptom]
                elif mode_for_SSE == 'reverse':  # just flip the sign to assign higher weight for higher worsening
                    estim_symp_improv_dict[symptom] = abs(estim_symp_improv_dict[symptom])

            estim_symp_weights_norm[symp_inx] = Rest_weight * estim_symp_improv_dict[symptom] / Impr_non_fixed

            #estim_symp_weights_norm[symp_inx] = Rest_weight * (1 - symp_distances[symp_inx] / sum_symp_nonfixed)
            total_estim_weighted_score += estim_symp_weights_norm[symp_inx] * symp_distances[symp_inx]
        else:
            estim_symp_weights_norm[symp_inx] = fixed_symptom_weights[symptom]
            total_estim_weighted_score += estim_symp_weights_norm[symp_inx] * symp_distances[symp_inx]

        estim_symp_weight_norm_dict[symptom] = estim_symp_weights_norm[symp_inx]

        symp_inx += 1


    # save estimated weights and total score to a file to trace progress
    if os.path.isfile(os.environ['STIMDIR'] + '/NB_' + str(side) + '/All_iters_estim_weights_and_total_score.csv'):
        iter_data = np.genfromtxt(
            os.environ['STIMDIR'] + '/NB_' + str(side) + '/All_iters_estim_weights_and_total_score.csv', dtype=float, delimiter=' ')
        if len(iter_data.shape) == 1:
            iter_now = 2
        else:
            iter_now = iter_data.shape[0] + 1

    iter_estim_weights_and_total_score = np.append(estim_symp_weights_norm, total_estim_weighted_score)
    with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/All_iters_estim_weights_and_total_score.csv', 'a') as f_handle:
        np.savetxt(f_handle, np.vstack((iter_estim_weights_and_total_score)).T)


    # check if the result improved. If yes, update the output files
    if os.path.isfile(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Estim_weights_and_total_score.csv'):
        estim_weights_and_total_score = np.genfromtxt(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Estim_weights_and_total_score.csv', dtype=float, delimiter=' ')

        if len(estim_weights_and_total_score.shape) == 1:
            estim_weights_and_total_score = np.array([estim_weights_and_total_score])

        if total_estim_weighted_score < estim_weights_and_total_score[-1,-1]:

            #print(estim_weights_and_total_score[-1,:])
            print("it: ", iter_now, ";  f_goal = ", total_estim_weighted_score)

            np.savetxt(os.environ['STIMDIR'] + '/NB_' + str(side) + '/NW_symptom_distances.csv', symp_distances, delimiter=" ")

            estim_weights_and_total_score = np.append(estim_symp_weights_norm, total_estim_weighted_score)
            with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Estim_weights_and_total_score.csv', 'a') as f_handle:
                np.savetxt(f_handle, np.vstack((estim_weights_and_total_score)).T)

            # save json
            if side == 0:
                with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Estim_weights_rh.json', 'w') as save_as_dict:
                    json.dump(estim_symp_weight_norm_dict, save_as_dict)
            else:
                with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Estim_weights_lh.json', 'w') as save_as_dict:
                    json.dump(estim_symp_weight_norm_dict, save_as_dict)

            #np.savetxt(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Estim_symp_weights.csv', estim_symp_weights, delimiter=" ")
    else:
        np.savetxt(os.environ['STIMDIR'] + '/NB_' + str(side) + '/NW_symptom_distances.csv', symp_distances, delimiter=" ")
        #np.savetxt(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Estim_symp_weights.csv', estim_symp_weights, delimiter=" ")

        estim_weights_and_total_score = np.append(estim_symp_weights_norm, total_estim_weighted_score)
        with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Estim_weights_and_total_score.csv', 'w') as f_handle:
            np.savetxt(f_handle, np.vstack((estim_weights_and_total_score)).T)

    #print(weighted_symp_scores, sum(weighted_symp_scores), sum(symp_distances), total_estim_weighted_score)
    #return np.sum(weighted_symp_scores)  # + np.sum(soft_SE_scores)

    return total_estim_weighted_score
