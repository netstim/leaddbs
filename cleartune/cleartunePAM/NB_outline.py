# this script is a sketch of a network blending algorithm

import numpy as np
import pandas
import json
import os
import h5py
from scipy.stats import qmc
import csv

def determine_el_type(el_model):

    concentric4 = ['PINS_Medical_L303','PINS_Medical_L302','PINS_Medical_L301', 'Medtronic_3389','Medtronic_3387','Medtronic_3391','St._Jude ActiveTip_(6142-6145)', 'St._Jude_ActiveTip_(6146-6149)']
    concentric8 = ['Boston_Scientific_Vercise']
    segmented8 = ['St._Jude_Directed_6172_(short)', 'St._Jude_Directed_6180','St._Jude_Directed_6173_(long)','Boston_Scientific_Vercise_Directed',]
    if el_model in concentric4:
        return 'concentric4'
    elif el_model in concentric8:
        return 'concentric8'
    elif el_model in segmented8:
        return 'segmented8'
    else:
        print('The electrode model was not recognized')
        return 'Not classified electrode model'


def load_AP_from_LeadDBS(side, inters_as_stim=False):
    """ if inters_as_stim == True, fibers inside encapsulation and/or outside of the domain will be treated as activated"""

    # the function will work only for a proper Lead-DBS import (connectome folder, oss-dbs_parameters.mat)
    from Pathways_Stats import get_simulated_pathways
    Pathways, axons_in_path = get_simulated_pathways(side)
    # we will take all pathways in this case

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




def launch_weight_optimizer(el_model, fixed_symptom_weights, side, approx_pathways):

    # Algorithm parameters
    similiarity_metric = 'Canberra'   # or Bray-Curtis, Euclidean, etc
    optim_alg = 'dual_annealing'      # or PSO
    num_iterations = 100 #  number of ANN iterations to optimize current at the given electrode position
    I_lim_conc = [-2.5, 2.5]  # current limits for concentric contacts
    I_lim_segm = [-1.0, 1.0]  # current limits for segmented contacts

    el_type = determine_el_type(el_model)  # 'concentric4', 'concentric8', 'segmented8', etc

    # load previously approved symptom-specific profiles
    with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/profile_dict.json', 'r') as fp:
        profile_dict = json.load(fp)
    fp.close()

    with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Soft_SE_dict.json', 'r') as fp:
        Soft_SE_dict = json.load(fp)
    fp.close()

    with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/SE_dict.json', 'r') as fp:
        SE_dict = json.load(fp)
    fp.close()

    if el_type == 'concentric4':
        min_bound_per_contact = 4 * [I_lim_conc[0]]
        max_bound_per_contact = 4 * [I_lim_conc[1]]
    elif el_type == 'concentric8':
        min_bound_per_contact = 8 * [I_lim_conc[0]]
        max_bound_per_contact = 8 * [I_lim_conc[1]]
    elif el_type == 'segmented8':
        min_bound_per_contact = [I_lim_conc[0], I_lim_segm[0], I_lim_segm[0], I_lim_segm[0], I_lim_segm[0],
                                 I_lim_segm[0], I_lim_segm[0], I_lim_conc[0]]
        max_bound_per_contact = [I_lim_conc[1], I_lim_segm[1], I_lim_segm[1], I_lim_segm[1], I_lim_segm[1],
                                 I_lim_segm[1], I_lim_segm[1], I_lim_conc[1]]
    else:
        print('Electrode type is not recognized')
        raise SystemExit

    if optim_alg == 'dual_annealing':

        # optimize in respect to all symptoms with some weights W fixed
        # solve min(W*distance(main_symptoms) + (1-W)*sum(distance(others)))
        # weights_others ~ 1/distance(others) - this will be stored in a separate file and read in

        from tensorflow.keras.models import load_model
        approx_model = load_model(os.environ['STIMDIR'] + '/NB_' + str(side) + '/ANN_approved_model')

        # this part is with ANN
        from scipy.optimize import dual_annealing
        from Optim_strategies import choose_weights_minimizer
        res = dual_annealing(choose_weights_minimizer,
                             bounds=list(zip(min_bound_per_contact, max_bound_per_contact)),
                             args=[approx_model, fixed_symptom_weights, similiarity_metric, profile_dict, Soft_SE_dict, SE_dict, side, approx_pathways], maxfun=num_iterations, seed=42, visit=2.62,
                             no_local_search=True)





    # ============================================ Ploting ============================================================#

    # should be moved to separate functions

    # estimated norm weights, predicted improvement

    optim_stim = np.reshape(np.array(res.x), (-1, len(res.x)))
    activation_profile = approx_model.predict(optim_stim, verbose=0)
    activation_profile = activation_profile[0]  # get the actual array

    estim_weights_and_total_score = np.genfromtxt(
        os.environ['STIMDIR'] + '/NB_' + str(side) + '/Estim_weights_and_total_score.csv', dtype=float,
        delimiter=' ')

    non_weighted_symptom_dist = np.genfromtxt(
        os.environ['STIMDIR'] + '/NB_' + str(side) + '/NW_symptom_distances.csv', dtype=float,
        delimiter=' ')

    from RoutinesForResults import get_activation_prediction
    get_activation_prediction(res.x, activation_profile, approx_pathways, non_weighted_symptom_dist, profile_dict,
                              Soft_SE_dict, side, score_symptom_metric='Canberra',
                              estim_weights_and_total_score=estim_weights_and_total_score, fixed_symptom_weights=fixed_symptom_weights)


    # # # estimate the improvement: canberra distance at null activation vs the optimized
    # # null_protocol = len(min_bound_per_contact) * [0.0]
    # # null_activation_profile = approx_model.predict(np.reshape(np.array(null_protocol), (-1, len(null_protocol))), verbose=0)
    # # null_activation_profile = null_activation_profile[0]  # get the actual array
    #
    # # or just assign directly. Note, the former might be more reasonable if damaged neurons are considered as activated!
    # null_activation_profile = np.zeros(activation_profile.shape[0])
    #
    # # here we can merge target profiles for symptoms and threshold profiles for soft side-effects
    # profile_dict.update(Soft_SE_dict)
    #
    # from Optim_strategies import get_symptom_distances
    # [__, null_symptom_diff] = get_symptom_distances(null_activation_profile, profile_dict, Soft_SE_dict, fixed_symptom_weights, approx_pathways, side, similiarity_metric)
    #
    #
    # # also get symptom distances for 100% activation to estimate worst case scenario for soft-side effects
    # max_activation_profile = 100.0 * np.ones(activation_profile.shape[0])
    # [__, max_symptom_diff] = get_symptom_distances(max_activation_profile, profile_dict, Soft_SE_dict, fixed_symptom_weights, approx_pathways, side,
    #                                                 similiarity_metric)
    #
    #
    #
    # # plot predicted improvement and corresponding weights in a combined bar plot
    # Impr_pred = np.zeros((len(profile_dict),2), float) # in the second dimension, we will store the estimated weight
    # symp_inx = 0
    #
    # symptom_labels_marked = []
    #
    # for symptom in profile_dict:
    #
    #     if symptom in Soft_SE_dict: # we assume there are no soft side-effects at null protocol
    #
    #         # IMPORTANT: for soft-side effect we calculate predicted worsening to the maximum (100% activation worsening
    #         if max_symptom_diff[symp_inx] == 0.0:
    #             Impr_pred[symp_inx,0] = 0.0
    #         else:
    #             # Here the value is always negative, i.e. worsening
    #             Impr_pred[symp_inx,0] = (max_symptom_diff[symp_inx] - non_weighted_symptom_dist[symp_inx]) / max_symptom_diff[symp_inx] - 1.0
    #     else:
    #         Impr_pred[symp_inx,0] = (null_symptom_diff[symp_inx] - non_weighted_symptom_dist[symp_inx]) / null_symptom_diff[symp_inx]
    #
    #     # estimated weight for the symptom, the order was preserved (we always iterate over symptom dictionary)
    #     Impr_pred[symp_inx, 1] = estim_weights_and_total_score[-1,symp_inx]
    #
    #     symp_inx += 1
    #
    #
    #     # mark fixed symptoms for plotting
    #     if symptom in fixed_symptom_weights:
    #         symptom_labels_marked.append(symptom + " (fixed)")
    #     else:
    #         symptom_labels_marked.append(symptom)
    #
    #
    # import matplotlib.pyplot as plt
    # import seaborn as sns
    # sns.set()
    #
    # pos = np.arange(len(symptom_labels_marked))  # the x locations for the groups
    # pos_adjusted = pos# - 0.5
    #
    # fig, ax = plt.subplots(figsize=(12, 3))
    # width = 0.25
    # colors_opt = ['C0', 'C1']
    # #labels_opt = ['Activated', 'Not Activated', 'Encapsulation', 'CSF', 'Other']
    # # bars_status =
    # # for i in range(Impr_pred.shape[1]):
    # #     ax.bar(pos_adjusted - 0.5 + i * width, Impr_pred[:, i] , width,
    # #            color=colors_opt[i])
    #
    #
    # ax.bar(pos_adjusted - 0.5 * width, Impr_pred[:, 0] * 100.0, width,
    #        color=colors_opt[0])
    # ax.set_ylabel("Improvement / worsening, %", color='C0')
    # #ax.axhline(0)
    # #ax.set_ylim(-100,100)
    #
    # ax2 = ax.twinx()
    # ax2.bar(pos_adjusted + 0.5 * width, Impr_pred[:, 1], width,
    #        color=colors_opt[1])
    # ax2.grid(False)
    #
    # # autolimit to align axes (potentially buggy!)
    # if np.any(Impr_pred[:, 0] < 0.0):
    #     print('here')
    #     #lower_w_lim = np.min(Impr_pred[:, 0]) / np.max(Impr_pred[:, 0])
    #     lower_w_lim = np.min(Impr_pred[:, 0]) / 1.0
    #
    #     ax.set_ylim(np.min(Impr_pred[:, 0])*100 , 100.0)
    #     ax2.set_ylim(lower_w_lim, 1)
    # else:
    #     ax.set_ylim(0, 100)
    #     ax2.set_ylim(0, 1)
    #
    #
    # print(Impr_pred)
    #
    # #ax2.set_ylim(0, 1)
    # #ax2.axhline(0)
    #
    # #ax.legend(loc='upper right', bbox_to_anchor=(0, 1.25, 1, 0), ncol=5, mode="expand", borderaxespad=0.)
    # ax.set_xticks(pos_adjusted)
    # ax.set_xticklabels(symptom_labels_marked)
    # ax2.set_ylabel("Suggested symptom weights", color='C1')
    # fig.tight_layout()
    # plt.savefig(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Symptom_profiles_' + str(side) + '.png',
    #             format='png',
    #             dpi=1000)
    #
    #
    #
    #  # activation profile for optimized res.x, current protocol,
    # fig, ax = plt.subplots(figsize=(12, 8))
    # width = 0.5
    #
    #
    # # only plot with non-zero activation
    # approx_pathways_activated = []
    # activation_nonzero = []
    # for i in range(activation_profile.shape[0]):
    #     if activation_profile[i] > 0.0:
    #         approx_pathways_activated.append(approx_pathways[i])
    #         activation_nonzero.append(activation_profile[i])
    # activation_nonzero = np.array(activation_nonzero)
    #
    # pos = np.arange(len(approx_pathways_activated))  # the x locations for the groups
    #
    # ax.bar(pos, activation_nonzero * 100.0, width,
    #        color=colors_opt[0])
    #
    # if len(min_bound_per_contact) == 4:
    #     textstr = '\n'.join((
    #         r'Optimized Currents (mA), Lead-DBS notation ',
    #         r'$k_{3}=%.2f$' % (res.x[3]),
    #         r'$k_{2}=%.2f$' % (res.x[2]),
    #         r'$k_{1}=%.2f$' % (res.x[1]),
    #         r'$k_{3}=%.2f$' % (res.x[0])))
    # elif len(min_bound_per_contact) == 8:
    #     textstr = '\n'.join((
    #         r'Optimized Currents (mA), Lead-DBS notation',
    #         r'$k_{3}=%.2f \ \ k_{7}=%.2f$' % (res.x[3], res.x[7]),
    #         r'$k_{2}=%.2f \ \ k_{6}=%.2f$' % (res.x[2], res.x[6]),
    #         r'$k_{1}=%.2f \ \ k_{5}=%.2f$' % (res.x[1], res.x[5]),
    #         r'$k_{0}=%.2f \ \ k_{4}=%.2f$' % (res.x[0], res.x[4])))
    # else:
    #     print("The electrode model was not recognized")
    #
    # #
    # # ax.legend(loc='upper right', bbox_to_anchor=(0, 1.25, 1, 0), ncol=5, mode="expand", borderaxespad=0.)
    # # ax.set_xticks(pos_adjusted)
    # ax.set_xticklabels(approx_pathways_activated)
    # ax.set_ylabel("Percent Activation, %")
    #
    # # these are matplotlib.patch.Patch properties
    # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    # # place a text box in upper left in axes coords
    # ax.text(0.0, 1.5, textstr, transform=ax.transAxes, fontsize=14,
    #         verticalalignment='top', bbox=props)
    #
    # ax.set_xticks(pos)
    # plt.xticks(rotation=45)
    # fig.tight_layout()
    # plt.savefig(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Activation_profile_' + str(side) + '.png',
    #             format='png',
    #             dpi=1000)



    # load Estim_symp_weights.csv to Lead-DBS
    # IMPORTANT: these are not predicted improvements, but the weights of importance for the symptom modulation
    # They should not sum up to 1.0 or 1.0 - main symptom, but should not exceed 1.0 - main symptom

    # after the weights are suggested, the actual optimization (minimization of all symptom scores with these weights)
    # will be conducted inside OSS-DBS (Launcher_OSS_NetBlend.py), res.x is the initial guess