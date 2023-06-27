'''
    By K. Butenko
    Estimate improvement for the given current protocol (either real or approximation model)
    and plot the results in the stimulation folder / NB_
'''

import numpy as np
import os
import json
import copy


def plot_results_with_weights(current_protocol, activation_profile, pathways, Impr_pred, symptom_labels_marked, side, NB_result=False):

    ''' call via get_activation_prediction() '''

    import matplotlib.pyplot as plt
    #import seaborn as sns
    #sns.set()

    #============= Plot predicted improvement and corresponding weights =============#
    pos = np.arange(len(symptom_labels_marked))  # the x locations for the groups
    pos_adjusted = pos  # - 0.5

    fig, ax = plt.subplots(figsize=(12, 8))
    width = 0.25
    colors_opt = ['C0', 'C1']  # different colors for predicted improvement and weights

    # first bars for the improvement
    ax.bar(pos_adjusted - 0.5 * width, Impr_pred[:, 0] * 100.0, width,
           color=colors_opt[0])
    ax.set_ylabel("Improvement / worsening, %", color='C0')
    # ax.axhline(0)
    # ax.set_ylim(-100,100)

    #seond bars for the weights
    ax2 = ax.twinx()
    ax2.bar(pos_adjusted + 0.5 * width, Impr_pred[:, 1], width,
            color=colors_opt[1])
    ax2.grid(False)

    # autolimit to align axes (potentially buggy!)
    if np.any(Impr_pred[:, 0] < 0.0):
        # lower_w_lim = np.min(Impr_pred[:, 0]) / np.max(Impr_pred[:, 0])
        lower_w_lim = np.min(Impr_pred[:, 0]) / 1.0

        ax.set_ylim(np.min(Impr_pred[:, 0]) * 100, 100.0)
        ax2.set_ylim(lower_w_lim, 1)
    else:
        ax.set_ylim(0, 100)
        ax2.set_ylim(0, 1)

    #print(Impr_pred)

    # ax2.set_ylim(0, 1)
    # ax2.axhline(0)

    # ax.legend(loc='upper right', bbox_to_anchor=(0, 1.25, 1, 0), ncol=5, mode="expand", borderaxespad=0.)
    ax2.set_ylabel("Suggested symptom weights", color='C1')
    ax.set_xticks(pos_adjusted)
    ax.set_xticklabels(symptom_labels_marked, rotation=45)
    fig.tight_layout()
    plt.savefig(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Symptom_profiles_' + str(side) + '.png',
                format='png',
                dpi=1000)


    # ============= Plot activation profile and current procotol (Lead-DBS notation) =============#

    fig, ax = plt.subplots(figsize=(12, 8))
    width = 0.5

    ## I think this is not necessary. We might have 0 activation pathways here if they had non zeros in training
    # if NB_result == True:
    #     # only plot with non-zero activation for Network Blending
    #     approx_pathways_activated = []
    #     activation_nonzero = []
    #     for i in range(activation_profile.shape[0]):
    #         if activation_profile[i] > 0.0:
    #             approx_pathways_activated.append(pathways[i])
    #             activation_nonzero.append(activation_profile[i])
    #     activation_nonzero = np.array(activation_nonzero)
    #
    #     pos = np.arange(len(approx_pathways_activated))  # the x locations for the groups
    #     ax.bar(pos, activation_nonzero * 100.0, width,
    #            color=colors_opt[0])
    # else:
    pos = np.arange(len(activation_profile))  # the x locations for the groups
    ax.bar(pos, activation_profile * 100.0, width,
           color=colors_opt[0])

    # convert to mA
    current_protocol = np.array(current_protocol) * 1000.0

    if len(current_protocol) == 4:
        textstr = '\n'.join((
            r'Optimized Currents (mA), Lead-DBS notation ',
            r'$k_{3}=%.2f$' % (current_protocol[3]),
            r'$k_{2}=%.2f$' % (current_protocol[2]),
            r'$k_{1}=%.2f$' % (current_protocol[1]),
            r'$k_{0}=%.2f$' % (current_protocol[0])))
    elif len(current_protocol) == 8:
        textstr = '\n'.join((
            r'Optimized Currents (mA), Lead-DBS notation',
            r'$k_{3}=%.2f \ \ k_{7}=%.2f$' % (current_protocol[3], current_protocol[7]),
            r'$k_{2}=%.2f \ \ k_{6}=%.2f$' % (current_protocol[2], current_protocol[6]),
            r'$k_{1}=%.2f \ \ k_{5}=%.2f$' % (current_protocol[1], current_protocol[5]),
            r'$k_{0}=%.2f \ \ k_{4}=%.2f$' % (current_protocol[0], current_protocol[4])))
    else:
        print("The electrode model was not recognized")

    #
    # ax.legend(loc='upper right', bbox_to_anchor=(0, 1.25, 1, 0), ncol=5, mode="expand", borderaxespad=0.)
    # ax.set_xticks(pos_adjusted)
    ax.set_xticklabels(pathways)
    ax.set_ylabel("Percent Activation, %")

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    # place a text box in upper left in axes coords
    ax.text(0.0, 1.5, textstr, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)

    ax.set_xticks(pos)
    plt.xticks(rotation=45)
    fig.tight_layout()
    plt.savefig(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Activation_profile_' + str(side) + '.png',
                format='png',
                dpi=1000)


def get_improvement_from_distance(profile_dict, Soft_SE_dict, side, symp_distances, max_symptom_diff, null_symptom_diff, estim_weights_and_total_score=0,fixed_symptom_weights=[]):

    # here we can merge target profiles for symptoms and threshold profiles for soft side-effects
    profile_dict_and_SE = copy.deepcopy(profile_dict)
    profile_dict_and_SE.update(Soft_SE_dict)

    # check how many symptoms / soft side-effects we have for that hemisphere
    N_symptoms_side = 0
    for key in profile_dict_and_SE:
        if side == 0 and "_rh" in key:
            N_symptoms_side += 1
        elif side == 1 and "_lh" in key:
            N_symptoms_side += 1
    Impr_pred = np.zeros((N_symptoms_side, 2), float)  # in the second dimension, we store the estimated weight

    # iterate over symptoms to estimate symptom improvement (from null activation)
    # or side-effect worsening (assuming the worst case at 100% activation)
    symptom_labels_marked = []
    symp_inx = 0
    estim_symp_improv_dict = {}

    for symptom in profile_dict_and_SE:

        if side == 0 and not ("_rh" in symptom):
            continue
        elif side == 1 and not ("_lh" in symptom):
            continue

        if symptom in Soft_SE_dict:  # we assume there are no soft side-effects at null protocol

            # IMPORTANT: for soft-side effect we calculate predicted worsening in comparison to the maximum worsening at 100% activation
            if max_symptom_diff[symp_inx] == 0.0:
                Impr_pred[symp_inx,0] = 0.0
            else:
                # Here the value is always negative, i.e. worsening
                Impr_pred[symp_inx,0] = (max_symptom_diff[symp_inx] - symp_distances[symp_inx]) / max_symptom_diff[symp_inx] - 1.0
        else:
            # we might have all pathways excluded for the symptom
            # in this case predict zero improvement
            if null_symptom_diff[symp_inx] == 0.0:
                Impr_pred[symp_inx, 0] = 0.0
            else:
                Impr_pred[symp_inx,0] = (null_symptom_diff[symp_inx] - symp_distances[symp_inx]) / null_symptom_diff[symp_inx]

        estim_symp_improv_dict[symptom] = Impr_pred[symp_inx,0]

        # add info for weights if Network Blending was conducted
        if np.any(estim_weights_and_total_score != 0):
            # estimated weight for the symptom, the order was preserved (we always iterate over the symptom dictionary)
            Impr_pred[symp_inx, 1] = estim_weights_and_total_score[-1,symp_inx]

            # mark fixed symptoms for plotting
            if symptom in fixed_symptom_weights:
                symptom_labels_marked.append(symptom + " (fixed)")
            else:
                symptom_labels_marked.append(symptom)

        else:
            Impr_pred[symp_inx, 1] = 1.0  # no weight optimization
            symptom_labels_marked.append(symptom + " (default)")

        symp_inx += 1

    return Impr_pred, estim_symp_improv_dict, symptom_labels_marked

def get_activation_prediction(current_protocol, activation_profile, pathways, symp_distances, profile_dict, Soft_SE_dict, side, plot_results = False, score_symptom_metric='Canberra', estim_weights_and_total_score=0,fixed_symptom_weights=[]):

    ''' call via Improvement4Protocol.py '''

    # # estimate the improvement: canberra distance at null activation vs the optimized
    # null_protocol = len(min_bound_per_contact) * [0.0]
    # null_activation_profile = approx_model.predict(np.reshape(np.array(null_protocol), (-1, len(null_protocol))), verbose=0)
    # null_activation_profile = null_activation_profile[0]  # get the actual array

    # or just assign directly. Note, the former might be more reasonable if damaged neurons are considered as activated!
    null_activation_profile = np.zeros(activation_profile.shape[0])

    from Optim_strategies import get_symptom_distances
    [__, null_symptom_diff, symptoms_list] = get_symptom_distances(null_activation_profile, profile_dict, Soft_SE_dict, [], pathways, side, score_symptom_metric)

    # also get symptom distances for 100% activation to estimate worst case scenario for soft-side effects
    max_activation_profile = 100.0 * np.ones(activation_profile.shape[0])
    [__, max_symptom_diff, symptoms_list] = get_symptom_distances(max_activation_profile, profile_dict, Soft_SE_dict, [], pathways, side, score_symptom_metric)

    Impr_pred, estim_symp_improv_dict, symptom_labels_marked = get_improvement_from_distance(profile_dict, Soft_SE_dict, side, symp_distances, max_symptom_diff, null_symptom_diff, estim_weights_and_total_score,fixed_symptom_weights)

    # save json
    if side == 0:
        with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Estim_symp_improv_rh.json', 'w') as save_as_dict:
            json.dump(estim_symp_improv_dict, save_as_dict)
    else:
        with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Estim_symp_improv_lh.json', 'w') as save_as_dict:
            json.dump(estim_symp_improv_dict, save_as_dict)

    if plot_results == True:
        plot_results_with_weights(current_protocol, activation_profile, pathways, Impr_pred, symptom_labels_marked, side)