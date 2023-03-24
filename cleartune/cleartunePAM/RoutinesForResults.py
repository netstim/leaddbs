import numpy as np
import os


def plot_results_with_weights(current_protocol, activation_profile, pathways, Impr_pred, symptom_labels_marked, side, NB_result=False):
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set()

    pos = np.arange(len(symptom_labels_marked))  # the x locations for the groups
    pos_adjusted = pos  # - 0.5

    fig, ax = plt.subplots(figsize=(12, 3))
    width = 0.25
    colors_opt = ['C0', 'C1']
    # labels_opt = ['Activated', 'Not Activated', 'Encapsulation', 'CSF', 'Other']
    # bars_status =
    # for i in range(Impr_pred.shape[1]):
    #     ax.bar(pos_adjusted - 0.5 + i * width, Impr_pred[:, i] , width,
    #            color=colors_opt[i])

    ax.bar(pos_adjusted - 0.5 * width, Impr_pred[:, 0] * 100.0, width,
           color=colors_opt[0])
    ax.set_ylabel("Improvement / worsening, %", color='C0')
    # ax.axhline(0)
    # ax.set_ylim(-100,100)

    ax2 = ax.twinx()
    ax2.bar(pos_adjusted + 0.5 * width, Impr_pred[:, 1], width,
            color=colors_opt[1])
    ax2.grid(False)

    # autolimit to align axes (potentially buggy!)
    if np.any(Impr_pred[:, 0] < 0.0):
        print('here')
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
    ax.set_xticks(pos_adjusted)
    ax.set_xticklabels(symptom_labels_marked)
    ax2.set_ylabel("Suggested symptom weights", color='C1')
    fig.tight_layout()
    plt.savefig(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Symptom_profiles_' + str(side) + '.png',
                format='png',
                dpi=1000)

    # activation profile for optimized res.x, current protocol,
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
            r'$k_{3}=%.2f$' % (current_protocol[0])))
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

def get_activation_prediction(current_protocol, activation_profile, pathways, symp_distances, profile_dict, Soft_SE_dict, side, score_symptom_metric='Canberra', estim_weights_and_total_score=0,fixed_symptom_weights=[]):

    # # estimate the improvement: canberra distance at null activation vs the optimized
    # null_protocol = len(min_bound_per_contact) * [0.0]
    # null_activation_profile = approx_model.predict(np.reshape(np.array(null_protocol), (-1, len(null_protocol))), verbose=0)
    # null_activation_profile = null_activation_profile[0]  # get the actual array

    # or just assign directly. Note, the former might be more reasonable if damaged neurons are considered as activated!
    null_activation_profile = np.zeros(activation_profile.shape[0])

    from Optim_strategies import get_symptom_distances
    [__, null_symptom_diff] = get_symptom_distances(null_activation_profile, profile_dict, Soft_SE_dict, [], pathways, side, score_symptom_metric)


    # also get symptom distances for 100% activation to estimate worst case scenario for soft-side effects
    max_activation_profile = 100.0 * np.ones(activation_profile.shape[0])
    [__, max_symptom_diff] = get_symptom_distances(max_activation_profile, profile_dict, Soft_SE_dict, [], pathways, side, score_symptom_metric)

    # here we can merge target profiles for symptoms and threshold profiles for soft side-effects
    profile_dict.update(Soft_SE_dict)

    # plot predicted improvement and corresponding weights in a combined bar plot

    N_symptoms_side = 0
    for key in profile_dict:
        if side == 0 and "_rh" in key:
            N_symptoms_side += 1
        elif side == 1 and "_lh" in key:
            N_symptoms_side += 1
    Impr_pred = np.zeros((N_symptoms_side,2), float) # in the second dimension, we store the estimated weight
    symp_inx = 0

    symptom_labels_marked = []

    for symptom in profile_dict:

        if side == 0 and not ("_rh" in symptom):
            continue
        elif side == 1 and not ("_lh" in symptom):
            continue

        if symptom in Soft_SE_dict: # we assume there are no soft side-effects at null protocol

            # IMPORTANT: for soft-side effect we calculate predicted worsening to the maximum (100% activation worsening
            if max_symptom_diff[symp_inx] == 0.0:
                Impr_pred[symp_inx,0] = 0.0
            else:
                # Here the value is always negative, i.e. worsening
                Impr_pred[symp_inx,0] = (max_symptom_diff[symp_inx] - symp_distances[symp_inx]) / max_symptom_diff[symp_inx] - 1.0
        else:
            Impr_pred[symp_inx,0] = (null_symptom_diff[symp_inx] - symp_distances[symp_inx]) / null_symptom_diff[symp_inx]

        # add info for weights if Network Blending was conducted
        if estim_weights_and_total_score != 0:
            # estimated weight for the symptom, the order was preserved (we always iterate over symptom dictionary)
            Impr_pred[symp_inx, 1] = estim_weights_and_total_score[-1,symp_inx]

            # mark fixed symptoms for plotting
            if symptom in fixed_symptom_weights:
                symptom_labels_marked.append(symptom + " (fixed)")
            else:
                symptom_labels_marked.append(symptom)

        else:
            Impr_pred[symp_inx, 1] = 1.0 # no weight optimization
            symptom_labels_marked.append(symptom + " (default)")

        symp_inx += 1

    plot_results_with_weights(current_protocol, activation_profile, pathways, Impr_pred, symptom_labels_marked, side)