'''
    By K. Butenko
    Functions for PathwayTune (see description in the headers)
'''


import json
import os
import sys
from scipy.spatial.distance import canberra, cityblock, euclidean, braycurtis, cosine
import numpy as np
import copy
import h5py
import shutil

# def create_netblend_dictionaries(side, ActivProfileDict, disease='spontaneous human combustion'):
#
#     ''' Either imports Activation Profile Dictionary from Fiber Filtering (ActivProfileDict)
#         or takes a pre-defined dictionary for the given disease from TractSymptomLibrary
#         and saves them in separate dictionaries for modulated symptoms, "soft" side-effects and "hard" side-effects
#         in the stimulation folder / NB_
#
#     Parameters
#     ----------
#     side: int, hemisphere index (0 - right, 1 - left)
#     ActivProfileDict: str, optional, path to Activation Profile Dictionary from Fiber Filtering, otherwise uses a pre-defined dictionary from TractSymptomLibrary
#
#     '''
#
#     # create an output folder in the stim folder of the patient
#     try:
#         os.makedirs(stim_dir + '/NB_' + str(side))
#     except:
#         print("NB folder already exists")
#
#     # retrieve activation profiles for the specific disease
#     # "side" in the name does not actually play the role, it will be checked "on site"
#     if ActivProfileDict != 'no dictionary':
#         with open(ActivProfileDict, 'r') as fp:
#             profiles = json.load(fp)
#         fp.close()
#
#         if 'profile_dict' in profiles.keys():
#             with open(stim_dir + '/NB_' + str(side) + '/profile_dict.json', 'w') as save_as_dict:
#                 json.dump(profiles['profile_dict'], save_as_dict)
#             profile_dict = profiles['profile_dict']
#         else:
#             print("profile_dict is missing, exiting...")
#             raise SystemExit
#
#         if 'Soft_SE_dict' in profiles.keys():
#             Soft_SE_dict = profiles['Soft_SE_dict']
#         else:
#             print("Soft_SE_dict is missing")
#             print("Will continue, but consider exporting negative tracts from FF")
#             Soft_SE_dict = {}
#
#         with open(stim_dir + '/NB_' + str(side) + '/Soft_SE_dict.json', 'w') as save_as_dict:
#             json.dump(Soft_SE_dict, save_as_dict)
#
#         # IMPORTANT: What about SE_dict? Take it from some prior studies?
#         if 'SE_dict' in profiles.keys():
#             SE_dict = profiles['SE_dict']
#         else:
#             print("SE_dict is missing")
#             print("Will continue, but consider exporting negative tracts from FF")
#             SE_dict = {}
#
#         with open(stim_dir + '/NB_' + str(side) + '/SE_dict.json', 'w') as save_as_dict:
#             json.dump(SE_dict, save_as_dict)
#
#     else:
#         from TractSymptomLibrary import get_disease_profiles
#         [profile_dict, Soft_SE_dict, SE_dict] = get_disease_profiles(disease)
#         # here we need to implement a selection mechanism (but no adjustment)
#
#         with open(stim_dir + '/NB_' + str(side) + '/profile_dict.json', 'w') as save_as_dict:
#             json.dump(profile_dict, save_as_dict)
#
#         with open(stim_dir + '/NB_' + str(side) + '/Soft_SE_dict.json', 'w') as save_as_dict:
#             json.dump(Soft_SE_dict, save_as_dict)
#
#         with open(stim_dir + '/NB_' + str(side) + '/SE_dict.json', 'w') as save_as_dict:
#             json.dump(SE_dict, save_as_dict)
#
#     return profile_dict, Soft_SE_dict, SE_dict


class ResultPAM:

    """ Pathway percent activations for a given current protocol that can be used to estimate stimulation outcome"""

    def __init__(self, side, stim_dir, current_protocol=None, inters_as_stim=False):

        # create an output folder in the stim folder of the patient
        self.target_profiles = None
        self.side = side
        self.stim_dir = stim_dir
        if current_protocol:
            self.current_protocol = current_protocol
        else:
            self.current_protocol = self.get_current_protocol()

        if self.side == 0:
            self.side_suffix = '_rh'
        else:
            self.side_suffix = '_lh'

        self.SE_dict = {}  # will be filled out if SE_dict is in self.target_profiles

        try:
            os.makedirs(self.stim_dir + '/NB' + self.side_suffix)
        except:
            print("NB folder already exists")

        # we can just check the distance for one activation profile across all simulated fibers
        self.activation_profile, self.sim_pathways = self.load_AP_from_OSSDBS(inters_as_stim)

        # activation_profile: Nx1 numpy.ndarray, percent activation for simulated pathways
        # sim_pathways: list, pathways simulated in OSS-DBS for this patient

    def get_target_profiles(self, ActivProfileDict=None, disease='spontaneous human combustion'):

        """

        Parameters
        ----------
        ActivProfileDict: str, optional, path to Activation Profile Dictionary from Fiber Filtering, otherwise uses a pre-defined dictionary from TractSymptomLibrary
        disease: str, optional, key to retrieve Activation Profile Dictionary from TractSymptomLibrary

        Returns
        ----------
        target_profiles: dict, Activation Profile Dictionary

        """

        if ActivProfileDict:
            with open(ActivProfileDict, 'r') as fp:
                self.target_profiles = json.load(fp)
            fp.close()
        else:
            from TractSymptomLibrary import get_disease_profiles
            self.target_profiles = get_disease_profiles(disease)

        # copy to NB/ in stim folder for the reference
        shutil.copyfile(ActivProfileDict,self.stim_dir + '/NB' + self.side_suffix + '/target_profiles.json')


    def load_AP_from_OSSDBS(self,inters_as_stim=False):

        """ Load activation profile and pathway labels from OSS-DBS results

        Parameters
        ----------
        inters_as_stim: bool, optional,  if true, fibers inside encapsulation and/or outside of the domain will be treated as activated

        Returns
        -------
        percent activation, Nx1 numpy.ndarray
        simulated pathways, list

        """

        # get all pathways that survived Kuncel(!) pre-filtering and original(!) number of fibers
        # this function will work only for a proper Lead-DBS import (connectome folder, oss-dbs_parameters.mat)
        from Pathways_Stats import get_simulated_pathways
        pathways, axons_in_path = get_simulated_pathways(self.side, self.stim_dir)

        res_folder = self.stim_dir + '/' + 'Results' + self.side_suffix

        perc_activation = np.zeros(len(pathways), float)
        pathway_index = 0
        for pathway in pathways:

            with open(res_folder + '/Pathway_status_' + pathway + '.json', 'r') as fp:
                pathway_results = json.load(fp)
            fp.close()

            if inters_as_stim == True:
                perc_activation[pathway_index] = (pathway_results['percent_activated'] + pathway_results['percent_damaged']) * 0.01 # rescale
            else:
                perc_activation[pathway_index] = pathway_results['percent_activated'] * 0.01 # rescale

            pathway_index += 1

        return perc_activation, pathways
    
    def get_symptom_distances(self, profile_to_check, fixed_symptom_weights,
                              score_symptom_metric='Canberra'):

        """ Compute distances in pathway activation space from target_profiles of symptoms and soft-side effects
            to the activation_profile

        Parameters
        ----------
        profile_to_check: Nx1 numpy.ndarray, percent activation for simulated pathways
        fixed_symptom_weights: dictionary with fixed weights in network blending
        score_symptom_metric: str, optional, metric to compute distances in symptom space

        Returns
        ----------
        symp_dist: Nx1 numpy.ndarray, distances to target profiles for simulated pathways
        symp_list: list, labels of the corresponding symptoms
        sum_symp_dist_nf: float, sum of distances for symptoms with non-fixed weighting

        """

        # "flatten" target profiles for symptoms and threshold profiles for soft side-effects
        Target_profiles_and_SE = copy.deepcopy(self.target_profiles['profile_dict'])
        if 'Soft_SE_dict' in self.target_profiles:
            Target_profiles_and_SE.update(self.target_profiles['Soft_SE_dict'])

        N_symptoms_side = 0
        for key in Target_profiles_and_SE:
            if self.side == 0 and "_rh" in key:
                N_symptoms_side += 1
            elif self.side == 1 and "_lh" in key:
                N_symptoms_side += 1

        symp_dist = np.zeros(N_symptoms_side, float)
        symp_inx = 0
        sum_symp_dist_nf = 0
        symp_list = []

        for key in Target_profiles_and_SE:
            if self.side == 0 and not ("_rh" in key):
                continue
            elif self.side == 1 and not ("_lh" in key):
                continue

            target_rates = []
            predicted_rates = []
            weights_for_pathways = []

            activ_target_profile = list(Target_profiles_and_SE[key].keys())

            for i in range(len(activ_target_profile)):

                target_rates.append(Target_profiles_and_SE[key][activ_target_profile[i]][0])
                weights_for_pathways.append(Target_profiles_and_SE[key][activ_target_profile[i]][2])

                if activ_target_profile[i] in self.sim_pathways:

                    inx = self.sim_pathways.index(activ_target_profile[i])

                    # if the activation is below the threshold, assign the threshold (so that the distance is 0)
                    if 'Soft_SE_dict' in self.target_profiles and key in self.target_profiles['Soft_SE_dict'] and \
                            Target_profiles_and_SE[key][activ_target_profile[i]][0] > profile_to_check[inx]:
                        predicted_rates.append(Target_profiles_and_SE[key][activ_target_profile[i]][0])
                    else:
                        predicted_rates.append(profile_to_check[inx])

                else:  # if not a part of the approx model, assign the threshold (so that the distance is 0)
                    predicted_rates.append(Target_profiles_and_SE[key][activ_target_profile[i]][0])
                    # print("Percent activation was not found for pathway ", activ_target_profile[i], "assigning null distance")

            # within the symptom, weights_for_pathways should sum up to 1
            # based on the acceptance of symptom-tract val equality
            weights_for_pathways = np.array(weights_for_pathways) / sum(weights_for_pathways)

            # distance in less important pathways is less penalized WITHIN the symptom
            if score_symptom_metric == 'Canberra':
                symp_dist[symp_inx] = canberra(predicted_rates, target_rates,
                                                    w=weights_for_pathways)  # / len(weights_for_pathways)
            elif score_symptom_metric == 'Manhattan':
                symp_dist[symp_inx] = cityblock(predicted_rates, target_rates,
                                                     w=weights_for_pathways)  # / len(weights_for_pathways)
            elif score_symptom_metric == 'Euclidean':
                symp_dist[symp_inx] = euclidean(predicted_rates, target_rates,
                                                     w=weights_for_pathways)  # / len(weights_for_pathways)
            elif score_symptom_metric == 'Cosine':
                symp_dist[symp_inx] = cosine(predicted_rates, target_rates,
                                                  w=weights_for_pathways)  # / len(weights_for_pathways)
            elif score_symptom_metric == 'Bray-Curtis':
                symp_dist[symp_inx] = braycurtis(predicted_rates, target_rates,
                                                      w=weights_for_pathways)  # / len(weights_for_pathways)
            else:
                print("Metric is not supported")
                raise SystemExit

            symp_list.append(key)  # symptoms / soft-side effects for the given hemisphere (defined as side)

            # also return symptom distances for non-fixed (adjusted) weights
            if key not in fixed_symptom_weights:
                sum_symp_dist_nf += symp_dist[symp_inx]

            symp_inx += 1

        return symp_dist, symp_list, sum_symp_dist_nf

    def check_for_side_effects(self, profile_to_check):

        """ Check if induced pathway activations are above threshold for critical side-effects.

        Parameters
        ----------
        profile_to_check: Nx1 numpy.ndarray, percent activation for simulated pathways

        Returns:
        ----------
        dict, side-effect dictionary with response status

        """

        SE_threshold_profile = copy.deepcopy(self.target_profiles['SE_dict'])
        SE_threshold_profile_side = {}
        for key in SE_threshold_profile:
            if self.side == 0 and not ("_rh" in key):
                continue
            elif self.side == 1 and not ("_lh" in key):
                continue
            else:
                SE_threshold_profile_side[key] = SE_threshold_profile[key]

            target_rates = []
            predicted_rates = []

            activ_target_profile = list(SE_threshold_profile[key].keys())

            for i in range(len(activ_target_profile)):

                target_rates.append(SE_threshold_profile[key][activ_target_profile[i]][0])

                if activ_target_profile[i] in self.sim_pathways:

                    inx = self.sim_pathways.index(activ_target_profile[i])
                    predicted_rates.append(profile_to_check[inx])

                else:  # if not a part of the approx model, assign 0 for side-effect pathways
                    predicted_rates.append(0.0)
                    # print("Percent activation was not found for pathway ", activ_target_profile[i], "assigning null distance")

                # check if above the threshold
                if predicted_rates[-1] >= target_rates[-1]:
                    SE_threshold_profile_side[key]["predicted"] = 1
                    break
                else:
                    SE_threshold_profile_side[key]["predicted"] = 0

        return SE_threshold_profile_side

    def get_current_protocol(self):

        """ Based on Lead-DBS files for OSS-DBS, get the simulated current protocol (in A)

        Returns:
        ----------
        Pulse_imp: list, currents across electrode contacts in A

        """

        file = h5py.File(self.stim_dir + '/oss-dbs_parameters.mat', mode='r')

        # if file['settings']['current_control'][0][index_side] != 1:
        #    print('The imported protocol is not current-controlled!')
        #    raise SystemExit

        Pulse_amp = file['settings']['Phi_vector'][:, self.side]
        Pulse_amp = Pulse_amp * 0.001  # because Lead-DBS uses mA as the input, switch to A here for consistency

        Pulse_amp = list(Pulse_amp)
        import math
        for i in range(len(Pulse_amp)):
            if math.isnan(Pulse_amp[i]):
                Pulse_amp[i] = 0.0  # IMPORTANT: not grounding here, but a 0A contact

        return Pulse_amp

    def get_improvement_from_distance(self, symp_dist, max_symptom_dist, null_symptom_dist,
                                      estim_weights_and_total_score=0, fixed_symptom_weights=[]):

        """ Estimate improvement/worsening from distances computed in symptom space

        Parameters
        ----------
        symp_dist: Nx1 numpy.ndarray, distances to target profiles for simulated pathways
        max_symptom_dist: Nx1 numpy.ndarray, distances to target profiles from 100% activation
        null_symptom_dist: dNx1 numpy.ndarray, distances to target profiles from 0% activation
        estim_weights_and_total_score: dNxM numpy.ndarray,
        fixed_symptom_weights: list, optional, labels for symptoms with fixed weighting

        Returns
        ----------

        I_hat: Nx2 numpy.ndarray, estimated improvement and weight for symptoms
        estim_symp_improv_dict: dict,
        symptom_labels_marked: list, shows whether a symptom weight was fixed

        """

        # "flatten" target profiles for symptoms and threshold profiles for soft side-effects
        Target_profiles_and_SE = copy.deepcopy(self.target_profiles['profile_dict'])
        if 'Soft_SE_dict' in self.target_profiles:
            Target_profiles_and_SE.update(self.target_profiles['Soft_SE_dict'])

        # check how many symptoms / soft side-effects we have for that hemisphere
        N_symptoms_side = 0
        for key in Target_profiles_and_SE:
            if self.side == 0 and "_rh" in key:
                N_symptoms_side += 1
            elif self.side == 1 and "_lh" in key:
                N_symptoms_side += 1
        I_hat = np.zeros((N_symptoms_side, 2), float)  # in the second dimension, we store the estimated weight

        # iterate over symptoms to estimate symptom improvement (from null activation)
        # or side-effect worsening (assuming the worst case at 100% activation)
        symptom_labels_marked = []
        symp_inx = 0
        estim_symp_improv_dict = {}

        for symptom in Target_profiles_and_SE:

            if self.side == 0 and not ("_rh" in symptom):
                continue
            elif self.side == 1 and not ("_lh" in symptom):
                continue


            if 'Soft_SE_dict' in self.target_profiles and symptom in self.target_profiles[
                'Soft_SE_dict']:  # we assume there are no soft side-effects at null protocol

                # IMPORTANT: for soft-side effect we calculate predicted worsening in comparison to the maximum worsening at 100% activation
                if max_symptom_dist[symp_inx] == 0.0:
                    I_hat[symp_inx, 0] = 0.0
                else:
                    # Here the value is always negative, i.e. worsening
                    I_hat[symp_inx, 0] = (max_symptom_dist[symp_inx] - symp_dist[symp_inx]) / max_symptom_dist[
                        symp_inx] - 1.0
            else:
                # we might have all pathways excluded for the symptom
                # in this case predict zero improvement
                if null_symptom_dist[symp_inx] == 0.0:
                    I_hat[symp_inx, 0] = 0.0
                else:
                    I_hat[symp_inx, 0] = (null_symptom_dist[symp_inx] - symp_dist[symp_inx]) / \
                                             null_symptom_dist[symp_inx]

            estim_symp_improv_dict[symptom] = I_hat[symp_inx, 0]

            # add info for weights if Network Blending was conducted
            if np.any(estim_weights_and_total_score != 0):
                # estimated weight for the symptom, the order was preserved (we always iterate over the symptom dictionary)
                I_hat[symp_inx, 1] = estim_weights_and_total_score[-1, symp_inx]

                # mark fixed symptoms for plotting
                if symptom in fixed_symptom_weights:
                    symptom_labels_marked.append(symptom + " (fixed)")
                else:
                    symptom_labels_marked.append(symptom)

            else:
                I_hat[symp_inx, 1] = 1.0  # no weight optimization
                symptom_labels_marked.append(symptom + " (default)")

            symp_inx += 1

        return I_hat, estim_symp_improv_dict, symptom_labels_marked

    def get_activation_prediction(self, symp_dist,
                                  plot_results=False, score_symptom_metric='Canberra', estim_weights_and_total_score=False,
                                  fixed_symptom_weights=[]):

        """ Predict and store symptom-profile improvement for given symptom distances, optimally generate relevant plots

        Parameters
        ----------
        symp_dist: Nx1 numpy.ndarray, distances to target profiles for simulated pathways
        plot_results: bool, optional
        score_symptom_metric: str, optional, metric to compute distances in symptom space
        estim_weights_and_total_score: NxM numpy.ndarray,
        fixed_symptom_weights: list, optional, labels for symptoms with fixed weighting

        """


        # # estimate the improvement: canberra distance at null activation vs the optimized
        # null_protocol = len(min_bound_per_contact) * [0.0]
        # null_activation_profile = approx_model.predict(np.reshape(np.array(null_protocol), (-1, len(null_protocol))), verbose=0)
        # null_activation_profile = null_activation_profile[0]  # get the actual array

        # or just assign directly. Note, the former might be more reasonable if damaged neurons are considered as activated!
        null_activation_profile = np.zeros(self.activation_profile.shape[0])

        [null_symp_dist, symptoms_list, __] = self.get_symptom_distances(null_activation_profile, [],
                                                                       score_symptom_metric)

        # also get symptom distances for 100% activation to estimate worst case scenario for soft-side effects
        max_activation_profile = 100.0 * np.ones(self.activation_profile.shape[0])
        [max_symp_dist, __, __] = self.get_symptom_distances(max_activation_profile, [],
                                                                      score_symptom_metric)

        I_hat, estim_symp_improv_dict, symptom_labels_marked = self.get_improvement_from_distance(symp_dist,
                                                                                                 max_symp_dist,
                                                                                                 null_symp_dist,
                                                                                                 estim_weights_and_total_score,
                                                                                                 fixed_symptom_weights)

        # save json
        with open(self.stim_dir + '/NB' + self.side_suffix + '/Estim_symp_improv' + self.side_suffix + '.json', 'w') as save_as_dict:
            json.dump(estim_symp_improv_dict, save_as_dict)

        if plot_results == True:
            self.plot_results_with_weights(I_hat, symptom_labels_marked)

    def plot_results_with_weights(self, I_hat, symptom_labels_marked):

        """ Plot activation profile and predicted improvements for the selected stimulation

        Parameters
        ----------
        I_hat: Nx2 numpy.ndarray, improvement and weight for each symptoms
        symptom_labels_marked: list, shows whether a symptom weight was fixed

        """

        import matplotlib.pyplot as plt
        # import seaborn as sns
        # sns.set()

        # ============= Plot predicted improvement and corresponding weights =============#
        pos = np.arange(len(symptom_labels_marked))  # the x locations for the groups
        pos_adjusted = pos  # - 0.5

        fig, ax = plt.subplots(figsize=(12, 8))
        width = 0.25
        colors_opt = ['C0', 'C1']  # different colors for predicted improvement and weights

        # first bars for the improvement
        ax.bar(pos_adjusted - 0.5 * width, I_hat[:, 0] * 100.0, width,
               color=colors_opt[0])
        ax.set_ylabel("Improvement / worsening, %", color='C0')
        # ax.axhline(0)
        # ax.set_ylim(-100,100)

        # seond bars for the weights
        ax2 = ax.twinx()
        ax2.bar(pos_adjusted + 0.5 * width, I_hat[:, 1], width,
                color=colors_opt[1])
        ax2.grid(False)

        # autolimit to align axes (potentially buggy!)
        if np.any(I_hat[:, 0] < 0.0):
            # lower_w_lim = np.min(I_hat[:, 0]) / np.max(I_hat[:, 0])
            lower_w_lim = np.min(I_hat[:, 0]) / 1.0

            ax.set_ylim(np.min(I_hat[:, 0]) * 100, 100.0)
            ax2.set_ylim(lower_w_lim, 1)
        else:
            ax.set_ylim(0, 100)
            ax2.set_ylim(0, 1)

        # ax2.set_ylim(0, 1)
        # ax2.axhline(0)

        # ax.legend(loc='upper right', bbox_to_anchor=(0, 1.25, 1, 0), ncol=5, mode="expand", borderaxespad=0.)
        ax2.set_ylabel("Suggested symptom weights", color='C1')
        ax.set_xticks(pos_adjusted)
        ax.set_xticklabels(symptom_labels_marked, rotation=45)
        fig.tight_layout()
        plt.savefig(self.stim_dir + '/NB' + self.side_suffix + '/Symptom_profiles' + self.side_suffix + '.png',
                    format='png',
                    dpi=500)

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
        pos = np.arange(len(self.activation_profile))  # the x locations for the groups
        ax.bar(pos, self.activation_profile * 100.0, width,
               color=colors_opt[0])

        # convert to mA
        current_protocol = np.array(self.current_protocol) * 1000.0

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

        # ax.legend(loc='upper right', bbox_to_anchor=(0, 1.25, 1, 0), ncol=5, mode="expand", borderaxespad=0.)
        # ax.set_xticks(pos_adjusted)
        ax.set_xticklabels(self.sim_pathways)
        ax.set_ylabel("Percent Activation, %")

        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0.0, 1.5, textstr, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)

        ax.set_xticks(pos)
        plt.xticks(rotation=45)
        fig.tight_layout()
        plt.savefig(self.stim_dir + '/NB' + self.side_suffix + '/Activation_profile' + self.side_suffix + '.png',
                    format='png',
                    dpi=500)

    def make_prediction(self, score_symptom_metric, ActivProfileDict=None, fixed_symptoms_dict=None, plot_results=True, disease='spontaneous human combustion'):

        """ Predict symptom-profile improvement for a given activation profile based on the target activation profiles
            fixed symptom dictionary is only needed for weight optimization

        Parameters
        ----------
        score_symptom_metric: str, metric to compute distances in symptom space
        ActivProfileDict: str, optional, path to Activation Profile Dictionary from Fiber Filtering, otherwise uses a pre-defined dictionary from TractSymptomLibrary
        fixed_symptoms_dict: str, optional, path to dictionary with fixed weights in network blending
        plot_results: bool, true to generate activation profile and symptom profile plots
        disease: str, optional, key to retrieve Activation Profile Dictionary from TractSymptomLibrary

        """

        # load fixed weights (they play role only for network blending, not simple prediction)
        if fixed_symptoms_dict == 0:
            fixed_symptom_weights = []  # placeholder
        else:
            with open(fixed_symptoms_dict, 'r') as fp:
                fixed_symptom_weights = json.load(fp)
            fp.close()

        self.get_target_profiles(ActivProfileDict,disease)

        # get symptom-wise difference between activation and target profiles
        # we do not need non-fixed symptom distances here
        [symptom_diff, self.symptom_list, __] = self.get_symptom_distances(self.activation_profile,
                                                   fixed_symptom_weights, score_symptom_metric=score_symptom_metric)
        # symptom_diff is in the symptom space, not pathway! So it might have a different dimensionality

        self.get_activation_prediction(symptom_diff, plot_results, score_symptom_metric=score_symptom_metric, estim_weights_and_total_score=False, fixed_symptom_weights=fixed_symptom_weights)

        # check critical side-effects if available
        if 'SE_dict' in self.target_profiles:
            self.SE_dict = self.check_for_side_effects(self.activation_profile)


if __name__ == '__main__':

    # called from MATLAB
    # sys.argv[1] - stim folder
    # sys.argv[2] - side
    # sys.argv[3] - Activation Profile Dictionary based on Fiber Filtering
    # sys.argv[4] - Fixed Symptoms Dictionary
    # sys.argv[5] - Score Symptom Metric

    stim_dir = sys.argv[1]

    stim_result = ResultPAM(int(sys.argv[2]),sys.argv[1])
    stim_result.make_prediction(sys.argv[5], sys.argv[3], sys.argv[4])
