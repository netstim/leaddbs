import pandas as pd
import numpy as np
import sys
import os
import json
import subprocess
from scipy.optimize import dual_annealing
from Improvement4Protocol import ResultPAM



class PamOptimizer:

    """ Class for current optimization based on pathway activation modeling for symptom tracts

    Parameters
    ----------
    side : int
        hemisphere ID (0 - right, 1 - left)
    stim_folder : str
        path to the stimulation folder
    optim_settings_file : str
        path to Activation Profile Dictionary from Fiber Filtering
    neuron_folder: str
        path to the folder containing NEURON models
    PAM_caller_script: str
        path to the high-level script for NEURON
    """

    def __init__(self, side, stim_folder, optim_settings_file, neuron_folder, PAM_caller_script):

        self.side = side
        self.stim_folder = stim_folder
        self.neuron_folder = neuron_folder
        self.PAM_caller_script = PAM_caller_script

        # import external function
        sys.path.insert(1, os.path.dirname(self.PAM_caller_script))

        if self.side == 0:
            self.side_suffix = '_rh'
        else:
            self.side_suffix = '_lh'

        self.results_folder = os.path.join(stim_folder, 'Results' + self.side_suffix)
        self.timeDomainSolution = self.results_folder + '/oss_time_result_PAM.h5'
        self.pathwayParameterFile = self.stim_folder + "/Allocated_axons_parameters.json"

        with open(optim_settings_file, 'r') as fp:
            optim_settings = json.load(fp)
        fp.close()
        self.optim_settings = optim_settings['netblendict']

        if self.optim_settings['optim_alg'] == 'Dual Annealing':

            if os.path.isfile(os.path.join(self.results_folder,'Best_scaling_yet.csv')):
                initial_scaling = np.genfromtxt(os.path.join(self.results_folder,'Best_scaling_yet.csv'), delimiter=' ')
                optimized_current = dual_annealing(self.compute_global_score,
                                     bounds=list(zip(self.optim_settings['min_bound_per_contact'], self.optim_settings['max_bound_per_contact'])),
                                     x0=initial_scaling, maxfun=self.optim_settings['num_iterations'], seed=42, visit=2.62,
                                     no_local_search=True)
            else:
                optimized_current = dual_annealing(self.compute_global_score,
                                     bounds=list(zip(self.optim_settings['min_bound_per_contact'], self.optim_settings['max_bound_per_contact'])),
                                     maxfun=self.optim_settings['num_iterations'], seed=42, visit=2.62,
                                     no_local_search=True)

        elif optim_settings['optim_alg'] == 'PSO':

            from pyswarms.single.global_best import GlobalBestPSO
            bounds = (self.optim_settings['min_bound_per_contact'], self.optim_settings['max_bound_per_contact'])
            options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9}
            optimizer = GlobalBestPSO(n_particles=20, dimensions=len(self.optim_settings['min_bound_per_contact']), options=options, bounds=bounds)

            cost, optimized_current = optimizer.optimize(self.prepare_swarm, iters=self.optim_settings['num_iterations_ANN'])


    def store_iteration_results(self, S_vector, global_score, symptoms_list, estim_Ihat, SE_dict={}):

        """ Store data on each iteration in a .csv file

        Inputs
        ------
        S_vector: numpy.ndarray, stimulation protocol tested by the optimizer
        global_score: float, summed up weighted symptom improvements
        symptoms_list: list, labels of the analyzed symptoms
        estim_Ihat: dict, estimated symptom improvements
        SE_dict: dict, info on critical side-effects

        """

        # store info about the iteration
        print(global_score)
        df = pd.DataFrame(
            {
                "weighted_total_score": [global_score],
            }
        )

        # stimulation current
        current_data = []
        for i in range(len(S_vector)):
            current_data.append(S_vector[i])
        current_data_array = np.vstack((current_data)).T

        # predicted improvement of symptoms
        for key in symptoms_list:
            df[key] = estim_Ihat[key]

        for i in range(len(S_vector)):
            df['Contact_' + str(i)] = S_vector[i]  # Lead-DBS notation!

        iter_file = self.stim_folder + '/NB' + self.side_suffix + '/optim_iterations.csv'
        df.to_csv(iter_file, mode='a', header=not os.path.exists(iter_file))

        # critical side-effect status if provided
        if SE_dict:
            for key in SE_dict:
                df[key] = SE_dict[key]["predicted"]

        iter_file = self.stim_folder + '/NB' + self.side_suffix + '/optim_iterations_CSE.csv'
        df.to_csv(iter_file, mode='a', header=not os.path.exists(iter_file))

    def prepare_swarm(self, x, args_to_pass=[]):

        """ Higher-level method to do forward_prop in the
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
        j = [self.compute_global_score(x[i], args_to_pass[:]) for i in range(n_particles)]
        return np.array(j)

    def compute_global_score(self, S_vector):

        """ Estimate the goal function (weighted symptom improvements) for a given stimulation protocol

        Inputs
        ------
        S_vector: numpy.ndarray, stimulation protocol tested by the optimizer

        Returns
        -------
        float, goal function to maximize - summed up weighted symptom improvements
        """

        # we use scaling, but later we will switch to scaling_vector
        from PAM_caller import launch_PAM
        launch_PAM(self.neuron_folder, self.results_folder, self.timeDomainSolution, self.pathwayParameterFile, S_vector[0] * 100)
        # the original solution for 10 mA
        # so scale by 100

        # make a prediction
        stim_result = ResultPAM(self.side, self.stim_folder)
        stim_result.make_prediction(self.optim_settings['similarity_metric'], self.optim_settings['ActivProfileDict'], self.optim_settings['symptom_weights_file'], plot_results=False)

        # put this in a separate function
        # load predicted symptom improvement and weights
        estim_Ihat_json = self.stim_folder + '/NB' + self.side_suffix + '/Estim_symp_improv' + self.side_suffix + '.json'

        with open(estim_Ihat_json, 'r') as fp:
            estim_Ihat = json.load(fp)
        fp.close()

        # load fixed symptom weights
        fp = open(self.optim_settings['symptom_weights_file'])
        symptom_weights = json.load(fp)
        symptom_weights = symptom_weights['fixed_symptom_weights']
        remaining_weights = 1.0
        N_fixed = 0
        # is it iterating only across the correct side?
        for key in symptom_weights:
            remaining_weights = remaining_weights - symptom_weights[key]
            N_fixed += 1

        # compute weight for non-fixed as the equal distribution of what remained
        if len(stim_result.symptom_list) != N_fixed:
            rem_weight = remaining_weights / (len(stim_result.symptom_list) - N_fixed)
        else:
            rem_weight = 0.0

        print(stim_result.symptom_list)

        # for now, the global score is a simple summation (exactly like in Optim_strategies.py)
        global_score = 0.0
        for key in stim_result.symptom_list:
            if self.side == 0 and not ("_rh" in key):
                continue
            elif self.side == 1 and not ("_lh" in key):
                continue

            # soft side-effect have negative estim_Ihat
            if key in symptom_weights:
                global_score = global_score + symptom_weights[key] * estim_Ihat[key]
            else:
                global_score = global_score + rem_weight * estim_Ihat[key]

        self.store_iteration_results(S_vector, global_score, stim_result.symptom_list, estim_Ihat, stim_result.SE_dict)

        # check if any side-effect responses were predicted
        if stim_result.SE_dict:
            for key in stim_result.SE_dict:
                if stim_result.SE_dict[key]["predicted"]:
                    return 1e9

        return -1 * global_score

if __name__ == '__main__':

    # called from MATLAB
    PAM_caller_script = sys.argv[1]
    neuron_folder = sys.argv[2]
    optim_settings_file = sys.argv[3]
    stim_folder = sys.argv[4]
    side = int(sys.argv[5])
    #proper_python = sys.argv[5]

    # side and stim_folder - from ea_optimizePAM_butenko
    optimization = PamOptimizer(side, stim_folder, optim_settings_file, neuron_folder, PAM_caller_script)


