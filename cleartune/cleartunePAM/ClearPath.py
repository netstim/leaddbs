'''
This script trains a MARS model to approximate pathway activation for a given position

Input: Preprocessed patient stimulation folder (run them once using Lead-DBS)
Steps:
    1) Solve for training set defined by LHS
    2) Solve for test set defined by LHS or other sampling
    3) Check the accuracy
    4) Pickle the model
'''

import os
import subprocess
import h5py
import csv
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set()
import sys
import pickle
import json

if __name__ == '__main__':


    print(sys.argv[1:])
    #raise SystemExit

    ''' User Input '''

    lead_dbs_folder = sys.argv[1]
    os.environ['OSSDIR'] = lead_dbs_folder + 'ext_libs/OSS-DBS'
    netblend_dict_file = sys.argv[2]
    side = sys.argv[3]
    os.environ['STIMDIR'] = sys.argv[4]


    # hardwired for now
    check_trivial = False  # if True, will check additional monopolar and bipolar protocols

    with open(netblend_dict_file, 'r') as fp:
        netblend_dict = json.load(fp)
    fp.close()
    netblend_dict['similarity_metric'] = 'Canberra'  # or Bray-Curtis, Euclidean, etc
    netblend_dict['optim_alg'] = 'dual_annealing'  # or PSO
    netblend_dict['num_iterations_ANN'] = 100  # number of ANN iterations to optimize current at the given electrode position




    # "hardwired" parameters
    one_pol_current_threshold = 8.0  # in mA
    total_current_threshold = 8.0
    FF_dict = True  # always called from FF
    if sys.platform == 'linux':
        docker_image = 'ningfei/oss-dbs:custom'
    elif sys.platform == 'darwin' or sys.platform == 'win32':
        docker_image = 'ningfei/oss-dbs:latest'
    else:
        print("The system's OS does not support the OSS-DBS docker image")


    from Improvement4Protocol import create_NB_dictionaries
    create_NB_dictionaries(os.environ['STIMDIR'], side, netblend_dict['FF_dictionary'])


    # run OSS-DBS as just for StimSets (the folder should be preproced in lead-dbs as "prepare for cluster", but in this case you need to call Axon_allocation and run as cluster)
    # call Axon_allocation.py, then Integrator (see how you run it on ERIS)

    # we can also just run it directly from Lead-DBS
    #output = subprocess.run(
    #    ['docker', 'run', '-e', 'PATIENTDIR', '-e', 'TZ', '--volume', os.environ['OSSDIR'] + ':/opt/OSS-DBS',
    #     '--volume', os.environ['STIMDIR'] + ':/opt/Patient',
    #      '--rm', docker_image, 'python3', 'Launcher_OSS_lite.py', '1'])  #


    # load StimSets_parameters (were created by Train_Test_Generator.py)
    with open(os.environ['STIMDIR'] + '/StimSets_info.json', 'r') as fp:
        StimSets_info = json.load(fp)
    fp.close()

    # load fixed symptoms
    try:
        with open(os.environ['STIMDIR'] + '/Fixed_symptoms.json', 'r') as fp:
            fixed_symptom_weights = json.load(fp)
        fp.close()
    except:
        print("Either no symptoms were fixed or the file is missing, proceeding...")
        fixed_symptom_weights = {}

    """ Train and check ANN model, also stores it in STIMDIR/NB_side/ """
    from ANN_module import train_test_ANN

    # those names will be actually hardcoded according to the OSS-DBS notation,
    # but we could also pass them with the json
    TrainTest_currents_file = os.environ['STIMDIR'] + "/Current_protocols.csv"
    TrainTest_activation_file = os.environ['STIMDIR'] + "/Activations_over_iterations.csv"

    approx_pathways = train_test_ANN(TrainTest_currents_file, TrainTest_activation_file, StimSets_info['trainSize_actual'], netblend_dict['Err_threshold'], netblend_dict['SE_err_threshold'], side, check_trivial)

    # at this point we need to make a decision whether pathways NOT in approx_pathways can be just set to 0

    """ Run the optimization """

    from NB_outline import launch_weight_optimizer
    launch_weight_optimizer(netblend_dict, fixed_symptom_weights, side, approx_pathways)

    """ BELOW IS THE MARS BLOCK, NOT RELEVANT ATM """
    #
    # ##  retrieve the error, the test currents and the activation
    # if side == 0:
    #     res_folder = 'Results_rh'
    # else:
    #     res_folder = 'Results_lh'
    #
    # hf = h5py.File(os.environ['STIMDIR'] + '/' + res_folder + '/Summary_status_0.h5', 'r')
    # Simulated_Pathways = list(hf.keys())
    # hf.close()
    # ActivationResults = np.genfromtxt(os.environ['STIMDIR'] + '/' + res_folder + '/Activations_over_iterations.csv',
    #                                   delimiter=' ')
    # ActivationResults_training = ActivationResults[:trainSize_actual, :]
    # ActivationResults_test = ActivationResults[trainSize_actual:, :]
    #
    # # reload current samples, because some could have been excluded
    # all_samples = np.genfromtxt(os.environ['STIMDIR'] + '/Current_protocols_' + str(side) + '.csv', delimiter=',',
    #                             skip_header=True)
    # training_samples = all_samples[:trainSize_actual]
    # test_samples = all_samples[trainSize_actual:, :]
    #
    # errors_on_test = np.genfromtxt(os.environ['STIMDIR'] + '/' + res_folder + '/errors_on_test.csv', delimiter=' ')
    #
    #
    # # call PAM_plotter to get axons_in_path
    # # checking from here excludes the downsampling issue
    # axons_in_path = np.zeros(len(Simulated_Pathways),int)
    #
    #
    # # plot the errors
    # plt.figure()
    # for i in range(len(Simulated_Pathways)):
    #
    #     #print(np.mean(errors_on_test[:, i]))
    #     sns.kdeplot(errors_on_test[:, i]/axons_in_path[i], bw_adjust=0.5, label=Simulated_Pathways[i])
    #     #err_dif = np.abs(errors_on_fit_rnd[:, i - 1] / axons_in_path[i - 1]) - np.abs(error_ANN[:, loc_ind])
    #     #print(err_dif.mean())
    #
    # #plt.title('Abs error difference between MARS and ANN for percent activation')
    # plt.title('Percent Activation error using MARS 2nd degree')
    # plt.legend()
    # plt.savefig(os.environ['STIMDIR'] + '/Images/Percent_activation_profile_' + str(side) + '.png', format='png',
    #                     dpi=1000)
    #
    #
    #
    #
    #
    # # run NB_outline in the container
    #
    # # this will work only within the container with PyEarth
    # # we actually have a model for each pathway!
    # with open(os.environ['STIMDIR'] + '/' + res_folder + '/MARS_model.file', "rb") as f:
    #     MARS_model = pickle.load(f)

    # if accepted, return True


    # run NB_outline
