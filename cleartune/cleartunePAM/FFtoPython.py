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


    # add more buttons to Currentune

    # 0) Button: "choose a patient folder" with textbox
        # a) Currentune creates its own stim folder in native
    # 1) Predict patient button with a text above "Run this patient in Lead-DBS with the exported connectome"
    # 2) Prepare training-test sets for ANN (calls create_Training_Test_sets(), stores in that folder)
    # 3) Launch PAM calculation button
        # a) Throws a message: It will take time! You can re-run training for your specific axon parameters via Lead-DBS
        # b) Call ea_genvat_butenko with some hardcoded parameters. Try to save varargin{1} and {2} to Currentune folder
        # c) Runs if ok is pressed. Upon completion, throws a "success" message
    # 4) Launch Network Blending Button
        # a) if all or all but one fixed, we should consider running current optimizer directly in OSS-DBS / Cleartune

    print(sys.argv[1:])
    raise SystemExit

    ''' User Input '''

    # Training - Test parameters
    sample_size = 12000  # half training, half testing (compare with 4 contacts systems)
    segm_threshold = 2.0  # in mA
    conc_threshold = 4.0
    one_pol_current_threshold = 8.0  # in mA
    total_current_threshold = 8.0
    SE_err_threshold = 0.10  # Side-effect thresh, refuse if any error > 10% or if 1% of errors > 10% / 2
    Err_threshold = 0.15     # refuse if any error > 15% or if 1% of errors > 15% / 2
    check_trivial = False  # if True, will check additional monopolar and bipolar protocols

    # pass directly from Lead-DBS (not oss-dbs_parameters.mat!)
    #os.environ['OSSDIR'] = '/home/cerebellum/Documents/leaddbs-develop/ext_libs/OSS-DBS'
    #os.environ['STIMDIR'] = '/home/cerebellum/Documents/Data/NetBlend'
    os.environ['STIMDIR'] = '/home/konstantin/Documents/example_pt/5/stimulations/MNI_ICBM_2009b_NLIN_ASYM/20230313191938'


    # "hardwired" parameters
    os.environ['OSSDIR'] = lead_dbs_folder + 'ext_libs/OSS-DBS'
    FF_dict = True
    if sys.platform == 'linux':
        docker_image = 'ningfei/oss-dbs:custom'
    elif sys.platform == 'darwin' or sys.platform == 'win32':
        docker_image = 'ningfei/oss-dbs:latest'
    else:
        print("The system's OS does not support the OSS-DBS docker image")


    Electrode_model = 'Boston Scientific Vercise Directed'
    side = 0  # 1 for lh
    FF_dict_path = '/home/konstantin/Documents/Nando.json'

    create_NB_dictionaries(os.environ['STIMDIR'], side, FF_dict_path, disease='PD')

    ## prepare Training - Test protocols
    #from TrainTest_Generator import create_Training_Test_sets
    #create_Training_Test_sets(os.environ['STIMDIR'], Electrode_model, conc_threshold, segm_threshold, side)

    """ call docker with Launcher_MARS, pass sample sizes for splitting  """
    # # stores activations, simulated pathways and errors (actual - predicted)
    # output = subprocess.run(
    #     ['docker', 'run', '-e', 'PATIENTDIR', '-e', 'TZ', '--volume', OSS_DBS_folder + ':/opt/OSS-DBS',
    #      '--volume', path_to_stim_folder + ':/opt/Patient',
    #      '-it', '--rm', docker_image, 'python3', 'Launcher_MARS.py', trainSize_actual, testSize_actual])  #


    # # alternatively, just compute activations

    #output = subprocess.run(
    #    ['docker', 'run', '-e', 'PATIENTDIR', '-e', 'TZ', '--volume', os.environ['OSSDIR'] + ':/opt/OSS-DBS',
    #     '--volume', os.environ['STIMDIR'] + ':/opt/Patient',
    #      '--rm', docker_image, 'python3', 'Launcher_OSS_lite.py', '1'])  #





    """ Train and check ANN model, also stores it in STIMDIR/NB_side/ """
    from ANN_module import train_test_ANN

    # those names will be actually hardcoded according to the OSS-DBS notation
    TrainTest_currents_file = "/home/cerebellum/Documents/Data/NetBlend/Current_protocols_8615.csv"
    TrainTest_activation_file = "/home/cerebellum/Documents/Data/NetBlend/Activations_over_iterations_8615.csv"
    trainSize_actual = 4898  # this is defined above
    #approx_pathways = ['ACC_cp_right', 'M1_cf_face_right', 'M1_cf_lowerex_right', 'M1_cf_upperex_right', 'PreMotor_cf_right', 'R_ACC_hdp_right', 'R_M1_hdp_face_right', 'R_M1_hdp_lowerex_right', 'R_M1_hdp_upperex_right', 'R_PreMotor_hdp_right', 'R_SMA_hdp_right', 'R_dlPFC_hdp_right', 'R_vmPFC_hdp_right', 'SMA_cf_right', 'ansa_lenticularis_right', 'cerebellothalamic_right', 'dlPFC_cp_right', 'dmPFC_cp_right', 'gpe2stn_ass_right', 'gpe2stn_sm_right', 'lenticular_fasciculus_right', 'vlPFC_cp_right', 'vmPFC_cp_right']

    approx_pathways = train_test_ANN(TrainTest_currents_file, TrainTest_activation_file, trainSize_actual, Err_threshold, SE_err_threshold, side, check_trivial)

    # at this point we need to make a decision whether pathways NOT in approx_pathways can be just set to 0


    """ Run the optimization """

    from NB_outline import launch_weight_optimizer
    launch_weight_optimizer(Electrode_model, fixed_symptom_weights, side, approx_pathways)

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
