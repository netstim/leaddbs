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

''' User Input '''

# Clinical Input
disease = "PD"  # not needed if exported from FF
fixed_symptom_weights = {     # includes both symptoms and soft-side effects, the sum should not exceed 1
    #'Rigidity_rh': 0.1,
    #'Tremor_rh': 0.25,

    'Motor_Mixed_Postop_Bradykinesia_perc_improvements_rh': 0.5,
    'Motor_Mixed_Postop_Rigidity_perc_improvements_rh': 0.25,
}


# Training - Test parameters
sample_size = 12000  # half training, half testing
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
os.environ['OSSDIR'] = '/home/konstantin/Documents/GitHub/leaddbs/ext_libs/OSS-DBS'
os.environ['STIMDIR'] = '/home/konstantin/Documents/example_pt/5/stimulations/MNI_ICBM_2009b_NLIN_ASYM/20230313191938'

docker_image = 'ningfei/oss-dbs:custom'
Electrode_model = 'Boston Scientific Vercise Directed'
side = 0  # 1 for lh
FF_dict = True
FF_dict_path = '/home/konstantin/Documents/Nando.json'


# create an output folder in the stim folder of the patient
try:
    os.makedirs(os.environ['STIMDIR'] + '/NB_' + str(side))
except:
    print("NB folder already exists")


# retrieve activation profiles for the specific disease
# "side" in the name does not actually play the role, it will be checked "on site"
if FF_dict == True:
    with open(FF_dict_path, 'r') as fp:
        profiles = json.load(fp)
    fp.close()

    if 'profile_dict' in profiles.keys():
        with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/profile_dict.json', 'w') as save_as_dict:
            json.dump(profiles['profile_dict'], save_as_dict)
        profile_dict = profiles['profile_dict']
    else:
        print("profile_dict is missing, exiting...")
        raise SystemExit

    if 'Soft_SE_dict' in profiles.keys():
        Soft_SE_dict = profiles['Soft_SE_dict']
    else:
        print("Soft_SE_dict is missing")
        print("Will continue, but consider exporting negative tracts from FF")
        Soft_SE_dict = {}

    with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Soft_SE_dict.json', 'w') as save_as_dict:
        json.dump(Soft_SE_dict, save_as_dict)

    # IMPORTANT: What about SE_dict? Take it from some prior studies?
    SE_dict = {}
    with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/SE_dict.json', 'w') as save_as_dict:
        json.dump(SE_dict, save_as_dict)

else:
    from TractSymptomLibrary import get_disease_profiles
    [profile_dict, Soft_SE_dict, SE_dict] = get_disease_profiles(disease)
    # here we need to implement a selection mechanism (but no adjustment)

    with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/profile_dict.json', 'w') as save_as_dict:
        json.dump(profile_dict, save_as_dict)

    with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Soft_SE_dict.json', 'w') as save_as_dict:
        json.dump(Soft_SE_dict, save_as_dict)

    with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/SE_dict.json', 'w') as save_as_dict:
        json.dump(SE_dict, save_as_dict)

## prepare Training - Test protocols
#from TrainTest_Generator import create_Training_Test_sets
#create_Training_Test_sets(Electrode_model, sample_size, conc_threshold, segm_threshold, total_current_threshold, one_pol_current_threshold, side)



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


# we can just check the distance for one activation profile across all simulated fibers
from NB_outline import load_AP_from_OSSDBS
activation_profile, Pathways = load_AP_from_OSSDBS(side, inters_as_stim=False)

# get symptom-wise difference between activation and target profiles
from Optim_strategies import get_symptom_distances
[__, symptom_diff, symptoms_list] = get_symptom_distances(activation_profile, profile_dict, Soft_SE_dict,
                                                fixed_symptom_weights, Pathways, side, score_symptom_metric='Canberra')
# symptom_diff is in the symptom space, not pathway! So it might have a smaller dimensionality

from Pathways_Stats import get_current_protocol
current_protocol = get_current_protocol(side)

from RoutinesForResults import get_activation_prediction
get_activation_prediction(current_protocol, activation_profile, Pathways, symptom_diff, profile_dict,
                          Soft_SE_dict, side, plot_results=True, score_symptom_metric='Canberra')


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
