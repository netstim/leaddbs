'''
    By K. Butenko
    This script trains and test an ANN model to approximate pathway activation for a given electrode position
    IMPORTANT: Based on tensorflow with a hard sigmoid swapped to abs(LeakyReLU) with alpha 1.25 (steeper slope for cathode)
'''



import matplotlib.pyplot as plt
import numpy as np
import os
import h5py
import seaborn as sns
sns.set()
import sys
import json

from sklearn.model_selection import train_test_split
import tensorflow as tf
from sklearn import preprocessing
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, BatchNormalization, Lambda
from sklearn.metrics import accuracy_score, confusion_matrix, precision_score, recall_score, f1_score, precision_recall_curve, auc
from tensorflow.keras import optimizers
from keras.layers import LeakyReLU

### Input

# some ANN parameters
learn_rate = 0.0025
N_epochs = 5000
min_activ_threshold = 0.05   # if less than 5% of fibers in the pathway were activated over all StimSets, ANN will not train on it

def train_test_ANN(TrainTest_currents_file, TrainTest_activation_file, trainSize, Err_threshold, SE_err_threshold, side, check_trivial, VAT_recruit = False):

    import os

    # load currents
    Currents = np.genfromtxt(TrainTest_currents_file, delimiter=',', skip_header=True)

    # load activation results
    ActivationResults = np.genfromtxt(TrainTest_activation_file, delimiter=' ')

    # we can also estimate errors for some simple pre-defined protocols
    check_trivial = False  # disable for now
    if check_trivial == True:
        ActivationResults_Bipolar1 = np.genfromtxt('Activations_over_iterationsBipolar.csv', delimiter=' ')
        Currents_Bipolar1 = np.genfromtxt('/home/cerebellum/Documents/Data/NetBlend/Current_protocols_0_Bipolar.csv', delimiter=',', skip_header=True)

        ActivationResults_Monopolar21 = np.genfromtxt('Activations_over_iterations_Monopolar21_79.csv', delimiter=' ')
        Currents_Monopolar21 = np.genfromtxt('Current_protocols_0_Monopolars21_79.csv', delimiter=',', skip_header=True)

    if VAT_recruit == True:

        from VAT_pathway_recruitment import remove_failed_protocols, get_VAT_pathways

        Currents, ActivationResults = remove_failed_protocols(Currents, ActivationResults)
        Pathways, axons_in_path = get_VAT_pathways(side)
        # note that the train size was not adjusted for the failed protocols!
    else:
        # the function will work only for a proper Lead-DBS import (connectome folder, oss-dbs_parameters.mat)
        # get all pathways that survived Kuncel(!) pre-filtering and original(!) number of fibers
        from Pathways_Stats import get_simulated_pathways
        Pathways, axons_in_path = get_simulated_pathways(side)

    #=============================================== Prepare the data =================================================#
    X_train = Currents[:trainSize,:] * 0.001  # convert to A
    X_test = Currents[trainSize:,:] * 0.001
    y_train_prelim = ActivationResults[:trainSize,1:] / axons_in_path  # from 1, because 0 is the index of the protocol
    y_test_prelim = ActivationResults[trainSize:,1:] / axons_in_path

    y_train = -100 * np.ones((y_train_prelim.shape), float)  # initialize with -100 to remove non-filled value later
    y_test = -100 * np.ones((y_test_prelim.shape), float)

    # optionally check some hardcoded trivial protocols
    if check_trivial == True:

        X_bipolar = Currents_Bipolar1 * 0.001
        X_monopolar = Currents_Monopolar21 * 0.001

        y_bipolar_prelim = ActivationResults_Bipolar1[:,1:] / axons_in_path
        y_monopolar_prelim = ActivationResults_Monopolar21[:,1:] / axons_in_path
        y_bipolar = -100 * np.ones((y_bipolar_prelim.shape), float)
        y_monopolar = -100 * np.ones((y_monopolar_prelim.shape), float)

    # only consider those pathways, where max activation >= min_activ_threshold
    pathway_filtered = []

    for i in range(y_train_prelim.shape[1]):
        # only compute for pathways with some percent activation and minimal number of fibers (10)
        if axons_in_path[i] > 9 and np.max(y_train_prelim[:, i]) >= min_activ_threshold and np.max(y_test_prelim[:, i]) >= min_activ_threshold:
            y_train[:,i] = y_train_prelim[:,i]
            y_test[:,i] = y_test_prelim[:,i]
            pathway_filtered.append(Pathways[i])

            if check_trivial == True:
                y_bipolar[:,i] = y_bipolar_prelim[:,i]
                y_monopolar[:, i] = y_monopolar_prelim[:, i]

    # remove entries for pathways with max activation < min_activ_threshold
    y_train = y_train[:, ~np.all(y_train == -100.0, axis=0)]
    y_test = y_test[:, ~np.all(y_test == -100.0, axis=0)]

    if check_trivial == True:
        y_bipolar = y_bipolar[:,~np.all(y_bipolar == -100.0, axis=0)]
        y_monopolar = y_monopolar[:,~np.all(y_monopolar == -100.0, axis=0)]

    # inject 10% of zero protocols to the training (poor solution)
    N_zero = int(0.1 * y_train.shape[0])
    y_train_exp = np.zeros((y_train.shape[0] + N_zero,y_train.shape[1]), float)
    X_train_exp = np.zeros((X_train.shape[0] + N_zero,X_train.shape[1]), float)
    y_train_exp[:y_train.shape[0],:] = y_train
    X_train_exp[:X_train.shape[0],:] = X_train
    y_train = y_train_exp       # already from 0 to 1
    X_train = X_train_exp       # normalization seems to be not necessary here

    #================================================== Train ANN =====================================================#

    model = Sequential(layers=None, name=None)
    model.add(Dense(128, input_shape=(X_train.shape[1],), activation='linear'))
    model.add(Dense(1024, activation=tf.keras.layers.LeakyReLU(alpha=-1.25)))  # alpha -1.25 to have a steeper slope for cathode
    model.add(Dense(np.sum(axons_in_path), activation='sigmoid'))  # following the percent activation curves
    #model.add(Dense(y_train.shape[1], activation='tanh'))
    model.add(Dense(y_train.shape[1], activation='sigmoid'))

    adam = optimizers.Adamax(lr=learn_rate)
    model.compile(optimizer=adam, loss='mean_squared_error', metrics=['accuracy'])
    model.fit(X_train, y_train, epochs=N_epochs, verbose=1)
    results = model.evaluate(X_test, y_test)

    # on Test
    y_predicted = model.predict(X_test)
    error_ANN = y_test - y_predicted

    if check_trivial == True:

        y_predicted_bi = model.predict(X_bipolar)
        error_ANN_bi = y_bipolar - y_predicted_bi
        results_bi = model.evaluate(X_bipolar, y_bipolar)

        y_predicted_mono = model.predict(X_monopolar)
        error_ANN_mono = y_monopolar - y_predicted_mono

    ## null check
    #X_null = np.array([[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[-0.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0]])
    #y_null_predict = model.predict(X_null) * 100.0

    MSEs = np.zeros(error_ANN.shape[1], float)

    from sklearn.metrics import mean_squared_error
    for i in range(error_ANN.shape[1]):
        MSEs[i] = mean_squared_error(y_test[:,i], y_predicted[:,i])
        #print(Pathways[i], ": ", MSEs[i])

    # from sklearn.metrics import mean_squared_error
    # for i in range(error_ANN_bi.shape[1]):
    #     MSEs[i] = mean_squared_error(y_bipolar[:,i], y_predicted_bi[:,i])
    #     print(Pathways[i], ": ", MSEs[i])

    # =========================================== Plot the errors =====================================================#

    import matplotlib
    matplotlib.rcParams['figure.dpi'] = 200
    plt.figure()

    pathways_max_errors = {}  # also store
    for i in range(len(pathway_filtered)):
        pathways_max_errors[pathway_filtered[i]] = np.max(abs(error_ANN[:, i]))
        if np.max(abs(error_ANN[:, i])) > 0.05:
            sns.kdeplot(error_ANN[:,i], bw_adjust=0.5, label=pathway_filtered[i])

    with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/ANN_abs_errors.json', 'w') as save_as_dict:
        json.dump(pathways_max_errors, save_as_dict)

    plt.legend()
    plt.title('Abs errors for ANN on Test')
    plt.xlim([-0.25,0.25])
    plt.savefig(os.environ['STIMDIR'] + '/NB_' + str(side) + '/ANN_abs_errors_on_Test_' + str(side) + '.png', format='png',
                dpi=1000)


    if check_trivial == True:
        plt.figure()
        for i in range(len(pathway_filtered)):

            if np.max(abs(error_ANN_bi[:,i])) > 0.05:
                sns.kdeplot(error_ANN_bi[:,i], bw_adjust=0.5, label=pathway_filtered[i])

        plt.legend()
        plt.title('Abs errors for ANN on Bipolar')
        #plt.xlim([-0.15,0.15])
        plt.savefig(os.environ['STIMDIR'] + '/NB_' + str(side) + '/ANN_abs_errors_on_Bipolar_' + str(side) + '.png',
                    format='png',
                    dpi=1000)


        plt.figure()
        for i in range(len(pathway_filtered)):

            if np.max(abs(error_ANN_mono[:,i])) > 0.05:
                sns.kdeplot(error_ANN_mono[:,i], bw_adjust=0.5, label=pathway_filtered[i])

        plt.legend()
        plt.title('Abs errors for ANN on Bipolar')
        # plt.xlim([-0.15,0.15])
        plt.savefig(os.environ['STIMDIR'] + '/NB_' + str(side) + '/ANN_abs_errors_on_Monopolar_' + str(side) + '.png',
                    format='png',
                    dpi=1000)

    # ====================================== Check if errors acceptable ===============================================#

    # iterate over each previously approved profile
    with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/profile_dict.json', 'r') as fp:
        profile_dict = json.load(fp)
    fp.close()

    with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/Soft_SE_dict.json', 'r') as fp:
        Soft_SE_dict = json.load(fp)
    fp.close()

    with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/SE_dict.json', 'r') as fp:
        SE_dict = json.load(fp)
    fp.close()

    # we can discard the error sign here
    error_ANN = abs(error_ANN)

    # first check side-effects
    for key in SE_dict:
        if side == 0 and not("_rh" in key):
            continue
        elif side == 1 and not("_lh" in key):
            continue

        activ_threshold_profile = list(SE_dict[key].keys())
        for i in range(len(activ_threshold_profile)):

            if activ_threshold_profile[i] in Pathways and not(activ_threshold_profile[i] in pathway_filtered):
                print(activ_threshold_profile[i], " had a low activation for training set, and was not added to ANN")
            elif not(activ_threshold_profile[i] in Pathways):
                print(activ_threshold_profile[i],
                      " was not in the training set. Perhaps, it is too far from the electrode")
            else:
                inx = pathway_filtered.index(activ_threshold_profile[i])
                if np.any(error_ANN[:, inx] > SE_err_threshold):
                    print('Error threshold for the side-effect implicated ', activ_threshold_profile[i],' was exceeded, the approximation model has to be revised')
                    print(activ_threshold_profile[i])
                    return False
                else:
                    # check number of errors > half of the threshold
                    N_half_errors = (error_ANN[:, inx] > SE_err_threshold / 2.0).sum()
                    # refuse if > 1%
                    if N_half_errors > 0.01 * error_ANN.shape[0]:
                        print(
                            '0.5 * error threshold for the side-effect implicated ', activ_threshold_profile[i],' was exceeded for more than 1% of tests, the approximation model has to be revised')
                        print(activ_threshold_profile[i])
                        print(max(error_ANN[:, inx]),SE_err_threshold)
                        return False

                if check_trivial == True:
                    if np.any(error_ANN_bi[:, inx] > SE_err_threshold) or np.any(error_ANN_mono[:, inx] > SE_err_threshold):
                        print(
                            'Error threshold for the side-effect implicated ', activ_threshold_profile[i],' was exceeded, the approximation model has to be revised')
                        return False


    # here we can merge target profiles for symptoms and threshold profiles for soft side-effects
    profile_dict.update(Soft_SE_dict)

    for key in profile_dict:
        if side == 0 and not ("_rh" in key):
            continue
        elif side == 1 and not ("_lh" in key):
            continue

        activ_threshold_profile = list(profile_dict[key].keys())
        for i in range(len(activ_threshold_profile)):

            if activ_threshold_profile[i] in Pathways and not (activ_threshold_profile[i] in pathway_filtered):
                print(activ_threshold_profile[i], " had a low activation for training set, and was not added to ANN")
            elif not (activ_threshold_profile[i] in Pathways):
                print(activ_threshold_profile[i],
                      " was not in the training set. Perhaps, it is too far from the electrode")
            else:
                inx = pathway_filtered.index(activ_threshold_profile[i])
                if np.any(error_ANN[:, inx] > Err_threshold):
                    print('Error threshold for ', activ_threshold_profile[i],' was exceeded, the model has to be revised')
                    return False
                else:
                    # check number of errors > half of the threshold
                    N_half_errors = (error_ANN[:, inx] > Err_threshold / 2.0).sum()
                    # refuse if > 1%
                    if N_half_errors > 0.01 * error_ANN.shape[0]:
                        print('0.5 * error threshold for ', activ_threshold_profile[i],' was exceeded for more than 1% of tests, the model has to be revised')
                        return False

                if check_trivial == True:
                    if np.any(error_ANN_bi[:, inx] > Err_threshold) or np.any(error_ANN_mono[:, inx] > Err_threshold):
                        print(
                            'Error threshold for ', activ_threshold_profile[i],' was exceeded, the approximation model has to be revised')
                        return False

    model.save(os.environ['STIMDIR'] + '/NB_' + str(side) + '/ANN_approved_model')
    return pathway_filtered

if __name__ == '__main__':

    # called from MATLAB
    # sys.argv[1] - stim folder
    # sys.argv[2] - side (0 - right hemisphere)

    os.environ['STIMDIR'] = sys.argv[1]
    side = int(sys.argv[2])

    if side == 0:
        res_folder = 'Results_rh/'
    else:
        res_folder = 'Results_lh/'

    TrainTest_currents_file = os.environ['STIMDIR'] + '/Current_protocols_' + str(side) + '.csv'  # current protocols
    if side == 0:
        TrainTest_activation_file = os.environ['STIMDIR'] + '/' + res_folder + 'Activations_over_StimSets_rh.csv'
    else:
        TrainTest_activation_file = os.environ['STIMDIR'] + '/' + res_folder + 'Activations_over_StimSets_lh.csv'

    # load parameters from .json folder generated in previous steps
    with open(os.environ['STIMDIR'] + '/netblend_dict_file.json', 'r') as fp:
        netblend_dict = json.load(fp)
    fp.close()
    netblend_dict = netblend_dict['netblendict']

    # load StimSets_parameters (were created by Train_Test_Generator.py)
    with open(os.environ['STIMDIR'] + '/NB_' + str(side) + '/StimSets_info.json', 'r') as fp:
        StimSets_info = json.load(fp)
    fp.close()

    # #better regenerate them
    # from Improvement4Protocol import create_NB_dictionaries
    # profile_dict, Soft_SE_dict, SE_dict = create_NB_dictionaries(side, FF_dictionary)

    approx_pathways = train_test_ANN(TrainTest_currents_file, TrainTest_activation_file, StimSets_info['trainSize_actual'],
                   netblend_dict['Err_threshold'], netblend_dict['SE_err_threshold'], side, check_trivial=False)

