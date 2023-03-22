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


### Input

# some ANN parameters
learn_rate = 0.0025
N_epochs = 2000
min_activ_threshold = 0.1

def train_test_ANN(TrainTest_currents_file, TrainTest_activation_file, trainSize, Err_threshold, SE_err_threshold, side, check_trivial):

    import os
    #print(os.environ['STIMDIR'])

    # load currents
    Currents = np.genfromtxt(TrainTest_currents_file, delimiter=',', skip_header=True)

    # load activation results
    ActivationResults = np.genfromtxt(TrainTest_activation_file, delimiter=' ')

    # we can also estimate errors for some simple pre-defined protocols
    if check_trivial == True:
        ActivationResults_Bipolar1 = np.genfromtxt('Activations_over_iterationsBipolar.csv', delimiter=' ')
        Currents_Bipolar1 = np.genfromtxt('/home/cerebellum/Documents/Data/NetBlend/Current_protocols_0_Bipolar.csv', delimiter=',', skip_header=True)

        ActivationResults_Monopolar21 = np.genfromtxt('Activations_over_iterations_Monopolar21_79.csv', delimiter=' ')
        Currents_Monopolar21 = np.genfromtxt('Current_protocols_0_Monopolars21_79.csv', delimiter=',', skip_header=True)

    """To be removed"""

    # #ActivationResults = np.genfromtxt('/home/konstantin/Downloads/Activations_over_iterations_TrainTest.csv', delimiter=' ')
    # #ActivationResults = np.genfromtxt('/home/cerebellum/Documents/Data/NetBlend/Activations_over_iterations7287.csv', delimiter=' ')
    # #Currents = np.genfromtxt('/home/konstantin/Downloads/Current_protocols_0.csv', delimiter=',', skip_header=True)
    # Currents = np.genfromtxt('/home/cerebellum/Documents/Data/NetBlend/Current_protocols_0.csv', delimiter=',', skip_header=True)
    # Currents = Currents[:7287,:]
    #
    # # additionally, test bipolar and monopolar settings
    # ActivationResults_Bipolar1 = np.genfromtxt('/home/cerebellum/Documents/Data/NetBlend/Activations_over_iterationsBipolar1.csv', delimiter=' ')
    # Currents_Bipolar1 = np.genfromtxt('/home/cerebellum/Documents/Data/NetBlend/Current_protocols_0_Bipolar1.csv', delimiter=',', skip_header=True)
    # ActivationResults_Monopolar21 = np.genfromtxt('/home/cerebellum/Documents/Data/NetBlend/Activations_over_iterations_Monopolars_2_1.80.csv', delimiter=' ')
    # Currents_Monopolar21 = np.genfromtxt('/home/cerebellum/Documents/Data/NetBlend/Current_protocols_0_Monopolars_2_1.csv', delimiter=',', skip_header=True)



    # this is defined by the PAM model
    # the function will work only for a proper Lead-DBS import (connectome folder, oss-dbs_parameters.mat)
    from Pathways_Stats import get_simulated_pathways
    Pathways, axons_in_path = get_simulated_pathways(side)

    # axons_in_path = np.array([250,250,250,250,250,100,100,100,100,100,100,100,100,250,500,500,250,250,250,250,500,500,250,250])
    # Pathways = ['ACC_cp_right',
    # 'M1_cf_face_right',
    # 'M1_cf_lowerex_right',
    # 'M1_cf_upperex_right',
    # 'PreMotor_cf_right',
    # 'R_ACC_hdp_right',
    # 'R_M1_hdp_face_right',
    # 'R_M1_hdp_lowerex_right',
    # 'R_M1_hdp_upperex_right',
    # 'R_PreMotor_hdp_right',
    # 'R_SMA_hdp_right',
    # 'R_dlPFC_hdp_right',
    # 'R_vmPFC_hdp_right',
    # 'SMA_cf_right',
    # 'ansa_lenticularis_right',
    # 'cerebellothalamic_right',
    # 'dlPFC_cp_right',
    # 'dmPFC_cp_right',
    # 'gpe2stn_ass_right',
    # 'gpe2stn_sm_right',
    # 'lenticular_fasciculus_right',
    # 'medial_lemniscus_right',
    # 'vlPFC_cp_right',
    # 'vmPFC_cp_right']


    ### Prepare the data
    X_train = Currents[:trainSize,:]
    X_test = Currents[trainSize:,:]
    y_train_prelim = ActivationResults[:trainSize,1:] / axons_in_path # from 1, because 0 is the index of the protocol
    y_test_prelim = ActivationResults[trainSize:,1:] / axons_in_path
    y_train = -100 * np.ones((y_train_prelim.shape), float)
    y_test = -100 * np.ones((y_test_prelim.shape), float)

    # optionally check some hardcoded trivial protocols
    if check_trivial == True:

        X_bipolar = Currents_Bipolar1
        X_monopolar = Currents_Monopolar21

        y_bipolar_prelim = ActivationResults_Bipolar1[:,1:] / axons_in_path
        y_monopolar_prelim = ActivationResults_Monopolar21[:,1:] / axons_in_path
        y_bipolar = -100 * np.ones((y_bipolar_prelim.shape), float)
        y_monopolar = -100 * np.ones((y_monopolar_prelim.shape), float)

    # only consider those pathways, where max activation > 5%
    pathway_filtered = []

    for i in range(y_train_prelim.shape[1]):
        # only compute for pathways with some percent activation
        if np.max(y_train_prelim[:, i]) >= min_activ_threshold and np.max(y_test_prelim[:, i]) >= min_activ_threshold:
            y_train[:,i] = y_train_prelim[:,i]
            y_test[:,i] = y_test_prelim[:,i]
            pathway_filtered.append(Pathways[i])

            if check_trivial == True:
                y_bipolar[:,i] = y_bipolar_prelim[:,i]
                y_monopolar[:, i] = y_monopolar_prelim[:, i]

    y_train = y_train[:, ~np.all(y_train == -100.0, axis=0)]
    y_test = y_test[:, ~np.all(y_test == -100.0, axis=0)]

    if check_trivial == True:
        y_bipolar = y_bipolar[:,~np.all(y_bipolar == -100.0, axis=0)]
        y_monopolar = y_monopolar[:,~np.all(y_monopolar == -100.0, axis=0)]

    # inject 10% of zero protocols
    N_zero = int(0.1 * y_train.shape[0])
    y_train_exp = np.zeros((y_train.shape[0] + N_zero,y_train.shape[1]), float)
    X_train_exp = np.zeros((X_train.shape[0] + N_zero,X_train.shape[1]), float)
    y_train_exp[:y_train.shape[0],:] = y_train
    X_train_exp[:X_train.shape[0],:] = X_train
    y_train = y_train_exp
    X_train = X_train_exp

    ### train ANN
    # normalization seems to be not necessary here

    model = Sequential(layers=None, name=None)
    model.add(Dense(128, input_shape = (X_train.shape[1],), activation='linear'))
    model.add(Dense(1024, activation='hard_sigmoid')) # actually abs(LeakyReLU) with alpha 1.25 (steeper slope for cathode)
    model.add(Dense(y_train.shape[1], activation='tanh')) # we need 0 -> 0 (negative vals are removed on the previous level)

    ## sigmoid produces a shift for monopolar and bipolar
    #model.add(Dense(y_train.shape[1], activation='sigmoid')) # we need 0 -> 0 (negative vals are removed on the previous level)
                                                                            #
    adam = optimizers.Adamax(lr = learn_rate)
    model.compile(optimizer = adam, loss = 'mean_squared_error', metrics=['accuracy'])
    model.fit(X_train, y_train, epochs = N_epochs, verbose = 1)
    results = model.evaluate(X_test, y_test)

    # on Test
    error_ANN = np.zeros(y_test.shape,float)
    y_predicted = model.predict(X_test)
    error_ANN = y_test - y_predicted

    if check_trivial == True:

        error_ANN_bi = np.zeros(y_bipolar.shape,float)
        y_predicted_bi = model.predict(X_bipolar)
        error_ANN_bi = y_bipolar - y_predicted_bi
        results_bi = model.evaluate(X_bipolar, y_bipolar)

        results_mono = model.evaluate(X_monopolar, y_monopolar)
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

    # ========================================= Plot the errors =======================================================#

    import matplotlib
    matplotlib.rcParams['figure.dpi'] = 200

    plt.figure()
    for i in range(len(pathway_filtered)):

        if np.max(abs(error_ANN[:, i])) > 0.05:
            sns.kdeplot(error_ANN[:,i], bw_adjust=0.5, label=pathway_filtered[i])

    plt.legend()
    plt.title('Abs errors for ANN on Test')
    #plt.xlim([-0.15,0.15])
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
                    print('Error threshold for the side-effect implicated pathway activation was exceeded, the approximation model has to be revised')
                    print(activ_threshold_profile[i])
                    return False
                else:
                    # check number of errors > half of the threshold
                    N_half_errors = (error_ANN[:, inx] > SE_err_threshold / 2.0).sum()
                    # refuse if > 1%
                    if N_half_errors > 0.01 * error_ANN.shape[0]:
                        print(
                            '0.5 * error threshold for the side-effect implicated pathway activation was exceeded for more than 1% of tests, the approximation model has to be revised')
                        print(activ_threshold_profile[i])
                        print(max(error_ANN[:, inx]),SE_err_threshold)
                        return False

                if check_trivial == True:
                    if np.any(error_ANN_bi[:, inx] > SE_err_threshold) or np.any(error_ANN_mono[:, inx] > SE_err_threshold):
                        print(
                            'Error threshold for the side-effect implicated pathway activation was exceeded, the approximation model has to be revised')
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
                    print('Error threshold for this pathway activation was exceeded, the model has to be revised')
                    return False
                else:
                    # check number of errors > half of the threshold
                    N_half_errors = (error_ANN[:, inx] > Err_threshold / 2.0).sum()
                    # refuse if > 1%
                    if N_half_errors > 0.01 * error_ANN.shape[0]:
                        print('0.5 * error threshold for this pathway activation was exceeded for more than 1% of tests, the model has to be revised')
                        return False

                if check_trivial == True:
                    if np.any(error_ANN_bi[:, inx] > Err_threshold) or np.any(error_ANN_mono[:, inx] > Err_threshold):
                        print(
                            'Error threshold for the side-effect implicated pathway activation was exceeded, the approximation model has to be revised')
                        return False

    model.save(os.environ['STIMDIR'] + '/NB_' + str(side) + '/ANN_approved_model')
    return pathway_filtered