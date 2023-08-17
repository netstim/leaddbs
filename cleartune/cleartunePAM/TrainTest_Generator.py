'''
    By K. Butenko
    Generate training and test datasets for OSS-DBS to train ANN
    Output: Current_protocols_ and StimSets_info.json in stim. folder
'''


import os
import sys
import numpy as np
from scipy.stats import qmc
import csv
import h5py
import json

one_pol_current_threshold = 8.0  # in mA
total_current_threshold = 8.0

def create_Training_Test_sets(stim_folder, Electrode_model, conc_threshold, segm_threshold, side):

    # check the electrode configuration
    from NB_outline import determine_el_type
    el_type = determine_el_type(Electrode_model)
    if el_type == 'concentric4':
        N_contacts = 4
        sample_size = 5000  # half training, half test
    else:
        N_contacts = 8
        sample_size = 10000

    # ## split to train and test and sample
    # # see "Test Set Sizing Via Random Matrix Theory" by A. Dubbs
    # # but this is for linear regression models!
    # n = N_contacts
    # m = sample_size
    # trainSize = round(m**(2/3) * (n*(2+n))**(1/3) - m**(1/3) * 2*n*(1+n) / (3*(n*(2+n))**(1/3)) + (1/3) * (6+n+n**2) - m**(-1/3)*2*n**(2)*(216 + 230*n + 87*n**(2) + 24*n**(3) + 5*n**(4))/(81*(n*(2+n))**(5/3)));
    # testSize = m - trainSize

    # Otherwise, just split in half
    trainSize = int(sample_size / 2)
    testSize = sample_size - trainSize

    # LHS sampling for training
    sampler = qmc.LatinHypercube(d=N_contacts)
    training_samples = sampler.random(n=trainSize)

    # Random sampling for test
    test_samples = np.random.rand(testSize, N_contacts)

    samples = np.concatenate((training_samples, test_samples), axis=0)

    # scale sample [0 1] samples to [threshold0, threshold1]
    if el_type == 'segmented8':
        samples[:, 0] = samples[:, 0] * (conc_threshold[1] - conc_threshold[0]) + conc_threshold[0]
        samples[:, 7] = samples[:, 7] * (conc_threshold[1] - conc_threshold[0]) + conc_threshold[0]
        samples[:, 1:7] = samples[:, 1:7] * (segm_threshold[1] - segm_threshold[0]) + segm_threshold[0]
    else:
        samples[:, :] = samples[:, :] * (conc_threshold[1] - conc_threshold[0]) + conc_threshold[0]

    # randomly nullify entries in a 25% of samples to marginalize
    import random
    for i in range(samples.shape[0]):
        if i % 4 == 0:
            if N_contacts == 4:
                N_null = int(round(random.uniform(1, 3)))
                C_list = [0, 1, 2, 3]

            elif N_contacts == 8:
                N_null = int(round(random.uniform(1, 6)))
                C_list = [0, 1, 2, 3, 4, 5, 6, 7]
            else:
                print("The electrode configuration was not recognized")
                raise SystemExit

            inx_null = random.sample(C_list, N_null)
            for j in inx_null:
                samples[i, j] = 0.0

            # increase currents below 0.5 mA
            if np.all(abs(samples[i, :]) < 0.5):
                samples[i, :] = samples[i, :] * 2

    # Create a .csv file native to OSS-DBS.
    with open(stim_folder + '/Current_protocols_' + str(side) + '.csv', 'w') as fd:
        writer = csv.writer(fd)
        if N_contacts == 8:
            writer.writerow(
                ['Contact0', 'Contact1', 'Contact2', 'Contact3', 'Contact4', 'Contact5', 'Contact6', 'Contact7'])
        elif N_contacts == 4:
            writer.writerow(
                ['Contact0', 'Contact1', 'Contact2', 'Contact3'])
        else:
            print('Check the electrode configuration')
            raise SystemExit

    # fill out the .csv file with current protocols (in mA!)
    # exclude samples that violate current sum thresholds
    trainSize_actual, testSize_actual = [0,0]
    for i in range(samples.shape[0]):
        if np.sum(samples[i,:]) < total_current_threshold and np.where(samples[i,:] < 0, samples[i,:], 0).sum(
                0) > -1*one_pol_current_threshold and np.where(samples[i,:] > 0, samples[i,:], 0).sum(0) < one_pol_current_threshold:

            stim_prot = samples[i].tolist()
            with open(stim_folder + '/Current_protocols_' + str(side) + '.csv', 'a') as fd:
                writer = csv.writer(fd)
                writer.writerow(stim_prot)

            if i < trainSize:
                trainSize_actual += 1
            else:
                testSize_actual += 1

    # create a json that describes Current protocols
    StimSets_info = {
        'trainSize_actual': trainSize_actual,
        'testSize_actual': testSize_actual,
        'el_type': el_type,
        'conc_threshold': conc_threshold,
        'segm_threshold': segm_threshold,
        }

    if not os.path.exists(stim_folder + '/NB_' + str(side)):
        os.mkdir(stim_folder + '/NB_' + str(side))

    with open(stim_folder + '/NB_' + str(side) + '/StimSets_info.json', 'w') as save_as_dict:
        json.dump(StimSets_info, save_as_dict)

    return trainSize_actual, testSize_actual


if __name__ == '__main__':

    # called from MATLAB

    # passed from Currentune
    # sys.argv[1] - stimfolder
    # sys.argv[2] - electrode model (-1 if not implanted)
    # sys.argv[3] - side
    # sys.argv[4:] - min cylind, max cylind, min segm, max_segm

    create_Training_Test_sets(sys.argv[1], sys.argv[2], [float(sys.argv[4]), float(sys.argv[5])],
                              [float(sys.argv[6]), float(sys.argv[7])], sys.argv[3])

    # if sys.argv[2] != '-1':
    #     create_Training_Test_sets(sys.argv[1],sys.argv[2], [float(sys.argv[4]),float(sys.argv[5])], [float(sys.argv[6]),float(sys.argv[7])], sys.argv[3])
    #
    # if sys.argv[3] != '-1':
    #     create_Training_Test_sets(sys.argv[1], sys.argv[2], [float(sys.argv[4]), float(sys.argv[5])],
    #                               [float(sys.argv[6]), float(sys.argv[7])], sys.argv[3])
