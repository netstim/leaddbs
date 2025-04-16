'''
    By K. Butenko
'''


import os
import sys
import numpy as np
from scipy.stats import qmc
import csv
import h5py
import json
from scipy.optimize import minimize

# hardwired: max total currents allowed
one_pol_current_threshold = 8.0  # in mA
total_current_threshold = 8.0
abs_current_threshold = 12.0


def l1_pos_neg_max(array):
    """
    Computes the L1 norm of positive and negative values for each row of a 2D NumPy array,
    and returns a vector containing the maximum of the two norms for each row.

    Args:
        array (numpy.ndarray): A 2D NumPy array.

    Returns:
        numpy.ndarray: A 1D NumPy array where each element is the maximum of the L1 norm
                       of the positive and negative values in the corresponding row of the input array.
    """
    positive_values = np.clip(array, 0, np.inf)  # Get positive values, set negative to 0
    negative_values = np.clip(array, -np.inf, 0) # Get negative values, set positive to 0

    l1_positive = np.sum(positive_values, axis=1)  # L1 norm of positive values for each row
    l1_negative = np.sum(np.abs(negative_values), axis=1)  # L1 norm of negative values for each row

    max_l1 = np.maximum(l1_positive, l1_negative)  # Get the maximum of the two L1 norms for each row
    return max_l1

def scale_array_to_l1_norm(array, target_norm, norm_type):
    """
    Scales an array of vectors (where each vector is a row) to have a specified L1 norm.

    Args:
        array (numpy.ndarray): The array to scale.  Each row represents a vector.
        target_norm (float): The desired L1 norm for each vector.
        norm_type (string): type of norm to check

    Returns:
        numpy.ndarray: The array with scaled vectors.
    """

    if norm_type == 'L1':
        original_norms = np.sum(np.abs(array), axis=1)  # L1 norms of each row
    elif norm_type == 'L1_sign':
        original_norms = np.sum(array, axis=1)  # L1 norms of each row
    elif norm_type == 'L1_polarity':
        original_norms = l1_pos_neg_max(array)
        
    for vector_i in range(array.shape[0]):
    
        if original_norms[vector_i] < target_norm:
            continue
        else:
            scaling_factor = target_norm / original_norms[vector_i]
            #if norm_type == 'L1_polarity':
            #    print(scaling_factor, array[vector_i,:],original_norms[vector_i])
            array[vector_i,:] = array[vector_i,:] * scaling_factor

    return array


def create_Training_Test_sets(stim_folder, Electrode_model, conc_threshold, segm_threshold, side):

    """ Generate current sets to solve for training (using LHS) and testing (random) of the approximation model

    Inputs
    ------
    stim_folder : str
        path to the stimulation folder

    Returns
    -------
    trainSize_actual, int, number of protocols for training
    testSize_actual, int, number of protocols for testing
    """

    if side == 0:
        side_suffix = '_rh'
    else:
        side_suffix = '_lh'

    # check the electrode configuration
    from NB_outline import determine_el_type
    el_type = determine_el_type(Electrode_model)
    if el_type == 'concentric4':
        N_contacts = 4
        sample_size = 8000  # half training, half test
    else:
        N_contacts = 8
        sample_size = 16000

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

    # downscale (if necessary to abide current bounds)
    samples = scale_array_to_l1_norm(samples, abs_current_threshold,'L1')
    samples = scale_array_to_l1_norm(samples, total_current_threshold,'L1_sign')
    samples = scale_array_to_l1_norm(samples, one_pol_current_threshold,'L1_polarity')

    # randomly nullify entries in a 25% of samples to marginalize
    import random
    for i in range(samples.shape[0]):
        if i % 4 == 0:
            # set a random number (1-3 or 1-6) of contacts to 0 mA
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

            # double if all currents below 0.5 mA
            if np.all(abs(samples[i, :]) < 0.5):
                samples[i, :] = samples[i, :] * 2


    if not os.path.exists(os.path.join(stim_folder,'NB' + side_suffix)):
        os.mkdir(os.path.join(stim_folder,'NB' + side_suffix))

    # in this case, store Current_protocols in NB folder
    current_protocols_path = os.path.join(stim_folder, 'NB' + side_suffix, 'Current_protocols_' + str(side) + '.csv')
    with open(current_protocols_path, 'w') as fd:
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
        if np.round(np.sum(np.abs(samples[i,:])),6) <= abs_current_threshold \
            and np.sum(samples[i,:]) <= total_current_threshold \
            and np.round(np.where(samples[i,:] < 0, samples[i,:], 0).sum(0),6) >= -1*one_pol_current_threshold \
            and np.round(np.where(samples[i,:] > 0, samples[i,:], 0).sum(0),6) <= one_pol_current_threshold:
                
        #if np.sum(samples[i,:]) < total_current_threshold and np.where(samples[i,:] < 0, samples[i,:], 0).sum(
        #        0) > -1*one_pol_current_threshold and np.where(samples[i,:] > 0, samples[i,:], 0).sum(0) < one_pol_current_threshold:

            stim_prot = samples[i].tolist()
            with open(current_protocols_path, 'a') as fd:
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

    with open(os.path.join(stim_folder,'NB' + side_suffix,'StimSets_info.json'), 'w') as save_as_dict:
        json.dump(StimSets_info, save_as_dict)

    return trainSize_actual, testSize_actual


if __name__ == '__main__':

    # called from MATLAB

    # passed from Currentune
    # sys.argv[1] - stimfolder
    # sys.argv[2] - electrode model (-1 if not implanted)
    # sys.argv[3] - side (0-rh)
    # sys.argv[4:] - min cylind, max cylind, min segm, max_segm

    create_Training_Test_sets(sys.argv[1], sys.argv[2], [float(sys.argv[4]), float(sys.argv[5])],
                              [float(sys.argv[6]), float(sys.argv[7])], int(sys.argv[3]))

    # if sys.argv[2] != '-1':
    #     create_Training_Test_sets(sys.argv[1],sys.argv[2], [float(sys.argv[4]),float(sys.argv[5])], [float(sys.argv[6]),float(sys.argv[7])], sys.argv[3])
    #
    # if sys.argv[3] != '-1':
    #     create_Training_Test_sets(sys.argv[1], sys.argv[2], [float(sys.argv[4]), float(sys.argv[5])],
    #                               [float(sys.argv[6]), float(sys.argv[7])], sys.argv[3])
