'''
    By K. Butenko
    This script trains and test an ANN model to approximate pathway activation for a given electrode position
'''

import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'   # to avoid any CUDA issues
import sys
import json
import copy
from typing import Tuple, List, Optional, Dict
from scipy.stats import pearsonr

import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, BatchNormalization, Lambda
from tensorflow.keras import optimizers


SIDE_SUFFIX = ['_rh','_lh']
SIDE_LABEL = {0: '_right', 1: '_left'}

# ANN parameters
learn_rate = 0.0025
N_epochs = 500
MIN_AXON_NUMBER = 10
MIN_ACTIV_THRESHOLD = 0.05  # Define a default value if not specified elsewhere

def matplotlib_kdeplot(data, ax=None, bw_method=None, **kwargs):
    """
    Replicates seaborn's kdeplot functionality using matplotlib and scipy.stats.gaussian_kde.

    Parameters:
    - data: 1D array-like data.
    - ax: matplotlib Axes object, if None, a new figure and axes will be created.
    - bw_method: The method used to calculate the estimator bandwidth. Can be 'scott', 'silverman', or a scalar.
    - **kwargs: Additional keyword arguments passed to matplotlib's plot function.

    Returns:
    - matplotlib Axes object.
    """

    if ax is None:
        fig, ax = plt.subplots()

    kde = gaussian_kde(data, bw_method=bw_method)
    x_vals = np.linspace(min(data), max(data), 1000)  # Adjust range and resolution as needed
    y_vals = kde(x_vals)

    ax.plot(x_vals, y_vals, **kwargs)

    return ax


def load_activation_results(res_folder: str, pathway: str, currents_shape: Tuple[int, int]) -> np.ndarray:
    """Loads activation results for each pathway from JSON files.

    Args:
        res_folder: Path to the results subfolder.
        pathway: Name of the pathway for the ANN training.
        currents_shape: Shape of the currents array (number of protocols, number of contacts).

    Returns:
        A NumPy array containing activation results (percentage activated) for each protocol and pathway.
    """
    activation_results = np.zeros((currents_shape[0],1), dtype=float)
    for prot_i in range(currents_shape[0]):
        pathway_file = os.path.join(res_folder, f'Pathway_status_{pathway}_{prot_i}.json')
        if os.path.isfile(pathway_file):
            with open(pathway_file, 'r') as fp:
                pam_res_dict = json.load(fp)
                activation_results[prot_i,0] = 0.01 * pam_res_dict['percent_activated']
                
    return activation_results


def filter_pathways(pathway: str, axons_in_path: np.ndarray, y_train_all: np.ndarray, y_test_all: np.ndarray, no_test=False) -> Tuple[np.ndarray, np.ndarray, List[str], Optional[np.ndarray], Optional[np.ndarray]]:
    """Filters pathways based on minimum axon number and activation threshold.

    Args:
        pathway: Name of the pathway for the ANN training.
        axons_in_path: NumPy array containing the number of axons in each pathway.
        y_train_all: All training activation results.
        y_test_all: All testing activation results.
    Returns:
        A tuple containing filtered training activation results, filtered testing activation results,
        a list of filtered pathway names, and optionally filtered bipolar and monopolar activation results.
    """
    y_train = -100 * np.ones((y_train_all.shape), dtype=float)
    y_test = -100 * np.ones((y_test_all.shape), dtype=float)

    pathway_filtered = []
    for i in range(y_train_all.shape[1]):
        if axons_in_path[i] >= MIN_AXON_NUMBER and np.max(y_train_all[:, i]) >= MIN_ACTIV_THRESHOLD and (no_test or np.max(y_test_all[:, i]) >= MIN_ACTIV_THRESHOLD):
            y_train[:, i] = y_train_all[:, i]
            y_test[:, i] = y_test_all[:, i]
            pathway_filtered.append(pathway)

    y_train_filtered = y_train[:, ~np.all(y_train == -100.0, axis=0)]
    if no_test:
        y_test_filtered = None
    else:
        y_test_filtered = y_test[:, ~np.all(y_test == -100.0, axis=0)]

    return y_train_filtered, y_test_filtered, pathway_filtered


def augment_training_data_with_zeros(X_train: np.ndarray, y_train: np.ndarray, zero_protocol_percentage: float = 0.1) -> Tuple[np.ndarray, np.ndarray]:
    """Injects a percentage of zero protocols into the training data.

    Args:
        X_train: Training currents.
        y_train: Training activation results.
        zero_protocol_percentage: Percentage of zero protocols to inject (default: 0.1).

    Returns:
        A tuple containing augmented training currents and augmented training activation results.
    """
    n_zero = int(zero_protocol_percentage * y_train.shape[0])
    y_train_ext = np.zeros((y_train.shape[0] + n_zero, y_train.shape[1]), dtype=float)
    X_train_ext = np.zeros((X_train.shape[0] + n_zero, X_train.shape[1]), dtype=float)
    y_train_ext[:y_train.shape[0], :] = y_train
    X_train_ext[:X_train.shape[0], :] = X_train
    return X_train_ext, y_train_ext

def analyze_and_plot_ann_errors(
    stim_dir: str,
    side: int,
    pathway_filtered: List[str],
    X_test: np.ndarray,
    error_ann: np.ndarray,
    error_ann_bi: Optional[np.ndarray] = None,
    error_ann_mono: Optional[np.ndarray] = None,
    check_trivial: bool = False
):
    """
    Analyzes and plots the absolute errors of an Artificial Neural Network (ANN)
    on test data, correlating them with total current and absolute total current.
    Also generates kernel density estimate (KDE) plots of the errors for each
    filtered pathway using scipy.stats.gaussian_kde. Optionally plots errors
    for bipolar and monopolar protocols.

    Args:
        stim_dir: Path to the stimulation directory.
        side: Hemisphere index (0 - right, 1 - left).
        pathway_filtered: List of filtered pathway names.
        X_test: NumPy array of the test current protocols.
        error_ann: NumPy array of the ANN errors on the test data.
        error_ann_bi: Optional NumPy array of ANN errors on bipolar protocols.
        error_ann_mono: Optional NumPy array of ANN errors on monopolar protocols.
        check_trivial: Boolean flag indicating whether to plot trivial protocol errors.
    """
    matplotlib.rcParams['figure.dpi'] = 200
    nb_side_folder = os.path.join(stim_dir, 'NB' + SIDE_SUFFIX[side])

    # --- Plotting errors vs. total current ---
    plt.figure()
    total_current = np.sum(X_test, axis=1)
    pearson_r, pearson_p = pearsonr(total_current, np.abs(error_ann[:, 0]))
    mae = np.mean(np.abs(error_ann[:, 0]))
    plt.scatter(total_current*1000.0, np.abs(error_ann[:, 0]))
    plt.text(0.05, 0.90, f"Pearson's R = {pearson_r:.2f}, p = {pearson_p:.5f}",
             transform=plt.gca().transAxes, fontsize=10)
    plt.text(0.05, 0.80, f"MAE = {mae:.2f}",
             transform=plt.gca().transAxes, fontsize=10)
    plt.xlabel("Total current, mA")
    plt.ylabel("Absolute Error of Percent Activation")
    filename_total_current = f"{pathway_filtered[0]}_ANN_abs_errors_on_Test_versus_total_current{SIDE_LABEL[side]}.png"
    plt.savefig(os.path.join(nb_side_folder, filename_total_current), format='png', dpi=500)

    # --- Plotting errors vs. absolute total current ---
    plt.figure()
    total_abs_current = np.sum(np.abs(X_test), axis=1)
    pearson_r2, pearson_p2 = pearsonr(total_abs_current, np.abs(error_ann[:, 0]))
    plt.scatter(total_abs_current*1000.0, np.abs(error_ann[:, 0]))
    plt.text(0.05, 0.85, f"Pearson's R = {pearson_r2:.2f}, p = {pearson_p2:.5f}",
             transform=plt.gca().transAxes, fontsize=10)
    plt.text(0.05, 0.75, f"MAE = {mae:.2f}",
             transform=plt.gca().transAxes, fontsize=10)
    plt.xlabel("Absolute current, mA")
    plt.ylabel("Absolute Error of Percent Activation")
    filename_abs_total_current = f"{pathway_filtered[0]}_ANN_abs_errors_on_Test_versus_abs_current{SIDE_LABEL[side]}.png"
    plt.savefig(os.path.join(nb_side_folder, filename_abs_total_current), format='png', dpi=500)

    # --- Plotting KDE of absolute errors for each pathway ---
    plt.figure()
    pathways_max_errors: Dict[str, any] = {}
    for i, pathway in enumerate(pathway_filtered):
        max_error = np.max(np.abs(error_ann[:, i]))
        pathways_max_errors[pathway] = max_error
        inx_max_error = np.argmax(np.abs(error_ann[:, i]))
        pathways_max_errors['stim_protocol'] = X_test[inx_max_error, :].tolist()  # Store as list for JSON
        if max_error > 0.05:
            kde = gaussian_kde(error_ann[:, i])
            x_vals = np.linspace(min(error_ann[:, i]), max(error_ann[:, i]), 200) # Generate x values for plotting
            plt.plot(x_vals, kde(x_vals), label=pathway)
            
            plt.text(0.05, 0.80, f"MAE = {mae:.2f}",
                     transform=plt.gca().transAxes, fontsize=10)

    error_filename_base = f"{pathway_filtered[0]}_ANN_abs_errors"
    with open(os.path.join(nb_side_folder, f"{error_filename_base}.json"), 'w') as save_as_dict:
        json.dump(pathways_max_errors, save_as_dict, indent=4)  # Added indent for readability

    plt.legend()
    plt.title('Absolute Errors for ANN on Test')
    plt.xlim([-25.0, 25.0])
    plt.xlabel("Absolute Error of Percent Activation")
    filename_kde_test = f"{error_filename_base}_on_Test{SIDE_LABEL[side]}.png"
    plt.savefig(os.path.join(nb_side_folder, filename_kde_test), format='png', dpi=1000)

    # --- Plotting KDE of absolute errors for bipolar protocols (if check_trivial is True) ---
    if check_trivial and error_ann_bi is not None:
        plt.figure()
        for i, pathway in enumerate(pathway_filtered):
            if np.max(np.abs(error_ann_bi[:, i])) > 0.01:
                
                mae_bi = np.mean(np.abs(error_ann_bi[:, 0]))
                
                kde = gaussian_kde(error_ann_bi[:, i])
                x_vals = np.linspace(min(error_ann_bi[:, i]), max(error_ann_bi[:, i]), 200)
                plt.plot(x_vals, kde(x_vals), label=pathway)
                
                plt.text(0.05, 0.80, f"MAE = {mae_bi:.2f}",
                         transform=plt.gca().transAxes, fontsize=10)
        plt.legend()
        plt.title('Absolute Errors for ANN on Bipolar')
        plt.xlim(min(error_ann_bi[:, i]), max(error_ann_bi[:, i]))
        plt.xlabel("Absolute Error of Percent Activation")
        filename_kde_bipolar = f"{error_filename_base}_on_Bipolar{SIDE_LABEL[side]}.png"
        plt.savefig(os.path.join(nb_side_folder, filename_kde_bipolar), format='png', dpi=500)

    # --- Plotting KDE of absolute errors for monopolar protocols (if check_trivial is True) ---
    if check_trivial and error_ann_mono is not None:
        plt.figure()
        for i, pathway in enumerate(pathway_filtered):
            if np.max(np.abs(error_ann_mono[:, i])) > 0.01:
                mae_mono = np.mean(np.abs(error_ann_mono[:, 0]))
                
                plt.text(0.05, 0.80, f"MAE = {mae_mono:.2f}",
                         transform=plt.gca().transAxes, fontsize=10)
                
                kde = gaussian_kde(error_ann_mono[:, i])
                x_vals = np.linspace(min(error_ann_mono[:, i]), max(error_ann_mono[:, i]), 200)
                plt.plot(x_vals, kde(x_vals), label=pathway)
        plt.legend()
        plt.title('Absolute Errors for ANN on Monopolar')
        plt.xlim(min(error_ann_mono[:, i]), max(error_ann_mono[:, i]))
        plt.xlabel("Absolute Error of Percent Activation")
        filename_kde_monopolar = f"{error_filename_base}_on_Monopolar{SIDE_LABEL[side]}.png"
        plt.savefig(os.path.join(nb_side_folder, filename_kde_monopolar), format='png', dpi=500)


def check_pathway_errors(pathway, side, profiles, error_array, error_threshold, pathway_filtered, check_trivial, error_array_bi=None, error_array_mono=None):
    """
    Check error thresholds for a given pathway.
    """
    
    error_threshold = error_threshold * 100.0
    
    if pathway not in pathways:
        print(f"{pathway} was not in the training set. Perhaps, it is too far from the electrode")
        return True  # Continue to the next pathway

    if pathway not in pathway_filtered:
        print(f"{pathway} had a low activation for training set, and was not added to ANN")
        return True  # Continue to the next pathway

    inx = pathway_filtered.index(pathway)

    if np.any(error_array[:, inx] > error_threshold):
        print(f'Error threshold for {pathway} was exceeded, the model has to be revised')
        return False

    n_half_errors = (error_array[:, inx] > error_threshold / 2.0).sum()
    if n_half_errors > 0.01 * error_array.shape[0]:
        print(f'0.5 * error threshold for {pathway} was exceeded for more than 1% of tests, the model has to be revised')
        print(f"{pathway}")
        print(f"{np.max(error_array[:, inx])}, {error_threshold}")
        return False

    if check_trivial and (error_array_bi is not None and np.any(error_array_bi[:, inx] > error_threshold)) or \
                       (error_array_mono is not None and np.any(error_array_mono[:, inx] > error_threshold)):
        print(f'Error threshold for {pathway} was exceeded (trivial check), the approximation model has to be revised')
        return False

    return True

def train_test_ANN(stim_dir: str, side: int, pathway: str, check_trivial: bool, vat_recruit: bool = False):
    """Trains an ANN model and tests it.

    Parameters:
        stim_dir: Path to the lead-dbs stimulation folder.
        side: Hemisphere index (0 - right, 1 - left).
        pathway: Name of the pathway for training.
        check_trivial: If True, also train and test on predefined monopolar and bipolar protocols.
        vat_recruit: If True, train ANN for VATs.
    """
    
    res_folder = os.path.join(stim_dir, 'OSS_sim_files' + SIDE_SUFFIX[side], 'Results')

    # load parameters from .json folder generated in previous steps
    with open(stim_dir + '/netblend_dict_file.json', 'r') as fp:
        netblend_dict = json.load(fp)
    fp.close()
    netblend_dict = netblend_dict['netblendict']
    err_threshold = netblend_dict['Err_threshold']
    SE_err_threshold = netblend_dict['SE_err_threshold']

    # load StimSets_parameters (were created by Train_Test_Generator.py)
    with open(stim_dir + '/NB' + SIDE_SUFFIX[side] + '/StimSets_info.json', 'r') as fp:
        StimSets_info = json.load(fp)
    fp.close()
    train_size = StimSets_info['trainSize_actual']
    no_test = StimSets_info['testSize_actual'] == 0

    #================================================== Import Training and Test Data =====================================================#

    # Load currents for training and test protocols
    train_test_currents_file = os.path.join(stim_dir, 'OSS_sim_files' + SIDE_SUFFIX[side], 'Current_protocols_' + str(side) + '.csv')
    currents = np.genfromtxt(train_test_currents_file, delimiter=',', skip_header=True)

    # Convert currents to A
    X_train = currents[:train_size, :] * 0.001
    X_test = currents[train_size:, :] * 0.001

    # Load PAM results
    if vat_recruit:
        print("TBD")
        raise SystemExit()
    else:
        # Determine pathways and load activation results
        from Pathways_Stats import get_simulated_pathways
        pathways, axons_in_path = get_simulated_pathways(side, stim_dir)
        activation_results = load_activation_results(res_folder, pathway, currents.shape)

    # plot activation across the protocols 
    import matplotlib
    matplotlib.rcParams['figure.dpi'] = 200
    
    plt.figure()
    plt.scatter(np.sum(currents,axis=1),activation_results[:,0]*100.0) # Add label for legend
    plt.xlabel("Absolute current, mA")
    plt.ylabel("Percent Activation")
    plt.savefig(os.path.join(stim_dir,'NB' + SIDE_SUFFIX[side],pathway+'_Percent_Activation_versus_total_current' + SIDE_LABEL[side] + '.png'), format='png',
                dpi=500)

    # Filter pathways based on activity and axon number    
    y_train_filtered, y_test_filtered, pathway_filtered = filter_pathways(
        pathway, axons_in_path, activation_results[:train_size, :], activation_results[train_size:, :], no_test=no_test
    )
    
    if not pathway_filtered:
        print("Low activation levels for ", pathway)
        return False
    if not no_test:
        y_test = y_test_filtered 

    # Augment training data with zero protocols
    X_train_augmented, y_train_augmented = augment_training_data_with_zeros(X_train, y_train_filtered)

    if y_train_augmented.shape[1] > 1:
        print("Multiple pathways detected after filtering. Ensure the logic for selecting a single pathway is correct.")
        raise SystemExit
        # Consider if this should be an error or if the code should handle multiple pathways

    # inject monopolar review
    if not no_test:
        mono_protocols = X_test[X_test.shape[0]-56:,:]
        mono_solutions = y_test[y_test.shape[0]-56:,:]
    else:
        mono_protocols = X_train[X_train.shape[0]-56:,:]
        mono_solutions = y_train_filtered[y_train_filtered.shape[0]-56:,:]
    
    X_train_augmented = np.concatenate((X_train_augmented,mono_protocols,mono_protocols,mono_protocols,mono_protocols,mono_protocols))
    y_train_augmented = np.concatenate((y_train_augmented,mono_solutions,mono_solutions,mono_solutions,mono_solutions,mono_solutions))

    print(f"Training on pathway: {pathway}")
    print(f"Number of training samples: {X_train.shape[0]}")
    print(f"Number of testing samples: {X_test.shape[0]}")
    print(f"Filtered pathways: {pathway_filtered}")

    #================================================== Train ANN =====================================================#

    model = Sequential(layers=None, name=None)
    model.add(Dense(128, input_shape=(X_train_augmented.shape[1],), activation='linear'))
    model.add(Dense(1024, activation=tf.keras.layers.LeakyReLU(alpha=-1.25)))  # alpha -1.25 to have a steeper slope for cathode
    #model.add(Dropout(0.2))
    model.add(Dense(np.sum(axons_in_path), activation='sigmoid'))  # following the percent activation curves
    #model.add(Dense(y_train.shape[1], activation='tanh'))
    model.add(Dense(y_train_augmented.shape[1], activation='sigmoid'))

    adam = optimizers.Adamax(lr=learn_rate)
    model.compile(optimizer=adam, loss='mean_squared_error', metrics=['mse'])
    model.fit(X_train_augmented, y_train_augmented, epochs=N_epochs, verbose=1)
    
    if no_test:
        pathway_to_save = pathway_filtered[-1] if pathway_filtered else "default_pathway" # Save based on the final filtered pathway?
        save_path = os.path.join(stim_dir, 'NB' + SIDE_SUFFIX[side], f'ANN_approved_model_{pathway_to_save}')
        model.save(save_path)
        return True
    
    results = model.evaluate(X_test, y_test)

    # on Test
    y_predicted = model.predict(X_test)
    error_ANN = (y_test - y_predicted) * 100.0
    error_ANN_bi = None
    error_ANN_mono = None

    if check_trivial == True:
        nonzero_counts = np.sum(X_test != 0, axis=1)
        if np.any(np.where(nonzero_counts == 2)):
            error_ANN_bi = error_ANN[np.where(nonzero_counts == 2)[0],:]
            
        if np.any(np.where(nonzero_counts == 1)):
            error_ANN_mono = error_ANN[np.where(nonzero_counts == 1)[0],:]
            
        if not np.any(error_ANN_mono) and not np.any(error_ANN_bi):
            check_trivial = False

    # Plot the errors
    analyze_and_plot_ann_errors(stim_dir,side,pathway_filtered,X_test,error_ANN,error_ANN_bi,error_ANN_mono,check_trivial)

    # ====================================== Check if errors acceptable ===============================================#

    # iterate over each previously approved profile
    with open(os.path.join(stim_dir, 'target_profiles.json'), 'r') as fp:
        target_profiles = json.load(fp)
    fp.close()

    # we can discard the error sign here
    error_ANN = abs(error_ANN)

    # Check critical side-effects
    if 'SE_dict' in target_profiles:
        for key, profile in target_profiles['SE_dict'].items():
            if side == 0 and "_rh" not in key:
                continue
            elif side == 1 and "_lh" not in key:
                continue

            for pathway in profile.keys():
                if not check_pathway_errors(pathway, side, profile, error_ANN, SE_err_threshold, pathway_filtered, check_trivial, error_ANN_bi, error_ANN_mono):
                    return False

    # Merge target profiles for symptoms and threshold profiles for soft side-effects
    target_profiles_and_sse = copy.deepcopy(target_profiles.get('profile_dict', {}))
    target_profiles_and_sse.update(target_profiles.get('Soft_SE_dict', {}))

    for key, profile in target_profiles_and_sse.items():
        if side == 0 and "_rh" not in key:
            continue
        elif side == 1 and "_lh" not in key:
            continue

        for pathway in profile.keys():
            if not check_pathway_errors(pathway, side, profile, error_ANN, err_threshold, pathway_filtered, check_trivial, error_ANN_bi, error_ANN_mono):
                return False


    pathway_to_save = pathway_filtered[-1] if pathway_filtered else "default_pathway" # Save based on the final filtered pathway?
    save_path = os.path.join(stim_dir, 'NB' + SIDE_SUFFIX[side], f'ANN_approved_model_{pathway_to_save}')
    model.save(save_path)

    return pathway_filtered


if __name__ == '__main__':

    # called from MATLAB
    # sys.argv[1] - stim folder
    # sys.argv[2] - side (0 - right hemisphere)

    stim_dir = sys.argv[1]
    side = int(sys.argv[2])

    # train a separate ANN for each pathway
    from Pathways_Stats import get_simulated_pathways
    pathways, axons_in_path = get_simulated_pathways(side, stim_dir)

    for pathway in pathways:
        approx_pathways = train_test_ANN(stim_dir, side, pathway, check_trivial=True)
