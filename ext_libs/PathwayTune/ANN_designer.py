'''
    By K. Butenko
    This script trains (over different configs), validates and tests an ANN model to approximate pathway activation for a given electrode position
'''

import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np
import os
import sys
import json
import csv
from typing import Tuple, List, Optional, Dict

# Set environment variable to avoid CUDA issues if not needed
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, BatchNormalization, Lambda
from tensorflow.keras import optimizers

from Pathways_Stats import get_simulated_pathways

# --- Global Constants and Utilities ---

HEMI_SUFFIX = ['_rh','_lh']
hemi_idx_LABEL = {0: '_right', 1: '_left'}

# ANN parameters
LEARN_RATE = 0.0025
N_EPOCHS = 1000
MIN_AXON_NUMBER = 10
MIN_ACTIV_THRESHOLD = 5.0   # at least one train and one test protocol should have percent activation above this threshold
ZERO_PROTOCOLS_PERC = 0.0  # To soft-enforce zero activation for zero current

# --- Configuration and Data Classes ---
class DataProcessor:
    """Handles data loading, filtering, and augmentation."""

    def __init__(self, stim_dir: str, hemi_idx: int):
        self.stim_dir = stim_dir
        self.hemi_idx = hemi_idx
        
        self.res_folder = os.path.join(stim_dir, 'OSS_sim_files' + HEMI_SUFFIX[hemi_idx], 'Results')
        self.stimsets_info: Dict = self._load_stimsets_info()
        self.train_size = self.stimsets_info.get('trainSize_actual', 0)
        self.no_test = self.stimsets_info.get('testSize_actual', 0) == 0

    def _load_stimsets_info(self) -> Dict:
        """Loads StimSets_info """
        with open(os.path.join(self.stim_dir, 'master_dict.json'), 'r') as fp:
            return json.load(fp)['stim_sets' + HEMI_SUFFIX[self.hemi_idx]]
                
    @staticmethod
    def _load_activation_results(res_folder: str, pathway: str, currents_shape: Tuple[int, int]) -> np.ndarray:
        """Loads activation results for a pathway from JSON files.

        Args:
            res_folder: Path to the results subfolder.
            pathway: Name of the pathway for the ANN training.
            currents_shape: Shape of the currents array (number of protocols, number of contacts).

        Returns:
            A NumPy array containing activation results (percentage activated) for each protocol and pathway.
        """
        activation_results = np.zeros((currents_shape[0], currents_shape[1]), dtype=float) # Corrected to handle multiple pathways if passed, although the original code implies one pathway at a time
        
        # NOTE: Original code only loads for a single named pathway, but is structured for multiple.
        # Assuming the external logic ensures `pathway` is correct for the single pathway loaded.
        # The inner loop will iterate over protocols, loading the 'percent_activated' for the named pathway.
        
        # The refactored version assumes it needs to load a 1-column array (single pathway)
        activation_results = np.zeros((currents_shape[0],1), dtype=float)

        for prot_i in range(currents_shape[0]):
            pathway_file = os.path.join(res_folder, f'Pathway_status_{pathway}_{prot_i}.json')
            if os.path.isfile(pathway_file):
                with open(pathway_file, 'r') as fp:
                    pam_res_dict = json.load(fp)
                    # Convert to fraction (0.0 to 1.0)
                    activation_results[prot_i, 0] = 0.01 * pam_res_dict['percent_activated']
        return activation_results


    def load_data(self, pathway: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str], np.ndarray]:
        """Loads currents and activation results, then filters and segments data.
        
        Args:
            pathway: Name of the pathway for the ANN training.
        """
        currents_file = os.path.join(self.stim_dir, 'OSS_sim_files' + HEMI_SUFFIX[self.hemi_idx],
                                     f'Current_protocols_{self.hemi_idx}.csv')
        currents = np.genfromtxt(currents_file, delimiter=',', skip_header=True)

        # Convert currents to A
        currents_A = currents * 0.001
        X_train = currents_A[:self.train_size, :]
        X_test = currents_A[self.train_size:, :]

        # Determine pathways and load activation results
        _, axons_in_path = get_simulated_pathways(self.hemi_idx, self.stim_dir)
        # Note: The original code loads only for the *current* pathway name, but the rest of the flow
        # seems to expect activation_results to be filtered later.
        # We'll load the single pathway's results as an N x 1 array.
        activation_results = self._load_activation_results(self.res_folder, pathway, currents_A.shape)

        y_train = activation_results[:self.train_size, :]
        y_test = activation_results[self.train_size:, :]

        # Filter pathways based on activity and axon number
        y_train_filtered, y_test_filtered, pathway_filtered = self._filter_pathways(
            pathway, axons_in_path, y_train, y_test
        )

        return X_train, X_test, y_train_filtered, y_test_filtered, pathway_filtered, axons_in_path

    def _filter_pathways(self, pathway: str, axons_in_path: np.ndarray, y_train: np.ndarray, y_test: np.ndarray) -> Tuple[np.ndarray, Optional[np.ndarray], List[str]]:
        """Filters pathways based on minimum axon number and activation threshold.
    
        Args:
            pathway: Name of the pathway for the ANN training.
            axons_in_path: NumPy array containing the number of axons in each pathway.
            y_train: All training activation results.
            y_test: All testing activation results.
        Returns:
            np.ndarray, filtered Training activation results
            np.ndarray, filtered testing activation results,
            list, filtered pathway names, and optionally filtered bipolar and monopolar activation results.
        """
        
        pathway_filtered = []

        # Check only the single loaded pathway (index 0 in the N x 1 array)
        if y_train.shape[1] > 0 and y_train.shape[1] == 1:
            
            # Find the corresponding axon count index (simplification for refactoring)
            all_pathways, all_axons = get_simulated_pathways(self.hemi_idx, self.stim_dir)
            try:
                pathway_idx = all_pathways.index(pathway)
                pathway_axon_count = all_axons[pathway_idx]
            except ValueError:
                print(f"Warning: Pathway {pathway} not found in simulated pathways list.")
                pathway_axon_count = 0

            # Check criteria for the single pathway
            if (pathway_axon_count >= MIN_AXON_NUMBER and 
                np.max(y_train[:, 0]) >= MIN_ACTIV_THRESHOLD*0.01 and 
                (self.no_test or np.max(y_test[:, 0]) >= MIN_ACTIV_THRESHOLD*0.01)):
                
                pathway_filtered.append(pathway)
                y_train_filtered = y_train
                y_test_filtered = y_test if not self.no_test else None
                return y_train_filtered, y_test_filtered, pathway_filtered
        
        # If the check failed or the array shape was unexpected (e.g., N x 0), return empty
        return np.zeros((0,0)), None, []

    def augment_data(self, X_train: np.ndarray, y_train: np.ndarray, extra_training_protocols: range=range(0), n_injections: int=1) -> Tuple[np.ndarray, np.ndarray]:
        """Injects a percentage of zero protocols into the training data.
        Optionally, also duplicates specified protocol range.
        
        Args:
            X_train: Training currents.
            y_train: Training activation results.
            extra_training_protocols: range of training protocols to be injected several times
            n_injections: number of injected protocols
        
        Returns:
            A tuple containing augmented training currents and augmented training activation results.
        """
        # 1. Inject zero protocols
        n_zero = int(ZERO_PROTOCOLS_PERC * 0.01 * y_train.shape[0])
        y_train_ext = np.zeros((y_train.shape[0] + n_zero, y_train.shape[1]), dtype=float)
        X_train_ext = np.zeros((X_train.shape[0] + n_zero, X_train.shape[1]), dtype=float)
        y_train_ext[:y_train.shape[0], :] = y_train
        X_train_ext[:X_train.shape[0], :] = X_train
        
        # 2. Inject monopolar review
        if extra_training_protocols.stop != 0:
             # Monopolar protocols were mixed into the training set (in no_test case)
            mono_protocols = X_train[extra_training_protocols,:]
            mono_solutions = y_train[extra_training_protocols,:]
            
            # Concatenate 
            X_train_augmented = np.concatenate([X_train_ext] + [mono_protocols] * n_injections)
            y_train_augmented = np.concatenate([y_train_ext] + [mono_solutions] * n_injections)
            
            return X_train_augmented, y_train_augmented
            
        else:
            return X_train_ext, y_train_ext

# class ANNModel:
#     """Handles the creation, compilation, and training of the ANN."""

#     def __init__(self, input_shape: int, output_shape: int, total_axons: int):
#         self.model: Optional[Sequential] = None
#         self.input_shape = input_shape
#         self.output_shape = output_shape
#         self.total_axons = total_axons

#     def create_and_compile(self):
#         """Defines and compiles the Keras model."""
#         model = Sequential(name="PathwayANN")
#         model.add(Dense(512, input_shape=(self.input_shape,), activation='linear', use_bias=True))
#         # Using LeakyReLU with alpha -1.25 as in original code
#         model.add(Dense(1024, activation=tf.keras.layers.LeakyReLU(alpha=-1.25), use_bias=True))
        
#         # The original code uses total axon couыnt here, which is unusual for a typical ANN, 
#         # but kept for fidelity.
#         model.add(Dense(self.total_axons, activation='sigmoid', use_bias=True))
        
#         # Final output layer
#         model.add(Dense(self.output_shape, activation='sigmoid', use_bias=False))

#         adam = optimizers.Adamax(learning_rate=LEARN_RATE)
#         model.compile(optimizer=adam, loss='mean_squared_error', metrics=['mse'])
#         self.model = model

#     def train(self, X_train: np.ndarray, y_train: np.ndarray):
#         """Fits the model to the training data."""
#         if self.model is None:
#             raise ValueError("Model has not been created and compiled.")
        
#         print(f"Starting training for {N_EPOCHS} epochs...")
#         self.model.fit(X_train, y_train, epochs=N_EPOCHS, verbose=0)

#     def predict(self, X: np.ndarray) -> np.ndarray:
#         """Makes predictions using the trained model."""
#         if self.model is None:
#             raise ValueError("Model has not been trained.")
#         return self.model.predict(X)

#     def save_model(self, save_path: str):
#         """Saves the trained model."""
#         if self.model is None:
#             raise ValueError("Model has not been trained.")
#         self.model.save(save_path)
#         print(f"Model saved to: {save_path}")


class PathwayApproximator:
    """The main class orchestrating the entire ANN training and testing process."""
    
    def __init__(self, stim_dir: str, hemi_idx: int, pathway: str, check_trivial: bool, vat_recruit: bool = False):
        
        """
        Parameters:
            stim_dir: Path to the lead-dbs stimulation folder.
            hemi_idx: Hemisphere index (0 - right, 1 - left).
            pathway: Name of the pathway for training.
            check_trivial: If True, also train and test on predefined monopolar and bipolar protocols.
            vat_recruit: If True, train ANN for VATs.
        """
        
        if vat_recruit:
            print("VAT recruitment is TBD and not implemented.")
            raise NotImplementedError()
        
        self.stim_dir = stim_dir
        self.hemi_idx = hemi_idx
        self.pathway = pathway
        self.check_trivial = check_trivial
        
        self.data_processor = DataProcessor(stim_dir, hemi_idx)
        from ANN_report import AnalysisReporter
        self.reporter = AnalysisReporter(stim_dir, hemi_idx)
        #self.model: Optional[ANNModel] = None
        _, self.axons_in_path = get_simulated_pathways(hemi_idx, stim_dir) # Needed for error checking

    def run(self) -> Optional[List[str]]:
        """Executes the training and testing workflow."""

        # #1. Load, filter, and segment data
        X_train, X_validation_test, y_train_filtered, y_validation_test_filtered, pathway_filtered, axons_in_path = self.data_processor.load_data(self.pathway)
        
        ## store in a single file
        #pathway_activation_file = os.path.join(self.stim_dir,self.pathway + '_percent_activations_Current_protocols'+ HEMI_SUFFIX[self.hemi_idx] + '_' + str(X_train.shape[0]) + '_' + str(X_test.shape[0]))
        #np.savez(pathway_activation_file, X_train=X_train, X_test=X_test, y_train_filtered=y_train_filtered, y_test_filtered=y_test_filtered, pathway_filtered=pathway_filtered, axons_in_path=axons_in_path)
        
        # # # alternatively, load from a single file
        # pathway_activation_file = os.path.join(self.stim_dir,self.pathway + '_percent_activations_Current_protocols'+ HEMI_SUFFIX[self.hemi_idx] + '_5000_32.npz')        
        # if os.path.isfile(pathway_activation_file):
        #     pathway_act_dict = np.load(pathway_activation_file)
        #     X_train, X_validation_test, y_train_filtered, y_validation_test_filtered, pathway_filtered, axons_in_path = pathway_act_dict['X_train'], pathway_act_dict['X_test'], pathway_act_dict['y_train_filtered'], pathway_act_dict['y_test_filtered'], pathway_act_dict['pathway_filtered'], pathway_act_dict['axons_in_path']
        # else:
        #     return None
        
        # split validation and test
        X_validation, X_test = X_validation_test[:int(X_validation_test.shape[0]/2),:], X_validation_test[int(X_validation_test.shape[0]/2):,:]
        y_validation, y_test_filtered = y_validation_test_filtered[:int(X_validation_test.shape[0]/2),:], y_validation_test_filtered[int(X_validation_test.shape[0]/2):,:]

        if not pathway_filtered:
            print(f"Low activation levels for {self.pathway}. Skipping ANN training.")
            return None

        # 2. Augment training data
        X_train_augmented, y_train_augmented = self.data_processor.augment_data(X_train, y_train_filtered)

        # Ensure only one pathway is being trained
        if y_train_augmented.shape[1] > 1:
            print("Multiple pathways detected after filtering. Refactored logic expects single pathway output.")
            return None

        print(f"Training on pathway: {self.pathway}")
        print(f"Number of training samples (original): {X_train.shape[0]}")
        print(f"Number of training samples (augmented): {X_train_augmented.shape[0]}")
        print(f"Number of validation samples: {X_validation.shape[0]}")
        print(f"Number of test samples: {X_test.shape[0]}")
        print(f"Filtered pathways: {pathway_filtered}")

        # Define a function to create and train the model
        def create_model(config):
            model = Sequential(name=config['name'])
            model.add(Dense(config['input_neurons'], input_shape=(X_train_augmented.shape[1],), activation='linear'))
        
            for neurons in config.get('hidden_layers', []):
                model.add(Dense(neurons, activation=tf.keras.layers.LeakyReLU(alpha=-1.25), use_bias=True))
                if config.get('dropout', True):
                    model.add(Dropout(config['dropout']))
        
            # "layer of axons"
            model.add(Dense(np.sum(axons_in_path), activation='sigmoid', use_bias=True))
            
            # "pathway activation layer"
            model.add(Dense(1, activation='sigmoid', use_bias=False))
    
            adam = optimizers.Adamax(learning_rate=config['learning_rate'])
            model.compile(optimizer=adam, loss='mean_squared_error', metrics=['mse'])
            
            return model
        
        # Define a function to evaluate the model based on the absolute error threshold
        def evaluate_performance(model, X_test, y_test, error_threshold=0.05):
            predictions = model.predict(X_test)
    
            absolute_errors = np.abs(predictions - y_test)
            relative_errors = np.abs((predictions - y_test) / y_test) * 100
            above_threshold_count = np.sum(absolute_errors > error_threshold)
            return above_threshold_count
    
        # Generate a list of epochs where the median is 250
        def generate_epochs_list(num_values, median_value):
            """Generates a list of integers with a specified median."""
            if num_values <= 0:
                return []
        
            epochs_list = []
            for i in range(num_values):
                epochs_list.append(median_value + (i - num_values // 2) * 10)  # Linear distribution
            return epochs_list
    
        # Define a set of parameter configurations to test
        learning_rates = [0.005,0.0025,0.001]
        num_layers = 1
        neurons_options = [64, 128, 256, 512, 1024]
        input_neurons_options = [128,256,512]
        epochs_options = [250, 500, 1000]
        dropout_options = [None,0.2]
        
        # Calculate the total number of configurations
        num_configs = len(input_neurons_options) * len(learning_rates) * len(epochs_options) * len(dropout_options) * len(neurons_options)
        
        # Create a grid of hidden layer configurations and epochs
        parameter_configs = []
        config_index = 0
        
        for dropout in dropout_options:
            for learning_rate in learning_rates:
                for input_neurons in input_neurons_options:
                    for neurons_per_layer in neurons_options:
                        for epochs in epochs_options:
                            layers = [neurons_per_layer] * num_layers
                            config = {
                                'name': f'config_{config_index + 1}',
                                'input_neurons': input_neurons,
                                'hidden_layers': layers,
                                'learning_rate': learning_rate,
                                'epochs': epochs,
                                'dropout': dropout
                            }
                            parameter_configs.append(config)
                            config_index += 1
        
        # Define the filename for the CSV file
        csv_filename = os.path.join(self.stim_dir,'NB' + HEMI_SUFFIX[self.hemi_idx], "ann_performance_summary_" + self.pathway + ".csv")
        
        # Prepare the header row for the CSV file
        header = ['model_name', 'input_neurons', 'hidden_layers', 'learning_rate', 'epochs', 'dropout', 'test_samples_above_5_percent_error']
        
        # just initialization for now
        performance_best = X_validation.shape[0] 
        
        # Train and evaluate each model configuration and write results to CSV
        for config in parameter_configs:
        
            # Open the CSV file in write mode        
            with open(csv_filename, 'a') as csvfile:
                writer = csv.writer(csvfile)
        
                if config['name'] == 'config_1':
                    # Write the header row
                    writer.writerow(header)
                    best_config = config
            
                print(f"Training and evaluating model: {config['name']}")
                model = create_model(config)
                model.fit(X_train, y_train_filtered, epochs=config['epochs'], verbose=0)
                performance = evaluate_performance(model, X_validation, y_validation, error_threshold=0.05)
                row = [
                    config['name'],
                    config['input_neurons'],
                    config.get('hidden_layers', []),
                    config['learning_rate'],
                    config['epochs'],
                    config.get('dropout', False),
                    performance
                ]
                writer.writerow(row)
                print(row)
                print(f"Model '{config['name']}' - Number of test samples with > 5% absolute error: {performance}\n")
    
                if performance < performance_best:
                    # the less the value, the better
                    best_config = config
    
        print("Best Configuration:")
        print(best_config)
        print(f"\n--- Summary of Results stored in '{csv_filename}' --- \n")
        print("--- Testing the Best Configuration ---")

        # validate the best configuration
        model = create_model(best_config)
        model.fit(X_train, y_train_filtered, epochs=best_config['epochs'], verbose=0)

        y_predicted = model.predict(X_test)
        
        error_ANN, error_ANN_bi, error_ANN_mono = self.reporter.calculate_errors(
            X_test, y_test_filtered, y_predicted, self.check_trivial
        )

        self.reporter.analyze_and_plot_ann_errors(
            pathway_filtered, X_test, y_test_filtered, error_ANN, error_ANN_bi, error_ANN_mono, self.check_trivial
        )

        # 5. Check Error Thresholds
        errors_acceptable = self.reporter.check_error_thresholds(
            pathway_filtered, error_ANN, error_ANN_bi, error_ANN_mono, self.check_trivial
        )

        if errors_acceptable:
            pathway_to_save = pathway_filtered[0]
            save_path = os.path.join(self.stim_dir, 'NB' + HEMI_SUFFIX[self.hemi_idx], f'ANN_approved_model_{pathway_to_save}')
            self.model.save_model(save_path)
            return pathway_filtered
        else:
            print(f"ANN model for {self.pathway} failed error threshold check.")
            return None

if __name__ == '__main__':
    # Called from MATLAB
    # sys.argv[1] - stim folder
    # sys.argv[2] - hemi_idx (0 - right hemisphere, 1 - left hemisphere)

    if len(sys.argv) < 3:
        print("Usage: python ANN_module.py <stim_folder> <hemi_idx>")
        sys.exit(1)
        
    stim_dir = sys.argv[1]
    try:
        hemi_idx = int(sys.argv[2])
    except ValueError:
        print("hemi_idx must be an integer (0 or 1).")
        sys.exit(1)

    try:
        all_pathways, _ = get_simulated_pathways(hemi_idx, stim_dir)
    except Exception as e:
        print(f"Could not load simulated pathways: {e}")
        sys.exit(1)
        
    print(f"Starting ANN approximation for hemi_idx {hemi_idx} in directory {stim_dir}.")
    
    for pathway in all_pathways:
        # it is more efficient to train separate ANNs for each pathway
        
        if pathway != 'tracts_capsule_lower_rh':
            continue

        # Create and run the approximator for each pathway
        approximator = PathwayApproximator(
            stim_dir=stim_dir, 
            hemi_idx=hemi_idx, 
            pathway=pathway, 
            check_trivial=True,
            vat_recruit=False
        )
        approved_pathways = approximator.run()

        if approved_pathways:
            print(f"Successfully approved and saved model for pathway: {pathway} \n")
        else:
            print(f"Failed or skipped pathway: {pathway} \n")

        # except NotImplementedError:
        #     print("VAT recruitment is not yet supported. Exiting.")
        #     break
        # except Exception as e:
        #     print(f"An unexpected error occurred during processing of pathway {pathway}: {e}")
        #     # Optionally, continue to the next pathway or re-raise
        #     # raise # Uncomment to stop on first unhandled error
        #     continue
