'''
    By K. Butenko
    This script trains and tests an ANN model to approximate pathway activation for a given electrode position
'''

import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np
import os
import sys
import json
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
N_EPOCHS = 500
MIN_AXON_NUMBER = 10
MIN_ACTIV_THRESHOLD = 5.0   # at least one train and one test protocol should have percent activation above this threshold
ZERO_PROTOCOLS_PERC = 10.0  # To soft-enforce zero activation for zero current

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

    def augment_data(self, X_train: np.ndarray, y_train: np.ndarray, extra_training_protocols: range=range(0), n_injections: int=5) -> Tuple[np.ndarray, np.ndarray]:
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
        
        # 2. Inject monopolar review (5x repetition)
        # Find monopolar protocols (assumed to be the last 56)
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

class ANNModel:
    """Handles the creation, compilation, and training of the ANN."""

    def __init__(self, input_shape: int, output_shape: int, total_axons: int):
        self.model: Optional[Sequential] = None
        self.input_shape = input_shape
        self.output_shape = output_shape
        self.total_axons = total_axons

    def create_and_compile(self):
        """Defines and compiles the Keras model."""
        model = Sequential(name="PathwayANN")
        model.add(Dense(128, input_shape=(self.input_shape,), activation='linear'))
        # Using LeakyReLU with alpha -1.25 as in original code
        model.add(Dense(1024, activation=tf.keras.layers.LeakyReLU(alpha=-1.25)))
        
        # The original code uses total axon count here, which is unusual for a typical ANN, 
        # but kept for fidelity.
        model.add(Dense(self.total_axons, activation='sigmoid'))
        
        # Final output layer
        model.add(Dense(self.output_shape, activation='sigmoid'))

        adam = optimizers.Adamax(learning_rate=LEARN_RATE)
        model.compile(optimizer=adam, loss='mean_squared_error', metrics=['mse'])
        self.model = model

    def train(self, X_train: np.ndarray, y_train: np.ndarray):
        """Fits the model to the training data."""
        if self.model is None:
            raise ValueError("Model has not been created and compiled.")
        
        print(f"Starting training for {N_EPOCHS} epochs...")
        self.model.fit(X_train, y_train, epochs=N_EPOCHS, verbose=1)

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Makes predictions using the trained model."""
        if self.model is None:
            raise ValueError("Model has not been trained.")
        return self.model.predict(X)

    def save_model(self, save_path: str):
        """Saves the trained model."""
        if self.model is None:
            raise ValueError("Model has not been trained.")
        self.model.save(save_path)
        print(f"Model saved to: {save_path}")


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

        # 1. Load, filter, and segment data
        X_train, X_test, y_train_filtered, y_test_filtered, pathway_filtered, axons_in_path = self.data_processor.load_data(self.pathway)

        if not pathway_filtered:
            print(f"Low activation levels for {self.pathway}. Skipping ANN training.")
            return None

        # 2. Augment training data
        X_train_augmented, y_train_augmented = self.data_processor.augment_data(X_train, y_train_filtered)
        
        # Ensure only one pathway is being trained (as per original code's check)
        if y_train_augmented.shape[1] > 1:
            print("Multiple pathways detected after filtering. Refactor logic expects single pathway output.")
            return None
        
        # Determine total axon count (assuming it's the sum of all simulated pathways)
        total_axons_count = np.sum(axons_in_path)

        print(f"Training on pathway: {self.pathway}")
        print(f"Number of training samples (original): {X_train.shape[0]}")
        print(f"Number of training samples (augmented): {X_train_augmented.shape[0]}")
        print(f"Number of testing samples: {X_test.shape[0]}")
        print(f"Filtered pathways: {pathway_filtered}")

        # 3. Train ANN
        self.model = ANNModel(
            X_train_augmented.shape[1], y_train_augmented.shape[1], total_axons_count
        )
        self.model.create_and_compile()
        self.model.train(X_train_augmented, y_train_augmented)

        # 4. Test and Analyze
        if self.data_processor.no_test:
            # Save model if no test data exists
            pathway_to_save = pathway_filtered[0]
            save_path = os.path.join(self.stim_dir, 'NB' + HEMI_SUFFIX[self.hemi_idx], f'ANN_approved_model_{pathway_to_save}')
            self.model.save_model(save_path)
            return pathway_filtered
        
        # Testing phase
        y_test = y_train_filtered if self.data_processor.no_test else y_test_filtered
        y_predicted = self.model.predict(X_test)
        
        # Evaluate model (optional, for metrics printing)
        # self.model.model.evaluate(X_test, y_test)
        
        error_ANN, error_ANN_bi, error_ANN_mono = self.reporter.calculate_errors(
            X_test, y_test, y_predicted, self.check_trivial
        )

        self.reporter.analyze_and_plot_ann_errors(
            pathway_filtered, X_test, error_ANN, error_ANN_bi, error_ANN_mono, self.check_trivial
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
        try:
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
                print(f"Successfully approved and saved model for pathway: {pathway}")
            else:
                print(f"Failed or skipped pathway: {pathway}")

        except NotImplementedError:
            print("VAT recruitment is not yet supported. Exiting.")
            break
        except Exception as e:
            print(f"An unexpected error occurred during processing of pathway {pathway}: {e}")
            # Optionally, continue to the next pathway or re-raise
            # raise # Uncomment to stop on first unhandled error
            continue