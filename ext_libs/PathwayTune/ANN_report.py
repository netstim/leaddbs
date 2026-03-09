#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 11:45:32 2025

@author: forel
"""

import os
import sys
import numpy as np
import json
from scipy.stats import pearsonr
import copy

import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from typing import Tuple, List, Optional, Dict, Union


HEMI_SUFFIX = ['_rh','_lh']
hemi_idx_LABEL = {0: '_right', 1: '_left'}

def find_key_recursively(data, target_key):
    """
    Recursively searches for a target_key within nested dictionaries and lists.
    
    Args:
        data: The current structure (dict, list, or other value) to search.
        target_key: The string key to find.
        
    Returns:
        True if the key is found, False otherwise.
    """
    
    # 1. Check if the current data structure is a dictionary
    if isinstance(data, dict):
        
        # Check if the key is in the current dictionary
        if target_key in data:
            return True
        
        # Recursively search in all values of the current dictionary
        for value in data.values():
            if find_key_recursively(value, target_key):
                return True
                
    # 2. Check if the current data structure is a list (or tuple)
    elif isinstance(data, (list, tuple)):
        
        # Recursively search in every item of the list
        for item in data:
            if find_key_recursively(item, target_key):
                return True
                
    # If not a dict or list, the key cannot be here (base case for recursion)
    return False


def matplotlib_kdeplot(data, ax=None, bw_method=None, **kwargs):
    """
    Replicates seaborn's kdeplot functionality using matplotlib and scipy.stats.gaussian_kde.
    """
    if ax is None:
        fig, ax = plt.subplots()

    # Handle cases where data might be sparse or non-existent
    if len(data) < 2 or np.ptp(data) == 0:
        print("Warning: Data too sparse or constant for KDE plot.")
        return ax

    kde = gaussian_kde(data, bw_method=bw_method)
    x_vals = np.linspace(min(data), max(data), 1000)
    y_vals = kde(x_vals)

    ax.plot(x_vals, y_vals, **kwargs)
    return ax

class AnalysisReporter:
    """Handles testing, error calculation, plotting, and error threshold checking."""

    def __init__(self, stim_dir: str, hemi_idx: int):
        self.stim_dir = stim_dir
        self.hemi_idx = hemi_idx
        self.nb_hemi_idx_folder = os.path.join(stim_dir, 'NB' + HEMI_SUFFIX[hemi_idx])

    def _load_error_thresholds(self) -> Tuple[float, float]:
        """Loads error thresholds from netblend_dict_file.json."""
        try:
            with open(os.path.join(self.stim_dir, 'master_dict.json'), 'r') as fp:
                netblend_dict = json.load(fp)['netblendict']
            err_threshold = netblend_dict['Err_threshold']
            SE_err_threshold = netblend_dict['SE_err_threshold']
            return err_threshold, SE_err_threshold
        except FileNotFoundError:
            print("Error: master_dict.json not found.")
            sys.exit(1)
        except KeyError:
            print("Error: 'Err_threshold' or 'SE_err_threshold' missing in master_dict.json.")
            sys.exit(1)

    def calculate_errors(self, X_test: np.ndarray, y_test: np.ndarray, y_predicted: np.ndarray, check_trivial: bool) -> Tuple[np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]:
        """Calculates errors and segments them for trivial protocols."""
        
        error_ANN = (y_test - y_predicted) * 100.0 
        error_ANN_bi = None
        error_ANN_mono = None

        if check_trivial:
            # Count non-zero elements in each protocol
            nonzero_counts = np.sum(X_test != 0, axis=1)
            
            # Bipolar (2 non-zero elements)
            bi_indices = np.where(nonzero_counts == 2)[0]
            if len(bi_indices) > 0:
                error_ANN_bi = error_ANN[bi_indices, :]

            # Monopolar (1 non-zero element)
            mono_indices = np.where(nonzero_counts == 1)[0]
            if len(mono_indices) > 0:
                error_ANN_mono = error_ANN[mono_indices, :]
        
        return error_ANN, error_ANN_bi, error_ANN_mono

    def analyze_and_plot_ann_errors(self, pathway_filtered: List[str], X_test: np.ndarray, y_test: np.ndarray,
                                    error_ann: np.ndarray, error_ann_bi: Optional[np.ndarray], 
                                    error_ann_mono: Optional[np.ndarray], check_trivial: bool):
        """Generates various plots for error analysis."""
        
        if not pathway_filtered:
            print("No pathways filtered, skipping error analysis plots.")
            return

        matplotlib.rcParams['figure.dpi'] = 200
        pathway_name = pathway_filtered[0]
        hemi_idx_label = hemi_idx_LABEL[self.hemi_idx]
        
        # --- Plotting errors vs. total current ---
        plt.figure()
        total_current = np.sum(X_test, axis=1)
        # Using only the first pathway's errors for correlation/MAE, as per original code's index [:, 0]
        abs_error_p0 = np.abs(error_ann[:, 0]) 
        
        if len(total_current) > 1 and np.std(total_current) > 1e-9:
            pearson_r, pearson_p = pearsonr(total_current, abs_error_p0)
        else:
            pearson_r, pearson_p = 0.0, 1.0

        mae = np.mean(abs_error_p0)
        plt.scatter(total_current * 1000.0, abs_error_p0)
        plt.text(0.05, 0.90, f"Pearson's R = {pearson_r:.2f}, p = {pearson_p:.5f}",
                  transform=plt.gca().transAxes, fontsize=10)
        plt.text(0.05, 0.80, f"MAE = {mae:.2f}",
                  transform=plt.gca().transAxes, fontsize=10)
        plt.xlabel("Total current, mA")
        plt.ylabel("Absolute Error of Percent Activation")
        filename_total_current = f"{pathway_name}_ANN_abs_errors_on_Test_versus_total_current{hemi_idx_label}.png"
        plt.savefig(os.path.join(self.nb_hemi_idx_folder, filename_total_current), format='png', dpi=500)
        plt.close() # Close figure to free memory

        # --- Plotting errors vs. absolute total current ---
        plt.figure()
        total_abs_current = np.sum(np.abs(X_test), axis=1)
        if len(total_abs_current) > 1 and np.std(total_abs_current) > 1e-9:
            pearson_r2, pearson_p2 = pearsonr(total_abs_current, abs_error_p0)
        else:
            pearson_r2, pearson_p2 = 0.0, 1.0

        plt.scatter(total_abs_current * 1000.0, abs_error_p0)
        plt.text(0.05, 0.85, f"Pearson's R = {pearson_r2:.2f}, p = {pearson_p2:.5f}",
                  transform=plt.gca().transAxes, fontsize=10)
        plt.text(0.05, 0.75, f"MAE = {mae:.2f}",
                  transform=plt.gca().transAxes, fontsize=10)
        plt.xlabel("Absolute current, mA")
        plt.ylabel("Absolute Error of Percent Activation")
        filename_abs_total_current = f"{pathway_name}_ANN_abs_errors_on_Test_versus_abs_current{hemi_idx_label}.png"
        plt.savefig(os.path.join(self.nb_hemi_idx_folder, filename_abs_total_current), format='png', dpi=500)
        plt.close()

        # --- Plotting KDE of absolute errors for each pathway ---
        plt.figure()
        pathways_max_errors: Dict[str, Union[float, List]] = {}
        error_filename_base = f"{pathway_name}_ANN_errors"

        for i, pathway in enumerate(pathway_filtered):
            error_data = error_ann[:, i]
            max_error = np.max(np.abs(error_data))
            pathways_max_errors[pathway] = max_error
            
            # Find the protocol that caused the max error for the first pathway
            if i == 0:
                inx_max_error = np.argmax(np.abs(error_data))
                pathways_max_errors['stim_protocol'] = X_test[inx_max_error, :].tolist()

            if max_error > 0.05: # Only plot if max error is substantial
                matplotlib_kdeplot(error_data, ax=plt.gca(), label=pathway)
        
        plt.text(0.05, 0.80, f"MAE = {mae:.2f}", transform=plt.gca().transAxes, fontsize=10)
        
        with open(os.path.join(self.nb_hemi_idx_folder, f"{error_filename_base}.json"), 'w') as save_as_dict:
            json.dump(pathways_max_errors, save_as_dict, indent=4)

        plt.legend()
        #plt.title('Errors for ANN on Test')
        # Setting xlim based on typical error range (0-25) in original code
        plt.xlim([-25.0, 25.0]) 
        plt.xlabel("Error of Percent Activation (Actual - Predicted)")
        filename_kde_test = f"{error_filename_base}_on_Test{hemi_idx_label}.png"
        plt.savefig(os.path.join(self.nb_hemi_idx_folder, filename_kde_test), format='png', dpi=500)
        plt.close()

        # --- Trivial protocol error plots (Bipolar and Monopolar) ---
        if check_trivial:
            self._plot_trivial_errors(pathway_filtered, error_filename_base, hemi_idx_label, error_ann_bi, 'Bipolar')
            self._plot_trivial_errors(pathway_filtered, error_filename_base, hemi_idx_label, error_ann_mono, 'Monopolar')
            
            
        # correlation of percent activation and error
        for i, pathway in enumerate(pathway_filtered):
        
            # correlation of percent activations and errors
            print(y_test.shape, error_ann[:, i].shape)
            pearson_r2, pearson_p2 = pearsonr(y_test[:,0], abs(error_ann[:, i]))
            print(pearson_r2, pearson_p2)
    
            # 4. Create the scatter plot
            plt.scatter(y_test*100, error_ann[:, i], color='blue', alpha=0.6, edgecolors='w')
            
            # Add a horizontal line at 0 to indicate zero error
            plt.axhline(0, color='red', linestyle='--', linewidth=1)
            
            # Add labels and title using LaTeX formatting
            plt.xlabel("Percent Activation")
            plt.ylabel("Error of Percent Activation (Actual - Predicted)")
            
            # Add grid for readability
            plt.grid(True, linestyle=':', alpha=0.7)
            
            plt.text(0.05, 0.85, f"Abs. Error vs PA R = {pearson_r2:.2f}, p = {pearson_p2:.5f}",
                     transform=plt.gca().transAxes, fontsize=10)
            
            # 5. Save the plot
            filename_PA_error = f"{error_filename_base}_versus_PA_on_Test{hemi_idx_label}.png"
            plt.savefig(os.path.join(self.nb_hemi_idx_folder, filename_PA_error), format='png', dpi=500)
            plt.close()

    def _plot_trivial_errors(self, pathway_filtered: List[str], error_filename_base: str, 
                             hemi_idx_label: str, error_array: Optional[np.ndarray], protocol_type: str):
        """Helper to plot KDE for bipolar or monopolar errors."""
        if error_array is None or error_array.size == 0:
            return

        plt.figure()
        mae = np.mean(np.abs(error_array[:, 0]))
        
        for i, pathway in enumerate(pathway_filtered):
            if np.max(np.abs(error_array[:, i])) > 0.01:
                matplotlib_kdeplot(error_array[:, i], ax=plt.gca(), label=pathway)

        plt.text(0.05, 0.80, f"MAE = {mae:.2f}", transform=plt.gca().transAxes, fontsize=10)
        plt.legend()
        #plt.title(f'Absolute Errors for ANN on {protocol_type}')
        # Dynamic xlim based on data range
        plt.xlim(np.min(error_array), np.max(error_array)) 
        plt.xlabel("Error of Percent Activation (Actual - Predicted)")
        filename = f"{error_filename_base}_on_{protocol_type}{hemi_idx_label}.png"
        plt.savefig(os.path.join(self.nb_hemi_idx_folder, filename), format='png', dpi=500)
        plt.close()

    def check_error_thresholds(self, pathway_filtered: List[str], 
                               error_array: np.ndarray, error_array_bi: Optional[np.ndarray], 
                               error_array_mono: Optional[np.ndarray], check_trivial: bool) -> bool:
        """Checks if the ANN errors are within acceptable limits."""

        err_threshold, SE_err_threshold = self._load_error_thresholds()
        
        # assuming a single pathway model
        pathway = pathway_filtered[0]
        
        # Load target profiles (symptoms, SSE, CSE)
        try:
            with open(os.path.join(self.stim_dir, 'master_dict.json'), 'r') as fp:
                target_profiles = json.load(fp)['target_profiles']
        except FileNotFoundError:
            print("Error: target_profiles.json not found.")
            return False

        # We use absolute error for checking
        abs_error_ANN = np.abs(error_array)

        # Helper to check a single pathway against a threshold
        def _check_single_pathway_error(pathway: str, threshold: float) -> bool:
        
            error_threshold_percent = threshold * 100.0

            if np.any(abs_error_ANN[:] > error_threshold_percent):
                print(f'Error threshold for {pathway} was exceeded (Test), model must be revised.')
                return False

            n_half_errors = (abs_error_ANN[:] > error_threshold_percent / 2.0).sum()
            if n_half_errors > 0.01 * abs_error_ANN.shape[0]:
                print(f'0.5 * error threshold for {pathway} was exceeded for >1% of tests, model must be revised.')
                return False

            if check_trivial:
                error_bi_check = error_array_bi is not None and np.any(np.abs(error_array_bi[:]) > error_threshold_percent)
                error_mono_check = error_array_mono is not None and np.any(np.abs(error_array_mono[:]) > error_threshold_percent)
                if error_bi_check or error_mono_check:
                    print(f'Error threshold for {pathway} was exceeded (Trivial check), model must be revised.')
                    return False
            
            return True
        
        if 'SE_dict' in target_profiles and find_key_recursively(target_profiles['SE_dict'], pathway):
            if not _check_single_pathway_error(pathway, SE_err_threshold):
                return False
            
        if 'Soft_SE_dict' in target_profiles and find_key_recursively(target_profiles['Soft_SE_dict'], pathway):
            if not _check_single_pathway_error(pathway, SE_err_threshold):
                return False

        if find_key_recursively(target_profiles['profile_dict'], pathway):
            return _check_single_pathway_error(pathway, err_threshold)
                                           
           

        # # 1. Check Critical hemi_idx-Effects (CSE)
        # for key, profile in target_profiles.get('SE_dict', {}).items():
        #     if (self.hemi_idx == 0 and "_rh" not in key) or (self.hemi_idx == 1 and "_lh" not in key):
        #         continue
        #     for pathway in profile.keys():
        #         if not _check_single_pathway_error(pathway, SE_err_threshold, profile):
        #             return False

        # # 2. Check Symptoms and Soft hemi_idx-Effects (SSE)
        # target_profiles_and_sse = copy.deepcopy(target_profiles.get('profile_dict', {}))
        # target_profiles_and_sse.update(target_profiles.get('Soft_SE_dict', {}))
        
        # for key, profile in target_profiles_and_sse.items():
        #     if (self.hemi_idx == 0 and "_rh" not in key) or (self.hemi_idx == 1 and "_lh" not in key):
        #         continue
        #     for pathway in profile.keys():
        #         if not _check_single_pathway_error(pathway, err_threshold, profile):
        #             return False

        return True