#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copy of OSS-DBS LeadSettings
"""

import os
import numpy as np
#from scipy.io import loadmat
import h5py
from dataclasses import asdict

AMPLITUDE = -1.0 # in mA
HEMI_SUFFIX = ['_rh','_lh']

from ossdbs.electrodes.defaults import default_electrode_parameters

class LeadSettings:
    """Lead-DBS settings in OSS-DBS format.

    Parameters
    ----------
    mat_file_path : str
        lead-dbs settings created in ea_genvat_butenko.m

    Notes
    -----
    The lead-dbs settings are usually stored in oss-dbs_parameters.mat

    """

    def __init__(self, mat_file_path):
        
        try:
            self._file = h5py.File(str(mat_file_path), "a")
            self._h5 = True
        except ValueError:
            print(
                "\n Please, save oss-dbs_parameters using"
                "'save(oss-dbs_parameters_path, 'settings', '-v7.3')'"
            )
            raise
            # TODO  Fix non-binary .mat import

        self._settings = self._file["settings"]

    def _is_h5(self):
        return self._h5

    def _get_str(self, field_name):
        entry = ""
        if self._is_h5():
            ascii_codes = self._settings[field_name][:, 0]
            entry = "".join(np.vectorize(chr)(ascii_codes))
        else:
            entry = self._settings[field_name][0][0][0]
        if entry == "no dti":
            entry = ""
        return entry
    
    def _get_arr(self, field_name):
        if self._is_h5():
            return self._settings[field_name][:, :].T
        else:
            if len(self._settings[field_name][0][0]) == 1:
                return self._settings[field_name][0][0][0]
            else:
                return self._settings[field_name][0][0].astype(float)
            
            # Worth considering making different subclasses for different file types?
        
    def _get_num(self, field_name):
        if self._is_h5():
            return self._settings[field_name][0][0]
        else:
            return self._settings[field_name][0][0][0][0]
    
    def _set_num(self, field_name, value):
        """Private setter for numerical fields (single float/int)."""
        if self._is_h5():
            # h5py dataset values are immutable once created.
            # We must delete and recreate the dataset.
            if field_name in self._settings:
                del self._settings[field_name]
            # Create a new dataset. The shape needs to match what _get_num expects.
            self._settings.create_dataset(field_name, data=np.array([[value]]))
        else:
            # For scipy.io (non-h5), the structure is complex.
            # This is a placeholder as the non-h5 import is a TODO in the original code.
            # Assuming the original structure: settings[field][0][0][0][0] = value
            self._settings[field_name][0][0][0][0] = value

    def _set_str(self, field_name, value: str):
        """Private setter for string fields."""
        if value == "":
            # Set to "no dti" for DTI or other empty strings as per _get_str logic
            matlab_str = "no dti" if field_name == "DTI_data_name" else ""
        else:
            matlab_str = value

        if self._is_h5():
            if field_name in self._settings:
                del self._settings[field_name]

            # Convert string to a NumPy array of ASCII codes (required for h5py/Matlab)
            ascii_codes = np.array([ord(c) for c in matlab_str], dtype="uint16")
            # Must be a column vector
            data = ascii_codes.reshape(-1, 1)

            # Create the dataset as a reference to the actual string data
            dt = h5py.special_dtype(vlen=bytes)
            self._settings.create_dataset(field_name, data=data)
        else:
            # Placeholder for scipy.io
            self._settings[field_name][0][0] = matlab_str

    def _set_arr(self, field_name, value: np.ndarray):
        """Private setter for array fields."""
        value = np.asarray(value) # Ensure it is a numpy array
        if self._is_h5():
            if field_name in self._settings:
                del self._settings[field_name]
            
            # The original _get_arr expects T.T (transpose of transpose).
            # So, we save the array transposed to match the original structure.
            data = value.T
            self._settings.create_dataset(field_name, data=data)
        else:
            # Placeholder for scipy.io
            self._settings[field_name][0][0] = value

    def stretch_electrode(self, oss_electrode_name: str, hemi_idx: int):
        """Stretch electrode geometry.

        Parameters
        ----------
        oss_electrode_name: str
            electrode name in OSS-DBS format
        hemi_idx: int
            Index of brain side
        """
        stretch_list = ["tip_length", "contact_length", "contact_spacing"]
        specs_array_length = self.get_specs_array_length(oss_electrode_name)
        stretch_factor = self.get_stretch_factor(specs_array_length, hemi_idx)
        if abs(stretch_factor - 1.0) < 0.01:  # 1% tolerance
            stretch_factor = 1.0
        default_parameters = default_electrode_parameters[oss_electrode_name]
        stretched_parameters = {}
        for parameter, value in zip(
            asdict(default_parameters), asdict(default_parameters).values()
        ):
            if parameter in stretch_list:
                stretched_parameters[parameter] = value * stretch_factor
            else:
                stretched_parameters[parameter] = value
        return stretched_parameters

    def get_grid_parameters(
        self, electrode_name, hemis_idx, unit_directions, specs_array_length
    ):
        """Center lattice w.r.t. etimated stimulation volume and set resolution.

        Parameters
        ----------
        electrode_name: str
            OSS-DBS notation
        unit_directions: numpy.ndarray
            implantation trajectory
        specs_array_length: float
            length of the contact span
        hemis_idx: int
            hemisphere ID (0 - right, 1 - left)

        Returns
        -------
        grid_center: numpy.ndarray, center of the lattice model
        grid_resolution: float, resolution of the lattice model
        """
        # get grid center for lattice / voxel lattice model
        if np.any(np.isnan(self.get_stim_center()[hemis_idx, :])):
            self.grid_center = (
                self.get_imp_coord()[hemis_idx, :]
                + unit_directions[hemis_idx, :] * specs_array_length / 2
            )
        else:
            grid_center = self.get_stim_center()[hemis_idx, :]

        # set resolution
        # coarser resolution for large span electrodes
        # and large amplitudes (>5 mA or 5 V)
        phi_vector = self.get_phi_vec()[hemis_idx, :]
        if (
            electrode_name == "BostonScientificVercise"
            or electrode_name == "BostonScientificVerciseCustom"
            or np.max(np.abs(phi_vector[~np.isnan(phi_vector)])) > 5.0
        ):
            grid_resolution = 0.4
        else:
            grid_resolution = 0.33

        return grid_center, grid_resolution

    # TODO check if calculation matches Lead-DBS (Konstantin)
    def get_specs_array_length(self, oss_elec_name):
        """TODO description."""
        elec_params = default_electrode_parameters[oss_elec_name]
        return elec_params.get_distance_l1_l4()

    def get_tip_position(self, oss_elec_name: str, hemi_idx: int):
        """Get tip, implantation trajectory from head
        (Implantation_coordinate) and tail (Second_coordinate),
        and length of the contact span.
    
        Parameters
        ----------
        oss_elec_name: str
            electrode name in OSS-DBS format
        hemi_idx: int
            Index of brain side
    
        Returns
        -------
        numpy.ndarray, numpy.ndarray, float
        """
        elec_params = default_electrode_parameters[oss_elec_name]
        imp_coords = np.array(self.get_imp_coord())
        sec_coords = np.array(self.get_sec_coord())
        directions = sec_coords - imp_coords
        unit_directions = directions / np.linalg.norm(directions, axis=1)[:, np.newaxis]
        specs_array_length = self.get_specs_array_length(oss_elec_name)
        stretch_factor = self.get_stretch_factor(specs_array_length, hemi_idx)
    
        if abs(stretch_factor - 1.0) < 0.01:  # 1% tolerance
            stretch_factor = 1.0
        offset = elec_params.get_center_first_contact() * stretch_factor
        tip_position = imp_coords - offset * unit_directions
    
        return unit_directions, tip_position, specs_array_length

    def get_stretch_factor(self, specs_array_length, hemis_idx=0) -> float:
        """Compute stretch/squeeze factor between the first and the last
        contact. Relevant for normative space computations since the
        electrodes are non-linearly warped.

        Parameters
        ----------
        hemis_idx: int
            hemisphere ID (0 - right, 1 - left)
        specs_array_length: float
            distance between the first and the last contact

        Returns
        -------
        stretch_factor: float
            describes stretching for the electrode
        """
        C1_coords = self.get_imp_coord()[hemis_idx, :]
        C_last_coords = self.get_sec_coord()[hemis_idx, :]
        el_array_length = np.sqrt(
            (C1_coords[0] - C_last_coords[0]) ** 2
            + (C1_coords[1] - C_last_coords[1]) ** 2
            + (C1_coords[2] - C_last_coords[2]) ** 2
        )
        return el_array_length / specs_array_length

    def set_cntct_loc(self, hemis_idx, new_coords: np.ndarray):
        """
        Overwrites the data in the HDF5 dataset referenced by 'contactLocation', 
        handling necessary resizing if the new data exceeds current dimensions.
    
        Parameters:
            hemis_idx (int): The index used to retrieve the object reference.
            new_coords (np.ndarray): The new NumPy array to write to the dataset.
        """
        
        # 1. Get the HDF5 object reference
        obj_ref = self._settings["contactLocation"][hemis_idx, 0]
        
        # 2. Dereference the object reference to get the actual Dataset object.
        # Assumes self._file is the open h5py.File handle (must be opened in 'r+' or 'a' mode).
        target_dset = self._file[obj_ref]
        
        current_shape = target_dset.shape
        new_shape = new_coords.shape
        
        # Check if the new data's dimensions exceed the current dataset's dimensions
        if new_shape != current_shape:
            print(f"Warning: Data shape mismatch. Current shape: {current_shape}, New shape: {new_shape}")
    
            # Check if the dataset is even resizable (has maxshape defined)
            if target_dset.maxshape is None:
                raise RuntimeError(
                    f"Cannot resize dataset '{target_dset.name}'. "
                    f"New shape {new_shape} exceeds current shape {current_shape}, "
                    "and the dataset was created with a fixed size (maxshape is None)."
                )
    
            # Check if the new shape fits within the maxshape limits
            for i, (new_dim, max_dim) in enumerate(zip(new_shape, target_dset.maxshape)):
                # If max_dim is None, it's unlimited in that dimension, so it's fine.
                if max_dim is not None and new_dim > max_dim:
                     raise ValueError(
                        f"New dimension {i} ({new_dim}) exceeds the maxshape limit of {max_dim} "
                        f"for dataset '{target_dset.name}'."
                     )
            
            # 3. Explicitly resize the dataset
            try:
                target_dset.resize(new_shape)
                print(f"Dataset successfully resized from {current_shape} to {new_shape}.")
            except Exception as e:
                # Catch potential HDF5 errors during resize (e.g., if new shape violates maxshape rules)
                raise RuntimeError(f"Failed to resize dataset '{target_dset.name}' to {new_shape}: {e}")
    
        # 4. Overwrite the entire dataset's contents with the new coordinates.
        # At this point, the shapes match or the dataset was successfully resized.
        target_dset[...] = new_coords
        
        # Optional: Confirmation message
        print(f"Successfully wrote new contactCoords to referenced dataset: {target_dset.name}")


    def add_pat_fold(self, value: str):
        """Set patient folder."""
        self._set_str("Patient_folder", value)

    def add_est_in_temp(self, value: bool):
        """Set Estimate_In_Template parameter."""
        self._set_num("Estimate_In_Template", float(value))

    def add_mri_name(self, value: str):
        """Set path to MRI file."""
        self._set_str("MRI_data_name", value)

    def add_dti_name(self, value: str):
        """Set path to DTI file."""
        self._set_str("DTI_data_name", value)

    def add_gm_idx(self, value: int):
        """Set gray matter index in MRI."""
        self._set_num("GM_index", float(value))

    def add_wm_idx(self, value: int):
        """Set white matter index in MRI."""
        self._set_num("WM_index", float(value))

    def add_csf_idx(self, value: int):
        """Set CSF index in MRI."""
        self._set_num("CSF_index", float(value))

    def add_def_mat(self, value: str):
        """Set default material."""
        self._set_str("default_material", value)

    def add_elec_type(self, value: str):
        """Set electrode type."""
        self._set_str("Electrode_type", value)

    def add_y_mark_nat(self, value: np.ndarray):
        """Set native y-marker. Must be a 2D numpy array (2 hemispheres x 3 coordinates)."""
        self._set_arr("yMarkerNative", value)

    def add_y_mark_mni(self, value: np.ndarray):
        """Set MNI y-marker. Must be a 2D numpy array (2 hemispheres x 3 coordinates)."""
        self._set_arr("yMarkerMNI", value)

    def add_head_nat(self, value: np.ndarray):
        """Set native head. Must be a 2D numpy array (2 hemispheres x 3 coordinates)."""
        self._set_arr("headNative", value)

    def add_head_mni(self, value: np.ndarray):
        """Set MNI head. Must be a 2D numpy array (2 hemispheres x 3 coordinates)."""
        self._set_arr("headMNI", value)

    def add_imp_coord(self, value: np.ndarray):
        """Set implantation coordinate. Must be a 2D numpy array (2 hemispheres x 3 coordinates)."""
        self._set_arr("Implantation_coordinate", value)

    def add_sec_coord(self, value: np.ndarray):
        """Set second coordinate. Must be a 2D numpy array (2 hemispheres x 3 coordinates)."""
        self._set_arr("Second_coordinate", value)

    def add_stim_center(self, value: np.ndarray):
        """Set stimulation center. Must be a 2D numpy array (2 hemispheres x 3 coordinates)."""
        self._set_arr("stim_center", value)

    def add_stim_set_mode(self, value: bool):
        """Set stimulation mode (stimSetMode)."""
        self._set_num("stimSetMode", float(value))

    def add_cur_ctrl(self, value: np.ndarray):
        """Set current-controlled flag. Must be a 1D numpy array (2 elements)."""
        # The get method returns .T[0] meaning the stored data is 1xN.
        # It's safest to set the data as a 1xN array (or a list of lists/1-element-arrays).
        self._set_arr("current_control", value.reshape(-1, 1))

    def add_remove_electrode(self, value: bool):
        """Set remove electrode (collapse VTA?)."""
        self._set_num("removeElectrode", float(value))

    def add_neuron_model(self, value: str):
        """Set neuron model."""
        self._set_str("neuronModel", value)

    def add_signal_type(self, value: str):
        """Set signal type."""
        self._set_str("signalType", value)

    def add_pulse_width(self, value: np.ndarray):
        """Set pulse width. Must be a 1D numpy array (2 elements)."""
        self._set_arr("pulseWidth", value.reshape(-1, 1))

    def add_biphasic(self, value: bool):
        """Set biphasic pulse flag."""
        self._set_num("biphasic", float(value))

    def add_adaptive_ref(self, value: bool):
        """Set adaptive mesh refinement flag."""
        self._set_num("AdaptiveRef", float(value))

    def add_encapsulation_type(self, value: str):
        """Set encapsulation type."""
        self._set_str("encapsulationType", value)

    def add_act_thresh_vta(self, value: np.ndarray):
        """Set activation threshold for VTA. Must be a 1D numpy array (2 elements)."""
        self._set_arr("Activation_threshold_VTA", value.reshape(-1, 1))

    def add_phi_vec(self, value: np.ndarray):
        """Set Phi_vector. Must be a 2D numpy array (2 hemispheres x N contacts)."""
        self._set_arr("Phi_vector", value)

    def add_case_grnd(self, value: np.ndarray):
        """Set case grounding. Must be a 1D numpy array (2 elements)."""
        self._set_arr("Case_grounding", value.reshape(-1, 1))

    def add_calc_axon_act(self, value: bool):
        """Set calculate axon activation flag."""
        self._set_num("calcAxonActivation", float(value))

    def add_out_of_core(self, value: bool):
        """Set out of core flag."""
        self._set_num("outOfCore", float(value))

    def add_connectome(self, value: str):
        """Set connectome name."""
        self._set_str("connectome", value)

    def add_axon_len(self, value: np.ndarray):
        """Set axon length. Must be a 1D numpy array."""
        self._set_arr("axonLength", value)

    def add_fib_diam(self, value: np.ndarray):
        """Set fibre diameter. Must be a 1D numpy array."""
        self._set_arr("fiberDiameter", value)

    def add_conectome_path(self, value: str):
        """Set connectome path."""
        self._set_str("connectomePath", value)

    def add_pathway_params_file(self, value: str):
        """Set pathway parameters file name."""
        self._set_str("pathwayParameterFile", value)

    def add_inter_mode(self, value: bool):
        """Set interactive mode flag."""
        self._set_num("interactiveMode", float(value))

    def get_num_elecs(self):
        """Number of electrodes."""
        return self.NUM_ELECS

    def get_pat_fold(self):
        """Patient folder."""
        return self._get_str("Patient_folder")

    def get_est_in_temp(self):
        """TODO description."""
        return bool(self._get_num("Estimate_In_Template"))

    def get_mri_name(self):
        """Path to MRI file."""
        mri_file = self._get_str("MRI_data_name")
        return os.path.abspath(mri_file)

    def get_dti_name(self):
        """Path to DTI file."""
        return self._get_str("DTI_data_name")

    def get_gm_idx(self):
        """Gray matter index in MRI."""
        return int(self._get_num("GM_index"))

    def get_wm_idx(self):
        """White matter index in MRI."""
        return int(self._get_num("WM_index"))

    def get_csf_idx(self):
        """CSF index in MRI."""
        return int(self._get_num("CSF_index"))

    def get_def_mat(self):
        """Default material."""
        return self._get_str("default_material")

    def get_elec_type(self):
        """Electrode type."""
        return self._get_str("Electrode_type")
    
    # Used to re-compute rot_z
    def get_y_mark_nat(self):
        """Native y-marker."""
        return self._get_arr("yMarkerNative")

    def get_y_mark_mni(self):
        """MNI y-marker."""
        return self._get_arr("yMarkerMNI")

    def get_head_nat(self):
        """Native head."""
        return self._get_arr("headNative")

    def get_head_mni(self):
        """MNI head."""
        return self._get_arr("headMNI")

    def get_imp_coord(self):
        """Implantation coordinate."""
        return self._get_arr("Implantation_coordinate")

    def get_sec_coord(self):
        """Second coordinate."""
        return self._get_arr("Second_coordinate")

    def get_stim_center(self):
        """Stimulation center."""
        return self._get_arr("stim_center")

    def get_stim_set_mode(self):
        """Stimulation mode."""
        return self._get_num("stimSetMode")

    def get_cur_ctrl(self):
        """Current-controlled."""
        return self._get_arr("current_control").T[0]

    def get_neuron_model(self):
        """Neuron model."""
        return self._get_str("neuronModel")

    def get_signal_type(self):
        """Signal type."""
        return self._get_str("signalType")

    def get_pulse_width(self):
        """Pulse width."""
        return self._get_arr("pulseWidth")

    def check_biphasic(self):
        """Biphasic pulse."""
        return self._get_num("biphasic")

    def do_adaptive_ref(self):
        """Adaptive mesh refinement."""
        return self._get_num("AdaptiveRef")

    def get_cntct_loc(self, hemis_idx):
        """Contact location."""
        contactCoords = np.asarray(
            self._file[self._settings["contactLocation"][hemis_idx, 0]][:, :]
        )
        return contactCoords

    def get_encapsulation_type(self):
        """Encapsulation type."""
        return self._get_str("encapsulationType")

    def get_act_thresh_vta(self):
        """Activation threshold for VTA."""
        return self._get_arr("Activation_threshold_VTA")

    def get_phi_vec(self):
        """TODO description."""
        return self._get_arr("Phi_vector")

    def get_case_grnd(self):
        """Case grounding."""
        return self._get_arr("Case_grounding")

    def get_calc_axon_act(self):
        """Calculate axon activation."""
        return self._get_num("calcAxonActivation")

    def get_out_of_core(self):
        """Check if intermediate solution is unloaded."""
        return self._get_num("outOfCore")

    def get_connectome(self):
        """TODO description."""
        return self._get_str("connectome")

    def get_axon_len(self):
        """Axon length."""
        return self._get_arr("axonLength")

    def get_fib_diam(self):
        """Fibre diameter."""
        return self._get_arr("fiberDiameter")

    def get_conectome_path(self):
        """Connectome path."""
        return self._get_str("connectomePath")

    def get_pathway_params_file(self):
        """File with the pathway parameters file."""
        return self._get_str("pathwayParameterFile")

    def get_rot_z(self, index_side: int):
        """TODO description.

        Parameters
        ----------
        index_side: int
            Side of brain

        Notes
        -----
         Always recalculated from the other settings
        *IMPORTANT*: it is actually not native but scrf!
        """
        if self.get_est_in_temp():
            head_MNI = self.get_head_mni()[index_side, :]
            y = self.get_y_mark_mni()[index_side, :] - head_MNI
        else:
            head_nat = self.get_head_nat()[index_side, :]
            y = self.get_y_mark_nat()[index_side, :] - head_nat
        y_postop = y / np.linalg.norm(y)
        phi = np.arctan2(-y_postop[0], y_postop[1])
        return phi * 180.0 / np.pi