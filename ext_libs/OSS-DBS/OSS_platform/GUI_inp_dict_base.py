d = {
    # Simulation statuses
    'Segm_MRI_processed': 0,
    'DTI_processed': 0,
    'Init_neuron_model_ready': 0,
    'Init_mesh_ready': 0,
    'Adjusted_neuron_model_ready': 0,
    'Signal_generated': 0,
    'CSF_mesh_ready': 0,
    'Adapted_mesh_ready': 0,
    'Parallel_comp_ready': 0,
    'Parallel_comp_interrupted': 0,
    'IFFT_ready': 0,  # set this to one if all steps above were done and only NEURON model should be rerun

    # Details of segmented MRI and DTI
    'MRI_data_name': 'segmask.nii',
    'MRI_in_m': 0,
    'DTI_data_name': 'IITMeanTensor_NormMapping.nii.gz',
    'DTI_in_m': 0,
    'CSF_index': 3.0,  # indexing as from segmask
    'WM_index': 2.0,
    'GM_index': 1.0,
    'default_material': 3,  # Which material to use for background, gaps in segmask (3 - GM, 2 - WM, 1 - CSF)

    # geometry parameters
    'Electrode_type': 'Medtronic3389',
    # For exact names check OSS_platform/OSS-DBS_LeadDBS_integrator.py
    'Brain_shape_name': '0',  # '0' to create a default approximation, otherwise provide a file (in the stim. folder)
    'Aprox_geometry_center': [10.929, -12.117, -7.697],
    # center of brain approximation in the MRI space
    'Approximating_Dimensions': [80.0, 80.0, 80.0],
    'Implantation_coordinate_X': 10.929,  # in the MRI space (usually the middle of the lowest contact)
    'Implantation_coordinate_Y': -12.117,
    'Implantation_coordinate_Z': -7.697,
    'Second_coordinate_X': 10.929,
    'Second_coordinate_Y': -9.437,
    'Second_coordinate_Z': 3.697,
    'Rotation_Z': 0.0,

    # electrode-tissue interface
    'encap_thickness': 0.1,  # in mm
    'encap_tissue_type': 2,  # (3 - GM, 2 - WM, 1 - CSF)
    'encap_scaling_cond': 1.0,
    'encap_scaling_perm': 1.0,
    'CPE_activ': 0,
    # Constant-phase element (only for voltage-controlled, you will need to provide additional parameters in the GUI)

    # neuron model parameters
    'Axon_Model_Type': 'McIntyre2002',
    # also supports classic McNeal's model (but contact konstantin.butenko@charite.de)
    'pattern_model_name': '0',  # if you have a custom neuron model, provide it here and save in the stim folder
    'diam_fib': [5.7],  # for McIntyre's model, these diameters are discrete (2.0, 3.0, 5.7, ...)
    'n_Ranvier': [20.0],  # number of nodes of Ranvier per computational axon
    'v_init': -80.0,  # resting state potential
    'Neuron_model_array_prepared': 0,
    # if 0, use OSS-DBS GUI. Set to 1 if neuron models are supplied externally (including those allocated by OSS-DBS on supplied fiber trajectories)
    'Name_prepared_neuron_array': 'Allocated_axons.h5',  # filename of the supplied neuron models
    'Global_rot': 1,
    'x_seed': 10.929,
    'y_seed': -12.117,
    'z_seed': -7.697,
    'x_steps': 2,
    'y_steps': 0,
    'z_steps': 2,
    'x_step': 0.5,
    'y_step': 0.5,
    'z_step': 0.5,
    'alpha_array_glob': [0.0, 0.0, 0.0, 0.0, 45, 90, 135],
    'beta_array_glob': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    'gamma_array_glob': [0.0, 45.0, 90.0, 135.0, 0.0, 0.0, 0.0],

    # simulation details
    'EQS_core': 'QS',  # QS - quasitstatic, EQS - electro-quasistatic (see OSS-DBS publication)
    'Skip_mesh_refinement': 1,  # if 1, some additional parmeters are required, use OSS-DBS GUI
    'el_order': 2,  # FEM basis function order (use 3 for current-controlled)
    'number_of_processors': 0,
    # Number of physical cores used in the simulation (usually restricted by RAM). If 0 will take half of them over all sockets on Linux.
    'FEniCS_MPI': 0,  # currently disabled
    'Solver_Type': 'GMRES',  # if not sure, use 'default'
    'VTA_approx': 0,  # !!! Here it refers to |E|-field based VTA approximation !!!
    'Activation_threshold_VTA' : 0.15,  # relevant if 'VTA_approx' == 1 (either V/mm or V/mm2). For values, check Astrom et al, IEEE, 2015.
    'spectrum_trunc_method': 'Octave Band Method',
    # 'No Truncation' or 'Octave Band Method' (recommended, see Butenko et al, EMBC proceedings, 2019)
    'trunc_param': 130.0,  # octave bands will be deployed after this frequency
    'Truncate_the_obtained_full_solution': 0,  # true only if we want to use truncation of already computed full power spectrum solution (useful for benchmarks)

    # DBS signal (you can set charge-balancing in OSS_platform/Launcher_OSS_lite.py)
    'current_control': 1,  # 0 for voltage-controlled
    'Pulse_amp': [None, -0.001, None, None],  # can be negative (i.e. cathode in VC)
    # DBS pulse amplitude (V or A!) on contacts (order accor. to manufacturers) "None" - disabled, "0.0" - grounded (0 V even for current-controlled)
    'freq': 130.0,  # repetition rate of DBS signal (usually 130 Hz or 184 Hz)
    'T': 60.0,  # pulse width in uS
    't_step': 1.0,  # time step in uS
    'phi': 0.0,  # signal shift in uS
    'Signal_type': 'Rectangle',  # 'Rectangle' (classic), 'Increasing Ramp', 'Decreasing Ramp' and 'Central Triangle'
    'Ampl_scale': 1.0,
    # Scale up or down the DBS signal and rerun only NEURON, i.e. set 'IFFT_ready': 1 (works only for voltage-controlled or monopolar current-controlled)
    'external_grounding': 1,  # 1 if case grounding is used

    'patient_folder': 'Jane_Doe',  # default name, usually overwritten by Lead-DBS

    #  parameters passed EXCLUSIVELY from Lead-DBS, otherwise edit them manually
    'StimSets' : 0,        # test multiple protocols (requires Current_protocols_)
    'stretch'  : 1.0,      # in atlas space, the electrode geometry might be distorted. This parameter defines the stretch/squeeze along the electrode array
    'Stim_side': 0,        # 0 - right hemisphere, 1 - left hemisphere
}
