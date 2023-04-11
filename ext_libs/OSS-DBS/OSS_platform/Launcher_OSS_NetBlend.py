#    Copyright (C) 2019 Konstantin Butenko, konstantin.butenko@uni-rostock.de
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.


'''
#@author: K.Butenko (konstantin.butenko@uni-rostock.de), routines for Salome prepared in cooperation with P.Trieu. Script "tissue dielectrics.py" and axon files are adapted from FanPy (https://arxiv.org/abs/1705.10478)
GUI tree is developed in cooperation with S.A.Udongwo
The platform is in beta test and should be used with caution and a proper amount of healthy suspicion.
'''

'''This is the main script that governs the simulation flow
It contains functions:
run_master_study(), which is useful when the platform is run interactively to solve UQ or optimization problems.
It will evaluate which truncation method and with how many frequencies we need to use for a particular setup.
It can also estimate the difference in potential for different profile salome scripts to estimate the initial meshing quality (profile scripts should be created manually).
Check out the function before you use it. By default, it tries to match absolutely the result obtained with the full spectrum, which is normally unnecessary if you have thousands of neurons.
run_full_model(master_dict) governs the simulation flow taking master_dict as the input that will modify the input dictionary (defined in GUI_inp_dict.py)
'''
print("\nOSS-DBS by K.Butenko --- version 0.5")
print("Butenko K, Bahls C, Schroeder M, Koehling R, van Rienen U (2020) 'OSS-DBS: Open-source simulation platform for deep brain stimulation with a comprehensive automated modeling.' PLoS Comput Biol 16(7): e1008023. https://doi.org/10.1371/journal.pcbi.1008023")
print("____________________________________\n")
print("Check out the progress in complete_log_*.log in your stimulation folder")

import importlib
import json
import numpy as np
import os
import pickle
import subprocess
import time
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category = FutureWarning)
    import h5py

import logging
import shutil

# a function to copy folders
def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


def run_full_model(master_dict, netblendict):        # master_dict can be used for customization

    start_simulation_run = time.time()

    from progress.bar import Bar
    bar = Bar('Processing', max=11)
    print("\n Setting up the simulation...")

    import os
    netblend_folder = '/opt/NetBlend'
    os.chdir('/opt/OSS-DBS/OSS_platform')

    import sys
    #sys.path.insert(1, netblend_folder)
    sys.path.append(netblend_folder)
    sys.path.insert(1, os.environ['PATIENTDIR'])
    
    #  comment this out if you run via on MATLAB inside Singularity (cluster)
    os.environ['PATIENTDIR'] = '/opt/Patient'  # Use fixed mount path for docker
    os.environ['STIMDIR'] = os.environ['PATIENTDIR']  # legacy issue
    os.environ['LGFDIR'] = '/opt/Patient' # either '/Patient' (same folder) or '/temp' (mounted specifically for large files)
    # on clusters, always use TMPDIR option
                                               # IMPORTANT: this is actually a stim folder within the patient folder

    #===================================== Load and update input dictionary ===========================================#

    # load parameters from Lead-DBS
    try:
        with open(os.environ['PATIENTDIR'] + '/Lead_DBS_input.json', 'r') as fp:
            lead_dict = json.load(fp)
        fp.close()
    except:
        logging.critical('No input from Lead-DBS')
        lead_dict = {}
        lead_dict['stretch'] = 1.0    #  parameters passed EXCLUSIVELY from Lead-DBS, otherwise edit them manually
        lead_dict['StimSets'] = 0
        lead_dict['Stim_side'] = 0


    if master_dict['cluster_run'] == False: #  in this case, we favor GUI_inp_dict.py over Lead_DBS_input.json
                                            #  because users have a direct access to OSS-DBS GUI
        #  should be loaded this way for iterative studies (some simulation state variables change during a run)
        import GUI_inp_dict
        importlib.reload(GUI_inp_dict)
        from GUI_inp_dict import d
        # IMPORTANT: do not just update the dictionary, some parameters are explicitly changed upstream
        d['stretch'] = lead_dict['stretch']
        logging.critical('Electrode array stretch: {}'.format(d['stretch']))
        d['StimSets'] = lead_dict['StimSets']
        d['Stim_side'] = lead_dict['Stim_side']
    else:
        import GUI_inp_dict_base            #  here we just load default OSS-DBS parameters that are not in Lead_DBS_input.json
                                            #  if you still need to adjust them, either modify OSS_platform/GUI_inp_dict_base.py (preferably via OSS-DBS GUI)
                                            #  or insert them via a custom dictionary (see below)
        importlib.reload(GUI_inp_dict_base)
        from GUI_inp_dict_base import d
        d.update(lead_dict)

    try:  #  if a custom dictionary was provided (e.g. via run_pam_simulation.py), it will be loaded on top
        with open(os.environ['PATIENTDIR'] + '/Custom_input.json', 'r') as fp:
            custom_dict = json.load(fp)
        d.update(custom_dict)     # update from axon allocation dictionary directly for cluster version
        logging.critical('Custom dictionary was imported')
    except:
        logging.critical('No custom dictionary detected')

    from Dict_corrector import rearrange_Inp_dict
    d = rearrange_Inp_dict(d)  # misc. transformation of parameters to the platform's format
    d.update(master_dict)  # modifies the user provided input dictionary (e.g. for UQ study), check run_master_study() function . Warning: this does not work update the encap. layer properties and the solver during adaptive mesh refiment, because these data are realoaded from the original dictionary

    # for rodent electrodes, another settings for "activation volumes"
    if d['Electrode_type'] in ["AA_rodent_monopolar", "SR_rodent", "SNEX100"]:
        d['rodent_electrode'] = True
    else:
        d['rodent_electrode'] = False

    d['Current_sets_simple'] = False

    # to save diskspace, we truncate DBS afterpulse
    # '17' is empirically defined number (i.e. we need 16*T after DBS pulse)
    d['t_steps_trunc'] = int(d['phi'] / d['t_step']) + int(d['T'] / d['t_step']) * 17

    d['number_of_processors'] = 1
    d['el_order'] = 2
    #================================================== Optimizer =====================================================#

    # this is the internal loop of optimization (current protocol for a given position)
    d['Optimizer'] = 1  # Always on in network blending

    from Improvement4Protocol import create_NB_dictionaries
    profile_dict, Soft_SE_dict, SE_dict = create_NB_dictionaries(d['Stim_side'], netblendict['FF_dictionary'],
                                                                 disease='spontaneous human combustion')



    if d['Optimizer'] == 1:  # run simulated annealing (no local search) from scipy on the unit current solutions, i.e., just recompute the NEURON model
        # but beware that the mesh is not specifically refined for all protocols

        # just an example how it will look like (will be imported and transformed from Lead-DBS)
        d['num_iterations'] = 10
        #d['min_bound_per_contact'] = [-0.003, -0.003, -0.003, -0.003]  # in A! same length as the electrode
        #d['max_bound_per_contact'] = [0.003, 0.003, 0.003, 0.003]  # same length as the electrode
        #d['optimal_profile'] = [0.5, 0.4, 0.3, 0.4, 0.5]  # activation rates (1.0 = 100%) as many as .mat files you have for the connectome
        # if one value, provide it still in list

        # this snippet is only relevant if you have a profile, not a single .mat

        d.update(netblendict)  # imported


        d['similarity_metric'] = 'Canberra'  # Criterion to evaluate how similar two activation profiles are. Supported: Bray-Curtis, Canberra, Manhattan, Euclidean, Cosine

        if len(d['min_bound_per_contact']) != len(d['Pulse_amp']):
            logging.critical('Error: the length of the bound vector does not match the number of contacts')
            raise SystemExit

        if d["VTA_approx"] == 1:
            logging.critical("Optimization is yet not supported for VTA from E-field/Rattay's function")
            raise SystemExit

        d["current_control"], cc_multicontact, d["Skip_mesh_refinement"], d["external_grounding"], d["EQS_core"] = (1, 1, 1, True, "QS")
        logging.critical(
            "At the moment, optimization is limited to QS formulation of current-controlled mode without mesh refinement")
        logging.critical(
            "When running optimization, grounding is fixed to the casing (but optimization might create 'pseudo-grounding' contacts)")

        d['Pulse_amp'] = [1.0] * len(d['Pulse_amp'])  # unit vector

        # do we need an initial guess?

    #===================================== Charge-balancing of DBS signal =============================================#

    d['Charge_balancing'] = True  # only required for symmetric balancing, but can be also used for low amplitude
    if d['Charge_balancing'] == True:
        d['Balancing_type'] = 'Low_amplitude'  # Symmetric - Counter pulse is identical to the DBS pulse and follows right after
        # Low_amplitude - always a rectangular counter pulse that has a pulse width = 1/(DBS_Freq) - DBS_pulse_width, and the amplitude is adjusted accordingly to have the same charge

    #========================================== Settings for StimSets =================================================#

    # StimSets - solutions are derived by superposition of single contact - ground solutions
    if d['StimSets'] == 1 and os.path.isfile(
            os.environ['PATIENTDIR'] + '/Current_protocols_' + str(d['Stim_side']) + '.csv'):
        stim_protocols = np.genfromtxt(os.environ['PATIENTDIR'] + '/Current_protocols_' + str(d['Stim_side']) + '.csv',
                                       dtype=float, delimiter=',', names=True)
        if stim_protocols.size:
            d['Current_sets'] = True
            d["Skip_mesh_refinement"] = 1
            logging.critical(
                "When testing different current set, adaptive refinement is unavailable, make sure the mesh is prerefined")
            d["EQS_core"] = "QS"
            logging.critical("When testing different current set, only QS formulation is currently available")
            cc_multicontact = True # for cc_multicontact = True, active contacts will be also simulated as volumes
            d["spectrum_trunc_method"] = "Octave Band Method"
            logging.critical("When testing different current set, only Octave Band Method is currently available")
            d['Pulse_amp'] = [1.0] * len(d['Pulse_amp'])  # unit vector
            logging.critical("When testing different current sets, grounding is fixed to the casing (but you can imitate grounding by assigning a value of -1.0*sum(Icontacts)) to one of the contacts")
            d["external_grounding"] = True
            d["current_control"] = 1

            if d["VTA_approx"] == 1:
                logging.critical("Field superposition is yet not supported for VTA from E-field/Rattay's function")
                raise SystemExit

            import math

            Currents_to_check = []
            for i in range(stim_protocols.shape[0]):
                stim_prot = list(stim_protocols[i])
                for j in range(len(stim_prot)):
                    if math.isnan(stim_prot[j]):
                        stim_prot[j] = None
                    elif stim_prot[j] == 0.0:
                        logging.critical(
                            '0.0 always refers to grounding in OSS-DBS. Please, type "passive" or "float" for contacts that do not deliver currents.')
                        raise SystemExit
                    else:
                        stim_prot[j] = stim_prot[j] * 0.001  # Lead-DBS stores in mA
                if len(d['Pulse_amp']) != len(stim_prot):
                    logging.critical("Current protocols do not match the number of contacts on the electrode, exiting")
                    raise SystemExit
                Currents_to_check.append(stim_prot)

                if len(d['Pulse_amp']) == 1:  #  For one contact electrodes, the simplified StimSets is used
                    #logging.critical('For one contact electrodes, the simplified StimSets is used')
                    d['Current_sets'] = False
                    d['Current_sets_simple'] = True
                    d['Pulse_amp'] = Currents_to_check[0]
        else:
            d['Current_sets'] = False
    else:
        d['Current_sets'] = False

    #====================================== Setup for Astrom VTA approximation ========================================#

    if d["VTA_approx"] == 1:  # for this case, we use a neuron array that matches the VTA array in dimensions
        d['Axon_Model_Type'] = 'Reilly2016'
        d['diam_fib'] = 5.0
        d['x_steps'], d['y_steps'], d['z_steps'] = (20, 0, 20)  # we assume that Z-axis is ventra-dorsal in the MRI
        d['Global_rot'] = 1
        d['alpha_array_glob'], d['beta_array_glob'], d['gamma_array_glob'] = ([0], [0], [0])
        d["Neuron_model_array_prepared"] = 0
        if d['rodent_electrode'] == True:  # rodent VTA
            d['x_seed'], d['y_seed'], d['z_seed'] = (d['Implantation_coordinate_X'], d['Implantation_coordinate_Y'], d['Implantation_coordinate_Z'])  # rodent electrodes have only 1 or 2 contacts
            d['n_Ranvier'] = 3
            d['x_step'], d['y_step'], d['z_step'] = (0.1, 0.1, 0.1)
        else:
            d['x_seed'], d['y_seed'], d['z_seed'] = (
            d['Implantation_coordinate_X'], d['Implantation_coordinate_Y'] + 3.0, d['Implantation_coordinate_Z'] + 5.0)     # it makes sense to shift it a bit from the tip
            d['n_Ranvier'] = 22
            d['x_step'], d['y_step'], d['z_step'] = (1.0, 1.0, 1.0)

    #============================ Check simulation setup and state, load the corresponding data =======================#

    if d["Stim_side"] == 0:
        logging.critical("Processing right hemisphere\n")
    else:
        logging.critical("Processing left hemisphere\n")

    if d["current_control"] == 1 and d["CPE_activ"] == 1:
        d["CPE_activ"] = 0
        logging.critical("Disabling CPE for current-controlled simulation")

    # check how many contacts have assigned current or non-zero voltage
    Pulse_amp_active_non_zero = [x for x in d["Pulse_amp"] if (x is not None) and (x != 0.0)]

    # if cc_multicontact = True, active contacts will be also simulated as volumes
    if d["current_control"] == 1 and len(Pulse_amp_active_non_zero) > 1:  # multicontact current-controlled case
        cc_multicontact = True
    else:
        cc_multicontact = False

    from Sim_state import check_state
    check_state(d)  # switches simulation state depending on flags in the input dictionary, also manages project folders

    bar.next()
    #========================================== Formatting MRI and DTI data ===========================================#

    if d["Segm_MRI_processed"] == 0 and d["DTI_processed"] == 1:
        logging.critical("MRI data are new, the DTI data will be reprocessed")
        d["DTI_processed"] = 0

    #d["Segm_MRI_processed"] = 0
    from MRI_DTI_processing import MRI_segm_data
    MRI_seg_param = MRI_segm_data(d)  # also creates 'MRI_DTI_derived_data/Tissue_array_MRI.csv' and meta data

    bar.next()

    if d["DTI_data_name"] != 0:  # 0 if we provide DTI data
        tensor_conduct = 1  # flag to use conductivity tensors
        from MRI_DTI_processing import DTI_meta_data
        DTI_param = DTI_meta_data(d, MRI_seg_param)  # also creates 'MRI_DTI_derived_data/Tensor_array_DTI.csv' and meta data
    else:
        DTI_param = 0       # this class won't be used
        tensor_conduct = 0

    bar.next()
    #=========================================== Geometry generation ==================================================#

    # tmp fix, otherwise docker throws an MPI error when running on mac
    from Neural_array_processing import Neuron_array
    N_array = Neuron_array(d, MRI_seg_param)    # initialize a class instance that describes neuron models


    if d["Init_mesh_ready"] == 1:
        with open(os.environ['PATIENTDIR'] + '/Meshes/Mesh_ind.file', "rb") as f:
            Domains = pickle.load(f)    # instance of Mesh_ind class containing indices for mesh entities (defined in CAD_Salome.py)
        bar.next()
        bar.next()
    else:
        print("\n Preparing neuron models and FEM mesh...")
        if d["Brain_shape_name"] == 0 or d["Brain_shape_name"] == '0' or d["Brain_shape_name"] == '':  # Creates a brain approximation (ellisploid)
            from CAD_Salome import build_brain_approx
            # builds an elliptic brain approx. either using provided or MRI dimensions (if d['Approximating_Dimensions'] = 0)
            approx_dimensions = build_brain_approx(d['Approximating_Dimensions'], d['Aprox_geometry_center'], MRI_param = MRI_seg_param)  # also creates 'Brain_substitute.brep'
            d["Brain_shape_name"] = 'Brain_substitute.brep'

        #======================================= Initial neuron array generation ======================================#

        from Neural_array_processing import Neuron_array
        N_array = Neuron_array(d, MRI_seg_param)    # initialize a class instance that describes neuron models

        if d["Init_neuron_model_ready"] == 0 and d["Neuron_model_array_prepared"] == 0:
            logging.critical("----- Creating initial neuron array -----")
            N_array.build_neuron_models()  #  builds a pattern model, if not provided, then builds a neuron array and stores in 'Neuron_model_arrays/All_neuron_models.csv'
            with open(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/Neuron_array_class.file', "wb") as f:
                pickle.dump(N_array, f, pickle.HIGHEST_PROTOCOL)

        elif d["Neuron_model_array_prepared"] == 1 and d["Init_neuron_model_ready"] == 0:
            logging.critical("----- Creating initial neuron array from a provided neuron array -----")
            N_array.process_external_array()  #  adjusts the prepared neuron array to the computational domain (only for brain substitutes!), then stores the neuron array in 'Neuron_model_arrays/All_neuron_models.csv'
            with open(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/Neuron_array_class.file', "wb") as f:
                pickle.dump(N_array, f, pickle.HIGHEST_PROTOCOL)

        else:
            logging.critical("--- Initial neuron array was loaded\n")
            with open(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/Neuron_array_class.file', "rb") as f:
                N_array = pickle.load(f)

        bar.next()

        # if a brain substitute is used, it will be enlarged to encompass previously defined neuron array (if necessary)
        if d["Brain_shape_name"] == 'Brain_substitute.brep':
            from CAD_Salome import check_approx_dimensions
            check_approx_dimensions(approx_dimensions, d['Aprox_geometry_center'], N_array.ROI_radius)


        #======================================== Final geometry generation ===========================================#
        print("\n Generating geometry and mesh...")

        from CAD_Salome import build_final_geometry
        # creates and discretizes the geometry with the implanted electrode, encapsulation layer and ROI, converts to the approp. format. The files are stored in Meshes/
        # Domains is instance of Mesh_ind class containing indices for mesh entities (defined in CAD_Salome.py)
        Domains = build_final_geometry(d, MRI_seg_param, N_array.ROI_radius, cc_multicontact)

        bar.next()
    #==========================================Adjusting neuron array==================================================#

    if d["Adjusted_neuron_model_ready"] == 0:
        print("\n Adjusting neuron models")
        if d["Init_mesh_ready"] == 1:
            with open(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/Neuron_array_class.file', "rb") as f:
                N_array = pickle.load(f)

        # subtracts neurons from the previously define All_neuron_models.csv, if the are non-physical (inside encap. layer or CSF, outside of the domain, intersect with the electrode geometry.)
        N_array.adjust_neuron_models(Domains, MRI_seg_param)
        with open(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/Neuron_array_class.file', "wb") as f:
            pickle.dump(N_array, f, pickle.HIGHEST_PROTOCOL)
    else:
        logging.critical("--- Adjusted neuron array was loaded\n")
        with open(os.environ['PATIENTDIR'] + '/Neuron_model_arrays/Neuron_array_class.file', "rb") as f:
            N_array = pickle.load(f)

    bar.next()

    number_of_points = int(np.sum(N_array.pattern['num_segments'] * N_array.N_models))

    # later should be added to GUI
    if d['Axon_Model_Type'] == 'McIntyre2002_ds':
        d['Axon_Model_Type'] = 'McIntyre2002'
        d['downsampled_axon'] = True
    else:
        d['downsampled_axon'] = False

    #  In this case we actually don't need steps above
    #  If VTA is computed from |E| or |div(E)|, sample probing points equidistantly (VTA array)
    VTA_edge = 0 # just initialization
    if d["VTA_approx"] == 1:

        from VTA_from_array import create_VTA_array, get_VTA
        VTA_edge, number_of_points = create_VTA_array([d['x_seed'], d['y_seed'], d['z_seed']],
                                                                   d['rodent_electrode'], MRI_seg_param)
        # also resaves the points as /Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv
        N_array.pattern['num_segments'] = [number_of_points]

    #===============================================Signal creation====================================================#
    # in case there was an update of the signal
    Pulse_amp = [x for x in d["Pulse_amp"] if x is not None]  # now we don't need None values (floating potentials), Contacts contain only the active ones
    Domains.Amp_vector = Pulse_amp   #  WARNING: if the choice of the stim contacts was changed, re-run geometry/mesh generation!
    with open(os.environ['PATIENTDIR'] + '/Meshes/Mesh_ind.file', "wb") as f:
        pickle.dump(Domains, f, pickle.HIGHEST_PROTOCOL)

    from Signal_generator import Stim_pulse  # description of a single DBS period
    DBS_pulse = Stim_pulse(d, cc_multicontact)  # initialize a class instance that describes neuron models

    if d["Signal_generated"] == 0:
        print("\n Generating DBS pulse...")
        logging.critical("----- Generating DBS signal -----")
        from Signal_generator import Stim_pulse  # description of a single DBS period
        DBS_pulse.create_pulse(d)
        with open(os.environ['PATIENTDIR'] + '/Stim_Signal/DBS_pulse.file', "wb") as f:
            pickle.dump(DBS_pulse, f, pickle.HIGHEST_PROTOCOL)
    else:
        logging.critical("--- DBS signal is taken from the previous simulation\n")
        with open(os.environ['PATIENTDIR'] + '/Stim_Signal/DBS_pulse.file', "rb") as f:
            DBS_pulse = pickle.load(f)

    bar.next()
    # =============================CSF_refinement==================================#
    if d["Parallel_comp_ready"] == 0:
        if d["external_grounding"] == 1:
            d["external_grounding"] == True
            from Ext_ground_preref import refine_external_ground
            refine_external_ground(Domains)

    '''if we want to skip CSF and adaptive mesh refinement'''
    if d["Skip_mesh_refinement"] == 1:
        from CSF_refinement_new import Dummy_CSF
        Dummy_CSF()  # will resave 'Meshes/Mesh_unref.xml' as 'Results_adaptive/mesh_adapt.xml.gz' and the same for subdomains and boundaries
        d["CSF_mesh_ready"] = 1
        d["Adapted_mesh_ready"] = 1
        logging.critical("CSF and adaptive mesh refinement was skipped\n")
    else:
        print("\n Conducting mesh refinement")
        # from Signal_generator import pick_refinement_freqs
        if d["refinement_frequency"][0] == -1:  # frequencies basing on the power spectrum distribution
            from Signal_generator import pick_refinement_freqs
            ref_freqs = pick_refinement_freqs(DBS_pulse.FR_vector_signal, DBS_pulse.Xs_signal_norm, d["num_ref_freqs"])  # A is needed to switch back to unit pulse. We use unit pulse in adaptive to avoid confusion (math modules always scale the signal to the required value)
        else:
            ref_freqs = d["refinement_frequency"]

    if d["CSF_mesh_ready"] == 0:
        # from CSF_refinement_MPI import launch_CSF_refinement
        from CSF_refinement_new import launch_CSF_refinement
        Scaling_CSF = launch_CSF_refinement(d, MRI_seg_param, DTI_param, Domains, tensor_conduct, cc_multicontact, ref_freqs)

        from Visualization_files.paraview_find_arrayname import get_Para_Array_name
        from Parameter_insertion import insert_f_name
    elif d["Skip_mesh_refinement"] == 0:  # should be changed
        Scaling_CSF = np.genfromtxt(os.environ['PATIENTDIR'] + '/CSF_ref/Scaling_CSF.csv', delimiter=' ')

    # ====================Mesh adaptive algorithm==================================#

    if d["Adapted_mesh_ready"] == 0:
        from Mesh_adaption_hybrid import mesh_adapter
        Ampl_on_vert = mesh_adapter(MRI_seg_param, DTI_param, Scaling_CSF, Domains, d, tensor_conduct, cc_multicontact,
                                    ref_freqs)  # also saves adapted mesh
        np.savetxt(os.environ['PATIENTDIR'] + '/Results_adaptive/Ampl_on_vert.csv', Ampl_on_vert,
                   delimiter=" ")  # just to check

    bar.next()

    # ===================Truncate the frequency spectrum===========================#
    if d["spectrum_trunc_method"] != 'No Truncation':

        from Truncation_of_spectrum import Truncated_spectrum  # description of a single DBS period
        Truncated_pulse = Truncated_spectrum(d, DBS_pulse)  # initialize a class instance that describes neuron models

        if os.path.isfile(
                os.environ['PATIENTDIR'] + '/Stim_Signal/Truncated_pulse_' + str(d['trunc_param'] * 1.0) + '.file'):
            logging.critical("--- Truncated DBS spectrum is taken from the previous simulation\n")
            with open(
                    os.environ['PATIENTDIR'] + '/Stim_Signal/Truncated_pulse_' + str(d['trunc_param'] * 1.0) + '.file',
                    "rb") as f:
                Truncated_pulse = pickle.load(f)
                # check for consistency
                if d["spectrum_trunc_method"] != Truncated_pulse.trunc_method:
                    logging.critical('The truncation method was changed, calculating a new truncated spectrum')
                    Truncated_pulse.get_freqs_for_calc()
                    with open(os.environ['PATIENTDIR'] + '/Stim_Signal/Truncated_pulse_' + str(
                            d['trunc_param'] * 1.0) + '.file', "wb") as f:
                        pickle.dump(Truncated_pulse, f, pickle.HIGHEST_PROTOCOL)
            logging.critical("Number of frequencies after truncation: {}".format(
                Truncated_pulse.FR_vector_signal_new.shape[0]))
        else:
            Truncated_pulse.get_freqs_for_calc()
            with open(
                    os.environ['PATIENTDIR'] + '/Stim_Signal/Truncated_pulse_' + str(d['trunc_param'] * 1.0) + '.file',
                    "wb") as f:
                pickle.dump(Truncated_pulse, f, pickle.HIGHEST_PROTOCOL)

        d["trunc_param"] = Truncated_pulse.trunc_param  # only changed when using one_sample_10kHz

    # ==========Calculate freq in parallel and rearrange field array===============#

    if d['Current_sets'] == True or d['Optimizer'] == 1:

        if d["Parallel_comp_ready"] == 0:
            print("\n Solving FFEM for unit currents")
            if ["Parallel_comp_interrupted"] == 1:
                import os
                if not (os.path.isfile(os.environ['PATIENTDIR'] + '/Field_solutions/Phi_real_scaled_' + str(
                        d["freq"]) + 'Hz.pvd') or os.path.isfile(
                        os.environ['PATIENTDIR'] + 'Field_solutions/Phi_real_unscaled_' + str(
                                d["freq"]) + 'Hz.pvd')):  # to make sure that there were interrupted computations
                    logging.critical("There were no previous computations, 'Parallel_comp_interrupted' is put to 0")
                    ["Parallel_comp_interrupted"] == 0

            from Parallel_unit_current_calc import calculate_in_parallel
            '''calculate_in_parallel will save a sorted_solution array if IFFT is pointwise or will save the whole field for each frequency in Field_solutions_functions/ if full_IFFT is requested'''

            if d["spectrum_trunc_method"] == 'No Truncation':
                logging.critical("----- Calculating electric field in the frequency spectrum -----")
                calculate_in_parallel(d, DBS_pulse.FR_vector_signal, Domains, MRI_seg_param, DTI_param, tensor_conduct, number_of_points,
                                      cc_multicontact)

            if d["spectrum_trunc_method"] == 'High Amplitude Method' or d["spectrum_trunc_method"] == 'Cutoff Method' or \
                    d["spectrum_trunc_method"] == 'Octave Band Method':
                logging.critical("----- Calculating electric field in the truncated frequency spectrum -----")
                calculate_in_parallel(d, Truncated_pulse.FR_vector_signal_new, Domains, MRI_seg_param, DTI_param, tensor_conduct,
                                      number_of_points, cc_multicontact)
        else:
            logging.critical("--- Results of calculations in the frequency spectrum were loaded\n")
        bar.next()

        if d["spectrum_trunc_method"] == 'No Truncation':
            name_sorted_solution = os.environ['PATIENTDIR'] + '/Field_solutions/sorted_solution_per_contact.csv'
        else:
            name_sorted_solution = os.environ['PATIENTDIR'] + '/Field_solutions/sorted_solution_per_contact_' + str(
                d["spectrum_trunc_method"]) + '_' + str(d["trunc_param"]) + '.csv'

        d['number_of_processors'] = d['number_of_processors'] * 2  # use multithreading
        logging.critical("Switching to multithreading")

        if d["IFFT_ready"] == 0 and d["VTA_approx"] != 1:
            print("\n Conducting IFFT for unit current solutions")
            if os.path.isdir(os.environ['PATIENTDIR'] + '/Axons_in_time'):  # we always re-run NEURON simulation
                os.system('rm -fr ' + os.environ['PATIENTDIR'] + '/Axons_in_time')
                os.makedirs(os.environ['PATIENTDIR'] + '/Axons_in_time')

            from Current_scaler import conduct_unit_IFFT
            conduct_unit_IFFT(d, DBS_pulse, N_array, name_sorted_solution,
                              Truncated_pulse)

        bar.next()

        if d['Current_sets'] == True:
            print("\n Probing action potentials for the current protocols")
            it_num = 0
            for current_comb in Currents_to_check:
                from Current_scaler import find_activation
                Activation = find_activation(current_comb, d, MRI_seg_param, DBS_pulse.Xs_signal_norm, N_array, DBS_pulse.FR_vector_signal, DBS_pulse.t_vector,
                                             name_sorted_solution, Truncated_pulse.inx_start_octv, it_num, VTA_edge=VTA_edge)
                it_num += 1
            bar.next()
        elif d['Optimizer'] == 1:
            print("\n Solving the current optimization problem")
            logging.critical('Running optimization of current protocol')
            from scipy.optimize import dual_annealing
            from Current_scaler_NetBlend import compute_global_score

            args_all = [d, MRI_seg_param, DBS_pulse.Xs_signal_norm, N_array, DBS_pulse.FR_vector_signal, DBS_pulse.t_vector, name_sorted_solution,
                        Truncated_pulse.inx_start_octv]

            d, profile_dict, Soft_SE_dict, SE_dict, MRI_param, Xs_signal_norm, N_array, FR_vector_signal, t_vector, name_sorted_solution, inx_start_octv


            if os.path.isfile(os.environ['PATIENTDIR'] + '/Best_scaling_yet.csv'):
                initial_scaling = np.genfromtxt(os.environ['PATIENTDIR'] + '/Best_scaling_yet.csv', delimiter=' ')
                res = dual_annealing(compute_global_score,
                                     bounds=list(zip(d['min_bound_per_contact'], d['max_bound_per_contact'])),
                                     args=args_all, x0=initial_scaling, maxfun=d['num_iterations'], seed=42, visit=2.62,
                                     no_local_search=True)
            else:
                res = dual_annealing(compute_global_score,
                                     bounds=list(zip(d['min_bound_per_contact'], d['max_bound_per_contact'])),
                                     args=args_all, maxfun=d['num_iterations'], seed=42, visit=2.62,
                                     no_local_search=True)

            # if you have outer loop for electrode placement optimization, you would want to compare res.fun against the previous best solution here
            np.savetxt(os.environ['PATIENTDIR'] + 'Best_scaling_yet.csv', res.x, delimiter=' ')
            opt_res = np.append(res.x, res.fun)
            logging.critical("Current optimization results: {}".format(' '.join(map(str, list(opt_res)))))
            bar.next()
        if d["Stim_side"] == 0:
            subprocess.call(['touch', os.environ['PATIENTDIR'] + '/success_rh.txt'])
        else:
            subprocess.call(['touch', os.environ['PATIENTDIR'] + '/success_lh.txt'])
        bar.finish()
        return True

    if d["Parallel_comp_ready"] == 0:
        print("\n Solving FFEM")
        if ["Parallel_comp_interrupted"] == 1:
            if not (os.path.isfile(os.environ['PATIENTDIR'] + '/Field_solutions/Phi_real_scaled_' + str(
                    d["freq"]) + 'Hz.pvd') or os.path.isfile(
                    os.environ['PATIENTDIR'] + '/Field_solutions/Phi_real_unscaled_' + str(
                            d["freq"]) + 'Hz.pvd')):  # to make sure that there were interrupted computations
                logging.critical("There were no previous computations, 'Parallel_comp_interrupted' is put to 0")
                ["Parallel_comp_interrupted"] == 0

        from Parallel_field_calc import calculate_in_parallel
        '''calculate_in_parallel will save a sorted_solution array if IFFT is pointwise or will save the whole field for each frequency in Field_solutions_functions/ if full_IFFT is requested'''

        if d["spectrum_trunc_method"] == 'No Truncation':
            logging.critical("----- Calculating electric field in the frequency spectrum -----")
            calculate_in_parallel(d, DBS_pulse.FR_vector_signal, Domains, MRI_seg_param, DTI_param, tensor_conduct, number_of_points,
                                  cc_multicontact)

        if d["spectrum_trunc_method"] == 'High Amplitude Method' or d["spectrum_trunc_method"] == 'Cutoff Method' or d[
            "spectrum_trunc_method"] == 'Octave Band Method':
            logging.critical("----- Calculating electric field in the truncated frequency spectrum -----")
            calculate_in_parallel(d, Truncated_pulse.FR_vector_signal_new, Domains, MRI_seg_param, DTI_param, tensor_conduct, number_of_points,cc_multicontact)
    else:
        logging.critical("--- Results of calculations in the frequency spectrum were loaded\n")
    bar.next()

    if d["spectrum_trunc_method"] == 'No Truncation':
        name_sorted_solution = os.environ['PATIENTDIR'] + '/Field_solutions/sorted_solution.csv'
    else:
        name_sorted_solution = os.environ['PATIENTDIR'] + '/Field_solutions/sorted_solution_' + str(
            d["spectrum_trunc_method"]) + '_' + str(d["trunc_param"]) + '.csv'

        '''sorted_solution structure is'''
        '''''''''x y z Phi_r Phi_im Freq'''
        '''pnt1'''                '''130'''
        '''...'''                 '''260'''
        '''pnt1'''                '''1300000'''
        '''pnt2'''                '''130'''
        '''...'''
        '''Points are defined by the neuron array and might be duplicated. The order of points depends on the neural models'''

    # =============================Full Field IFFT=================================#
    if d["IFFT_ready"] == 1:
        logging.critical("--- Results of IFFT (FFEM) were loaded\n")

    # if IFFT_on_VTA_array == 1:
    if d["VTA_approx"] == 1:
        print("\n Creating VAT")
        # Max_signal_for_point=ifft_on_VTA_array(name_sorted_solution,d,FR_vector_signal,Xs_signal_norm,t_vector,d["T"],inx_start_octv,arrays_shape)

        if d["IFFT_ready"] == 0:
            from Parallel_IFFT_on_VTA_array import get_IFFT_on_VTA_array
            Max_signal_for_point = get_IFFT_on_VTA_array(DBS_pulse.FR_vector_signal, DBS_pulse.Xs_signal_norm, DBS_pulse.t_vector, d["T"], Truncated_pulse.inx_start_octv,number_of_points)
        else:
            Max_signal_for_point = np.load(os.environ['PATIENTDIR']+'/Field_solutions/Max_voltage_resp.npy')

        get_VTA(d, VTA_edge, Max_signal_for_point, np.array(MRI_seg_param.first_vox_coords))
        # home_dir=os.path.expanduser("~")
        if d["Stim_side"] == 0:
            subprocess.call(['touch', os.environ['PATIENTDIR'] + '/success_rh.txt'])
        else:
            subprocess.call(['touch', os.environ['PATIENTDIR'] + '/success_lh.txt'])
        bar.next()
        return True

        d["IFFT_ready"] = 1  # modification of dictionary

    '''Whole domain IFFT: will be enabled later'''
    # if d["IFFT_ready"] == 0 and d["Full_Field_IFFT"] == 1:
    #     from Full_IFFT_field_function import get_field_in_time

    #     VTA_size = get_field_in_time(d,FR_vector_signal,Xs_signal_norm,t_vector)        # also uses data from Field_solutions_functions/ and if spectrum truncation is applied, than data from Stim_Signal/
    #
    #     d["IFFT_ready"] = 1               #modification of dictionary

    #     if d["VTA_from_divE"]==True or d["VTA_from_E"]==True:           #else it will just jump to NEURON_direct_run.py
    #         minutes=int((time.time() - start_simulation_run)/60)
    #         secnds=int(time.time() - start_simulation_run)-minutes*60
    #         total_seconds=time.time() - start_simulation_run
    #         logging.critical("---Simulation run took ",minutes," min ",secnds," s ")

    #         return None,VTA_size

    # ==============Truncation of the already computed spectrum====================#

    if d["Truncate_the_obtained_full_solution"] == 1 and d["IFFT_ready"] == 0:
        if d["spectrum_trunc_method"] == 'High Amplitude Method' or d["spectrum_trunc_method"] == 'Cutoff Method':
            logging.critical("----- Conducting IFFT truncating already obtained full solution -----")
            from Field_IFFT_on_different_axons import convolute_signal_with_field_and_compute_ifft
            if isinstance(d["n_Ranvier"], list):  # if different populations
                last_point = 0
                hf = h5py.File(os.environ['PATIENTDIR'] + '/' + d["Name_prepared_neuron_array"], 'r')
                lst_population_names = list(hf.keys())
                hf.close()
                for i in range(len(d["n_Ranvier"])):
                    logging.critical("in {} population".format(lst_population_names[i]))
                    last_point = convolute_signal_with_field_and_compute_ifft(d, DBS_pulse, N_array.N_models[i],
                                                                              N_array.pattern['num_segments'][i],
                                                                              name_sorted_solution,
                                                                              dif_axons=True, last_point=last_point)
            else:
                convolute_signal_with_field_and_compute_ifft(d, DBS_pulse, N_array.N_models,
                                                             N_array.pattern['num_segments'], name_sorted_solution)
        else:
            logging.critical("Truncation of the obtained full solution is only for high. ampl and cutoff methods")

    # ==============================Parall. IFFT===================================#

    if d["IFFT_ready"] == 0 and d["Truncate_the_obtained_full_solution"] != 1:
        print("\n Conducting IFFT")
        logging.critical("----- Conducting signal scaling and IFFT -----")
        from Field_IFFT_on_different_axons import convolute_signal_with_field_and_compute_ifft

        if d["spectrum_trunc_method"] == 'No Truncation':
            if isinstance(d["n_Ranvier"], list):  # if different populations
                last_point = 0
                hf = h5py.File(os.environ['PATIENTDIR'] + '/' + d["Name_prepared_neuron_array"], 'r')
                lst_population_names = list(hf.keys())
                hf.close()
                for i in range(len(d["n_Ranvier"])):
                    logging.critical("in {} population".format(lst_population_names[i]))
                    last_point = convolute_signal_with_field_and_compute_ifft(d, DBS_pulse.Xs_signal_norm, N_array.N_models[i],
                                                                              N_array.pattern['num_segments'][i],
                                                                              DBS_pulse.FR_vector_signal, DBS_pulse.t_vector, DBS_pulse.A,
                                                                              name_sorted_solution, dif_axons=True,
                                                                              last_point=last_point)
            else:
                convolute_signal_with_field_and_compute_ifft(d, DBS_pulse.Xs_signal_norm, N_array.N_models,
                                                             N_array.pattern['num_segments'], DBS_pulse.FR_vector_signal,
                                                             DBS_pulse.t_vector, DBS_pulse.A, name_sorted_solution)

        if (d["spectrum_trunc_method"] == 'High Amplitude Method' or d["spectrum_trunc_method"] == 'Cutoff Method') and \
                d["Truncate_the_obtained_full_solution"] == 0:
            if isinstance(d["n_Ranvier"], list):  # if different populations
                last_point = 0
                hf = h5py.File(os.environ['PATIENTDIR'] + '/' + d["Name_prepared_neuron_array"], 'r')
                lst_population_names = list(hf.keys())
                hf.close()
                for i in range(len(d["n_Ranvier"])):
                    logging.critical("in {} population".format(lst_population_names[i]))
                    last_point = convolute_signal_with_field_and_compute_ifft(d, Xs_signal_norm_new,
                                                                              N_array.N_models[i],
                                                                              N_array.pattern['num_segments'][i],
                                                                              FR_vector_signal_new, DBS_pulse.t_vector, DBS_pulse.A,
                                                                              name_sorted_solution, dif_axons=True,
                                                                              last_point=last_point)
            else:
                convolute_signal_with_field_and_compute_ifft(d, Xs_signal_norm_new, N_array.N_models,
                                                             N_array.pattern['num_segments'], FR_vector_signal_new,
                                                             DBS_pulse.t_vector, DBS_pulse.A, name_sorted_solution)

        if d["spectrum_trunc_method"] == 'Octave Band Method':
            if isinstance(d["n_Ranvier"], list):  # if different populations
                last_point = 0
                hf = h5py.File(os.environ['PATIENTDIR'] + '/' + d["Name_prepared_neuron_array"], 'r')
                lst_population_names = list(hf.keys())
                hf.close()
                for i in range(len(d["n_Ranvier"])):
                    logging.critical("in {} population".format(lst_population_names[i]))
                    last_point = convolute_signal_with_field_and_compute_ifft(d, DBS_pulse.Xs_signal_norm, N_array.N_models[i],
                                                                              N_array.pattern['num_segments'][i],
                                                                              DBS_pulse.FR_vector_signal, DBS_pulse.t_vector, DBS_pulse.A,
                                                                              name_sorted_solution,
                                                                              inx_st_oct=inx_start_octv, dif_axons=True,
                                                                              last_point=last_point)
            else:
                convolute_signal_with_field_and_compute_ifft(d, DBS_pulse.Xs_signal_norm, N_array.N_models,
                                                             N_array.pattern['num_segments'], DBS_pulse.FR_vector_signal,
                                                             DBS_pulse.t_vector, DBS_pulse.A, name_sorted_solution,
                                                             inx_st_oct=inx_start_octv, dif_axons=False, last_point=0)
    bar.next()
    '''Just to compute impedance in time, on if CPE is added or current-control mode'''

    if (d["CPE_activ"] == 1 or d[
        "current_control"] == 1) and cc_multicontact == False:  # and d["IFFT_ready"]==0:        #modify later
        from Field_IFFT_on_different_axons import compute_Z_ifft
        logging.critical("----- Calculating impedance -----\n")
        if d["spectrum_trunc_method"] == 'No Truncation' or d["Truncate_the_obtained_full_solution"] == 1:
            Imp_in_time = compute_Z_ifft(d, DBS_pulse.Xs_signal_norm, DBS_pulse.t_vector, DBS_pulse.A)
        elif d["spectrum_trunc_method"] == 'Octave Band Method':
            Imp_in_time = compute_Z_ifft(d, DBS_pulse.Xs_signal_norm, DBS_pulse.t_vector, DBS_pulse.A, i_start_octv=inx_start_octv)
        else:
            Imp_in_time = compute_Z_ifft(d, Xs_signal_norm_new, DBS_pulse.t_vector, A)
    # ===========================NEURON model simulation===========================#

    oss_plat_cont = os.getcwd()

    # we need to copy Axon_files to the patient folder to ensure stability
    src = oss_plat_cont + "/Axon_files"
    dst = os.environ['PATIENTDIR'] + "/Axon_files"

    if os.path.isdir(dst):
        os.chdir(os.environ['PATIENTDIR'])
        os.system('rm -fr Axon_files')
    os.makedirs(dst)

    copytree(src, dst, symlinks=False, ignore=None)
    os.chdir(dst)  # we now operate in Axon_files/ in the stim folder

    logging.critical("----- Estimating neuron activity -----")
    print("\n Probing action potentials")
    start_neuron = time.time()

    if d["Axon_Model_Type"] == 'McIntyre2002':
        with open(os.devnull, 'w') as FNULL:
            subprocess.call('nocmodl axnode.mod', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        with open(os.devnull, 'w') as FNULL:
            subprocess.call('nrnivmodl', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    elif d["Axon_Model_Type"] == 'Reilly2016':
        os.chdir("Reilly2016/")
        with open(os.devnull, 'w') as FNULL:
            subprocess.call('nrnivmodl', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    from Axon_files.NEURON_run import run_simulation_with_NEURON

    if d['Current_sets_simple'] == True:
        for k in range(len(Currents_to_check)):
            d["Ampl_scale"] = Currents_to_check[k][0]/d['Pulse_amp'][0]
            if isinstance(d["n_Ranvier"], list):
                Number_of_activated = 0
                last_point = 0
                for i in range(len(d["n_Ranvier"])):
                    Number_of_activated_population = run_simulation_with_NEURON(d, N_array, np.array(MRI_seg_param.MRI_shift), population_index=i, last_point=last_point, scaling_index=k)
                    Number_of_activated = Number_of_activated + Number_of_activated_population

                    last_point = N_array.pattern['num_segments'][i] * N_array.N_models[i] + last_point
            else:
                Number_of_activated = run_simulation_with_NEURON(d, N_array, np.array(MRI_seg_param.MRI_shift), scaling_index=k)
    else:
        if isinstance(d["n_Ranvier"], list):
            Number_of_activated = 0
            last_point = 0
            for i in range(len(d["n_Ranvier"])):
                Number_of_activated_population = run_simulation_with_NEURON(d, N_array, np.array(MRI_seg_param.MRI_shift), population_index=i, last_point=last_point)
                Number_of_activated = Number_of_activated + Number_of_activated_population

                last_point = N_array.pattern['num_segments'][i] * N_array.N_models[i] + last_point
        else:
            Number_of_activated = run_simulation_with_NEURON(d, N_array, np.array(MRI_seg_param.MRI_shift))
    bar.next()
    os.chdir(oss_plat_cont)

    # if isinstance(d["n_Ranvier"],list) and len(d["n_Ranvier"])>1:
    # with open(os.devnull, 'w') as FNULL: subprocess.call('python Visualization_files/Paraview_connections_activation.py', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    # if d['Show_paraview_screenshots']==1:
    #    subprocess.call('xdg-open "Images/Axon_activation.png"',shell=True)
    # else:
    # subprocess.call('python Visualization_files/Paraview_csv_activation.py', shell=True)
    # if d['Show_paraview_screenshots']==1:
    #   subprocess.call('xdg-open os.environ['PATIENTDIR']+"/Images/Activated_neurons.png"',shell=True)

    minutes = int((time.time() - start_neuron) / 60)
    secnds = int(time.time() - start_neuron) - minutes * 60
    logging.critical("----- NEURON calculations took {} min {} sec ----- \n".format(minutes, secnds))
    # logging.critical("----- NEURON calculations took ",minutes," min ",secnds," s -----\n")

    minutes = int((time.time() - start_simulation_run) / 60)
    secnds = int(time.time() - start_simulation_run) - minutes * 60
    total_seconds = time.time() - start_simulation_run
    logging.critical("----- Simulation run took {} min {} sec -----".format(minutes, secnds))
    # logging.critical("---Simulation run took ",minutes," min ",secnds," s ")

    # home_dir=os.path.expanduser("~")
    if d["Stim_side"] == 0:
        subprocess.call(['touch', os.environ['PATIENTDIR'] + '/success_rh.txt'])
    else:
        subprocess.call(['touch', os.environ['PATIENTDIR'] + '/success_lh.txt'])
    bar.finish()
    return True


# ===========================Master Study======================================#


'''Master study is useful for optimization and UQ studies'''
'''It will evaluate which truncation method and with how many frequencies we need to use for a particular setup'''
'''If adaptive mesh refinement is disabled, it will also estimate the difference in potential for 2 profile salome scripts (profile scripts should be created manually)'''


def run_master_study():
    '''Computations on initial mesh generated in SALOME'''

    Salome_profile_list = ['SNEX100', 'SNEX100_UQ2',
                           'SNEX100_UQ3']  # create 'SNEX100_UQ2_profile','SNEX100_UQ3_profile' and 'SNEX100_UQ4_profile' from 'SNEX100_profile' (in Electrode_files/) by varying mesh parameters
    Salome_best_profile = 'SNEX100_UQ4'  # name of the script with a highest mesh requirements (which might be too strict for iterative studies). Will be used as a benchmark.

    master_dict = {'Electrode_type': Salome_best_profile}
    total_seconds_bench, Number_of_activated_bench = run_full_model(master_dict)

    Result_diff_profiles = np.zeros((len(Salome_profile_list), 2), float)

    for i in range(len(Salome_profile_list)):
        master_dict = {'Electrode_type': Salome_profile_list[i], 'Segm_MRI_processed': 1, 'DTI_processed': 1,
                       'Init_neuron_model_ready': 1}  # some steps were already done
        Result_diff_profiles[i, :] = run_full_model(master_dict)

    logging.critical("Number of activated neurons, the benchmark and the profiles : ", Number_of_activated_bench,
                     Result_diff_profiles[:, 1])
    logging.critical("Time in seconds, the benchmark and the profiles : ", total_seconds_bench,
                     Result_diff_profiles[:, 0])
    logging.critical(
        "Choose the best option for 'Electrode_type'. Comment out '''Computations on initial mesh generated in SALOME''' section")

    raise SystemExit

    '''Spectrum truncation'''

    Number_of_activated_high_ampl = 0  # initizialization
    Number_of_activated_octaves = 0  # initizialization

    # first, run for the full spectrum
    master_dict = {"spectrum_trunc_method": 'No Truncation', "Truncate_the_obtained_full_solution": 0}
    Run_time, Number_of_activated_benchmark = run_full_model(master_dict)

    # secondly, extract results from the full solution on freqs. with high ampl. of FFT
    trunc_param_N_freqs = 25  # we will start with truncation to 25 frequencies
    while Number_of_activated_high_ampl != Number_of_activated_benchmark:  # here in the next version we could give a percentage of deviation
        master_dict = {"spectrum_trunc_method": 'High Amplitude Method', "Truncate_the_obtained_full_solution": 1,
                       "trunc_param": trunc_param_N_freqs}
        Run_time_high_ampl, Number_of_activated_high_ampl = run_full_model(master_dict)
        trunc_param_N_freqs = trunc_param_N_freqs + 25

    # now run the same study, but without extracting from the full solution
    master_dict = {"spectrum_trunc_method": 'High Amplitude Method', "Truncate_the_obtained_full_solution": 0,
                   "trunc_param": trunc_param_N_freqs - 25}
    Run_time_high_ampl, Number_of_activated_high_ampl = run_full_model(master_dict)

    def func_cutoff_oct(cutoff_oct,
                        *data_pass):  # defines the frequency after which octaves will be placed depending on the total number of freqs (the final total number can deviate by 1 frequency)
        from Inp_dict import d
        Sim_time = 1.0 / d["freq"]  # one pulse per simulation
        freq_max = d["freq"] * Sim_time / d["t_step"]

        return int(
            np.log2((freq_max - cutoff_oct) * np.sqrt(2.0) / d["freq"]) + cutoff_oct / d["freq"]) + 1 - N_freqs_oct

    from scipy.optimize import fsolve

    N_freqs_oct = 10  # we will start with truncation to 10 frequencies
    while Number_of_activated_octaves != Number_of_activated_benchmark:  # here in the next version we could give a percentage of deviation
        data_pass = (N_freqs_oct)
        cutoff_octaves = fsolve(func_cutoff_oct, 10, args=data_pass)

        master_dict = {"spectrum_trunc_method": 'Octave Band Method', "Truncate_the_obtained_full_solution": 0,
                       "trunc_param": cutoff_octaves}
        Run_time_oct, Number_of_activated_octaves = run_full_model(master_dict)
        N_freqs_oct = N_freqs_oct + 5
        if Run_time_oct > Run_time_high_ampl:
            logging.critical("Octave method is already slower than high amplitude method., no need to continue.")
            break

    if Run_time_oct < Run_time_high_ampl and Run_time_oct < Run_time:
        logging.critical("Octave truncation is the fastest")
    elif Run_time_oct > Run_time_high_ampl and Run_time_high_ampl < Run_time:
        logging.critical("High amplitude truncation is the fastest")
    else:
        logging.critical("Truncation did not accelerate the computations")

    logging.critical("Time for Full_run,High_ampl,Octaves: ", Run_time, Run_time_high_ampl, Run_time_oct)
    logging.critical("Activation for Full_run,High_ampl,Octaves: ", Number_of_activated_benchmark,
                     Number_of_activated_high_ampl, Number_of_activated_octaves)


# master_dict = {}
# run_full_model(master_dict)

if __name__ == '__main__':
    from datetime import datetime

    # datetime object containing current date and time
    now = datetime.now()
    date_and_time = now.strftime("%d-%m-%Y___%H-%M-%S")  # EUROPEAN format

    import logging
    import sys
    if len(sys.argv) > 1:
        cluster_run = bool(int(sys.argv[1]))
    else:
        cluster_run = False

    if len(sys.argv) > 2:
        #FF_dictionary = bool(int(sys.argv[2]))
        with open(os.environ['PATIENTDIR'] + '/netblend_dict_file.json', 'r') as fp:
            netblend_dictionary = json.load(fp)
            netblend_dictionary = netblend_dictionary['netblendict']
        fp.close()
    else:
        netblend_dictionary = False

    # if running in from MATLAB in singularity, you might need to change '/opt/Patient' -> os.environ['PATIENTDIR']
    logging.basicConfig(filename='/opt/Patient' + '/complete_log_' + date_and_time + '.log',
                       format='[%(asctime)s]:%(message)s', level=logging.ERROR)

    try:
        master_dict = {'cluster_run': cluster_run}
        run_full_model(master_dict, netblend_dictionary)
    except Exception:
        logging.exception("Fatal error in main loop")
