# -*- coding: utf-8 -*-
"""

@author: Butenko K.

rearrange_Inp_dict converts from user-friendly units to the platform's format

"""


def rearrange_Inp_dict(dict_from_GUI):
    # first, change empty lines to 0
    for key in dict_from_GUI:
        if dict_from_GUI[key] == '' or dict_from_GUI[key] == '0':
            dict_from_GUI[key] = 0

    import numpy as np
    # second, some angles will be changed to rads
    if 'Global_rot' in dict_from_GUI.keys():
        if dict_from_GUI['Global_rot'] == 1:
            dict_from_GUI['alpha_array_glob'] = [i * np.pi / (180.0) for i in dict_from_GUI['alpha_array_glob']]
            dict_from_GUI['beta_array_glob'] = [i * np.pi / (180.0) for i in dict_from_GUI['beta_array_glob']]
            dict_from_GUI['gamma_array_glob'] = [i * np.pi / (180.0) for i in dict_from_GUI['gamma_array_glob']]
            dict_from_GUI['YZ_angles'], dict_from_GUI['XY_angles'], dict_from_GUI['ZX_angles'] = (0, 0, 0)
        else:
            dict_from_GUI['YZ_angles'] = [i * np.pi / (180.0) for i in dict_from_GUI['YZ_angles']]
            dict_from_GUI['XY_angles'] = [i * np.pi / (180.0) for i in dict_from_GUI['XY_angles']]
            dict_from_GUI['ZX_angles'] = [i * np.pi / (180.0) for i in dict_from_GUI['ZX_angles']]
            dict_from_GUI['alpha_array_glob'], dict_from_GUI['beta_array_glob'], dict_from_GUI['gamma_array_glob'] = (0, 0, 0)

    # put some variables to standard units (mm->m, ms->s and so on)
    dict_from_GUI['T'] = dict_from_GUI['T'] / 1000000.0
    dict_from_GUI['t_step'] = dict_from_GUI['t_step'] / 1000000.0
    dict_from_GUI['phi'] = dict_from_GUI['phi'] / 1000000.0

    # forcing to integer
    if dict_from_GUI['spectrum_trunc_method'] == 'Octave Band Method':
        dict_from_GUI['trunc_param'] = float(dict_from_GUI['trunc_param'])
    else:
        dict_from_GUI['trunc_param'] = int(dict_from_GUI['trunc_param'])

    # one value list to int
    if isinstance(dict_from_GUI['n_Ranvier'], list):
        if len(dict_from_GUI['n_Ranvier']) == 1:
            dict_from_GUI['n_Ranvier'] = int(dict_from_GUI['n_Ranvier'][0])

    if isinstance(dict_from_GUI['diam_fib'], list):
        if len(dict_from_GUI['diam_fib']) == 1:
            dict_from_GUI['diam_fib'] = float(dict_from_GUI['diam_fib'][0])

    if isinstance(dict_from_GUI['Aprox_geometry_center'], list):
        if len(dict_from_GUI['Aprox_geometry_center']) == 1:
            dict_from_GUI['Aprox_geometry_center'] = 0

    if not (isinstance(dict_from_GUI['Approximating_Dimensions'], list)):
        dict_from_GUI['Approximating_Dimensions'] = [0]

    if dict_from_GUI['Neuron_model_array_prepared'] == 0:
        dict_from_GUI['Name_prepared_neuron_array'] = 0

    # switch from percents
    if dict_from_GUI['Skip_mesh_refinement'] == 0:
        dict_from_GUI['rel_div_current'] = dict_from_GUI['rel_div_current'] / 100.0
        dict_from_GUI['rel_div_CSF'] = dict_from_GUI['rel_div_CSF'] / 100.0
        dict_from_GUI['rel_div'] = dict_from_GUI['rel_div'] / 100.0
        dict_from_GUI['Adaptive_frac_div'] = dict_from_GUI['Adaptive_frac_div'] / 100
        # dict_from_GUI['CSF_frac_div']=dict_from_GUI['CSF_frac_div']/100

    return dict_from_GUI
