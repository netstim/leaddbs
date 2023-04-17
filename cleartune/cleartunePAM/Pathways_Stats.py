'''
    By K. Butenko
    These functions allow to retrieve data about simulated pathways
    (after Kuncel pre-filtering, but keeping original number of tracts!) and stimulation protocol
'''
import matplotlib.pyplot as plt
import numpy as np
import os
import h5py
#import seaborn as sns
#sns.set()
import sys

os.environ['PATIENTDIR'] = '/opt/Patient'
sys.path.insert(1, os.environ['PATIENTDIR'])


def get_simulated_pathways(side):

    ''' based on Lead-DBS files for OSS-DBS, get names of Kuncel filtered pathways and original number of fibers per pathway'''

    # later we should store connectome name in GUI_inp_dict.py
    # load parameters from the file prepared by Lead-DBS
    #file_inp = h5py.File(os.environ['PATIENTDIR'] + '/oss-dbs_parameters.mat', mode='r')

    file_inp = h5py.File(os.environ['STIMDIR'] + '/oss-dbs_parameters.mat', mode='r')
    array_ascii = file_inp['settings']['connectome'][:]
    list_ascii = []
    for i in range(array_ascii.shape[0]):
        list_ascii.append(array_ascii[i][0])
    # list_ascii = map(lambda s: s.strip(), list_ascii)
    Connectome_name = ''.join(chr(i) for i in list_ascii)

    Projections = []
    for i in range(len(file_inp['settings']['connectomeTractNames'][0])):
        ext_string = file_inp[file_inp['settings']['connectomeTractNames'][0][i]]
        list_ascii = []
        for i in range(ext_string.shape[0]):
            list_ascii.append(ext_string[i][0])
        # list_ascii = map(lambda s: s.strip(), list_ascii)
        projection_name = ''.join(chr(i) for i in list_ascii)
        # print(projection_name)
        Projections.append(projection_name)

    # It will always be Multi-Tract
    # 'Multi-tract' connectomes contain multiple pathways in separate .mat files
    if 'Multi-Tract' in Connectome_name:
        #Full_paths = [
        #    os.environ['PATIENTDIR'] + '/' + Connectome_name.rsplit(' ', 1)[1] + '/data' + str(side + 1) + '.mat']
        Full_paths = [
            os.environ['STIMDIR'] + '/' + Connectome_name.rsplit(' ', 1)[1] + '/data' + str(side + 1) + '.mat']

        file = h5py.File(Full_paths[0], mode='r')
        number_original = []
        name_original = []
        # run this in loop, later you will check if you actually have this pathways in Summary_status.h5
        for projection_name in Projections:
            if 'origNum' in file[projection_name]:
                number_original.append(file[projection_name]['origNum'][:][0][0])
                name_original.append(projection_name)
            #else: # that means no fibers of the pathway were passed from Lead-DBS
            #    number_original.append(None)

    else:
        print('This connectome has only one pathway, exiting')
        raise SystemExit

    return name_original, number_original

def get_current_protocol(index_side):

    ''' based on Lead-DBS files for OSS-DBS, get the simulated current protocol (in A)'''

    file = h5py.File(os.environ['STIMDIR'] + '/oss-dbs_parameters.mat', mode='r')

    # if file['settings']['current_control'][0][index_side] != 1:
    #    print('The imported protocol is not current-controlled!')
    #    raise SystemExit

    Pulse_amp = file['settings']['Phi_vector'][:, index_side]
    Pulse_amp = Pulse_amp * 0.001  # because Lead-DBS uses mA as the input, switch to A here for consistency

    Pulse_amp = list(Pulse_amp)
    import math
    for i in range(len(Pulse_amp)):
        if math.isnan(Pulse_amp[i]):
            Pulse_amp[i] = 0.0  # IMPORTANT: not grounding here, but a 0A contact

    return Pulse_amp