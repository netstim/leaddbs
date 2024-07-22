import os
import sys
import numpy as np

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
from tensorflow.keras.models import load_model

SIDE_SUFFIX = ['_rh','_lh']

'''
    By K. Butenko
    MATLAB interfaced function to get pathway percent activation using a pre-trained ANN model
'''
def call_ANN(stim_folder, stim_vector, pathway, side):

    approx_model = load_model(os.path.join(stim_folder, 'NB' + str(SIDE_SUFFIX[side]), 'ANN_approved_model_' + pathway))
    stim_array = np.reshape(np.array(stim_vector), (-1, len(stim_vector)))
    activation_profile = approx_model.predict(stim_array, verbose=0)
    activation_profile = activation_profile[0]  # get the actual array

    print(activation_profile)  # this will be captured by MATLAB

if __name__ == '__main__':

    # called from MATLAB
    stim_dir = sys.argv[1]
    pathway = sys.argv[2]
    side = int(sys.argv[3])
    # the rest is assumed to be tested currents (should be in A!)
    stim_vector = [float(i) * 0.001 for i in sys.argv[4:]]

    call_ANN(stim_dir,stim_vector,pathway,side)