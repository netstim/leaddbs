# -*- coding: utf-8 -*-
"""

@author: Konstantin

These methods are used to scale the tissue-voltage response with the DBS power spectrum
and convert the solution to the time domain by means of Inverse Fourier Transformation

convolute_signal_with_field_and_compute_ifft - scaling and IFFT for the particular DBS protocol
convolute_signal_with_unit_field_and_compute_ifft - scaling and IFFT of unit solutions for each electrode contact
                                                  - The latter is used to get new field distributions with superposition
                                                  - in optimization and multiprotocol studies.

"""

import logging
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import os
import time as time_lib
from pandas import read_csv
import h5py

solution_sort = 'none'
Xs_signal_normalized = 'none'
t_vect = 'none'
FREQ_vector_signal = 'none'

def compute_Z_ifft(d, DBS_pulse, trunc_pulse = None):

    """ Computes impedance in time domain (only for 2-contact or contact-casing systems)"""

    # load computed impedance in the frequency domain
    Impedance_fr_get = read_csv(os.environ['PATIENTDIR']+'/Field_solutions/sorted_impedance.csv', delimiter=' ', header=None)
    Impedance_fr = Impedance_fr_get.values
    Z_Tr = np.vectorize(complex)(Impedance_fr[:,0], Impedance_fr[:,1])

    if d["spectrum_trunc_method"] == 'Octave Band Method':
        ## load the octave band frequencies where the impedance was computed
        #FR_vec_sign_octv = np.genfromtxt(os.environ['PATIENTDIR'] + '/Stim_Signal/FR_vector_signal_octaves' + str(d["trunc_param"] * 1.0)+'.csv', delimiter=' ')
        FR_vec_sign_octv = np.round(trunc_pulse.FR_vector_signal_new, 6)

        ## load the array that matches power spectrum frequencies with the octave band frequencies
        #Fr_corresp_ar = np.genfromtxt(os.environ['PATIENTDIR'] + '/Stim_Signal/Fr_corresp_array'+str(d["trunc_param"] * 1.0)+'.csv', delimiter=' ')
        Fr_corresp_ar = np.round(trunc_pulse.Fr_corresp_array, 6)

        # create the full size vectors in the frequency domain of the DBS power spectrum
        cutoff = int(np.ceil((DBS_pulse.t_vector.shape[0]+1)/2.))
        Z_Tr_full_real = np.zeros(cutoff, float)
        Z_Tr_full_imag = np.zeros(cutoff, float)

        # extrapolate the solution at octave band center frequencies to the whole DBS power spectrum
        stepper = 0
        for i_inx in range(Z_Tr.shape[0]):
            if i_inx >= trunc_pulse.inx_start_octv:
                # check which power spectrum frequencies belong to which octave
                rslt = np.where(Fr_corresp_ar[:,0] == np.round(FR_vec_sign_octv[i_inx], 6))
                step_octv = rslt[0].shape[0]   #size of the freq. pack in the octave

                Z_Tr_full_real[stepper:stepper+step_octv] = Z_Tr[i_inx].real
                Z_Tr_full_imag[stepper:stepper+step_octv] = Z_Tr[i_inx].imag
                stepper = stepper + step_octv
            else: # this case is for frequencies before the octave bands
                Z_Tr_full_real[stepper]=(Z_Tr[i_inx].real)
                Z_Tr_full_imag[stepper]=(Z_Tr[i_inx].imag)
                stepper += 1

        Z_Tr_full_complex = np.vectorize(complex)(Z_Tr_full_real,Z_Tr_full_imag)
        # Ohm's law (devided by A because we need a unit Xs_signal_normalized, not scaled with the signal ampl.)
        Z_conv = DBS_pulse.Xs_signal_norm * Z_Tr_full_complex / DBS_pulse.A
    elif d["spectrum_trunc_method"] == 'No Truncation':
        Z_conv = DBS_pulse.Xs_signal_norm * Z_Tr / DBS_pulse.A       # devided by A because we need a unit Xs_signal_normalized, not scaled with the signal ampl.
    else:
        logging.critical('For the selectred spetrum truncation method, the impedance might be highly underestimated')
        logging.critical('Please, consider using Octave Band Method')
        return None

    # now we do the inverse Fourier transformation to get the solution in time domain
    # following https://numpy.org/doc/stable/reference/generated/numpy.fft.ifft.html
    if np.mod(DBS_pulse.t_vector.shape[0], 2):  # if the FT vector is odd
        fv_conj = np.conjugate(Z_conv[-1:0:-1])
    else:  # if the FT vector is even
        fv_conj = np.conjugate(Z_conv[-2:0:-1])

    Y = np.concatenate((Z_conv, fv_conj), axis=0)
    Signal_t_Zconv = np.fft.ifft(Y).real
    logging.critical("Max impedance: {}".format(Signal_t_Zconv.max()))

    plt.figure(111122112)
    plt.plot(DBS_pulse.t_vector, Signal_t_Zconv.real)
    plt.xlim(0.000,d["T"]*5)
    plt.grid(True)
    plt.xlabel('t, sec')
    plt.ylabel('Zreal, Ohm')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.savefig(os.environ['PATIENTDIR']+'/Images/Z_convoluted.png', format='png', dpi=1000)
    np.savetxt(os.environ['PATIENTDIR']+'/Field_solutions/Z_R_TimeDomain.csv', Signal_t_Zconv.real, delimiter=" ")
    np.savetxt(os.environ['PATIENTDIR']+'/Field_solutions/Z_Im_TimeDomain.csv', Signal_t_Zconv.imag, delimiter=" ")

    return Signal_t_Zconv

#get electric potential in time over the axon
def convolute_and_ifft(d,last_point,i_axon,num_segments,ind_trunc, contact_i,output):

    # d['t_steps_trunc'] - number of steps after which the DBS signal is truncated (defined autom. in Launcher_OSS)
    # d['phi'] - signal phase shift in s
    # d['T'] - pulse width in s
    # d["spectrum_trunc_method"]
    # d["Truncate_the_obtained_full_solution"] - true only if we want to use truncation of already computed full power spectrum solution (useful for benchmarks)

    global solution_sort
    global Xs_signal_normalized
    global t_vect
    global FREQ_vector_signal

    N_freq = FREQ_vector_signal.shape[0]
    Signal_t_conv = np.zeros((int(num_segments), t_vect.shape[0]), float)

    for i_point in range(int(num_segments)):
        global_i_point=int(num_segments * i_axon + i_point) # index of the neuron compartment in Vert_of_Neural_model_NEURON.csv
        # get the real and imag voltage responses at the neuron compartment defined by 'global_i_point'
        Xs_Tr = np.vectorize(complex)(solution_sort[global_i_point * N_freq:(global_i_point * N_freq + N_freq), 3], solution_sort[global_i_point * N_freq:(global_i_point * N_freq + N_freq), 4])         #real and im parts for the first point in VTA

        if d["spectrum_trunc_method"] == 'No Truncation':
            # scale the voltage-tissue response by the DBS power spectrum (Xs_signal_normalized)
            Xs_conv = Xs_signal_normalized * Xs_Tr

        elif d["spectrum_trunc_method"] == 'Cutoff Method':
            # restore the full power spectrum vector shape and fill accordingly
            cutoff = int(np.ceil((t_vect.shape[0] + 1)/2.))
            Xs_signal_full = np.complex(1.0, 0.0) * np.zeros(cutoff, float)
            Xs_Tr_full = np.complex(1.0, 0.0) * np.ones(cutoff, float)

            if d["Truncate_the_obtained_full_solution"] == 0:
                np.put(Xs_Tr_full, np.arange(int(ind_trunc)), Xs_Tr)
                np.put(Xs_signal_full, np.arange(int(ind_trunc)), Xs_signal_normalized)

            if d["Truncate_the_obtained_full_solution"] == 1:
                Xs_signal_normalized_trunc = Xs_signal_normalized[np.arange(int(ind_trunc))]
                Xs_Tr_trunc = Xs_Tr[np.arange(int(ind_trunc))]

                np.put(Xs_Tr_full, np.arange(int(ind_trunc)), Xs_Tr_trunc)
                np.put(Xs_signal_full, np.arange(int(ind_trunc)), Xs_signal_normalized_trunc)

            Xs_conv = Xs_signal_full * Xs_Tr_full

        elif d["spectrum_trunc_method"] == 'High Amplitude Method':

            # restore the full power spectrum vector shape and fill accordingly
            cutoff = int(np.ceil((t_vect.shape[0]+1)/2.))
            Xs_signal_full = np.complex(1.0,0.0)*np.zeros(cutoff,float)
            Xs_Tr_full = np.complex(1.0,0.0)*np.ones(cutoff,float)

            if d["Truncate_the_obtained_full_solution"] == 0:
                np.put(Xs_Tr_full,ind_trunc.astype(int),Xs_Tr)
                np.put(Xs_signal_full,ind_trunc.astype(int),Xs_signal_normalized)

            if d["Truncate_the_obtained_full_solution"] == 1:
                Xs_signal_normalized_trunc=Xs_signal_normalized[ind_trunc.astype(int)]
                Xs_Tr_trunc=Xs_Tr[ind_trunc.astype(int)]

                np.put(Xs_Tr_full,ind_trunc.astype(int),Xs_Tr_trunc)
                np.put(Xs_signal_full,ind_trunc.astype(int),Xs_signal_normalized_trunc)

            Xs_conv = Xs_signal_full * Xs_Tr_full

        # now we do the inverse Fourier transformation to get the solution in time domain
        # following https://numpy.org/doc/stable/reference/generated/numpy.fft.ifft.html
        if np.mod(t_vect.shape[0], 2):  # if the FT vector is odd
            fv_conj = np.conjugate(Xs_conv[-1:0:-1])
        else:  # if the FT vector is even
            fv_conj = np.conjugate(Xs_conv[-2:0:-1])

        Y = np.concatenate((Xs_conv, fv_conj), axis=0)
        Signal_t_conv[i_point,:] = np.fft.ifft(Y).real

        if global_i_point+last_point == 0:
            plt.figure(11111234)
            plt.plot(t_vect, Signal_t_conv[i_point, :])
            plt.xlim(0.000, d['T'] * 5 + d['phi'])
            plt.grid(True)
            plt.xlabel('t, sec')
            plt.ylabel('Potential, V')
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            if contact_i == -1:
                plt.savefig(os.environ['PATIENTDIR']+'/Images/Signal_convoluted_1st_point.png', format='png', dpi=500)
            else:
                plt.savefig(os.environ['PATIENTDIR'] + '/Images/1A_scaled_convoluted_1st_point_' + str(contact_i) + '.png', format='png', dpi=500)

    # IMPORTANT: the file is indexed by the global point index of the last compartment of the neuron (axon)
    if contact_i == -1:
        np.save(os.environ['LGFDIR']+'/Axons_in_time/Signal_t_conv'+str(global_i_point + last_point), Signal_t_conv[:d['t_steps_trunc']])
    else:
        np.save(os.environ['LGFDIR'] + '/Axons_in_time/Signal_t_conv' + str(global_i_point + last_point) + "_" + str(contact_i), Signal_t_conv[:d['t_steps_trunc']])


    output.put(i_axon)

#get electric potential in time over the axon when using octave band method for the frequency spectrum approximation
def convolute_and_ifft_octaves(d,last_point,i_axon,num_segments,Fr_corresp_ar,i_start_octv,contact_i,output):

    # d['t_steps_trunc'] - number of steps after which the DBS signal is truncated (defined autom. in Launcher_OSS)
    # d['phi'] - signal phase shift in s
    # d['T'] - pulse width in s

    global solution_sort_octv
    global Xs_signal_normalized
    global t_vect
    global FREQ_vector_signal
    global FR_vec_sign_octv

    N_freq_octv = FR_vec_sign_octv.shape[0]
    Signal_t_conv = np.zeros((int(num_segments), t_vect.shape[0]),float)

    for i_point in range(int(num_segments)):
        global_i_point = int(num_segments * i_axon + i_point) # index of the neuron compartment in Vert_of_Neural_model_NEURON.csv
        # get the real and imag voltage responses at the neuron compartment defined by 'global_i_point'
        if contact_i == -1:
            Xs_Tr = np.vectorize(complex)(solution_sort_octv[global_i_point * N_freq_octv:(global_i_point * N_freq_octv + N_freq_octv), 3], solution_sort_octv[global_i_point*N_freq_octv:(global_i_point*N_freq_octv+N_freq_octv), 4])
        else:
            Xs_Tr = np.vectorize(complex)(
                solution_sort_octv[global_i_point * N_freq_octv:(global_i_point * N_freq_octv + N_freq_octv), 0],
                solution_sort_octv[global_i_point * N_freq_octv:(global_i_point * N_freq_octv + N_freq_octv), 1])

        # now we need to extrapolate the solution at octave band center frequencies to the whole DBS power spectrum
        Xs_Tr_full_real = np.zeros(FREQ_vector_signal.shape[0], float)
        Xs_Tr_full_imag = np.zeros(FREQ_vector_signal.shape[0], float)
        stepper = 0
        for i_inx in range(Xs_Tr.shape[0]):
            if i_inx >= i_start_octv:
                # check which power spectrum frequencies belong to which octave
                rslt = np.where(Fr_corresp_ar[:,0] == np.round(FR_vec_sign_octv[i_inx], 6))
                step_octv = rslt[0].shape[0]   # size of the freq. pack in the octave

                Xs_Tr_full_real[stepper:stepper + step_octv] = Xs_Tr[i_inx].real
                Xs_Tr_full_imag[stepper:stepper + step_octv] = Xs_Tr[i_inx].imag
                stepper = stepper + step_octv
            else:   # this case is for frequencies before the octave bands
                Xs_Tr_full_real[stepper] = Xs_Tr[i_inx].real
                Xs_Tr_full_imag[stepper] = Xs_Tr[i_inx].imag
                stepper += 1

        # store one solution for debugging
        if global_i_point == 0:
            np.savetxt(os.environ['PATIENTDIR']+'/Field_solutions/Xs_Tr_full_real.csv', Xs_Tr_full_real, delimiter=" ")
            np.savetxt(os.environ['PATIENTDIR']+'/Field_solutions/Xs_Tr_full_imag.csv', Xs_Tr_full_imag, delimiter=" ")

        # scale the voltage-tissue response by the DBS power spectrum (Xs_signal_normalized)
        Xs_conv = Xs_signal_normalized * np.vectorize(complex)(Xs_Tr_full_real,Xs_Tr_full_imag)

        # now we do the inverse Fourier transformation to get the solution in time domain
        # following https://numpy.org/doc/stable/reference/generated/numpy.fft.ifft.html
        if np.mod(t_vect.shape[0], 2):  # if the FT vector is odd
            fv_conj = np.conjugate(Xs_conv[-1:0:-1])
        else:  # if the FT vector is even
            fv_conj = np.conjugate(Xs_conv[-2:0:-1])

        Y = np.concatenate((Xs_conv, fv_conj), axis=0)
        Signal_t_conv[i_point,:] = np.fft.ifft(Y).real

        if global_i_point+last_point == 0:
            plt.figure(11111234)
            plt.plot(t_vect, Signal_t_conv[i_point, :])
            plt.xlim(0.000, d['T'] * 5 + d['phi'])
            plt.grid(True)
            plt.xlabel('t, sec')
            plt.ylabel('Potential, V')
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            if contact_i == -1:
                plt.savefig(os.environ['PATIENTDIR']+'/Images/Signal_convoluted_1st_point.png', format='png', dpi=500)
            else:
                plt.savefig(os.environ['PATIENTDIR'] + '/Images/1A_scaled_convoluted_1st_point_' + str(contact_i) + '.png', format='png', dpi=500)

    # IMPORTANT: the file is indexed by the global point index of the last compartment of the neuron (axon)
    if contact_i == -1:
        np.save(os.environ['LGFDIR']+'/Axons_in_time/Signal_t_conv'+str(global_i_point + last_point), Signal_t_conv[:d['t_steps_trunc']])
    else:
        np.save(os.environ['LGFDIR'] + '/Axons_in_time/Signal_t_conv' + str(global_i_point + last_point) + "_" + str(contact_i), Signal_t_conv[:d['t_steps_trunc']])

    output.put(i_axon)


def convolute_signal_with_field_and_compute_ifft(d, DBS_pulse, models_in_population, number_of_segments, name_sol, trunc_pulse=None, dif_axons=False, last_point=0):

    start_ifft = time_lib.time()

    # global variables to process arrays in parallel
    global FREQ_vector_signal
    global t_vect
    global Xs_signal_normalized

    if d["spectrum_trunc_method"] == 'Octave Band Method' or d["spectrum_trunc_method"] == 'No Truncation':
        Xs_signal_normalized = DBS_pulse.Xs_signal_norm
        FREQ_vector_signal = DBS_pulse.FR_vector_signal  # the IFFT results will be returned for the whole spectrum
    else:
        Xs_signal_normalized = trunc_pulse.Xs_signal_norm_new
        FREQ_vector_signal = trunc_pulse.FR_vector_signal_new  # the IFFT results will be returned at selected freqs.

    t_vect = DBS_pulse.t_vector

    if d["spectrum_trunc_method"] == 'Octave Band Method':

        # array that points which power spectrum freqs belong to which octave bands
        Fr_corresp_arr = np.round(trunc_pulse.Fr_corresp_array, 6)

        global FR_vec_sign_octv
        FR_vec_sign_octv = trunc_pulse.FR_vector_signal_new
        global solution_sort_octv  # stores V-distribution on neurons at octave band center frequencies

        if dif_axons == True:
            hf = h5py.File(name_sol[:-4]+'.h5', 'r')
            solution_sort_octv_full = hf.get('dataset_1')
            solution_sort_octv_full = np.array(solution_sort_octv_full)
            hf.close()
            solution_sort_octv = solution_sort_octv_full[last_point * FR_vec_sign_octv.shape[0]:(last_point + models_in_population * number_of_segments)*FR_vec_sign_octv.shape[0], :]
        else:
            hf = h5py.File(name_sol[:-4]+'.h5', 'r')
            solution_sort_octv = hf.get('dataset_1')
            solution_sort_octv = np.array(solution_sort_octv)
            hf.close()
    else:
        global solution_sort  # stores V-distribution on neurons at (some) frequencies of the DBS power spectrum

        if dif_axons == True:
            hf = h5py.File(name_sol[:-4]+'.h5', 'r')
            solution_sort_full = hf.get('dataset_1')
            solution_sort_full = np.array(solution_sort_full)
            hf.close()
            solution_sort = solution_sort_full[last_point*FREQ_vector_signal.shape[0]:(last_point+models_in_population*number_of_segments)*FREQ_vector_signal.shape[0],:]
        else:
            hf = h5py.File(name_sol[:-4]+'.h5', 'r')
            solution_sort = hf.get('dataset_1')
            solution_sort = np.array(solution_sort)
            hf.close()

        if d["spectrum_trunc_method"] == 'No Truncation':
            ind_trunc = 0  # not needed if the whole spectrum was computed
        else:
            ind_trunc = trunc_pulse.trunc_ind  # int if cutoff method, numpy array if 'High Amplitude Method'


    # to check the progress
    axons_quart = [int(models_in_population/4.0), int(2*models_in_population/4.0), int(3*models_in_population/4.0)]

    j = 0  # index of the process
    i = 0  # index of the population
    contact_i = -1  # we will not need to index the solutions by contact

    # conduct IFFT in parallel for multiple neurons
    # the output is saved directly by the process
    while i < models_in_population:
        output = mp.Queue()
        proc = []
        while j < d["number_of_processors"] and i < models_in_population:

            if d["spectrum_trunc_method"] == 'Octave Band Method':
                processes = mp.Process(target=convolute_and_ifft_octaves, args=(d, last_point, i, number_of_segments, Fr_corresp_arr, trunc_pulse.inx_start_octv, contact_i, output))
            else:
                processes = mp.Process(target=convolute_and_ifft, args=(d, last_point, i, number_of_segments, ind_trunc, contact_i, output))
            proc.append(processes)
            j += 1
            i += 1
            if i in axons_quart:
                logging.critical("{}% of neuron models were processed".format(int(i*100/models_in_population)+1))
        j = 0  # new pack
        for p in proc:
            p.start()
        for p in proc:
            p.join()

    if d["spectrum_trunc_method"] == 'Octave Band Method':
        del solution_sort_octv
        if dif_axons == True:
            del solution_sort_octv_full
    else:
        del solution_sort
        if dif_axons == True:
            del solution_sort_full

    last_point = last_point + i * number_of_segments

    minutes = int((time_lib.time() - start_ifft)/60)
    secnds = int(time_lib.time() - start_ifft)-minutes*60
    logging.critical("----- Signal scaling and IFFT took {} min {} sec -----\n".format(minutes, secnds))

    return last_point


def convolute_signal_with_unit_field_and_compute_ifft(d, DBS_pulse, models_in_population, number_of_segments, name_sol, trunc_pulse=None, dif_axons=False, last_point=0):
    # here we don't want to pass large arrays, so we use global variables
    start_ifft = time_lib.time()
    if d["spectrum_trunc_method"] != 'Octave Band Method':
        logging.critical('Field superposition and scaling is implemented only for the Octave Band FFEM')
        logging.critical('Consider it as the truncation method')
        raise SystemExit

    # global variables to process arrays in parallel
    global FREQ_vector_signal
    global t_vect
    global Xs_signal_normalized
    global FR_vec_sign_octv
    global solution_sort_octv

    FR_vec_sign_octv = trunc_pulse.FR_vector_signal_new
    Xs_signal_normalized = DBS_pulse.Xs_signal_norm
    FREQ_vector_signal = DBS_pulse.FR_vector_signal  # the IFFT results will be returned for the whole spectrum
    t_vect = DBS_pulse.t_vector

    # array that points which power spectrum freqs belong to which octave bands
    Fr_corresp_arr = np.round(trunc_pulse.Fr_corresp_array, 6)

    if dif_axons == True:
        hf = h5py.File(name_sol[:-4] + '.h5', 'r')
        solution_sort_octv_full = hf.get('dataset_1')
        solution_sort_octv_full = np.array(solution_sort_octv_full)
        hf.close()
        # solution_sort_octv=solution_sort_octv_full[last_point*FR_vec_sign_octv.shape[0]:(last_point+models_in_population*number_of_segments)*FR_vec_sign_octv.shape[0],:]
        solution_over_contacts = solution_sort_octv_full[last_point * FR_vec_sign_octv.shape[0]:(last_point + models_in_population * number_of_segments) *
                                                                                                FR_vec_sign_octv.shape[0], :]
    else:
        hf = h5py.File(name_sol[:-4] + '.h5', 'r')
        solution_over_contacts = hf.get('dataset_1')
        solution_over_contacts = np.array(solution_over_contacts)
        hf.close()

    solution_sort_octv = np.zeros((solution_over_contacts.shape[0], 2), float)
    N_contacts = solution_over_contacts.shape[1] - 1  # all (the last is the frequency)

    # apply IFFT for each contact-solution. Scaling will be done when solving the cable equations
    for contact_i in range(N_contacts):

        for point_i in np.arange(0, (models_in_population * number_of_segments) * FR_vec_sign_octv.shape[0],
                                 FR_vec_sign_octv.shape[0]):
            solution_sort_octv[point_i:(point_i + FR_vec_sign_octv.shape[0]), 0] = solution_over_contacts[point_i:(
                    point_i + FR_vec_sign_octv.shape[0]), contact_i]

        logging.critical("IFFT for contact: {}".format(contact_i))
        j = 0  # index of the process
        i = 0  # index of the population

        # just to check the progress
        axons_quart = [int(models_in_population / 4.0), int(2 * models_in_population / 4.0),
                       int(3 * models_in_population / 4.0)]

        # conduct IFFT in parallel for multiple neurons
        # the output is saved directly by the process
        while i < models_in_population:
            output = mp.Queue()  # defines an output queue
            proc = []
            while j < d["number_of_processors"] and i < (models_in_population):

                processes = mp.Process(target=convolute_and_ifft_octaves, args=(
                    d, last_point, i, number_of_segments, Fr_corresp_arr, trunc_pulse.inx_start_octv, contact_i, output))

                proc.append(processes)
                j += 1
                i += 1
                # if i in axons_quart:
                # print(int(i*100/models_in_population)+1,"% of neuron models were processed")
            j = 0  # new pack
            for p in proc:
                p.start()
            for p in proc:
                p.join()
            # convoluted_lists = [output.get() for p in proc]  # the files are already created

    del solution_sort_octv
    if dif_axons == True:
        del solution_sort_octv_full

    last_point = last_point + i * number_of_segments

    minutes = int((time_lib.time() - start_ifft) / 60)
    secnds = int(time_lib.time() - start_ifft) - minutes * 60
    logging.critical("----- Signal scaling and IFFT took {} min {} sec -----\n".format(minutes, secnds))

    return last_point



