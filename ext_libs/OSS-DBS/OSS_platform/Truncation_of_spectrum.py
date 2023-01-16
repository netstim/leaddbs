# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 15:53:14 2019

@author: Konstantin

Truncated_spectrum class contains methods necessary to describe a truncated DBS signal in the frequency domain.
It also stores info about the original signal and data necessary to 'map back' to the full DBS power spectrum.

Three different truncation methods are offered ('Octave Band Method','Cutoff Method','High Amplitude Method')
Preferred is 'Octave Band Method'. For detail, see Butenko et al, EMBC, 2019

"""
import numpy as np
import matplotlib.pyplot as plt
import os
import logging

class Truncated_spectrum(object):
    def __init__(self, inp_dict, DBS_pulse):

        self.trunc_method = inp_dict["spectrum_trunc_method"] #  'Octave Band Method' is recommended, see Butenko et al, EMBC proceedings, 2019
        self.trunc_param = inp_dict["trunc_param"] #  freq (Hz) after which octave bands are deployed,otherwise the total number of freqs
        self.rep_rate = inp_dict["freq"]         #  repetition rate of the DBS signal (usually 130-184 Hz)
        self.FR_vector_signal = DBS_pulse.FR_vector_signal  # vector of frequencies in the DBW power spectrum
        self.Xs_signal_norm = DBS_pulse.Xs_signal_norm      # vector of magnitude of the power
        self.t_vector = DBS_pulse.t_vector                  # time vector in s
        self.pw = DBS_pulse.pw                              # pulse width in s

    def get_specific_freq(self):

        """ Simple truncation methods, outdated """

        logging.critical("----- Truncating the frequency spectrum -----")
        Xs_signal_full = np.complex(1.0, 0.0) * np.zeros(self.Xs_signal_norm.shape[0], float)

        if self.trunc_method == 'High Amplitude Method':     #choosing freq. where FT has the highest amplitudes, here self.trunc_param is the number of frequencies in the pack
            ind_max = np.sort(np.argpartition(np.absolute(self.Xs_signal_norm), -1*self.trunc_param)[-1*self.trunc_param:])
            self.FR_vector_signal_new = self.FR_vector_signal[ind_max]
            self.Xs_signal_norm_new = self.Xs_signal_norm[ind_max]
            np.savetxt(os.environ['PATIENTDIR']+'/Stim_Signal/Indices_high_ampl.csv', ind_max, delimiter=" ")
            np.put(Xs_signal_full, ind_max.astype(int), self.Xs_signal_norm)
            self.trunc_ind = ind_max # numpy array!
        elif self.trunc_method == 'Cutoff Method':
            inx_cutoff = self.trunc_param
            self.FR_vector_signal_new = self.FR_vector_signal[:inx_cutoff]
            self.Xs_signal_norm_new = self.Xs_signal_norm[:inx_cutoff]
            np.savetxt(os.environ['PATIENTDIR'] + '/Stim_Signal/Indices_cutoff.csv',
                       np.array([self.Xs_signal_norm.shape[0], inx_cutoff]), delimiter=" ")
            self.trunc_ind = int(inx_cutoff)
            np.put(Xs_signal_full, np.arange(int(inx_cutoff)), self.Xs_signal_norm)

        '''To create the truncated signal'''
        Xs_ones = np.complex(1.0, 0.0) * np.ones(self.Xs_signal_norm.shape[0], float)
        Xs_signal_unit = Xs_signal_full * Xs_ones

        if np.mod(self.t_vector.shape[0], 2):  # if the FT vector is odd
            fv_conj = np.conjugate(Xs_signal_unit[-1:0:-1])
        else:  # if the FT vector is even
            fv_conj = np.conjugate(Xs_signal_unit[-2:0:-1])

        Y = np.concatenate((Xs_signal_unit, fv_conj), axis=0)
        Signal_t_conv = np.fft.ifft(Y).real

        plt.figure(111112341)
        plt.plot(self.t_vector, Signal_t_conv)
        plt.xlim(0.000, self.pw * 5)
        plt.grid(True)
        plt.xlabel('t, sec')
        plt.ylabel('Signal amplitude (A or V)')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.savefig(os.environ['PATIENTDIR']+'/Images/Signal_convoluted_truncated.png', format='png', dpi=1000)

    def get_octave_band_freqs(self):

        """ Compute center frequencies of octave bands based on the repetition rate (self.rep_rate)
         and the specified frequency, after which octave bands are deployed (self.trunc_param) """

        self.FR_vector_signal_new = -1 * np.ones((self.FR_vector_signal.shape[0]), float)  #  this vector will contain frequencies for the FEM problem
        self.Fr_octave_vector = []  #  this vector will contain octave bands defined by DBS repetition rate
        inx = 0

        one_sample_10kHz = False
        if one_sample_10kHz == True:
            self.trunc_param = 1.0

        for fr in self.FR_vector_signal:

            if fr <= self.trunc_param:  # include all frequencies before self.trunc_param (defines the start of the octaves)
                self.FR_vector_signal_new[inx] = fr
                inx += 1
            else:
                octave_scale = 1
                base_octaves = self.FR_vector_signal_new[inx-1]
                while (base_octaves + (self.rep_rate * octave_scale) / np.sqrt(2.0)) < self.FR_vector_signal[-1]:

                    Octave_freq = base_octaves + (self.rep_rate * octave_scale / np.sqrt(2.0))  # compute the octave center frequency
                    Fr_in_octave = []

                    # dirty implementation
                    if one_sample_10kHz == True and Octave_freq < 10000.0:
                        Octave_freq = 5000.0
                        for fr_init in self.FR_vector_signal:
                            if fr_init < 10000.0:
                                Fr_in_octave.append(fr_init)
                                if Octave_freq not in self.Fr_octave_vector:
                                    self.Fr_octave_vector.append(Octave_freq)

                        if ~np.isin(5000.0,self.FR_vector_signal_new):
                            self.FR_vector_signal_new[inx] = Octave_freq
                    else:
                        for fr_init in self.FR_vector_signal:
                            # check which/if the original frequencies are in the octave
                            if fr_init >= base_octaves + np.round((Octave_freq - base_octaves) / np.sqrt(2), 2) and fr_init < base_octaves + np.round((Octave_freq - base_octaves) * np.sqrt(2),2):
                                Fr_in_octave.append(fr_init)
                                self.FR_vector_signal_new[inx] = Octave_freq  # only use those octave frequencies, which contain the original ones
                                if Octave_freq not in self.Fr_octave_vector:
                                    self.Fr_octave_vector.append(Octave_freq)

                    inx += 1

                    # store the correspondence between the original and the octave frequencies
                    self.Fr_corresp_array = np.zeros((len(Fr_in_octave), 2), float)
                    for i in range(self.Fr_corresp_array.shape[0]):
                        self.Fr_corresp_array[i,0] = Octave_freq
                        self.Fr_corresp_array[i,1] = Fr_in_octave[i]

                    octave_scale = octave_scale * 2

                    f = open(os.environ['PATIENTDIR']+'/Stim_Signal/Fr_corresp_array'+str(self.trunc_param*1.0)+'.csv','ab')
                    np.savetxt(f, self.Fr_corresp_array, delimiter=" ")
                    f.close()

                break

        # load the complete correspondence array and assign it to  self.Fr_corresp_array
        Fr_corresp_arr_full = np.genfromtxt(os.environ['PATIENTDIR']+'/Stim_Signal/Fr_corresp_array'+str(self.trunc_param * 1.0)+'.csv', delimiter=' ')
        self.Fr_corresp_array = np.round(Fr_corresp_arr_full, 6)

        self.FR_vector_signal_new = self.FR_vector_signal_new[self.FR_vector_signal_new != -1.0] # remove not filled entries
        logging.critical("Number of frequencies after truncation with the octave method: {}".format(self.FR_vector_signal_new.shape[0]))

        self.Fr_octave_vector = np.asarray(self.Fr_octave_vector)
        np.savetxt(os.environ['PATIENTDIR']+'/Stim_Signal/Fr_octave_vector_'+str(self.trunc_param*1.0)+'.csv', self.Fr_octave_vector, delimiter=" ")

        # find the index of the first frequency from the octave bands
        inx_start_octv_rslt = np.where(np.round(self.FR_vector_signal_new, 6) == np.round(self.Fr_octave_vector[0],6))
        self.inx_start_octv = inx_start_octv_rslt[0][0]

        # plot the octave bands
        Xs_unit = 3 * (np.absolute(self.Xs_signal_norm)).max() * np.ones(self.FR_vector_signal_new.shape[0], float)
        plt.figure(161)
        plt.stem(self.FR_vector_signal, np.absolute(self.Xs_signal_norm), markerfmt=" ", linefmt='C0', basefmt="C0-", use_line_collection=True)
        plt.stem(self.FR_vector_signal_new, Xs_unit, markerfmt=" ", linefmt="C1--",
                 basefmt="C0-", use_line_collection=True)  # we need to scale appropriately. Think about this image
        plt.xscale("log")
        plt.xlabel('Frequency, Hz')
        plt.xlim(10e0, 10e6)
        # plt.ylabel('Magnitude')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.savefig(os.environ['PATIENTDIR'] + '/Images/FT_octave' + str(self.trunc_param * 1.0) + '.eps', format='eps',
                    dpi=1000)
        np.savetxt(
            os.environ['PATIENTDIR'] + '/Stim_Signal/FR_vector_signal_octaves' + str(self.trunc_param * 1.0) + '.csv',
            self.FR_vector_signal_new, delimiter=" ")
        logging.critical('New frequency vector can be found in Stim_Signal/FR_vector_signal_octaves{}\n.csv'.format(
            str(self.trunc_param * 1.0)))


    def get_freqs_for_calc(self):

        if self.trunc_method == 'Octave Band Method':
            self.get_octave_band_freqs()
        elif self.trunc_method == 'High Amplitude Method' or self.trunc_method == 'Cutoff Method':
            self.get_specific_freq()
            # this block could be outdated
            Xs_storage_trunc = np.zeros((np.real(self.Xs_signal_norm_new).shape[0], 2), float)
            Xs_storage_trunc[:,0] = np.real(self.Xs_signal_norm_new)
            Xs_storage_trunc[:,1] = np.imag(self.Xs_signal_norm_new)
            np.savetxt(os.environ['PATIENTDIR']+'/Stim_Signal/Xs_storage_high_ampl'+str(self.trunc_param)+'.csv', Xs_storage_trunc, delimiter=" ")
        else:
            logging.critical('The truncation method {} was not implemented!'.format(self.trunc_method))
            raise SystemExit

        return True
