# -*- coding: utf-8 -*-
"""
Created on Sat Sep 29 22:50:08 2018

@author: Konstantin and Trieu

Stim_pulse class contains methods necessary to describe a DBS signal in both time and frequency domains
Charged balancing is supported
To include new signals, you need to add their Fourier series expansion to define_signals_in_freq_domain(...)
"""
import matplotlib.pyplot as plt
import numpy as np
import time as time_lib
import os
import pickle
import logging

from multiprocessing import Pool #  Process pool
from multiprocessing import sharedctypes
from functools import partial

def get_vector_in_time(Hf_zero, Hf_signal, w0, Num_freqs, signal_shift, t_step, num_tsteps, t_ind):

    """ parallelized transformation of the signal from the frequency to the time domain """

    tmp = np.ctypeslib.as_array(shared_array)
    t = [t_step*x for x in range(num_tsteps)]

    Ht = Hf_zero + np.sum(Hf_signal[k]*np.exp(-1j*w0*(k+1)*(t[t_ind]-signal_shift)) for k in range(0, Num_freqs - 2))
    tmp[t_ind] = np.real(Ht)

class Stim_pulse(object):
    def __init__(self, inp_dict, cc_multi):

        self.Sim_time = 1.0/inp_dict["freq"]                             #  always one pulse per simulation
        self.num_tsteps = int(self.Sim_time/inp_dict["t_step"])          #  number of time steps
        self.t_vector = [inp_dict["t_step"]*x for x in range(self.num_tsteps)]  #  time vector, in seconds

        self.freq_max = inp_dict["freq"] * self.Sim_time / inp_dict["t_step"]  # the last frequency in the Fourier Transform
        self.FR_vector_signal = np.arange(0.0, self.freq_max, inp_dict["freq"])
        self.num_freqs = self.FR_vector_signal.shape[0]
        # logging.critical("Max frequency in the spectrum: {}".format(freq_max/2.0))

        # to generate the signal using analytical formulation (by Trieu)
        self.signal_shift = inp_dict["phi"]      #  in sec
        self.w0 = 2 * np.pi * inp_dict["freq"]   #  angular frequency
        self.pw = inp_dict["T"]                  #  pulse width in sec
        self.rep_rate = inp_dict["freq"]         #  repetition rate of the DBS signal (usually 130-184 Hz)

        # now we don't need None values (floating potentials), Contacts contain only the active ones
        pulse_amp = [x for x in inp_dict["Pulse_amp"] if x is not None]

        self.cc_multi = cc_multi
        # For monopolar current-controlled stimulation, the DBS pulse is assigned the actual amplitude!
        # (The field will be initially computed for the unit current, and then scaled during cross-correlation and IFFT)
        if inp_dict["current_control"] == 1 and self.cc_multi == False:
            self.A = max(pulse_amp[:], key=abs)  # Here we should have only one non-zero value. Negative if cathodes
        else:
            self.A = 1.0  #  Use the unit signal for voltage-controlled and miltipolar current-controlled
                     #  Because the electric field will be scaled before the cross-correlation and IFFT
            self.Amp_max = max(pulse_amp[:], key=abs)    #  This variable is used for plotting


    def create_pulse(self, d):

        start_signal_generation = time_lib.time()
        logging.critical("{} with repetition rate {} Hz and {} ms pulse width".format(d["Signal_type"], d["freq"], np.round(d["T"]*1000, 8)))

        # analyt. definition in the frequency domain
        Hf_zero, Hf_signal = self.define_signals_in_freq_domain(d["Signal_type"])

        #  for parallel array processing in Pool
        signal_out = np.ctypeslib.as_ctypes(np.zeros(len(self.t_vector), float))
        global shared_array
        shared_array = sharedctypes.RawArray(signal_out._type_, signal_out)
        p = Pool()

        # convert signal to the time domain using parallel processing of time steps
        time_ind = np.arange(len(self.t_vector))
        res = p.map(partial(get_vector_in_time, Hf_zero, Hf_signal, self.w0, self.num_freqs, self.signal_shift, d["t_step"], self.num_tsteps), time_ind)
        signal_out = np.ctypeslib.as_array(shared_array)
        p.terminate()

        del Hf_signal

        if d['Charge_balancing'] == True:   #  in this case we add another signal
            true_parameters = [self.A, self.pw, self.signal_shift]

            self.signal_shift2 = self.signal_shift + self.pw  # signal shift (in sec) + shift by pulse width
            if d['Balancing_type'] == 'Symmetric':
                self.A2 = -1.0 * self.A
                self.pw2 = self.pw
                Signal_type_counter = d["Signal_type"]
            elif d['Balancing_type'] == 'Low_amplitude':
                self.pw2 = self.Sim_time - self.pw
                if d["Signal_type"] == 'Rectangle':
                    self.A2 = -1.0 * self.A * self.pw / self.pw2
                else:
                    self.A2 = -0.5 * self.A * self.pw / self.pw2
                Signal_type_counter = 'Rectangle'

            self.A, self.pw, self.signal_shift = (self.A2, self.pw2, self.signal_shift2)   #  just temporary

            # analyt. definition in the frequency domain
            Hf_zero2, Hf_signal2 = self.define_signals_in_freq_domain(Signal_type_counter)

            self.A, self.pw, self.signal_shift = true_parameters[:]  #  restore true parameters

            signal_out2 = np.ctypeslib.as_ctypes(np.zeros(len(self.t_vector), float))
            shared_array = sharedctypes.RawArray(signal_out2._type_, signal_out2)

            p = Pool()
            time_ind = np.arange(len(self.t_vector))
            t = d["t_step"] * time_ind
            res = p.map(partial(get_vector_in_time, Hf_zero2, Hf_signal2, self.w0, self.num_freqs, self.signal_shift2, d["t_step"], self.num_tsteps), time_ind)
            signal_out2 = np.ctypeslib.as_array(shared_array)
            p.terminate()

            del Hf_signal2

            signal_out_real_sum = signal_out2.real + signal_out.real
        else:
            signal_out_real_sum = signal_out.real


        plt.figure(11111231)
        if d["current_control"]==1 and self.cc_multi == False:
            plt.plot(self.t_vector, signal_out_real_sum)
        else:
            signal_out_scaled=[i * self.Amp_max for i in signal_out_real_sum]
            plt.plot(self.t_vector, signal_out_scaled)
        plt.xlim(0.000, self.pw * 5) # don't show the whole signal (1.0/freq), just 5 pulse widths
        plt.grid(True)
        plt.xlabel('t, sec')
        plt.ylabel('Signal amplitude (A or V)')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.savefig(os.environ['PATIENTDIR']+'/Images/Signal.png', format='png', dpi=750)
        self.t_vector = np.asarray(self.t_vector)

        # get a Fourier transformation of the signal with np.fft.fft and recover to check with np.fft.ifft
        self.numpy_analog_digit_converter(signal_out_real_sum)

        #=============================== Plots =========================================#

        # these take time to generate

        # plt.figure(11)
        # plt.stem(self.FR_vector_signal, np.real(self.Xs_signal_norm), markerfmt=" ")
        # plt.xscale("log")
        # plt.xlabel('Frequency, Hz')
        # plt.ylabel('Real part')
        # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        # plt.savefig(os.environ['PATIENTDIR']+'/Images/FT_real.png', format='png', dpi=1000)

        # plt.figure(12)
        # plt.stem(self.FR_vector_signal, np.imag(self.Xs_signal_norm, markerfmt=" ")
        # plt.xscale("log")
        # plt.xlabel('Frequency, Hz')
        # plt.ylabel('Imaginary part')
        # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        # plt.savefig(os.environ['PATIENTDIR']+'/Images/FT_imag.png', format='png', dpi=1000)

        # plt.figure(111342)
        # plt.stem(self.FR_vector_signal, np.absolute(self.Xs_signal_norm), markerfmt=" ",linefmt='C0',basefmt="C0-")
        # plt.xscale("log")
        # plt.xlabel('Frequency, Hz')
        # plt.ylabel('Amplitude')
        # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        # plt.savefig(os.environ['PATIENTDIR']+'/Images/FT_full_ampl.png', format='png', dpi=750)

        minutes=int((time_lib.time() - start_signal_generation)/60)
        secnds=int(time_lib.time() - start_signal_generation)-minutes*60
        logging.critical("----- Signal generation took {} min {} sec -----\n".format(minutes, secnds))

        return True


    def define_signals_in_freq_domain(self, Signal_type):
        """ signals definitions should be added in this function """

        Hf_signal = np.complex(0, 0) * np.zeros(self.num_freqs - 1, float)
        k = np.arange(1, self.num_freqs)

        if (Signal_type== 'Increasing Ramp'):  # Ascending Ramp
            Hf_zero = self.A * self.pw / 2 / self.Sim_time
            Hf_signal[k - 1] = 2 * self.A / (self.Sim_time * self.pw) * (
                    self.pw * np.exp(1j * self.w0 * k * self.pw) / (1j * self.w0 * k) + (np.exp(1j * self.w0 * k * self.pw) - 1) / (
                    self.w0 * k) ** 2)

        elif (Signal_type == 'Decreasing Ramp'):  # Descending Ramp
            Hf_zero = self.A * self.pw / 2 / self.Sim_time
            Hf_signal[k - 1] = -2 * self.A / (self.Sim_time * self.pw) * (
                    self.pw / (1j * self.w0 * k) + (np.exp(1j * self.w0 * k * self.pw) - 1) / (self.w0 * k) ** 2)

        elif (Signal_type == 'Central Triangle'):  # Central Triangular
            Hf_zero = self.A * self.pw / 2 / self.Sim_time
            Hf_signal[k - 1] = 4 * self.A / (self.Sim_time * self.pw) * (
                    (2 * np.exp(1j * self.w0 * k * self.pw / 2) - np.exp(1j * self.w0 * k * self.pw) - 1) / (self.w0 * k) ** 2)

        elif (Signal_type == 'Rectangle'):  # Rectangular
            Hf_zero = self.A * self.pw / self.Sim_time
            Hf_signal[k - 1] = 2 * self.A / (self.Sim_time * 1j * self.w0 * k) * (np.exp(1j * self.w0 * k * self.pw) - 1)

        else:
            logging.critical("\nThe signal type is not implemented")
            logging.critical("Check available signals in Signal_generator.py")
            raise SystemExit

        return Hf_zero, Hf_signal

    def numpy_analog_digit_converter(self, signal_t_vect):

        # as in the numpy example for np.fft.fft
        cutoff = int(np.ceil((self.t_vector.shape[0] + 1) / 2.))
        sp = np.fft.fft(signal_t_vect)
        self.Xs_signal_norm = sp[:cutoff]
        self.FR_vector_signal = self.rep_rate*np.arange(0, cutoff, 1)

        '''if we want to rebuild the signal (just shape, the amplitude may differ)'''
        if np.mod(self.t_vector.shape[0], 2):  # if the FT vector is odd
            fv_conj = np.conjugate(self.Xs_signal_norm[-1:0:-1])
        else:  # if the FT vector is even
            fv_conj = np.conjugate(self.Xs_signal_norm[-2:0:-1])

        Y = np.concatenate((self.Xs_signal_norm, fv_conj), axis=0)
        signal_t_rebuilt = np.fft.ifft(Y).real
        #plt.show()
        plt.figure(11111232111)
        plt.plot(self.t_vector, signal_t_rebuilt)
        plt.xlim(0.000, self.pw*5)
        plt.grid(True)
        plt.xlabel('t, sec')
        plt.ylabel('Signal amplitude (A or V)')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.savefig(os.environ['PATIENTDIR']+'/Images/Signal_Shape_recovered.png', format='png', dpi=750)

        return True



def pick_refinement_freqs(Fr_vect, Xs_vect, num_freqs):
    """ This function chooses frequencies (as many as num_freqs) for mesh refinement
        depending on the DBS power spectrum """

    Xs_vect_magn = np.absolute(Xs_vect)

    ref_freqs = np.zeros(num_freqs, float)
    ref_Xs = (1+1j)*np.zeros(num_freqs, float)

    freq_i = 1
    last_index = 0
    loc_index = 0

    while (freq_i - 1) < num_freqs:
        Xs_sum_fraction = freq_i * np.sum(Xs_vect_magn) / (num_freqs + 1)
        el_sum = np.sum(np.absolute(Xs_vect[:last_index]))
        while el_sum < Xs_sum_fraction:
            el_sum += abs(Xs_vect[loc_index])
            loc_index += 1

        ref_freqs[freq_i-1] = Fr_vect[loc_index-1]
        ref_Xs[freq_i-1] = Xs_vect[loc_index-1]

        freq_i += 1
        last_index = loc_index

    logging.critical("Refinement frequencies: {} \n".format(' '.join(map(str, ref_freqs))))

    return ref_freqs








