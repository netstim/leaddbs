# -*- coding: utf-8 -*-
"""
Created on Sat Sep 29 22:50:08 2018

@author: trieu, butenko
"""
import matplotlib.pyplot as plt
import numpy as np
import time as time_lib

from multiprocessing import Pool #  Process pool
from multiprocessing import sharedctypes
from functools import partial

#generate_signal is the manager function (called in Launcher)

def pick_refinement_freqs(Fr_vect,Xs_vect_unscaled,num_freqs,A):
    
    Xs_vect=Xs_vect_unscaled[:]/A    #scale back to unit if necessary
    
    Xs_vect_magn=np.absolute(Xs_vect)
    
    ref_freqs=np.zeros((num_freqs),float)
    ref_Xs=(1+1j)*np.zeros((num_freqs),float)
    
    freq_i=1
    last_index=0
    loc_index=0
    
    
    while (freq_i-1) < num_freqs:
        Xs_sum_fraction=freq_i*np.sum(Xs_vect_magn)/(num_freqs+1)
        el_sum=np.sum(np.absolute(Xs_vect[:last_index]))
        while el_sum<Xs_sum_fraction:
            el_sum += abs(Xs_vect[loc_index])
            loc_index=loc_index+1
            
        ref_freqs[freq_i-1]=Fr_vect[loc_index-1]
        ref_Xs[freq_i-1]=Xs_vect[loc_index-1]

        freq_i=freq_i+1
        last_index=loc_index
        
    print("Refinement frequencies: ",ref_freqs,"\n")
    #print(ref_Xs)
    
    return ref_freqs,ref_Xs

def my_generate_signal(T,A,signal_freq,phi,freq_max,t_step,simulation_time):        #only positive ampl?           
    #import numpy as np
    FR_vector_signal=np.arange(signal_freq,freq_max,signal_freq)
    #print FR_vector_signal.shape[0]
    #t_step=10e-7
    t_vector=np.arange(0,simulation_time,t_step)
            
    #T=6.0e-5
    #A=0.0002
    #phi=T/2
    
    Xs_signal=[]
    for i in FR_vector_signal:
        res=A*T*np.sinc(T*i)*np.exp(-2j*np.pi*i*phi)
        Xs_signal.append(res)
    
    Xs_signal=np.array(Xs_signal)
    
    #        Out_arr=zeros((FR_vector_signal.shape[0],1),float)
    Signal_t=(1/(2*np.pi))*np.sum([Xs_signal[int((i-signal_freq)/signal_freq)]*np.exp(2j*np.pi*i*t_vector) for i in FR_vector_signal],axis=0)
    
    #plt.figure(0)
    #plt.plot(t_vector,Signal_t)
    #plt.xlim(0.0,0.0002)
    
    koefff=max(Signal_t.real)
    
    Xs_signal_normalized=(A/koefff)*np.array(Xs_signal)
    #print Xs_signal_normalized
    #print type(Xs_signal_normalized)
    Signal_t_normalized=(1/(2*np.pi))*np.sum([Xs_signal_normalized[int((i-signal_freq)/signal_freq)]*np.exp(2j*np.pi*i*t_vector) for i in FR_vector_signal],axis=0)        
    plt.figure(1123123123)
    plt.plot(t_vector,Signal_t_normalized)
    
    plt.xlim(0.000,T*5)
    plt.grid(True)
    plt.xlabel('t, sec')
    plt.ylabel('Signal amplitude (A or V)')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.savefig('Images/Signal_regenerated.png', format='png', dpi=750)

    return t_vector,Signal_t_normalized,Xs_signal_normalized,FR_vector_signal

def numpy_analog_digit_converter(t_vect,signal_t_vect,freq,freq_spect_size,T):

    # as in the numpy example for np.fft.fft
    cutoff = int(np.ceil((t_vect.shape[0]+1)/2.))  
    
    sp = np.fft.fft(signal_t_vect)
    sp_cut=sp[:cutoff]
        
    fr_vector=freq*np.arange(0, cutoff, 1)

    '''if we want to rebuild the signal'''   
    if np.mod(t_vect.shape[0], 2):  # if the FT vector is odd
        fv_conj = np.conjugate(sp_cut[-1:0:-1])
    else:  # if the FT vector is even
        fv_conj = np.conjugate(sp_cut[-2:0:-1])
    
    Y = np.concatenate((sp_cut, fv_conj), axis=0)
    
    signal_t_rebuilt=np.fft.ifft(Y).real
    #plt.show()
    plt.figure(11111232111)
    plt.plot(t_vect,signal_t_rebuilt)
    plt.xlim(0.000,T*5)
    plt.grid(True)
    plt.xlabel('t, sec')
    plt.ylabel('Signal amplitude (A or V)')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.savefig('Images/Signal_recovered_shape.png', format='png', dpi=750)

    return fr_vector,sp_cut


def manual_signal_out_generator(d,t,A):
    signal_out=np.zeros(len(t),float)
    num_steps_signal=int(d["T"]/d["t_step"])
    num_steps_shift=int(d["phi"]/d["t_step"])

    for signal_step in np.arange(num_steps_shift,num_steps_signal+num_steps_shift):
        if (d["Signal_type"] == 'Rectangle'):# Increase Ram
            signal_out[signal_step]=A
        elif (d["Signal_type"] == 'Increasing Ramp'):# Increase Ram
            signal_out[signal_step]=(A/num_steps_signal)*(signal_step-num_steps_shift)           
        elif (d["Signal_type"] == 'Decreasing Ramp'):# Increase Ram
            signal_out[signal_step]=-1*(A/num_steps_signal)*(signal_step-num_steps_shift-num_steps_signal)
        elif(d["Signal_type"] == 'Central Triangle'):# Central Triangular
            if signal_step>int(num_steps_signal/2)+num_steps_shift:
                signal_out[signal_step]=-2*(A/num_steps_signal)*(signal_step-num_steps_shift-num_steps_signal)
            else:
                signal_out[signal_step]=2*(A/num_steps_signal)*(signal_step-num_steps_shift)
    
    return signal_out

def get_vector_in_time(Hf_zero,Hf_signal,w0,Nmax,phi,t_step,n_time_max,t_ind):
    
    tmp = np.ctypeslib.as_array(shared_array)
    t = [t_step*x for x in range(n_time_max)]
    
    Ht = Hf_zero + np.sum (Hf_signal[k]*np.exp(-1j*w0*(k+1)*(t[t_ind]-phi)) for k in range (0,Nmax-2))  #Nmax-2
    tmp[t_ind]=np.real(Ht)
      


def generate_signal(d,A,amp_max,cc_multi):
    
    start_signal_generation=time_lib.clock() 
    
    print(d["Signal_type"]," with repetition rate ",d["freq"]," Hz and ",np.round(d["T"]*1000,8)," ms pulse width")            
    Sim_time=1.0/d["freq"]      # always one pulse per simulation
    freq_max=d["freq"]*Sim_time/d["t_step"]
    print("Max frequency in the spectrum: ", freq_max/2.0,"\n")
    
    FR_vector_signal=np.arange(0.0,freq_max,d["freq"])
    n_time_max = int (Sim_time/d["t_step"])
    t = [d["t_step"]*x for x in range(n_time_max)]

    # to generate the signal with its analytical formulation (by Trieu)  
    II = np.pi
    phi = d["phi"] # signal shift in sec
    w0 = 2*II*d["freq"]
    pw = d["T"] #pulse width
    Nmax=FR_vector_signal.shape[0]    
    
    #signal_out = []
    #Hf_signal = []
    
    signal_out=np.ctypeslib.as_ctypes(np.zeros(len(t),float))
    global shared_array
    shared_array = sharedctypes.RawArray(signal_out._type_, signal_out)
    
    Hf_signal=np.complex(0,0)*np.zeros(Nmax-1,float)
    
    for k in range (1,Nmax):
        if (d["Signal_type"] == 'Increasing Ramp'):# Ascending Ramp
            Hf_zero = A*pw/2/Sim_time # Hf1 at k=0
            #Hf1 = 2*A/(Sim_time*pw)*(pw*np.exp(1j*w0*k*pw)/(1j*w0*k) + (np.exp(1j*w0*k*pw)-1)/(w0*k)**2)
            #Hf_signal.append(Hf1)
            Hf_signal[k-1]=2*A/(Sim_time*pw)*(pw*np.exp(1j*w0*k*pw)/(1j*w0*k) + (np.exp(1j*w0*k*pw)-1)/(w0*k)**2)
        elif(d["Signal_type"] == 'Decreasing Ramp'): # Descending Ramp
            Hf_zero = A*pw/2/Sim_time # Hf2 at k=0
            #Hf2 = -2*A/(Sim_time*pw)*(pw/(1j*w0*k) + (np.exp(1j*w0*k*pw)-1)/(w0*k)**2)
            #Hf_signal.append(Hf2)
            Hf_signal[k-1]=-2*A/(Sim_time*pw)*(pw/(1j*w0*k) + (np.exp(1j*w0*k*pw)-1)/(w0*k)**2)
        elif(d["Signal_type"] == 'Central Triangle'):# Central Triangular
            Hf_zero = A*pw/2/Sim_time # Hf3 at k=0
            #Hf3 = 4*A/(Sim_time*pw)*(( 2*np.exp(1j*w0*k*pw/2) - np.exp(1j*w0*k*pw)-1 )/(w0*k)**2)
            #Hf_signal.append(Hf3)
            Hf_signal[k-1]=4*A/(Sim_time*pw)*(( 2*np.exp(1j*w0*k*pw/2) - np.exp(1j*w0*k*pw)-1 )/(w0*k)**2)
        elif(d["Signal_type"] == 'Rectangle'):# Rectangular
            Hf_zero = A*pw/Sim_time # Hf4 at k=0
            #Hf4 = 2*A/(Sim_time*1j*w0*k)*( np.exp(1j*w0*k*pw) -1 )
            #Hf_signal.append(Hf4)
            Hf_signal[k-1]=2*A/(Sim_time*1j*w0*k)*( np.exp(1j*w0*k*pw) -1 )


    p = Pool()
    time_ind=np.arange(len(t))
    res = p.map(partial(get_vector_in_time, Hf_zero,Hf_signal,w0,Nmax,phi,d["t_step"],n_time_max),time_ind)
    signal_out = np.ctypeslib.as_array(shared_array)
    p.terminate()
    
    #signal_out=np.asarray(signal_out)
    del Hf_signal
    signal_out_real=signal_out.real
    
    ## to construct the signal manually (quick approach but the signal is almost "untruncatable")
    #signal_out_real=manual_signal_out_generator(d,t,A) 

    plt.figure(11111231)
    if d["current_control"]==1 and cc_multi==False:
        plt.plot(t,signal_out_real)
    else:
        signal_out_scaled=[i * amp_max for i in signal_out_real]
        plt.plot(t,signal_out_scaled)
    plt.xlim(0.000,d["T"]*5)
    plt.grid(True)
    plt.xlabel('t, sec')
    plt.ylabel('Signal amplitude (A or V)')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.savefig('Images/Signal.png', format='png', dpi=750)
    t=np.asarray(t)

    # get a Fourier transformation of the signal with np.fft.fft and recover to check with np.fft.ifft
    Fr_vect,Xs_vect=numpy_analog_digit_converter(t,signal_out_real,d["freq"],FR_vector_signal.shape[0],d["T"]) 

    #==========Plots==========================================================#
           
    # these take time to generate 
    
    # plt.figure(11)
    # plt.stem(Fr_vect, np.real(Xs_vect), markerfmt=" ")
    # plt.xscale("log")    
    # plt.xlabel('Frequency, Hz')
    # plt.ylabel('Real part')
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    # plt.savefig('Images/FT_real.png', format='png', dpi=1000)
    
    # plt.figure(12)
    # plt.stem(Fr_vect, np.imag(Xs_vect), markerfmt=" ")
    # plt.xscale("log")
    # plt.xlabel('Frequency, Hz')
    # plt.ylabel('Imaginary part')
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    # plt.savefig('Images/FT_imag.png', format='png', dpi=1000)

    # '''Maybe we want to use less freq'''
    # plt.figure(111342)
    # plt.stem(Fr_vect, np.absolute(Xs_vect), markerfmt=" ",linefmt='C0',basefmt="C0-")
    # plt.xscale("log")
    # plt.xlabel('Frequency, Hz')
    # plt.ylabel('Amplitude')
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    # plt.savefig('Images/FT_full_ampl.eps', format='eps', dpi=1000)   

    minutes=int((time_lib.clock() - start_signal_generation)/60)
    secnds=int(time_lib.clock() - start_signal_generation)-minutes*60
    print("----- Signal generation took ",minutes," min ",secnds," s -----")     

    return t,signal_out_real,Xs_vect,Fr_vect
