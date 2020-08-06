#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 13:29:22 2020

@author: butenko
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 14:09:42 2020

@author: butenko
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt

from multiprocessing import Pool #  Process pool
from multiprocessing import sharedctypes
from functools import partial

import time as time_lib

# parallelized IFFT on VTA array


#d,FR_vector_signal,Xs_signal_norm,t_vector
def ifft_on_VTA_array(Xs_signal_normalized,num_freqs,N_freq_octv,FR_vec_sign_octv,Fr_corresp_ar,t_vect,T,i_start_octv,i_point):    #in mm, in MRI space

    tmp = np.ctypeslib.as_array(shared_array)

    Xs_Tr=np.vectorize(complex)(solution_sort_octv[(i_point)*N_freq_octv:(i_point*N_freq_octv+N_freq_octv),3],solution_sort_octv[i_point*N_freq_octv:(i_point*N_freq_octv+N_freq_octv),4])         #real and im parts for the first point in VTA
    Xs_Tr_full_real=np.zeros(num_freqs,float)
    Xs_Tr_full_imag=np.zeros(num_freqs,float)
    stepper=0
    for i_inx in range(Xs_Tr.shape[0]):
        if i_inx>=i_start_octv:
            rslt=np.where(Fr_corresp_ar[:,0]==np.round(FR_vec_sign_octv[i_inx],6))
            step_octv=rslt[0].shape[0]   #size of the freq. pack in the octave
            
            Xs_Tr_full_real[stepper:stepper+step_octv]=(Xs_Tr[i_inx].real) 
            Xs_Tr_full_imag[stepper:stepper+step_octv]=(Xs_Tr[i_inx].imag) 
            stepper=stepper+step_octv
        else:
            Xs_Tr_full_real[stepper]=(Xs_Tr[i_inx].real) 
            Xs_Tr_full_imag[stepper]=(Xs_Tr[i_inx].imag) 
            stepper=stepper+1
            
    if i_point==0:
        np.savetxt('/opt/Patient/Field_solutions/Xs_Tr_full_real.csv', Xs_Tr_full_real, delimiter=" ")
        np.savetxt('/opt/Patient/Field_solutions/Xs_Tr_full_imag.csv', Xs_Tr_full_imag, delimiter=" ")
        
    Xs_Tr_full_complex=np.vectorize(complex)(Xs_Tr_full_real,Xs_Tr_full_imag)

    Xs_conv=Xs_signal_normalized*Xs_Tr_full_complex

    if np.mod(t_vect.shape[0], 2):  # if the FT vector is odd
        fv_conj = np.conjugate(Xs_conv[-1:0:-1])
    else:  # if the FT vector is even
        fv_conj = np.conjugate(Xs_conv[-2:0:-1])
    
    Y = np.concatenate((Xs_conv, fv_conj), axis=0)
    
    Signal_t_conv=np.fft.ifft(Y).real
    
    tmp[i_point]=abs(max(Signal_t_conv[:], key=abs))
    
    if i_point==1:
         plt.figure(11111234)
         plt.plot(t_vect,Signal_t_conv)
         plt.xlim(0.000,T*5)
         plt.grid(True)
         plt.xlabel('t, sec')
         plt.ylabel('Potential, V')
         plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
         plt.savefig('/opt/Patient/Images/Signal_convoluted_1st_point.png', format='png', dpi=500)
   
    #np.save('Points_in_time/Signal_t_conv'+str(i_point), Signal_t_conv.real)
        
    #return(Max_signal_for_point)

def get_IFFT_on_VTA_array(num_of_proc,name_sol,d,FREQ_vector_signal,Xs_signal_normalized,t_vect,T,i_start_octv,arrays_shape):

    start_IFFT=time_lib.time()    

    global solution_sort_octv
    
    hf = h5py.File(name_sol[:-4]+'.h5', 'r')
    solution_sort_octv = hf.get('dataset_1')
    solution_sort_octv = np.array(solution_sort_octv)
    hf.close()    
        
    num_segments=sum(arrays_shape)

    Max_field_on_VTA_array=np.ctypeslib.as_ctypes(np.zeros(num_segments,float))
    global shared_array
    shared_array = sharedctypes.RawArray(Max_field_on_VTA_array._type_, Max_field_on_VTA_array)
    
    
    Fr_corresp_ar = np.genfromtxt('/opt/Patient/Stim_Signal/Fr_corresp_array'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')
    FR_vec_sign_octv = np.genfromtxt('/opt/Patient/Stim_Signal/FR_vector_signal_octaves'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')
        
    Fr_corresp_ar=np.round(Fr_corresp_ar,6)
    N_freq_octv=(FR_vec_sign_octv.shape[0])
    
    p = Pool(num_of_proc)
    res = p.map(partial(ifft_on_VTA_array, Xs_signal_normalized,FREQ_vector_signal.shape[0],N_freq_octv,FR_vec_sign_octv,Fr_corresp_ar,t_vect,T,i_start_octv),np.arange(num_segments))
    Max_field_on_VTA_array = np.ctypeslib.as_array(shared_array)
    p.terminate()
    
    minutes=int((time_lib.time() - start_IFFT)/60)
    secnds=int(time_lib.time() - start_IFFT)-minutes*60
    print("----- IFFT took ",minutes," min ",secnds," s -----")   
    
    return Max_field_on_VTA_array
    
