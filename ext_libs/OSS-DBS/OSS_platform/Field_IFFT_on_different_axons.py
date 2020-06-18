#convolute_signal_with_field_and_compute_ifft is the manager function (called in Launcher)
#compute_Z_ifft is called directly in Launcher

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import os
import time as time_lib
from pandas import read_csv
import h5py

solution_sort='none'
Xs_signal_normalized='none'
t_vect='none'
FREQ_vector_signal='none'

def compute_Z_ifft(d,Xs_signal_normalized,t_vector,A):      # get impedance in Time (only for 2-contact systems)
          
    Impedance_fr_get=read_csv('Field_solutions/sorted_impedance.csv', delimiter=' ', header=None)
    Impedance_fr=Impedance_fr_get.values
    Z_Tr=np.vectorize(complex)(Impedance_fr[:,0],Impedance_fr[:,1])
    
    if d["spectrum_trunc_method"]=='No Truncation':
        Z_conv=Xs_signal_normalized*Z_Tr/A       #devided by A because we need a normilized Xs_signal_normalized, not scaled with the signal ampl.
        
    if d["spectrum_trunc_method"]=='High Amplitude Method':
        Ind_trunc1 = np.genfromtxt('Stim_Signal/Indices_high_ampl.csv', delimiter=' ')
        cutoff = int(np.ceil((t_vector.shape[0]+1)/2.))             
        Xs_signal_full=np.complex(1.0,0.0)*np.zeros(cutoff,float)            
        Z_Tr_full=np.complex(1.0,0.0)*np.ones(cutoff,float)
        
        if d["Truncate_the_obtained_full_solution"]==0:            
            np.put(Z_Tr_full,Ind_trunc1.astype(int),Z_Tr)            
            np.put(Xs_signal_full,Ind_trunc1.astype(int),Xs_signal_normalized)
        
        if d["Truncate_the_obtained_full_solution"]==1:
            Xs_signal_normalized_trunc=Xs_signal_normalized[Ind_trunc1.astype(int)]
            Z_Tr_trunc=Z_Tr[Ind_trunc1.astype(int)]
            
            np.put(Z_Tr_full,Ind_trunc1.astype(int),Z_Tr_trunc)            
            np.put(Xs_signal_full,Ind_trunc1.astype(int),Xs_signal_normalized_trunc)
                
        Z_conv=Xs_signal_full*Z_Tr_full/A
        
    if d["spectrum_trunc_method"]=='Cutoff Method':
        Ind_trunc1 = np.genfromtxt('Stim_Signal/Indices_cutoff.csv', delimiter=' ')
        cutoff = int(np.ceil((t_vector.shape[0]+1)/2.))             
        Xs_signal_full=np.complex(1.0,0.0)*np.zeros(cutoff,float)            
        Z_Tr_full=np.complex(1.0,0.0)*np.ones(cutoff,float)
        
        if d["Truncate_the_obtained_full_solution"]==0:
            np.put(Z_Tr_full,np.arange(int(Ind_trunc1[1])),Z_Tr)            
            np.put(Xs_signal_full,np.arange(int(Ind_trunc1[1])),Xs_signal_normalized)
            
        if d["Truncate_the_obtained_full_solution"]==1:
            Xs_signal_normalized_trunc=Xs_signal_normalized[np.arange(int(Ind_trunc1[1]))]
            Z_Tr_trunc=Z_Tr[np.arange(int(Ind_trunc1[1]))]
            
            np.put(Z_Tr_full,np.arange(int(Ind_trunc1[1])),Z_Tr_trunc)            
            np.put(Xs_signal_full,np.arange(int(Ind_trunc1[1])),Xs_signal_normalized_trunc)

        Z_conv=Xs_signal_full*Z_Tr_full/A

    if np.mod(t_vector.shape[0], 2):  # if the FT vector is odd
        fv_conj = np.conjugate(Z_conv[-1:0:-1])
    else:  # if the FT vector is even
        fv_conj = np.conjugate(Z_conv[-2:0:-1])
      
    Y = np.concatenate((Z_conv, fv_conj), axis=0)
    Signal_t_Zconv=np.fft.ifft(Y).real
    
    plt.figure(111122112)
    plt.plot(t_vector,Signal_t_Zconv.real)
    plt.xlim(0.000,d["T"]*5)
    plt.grid(True)
    plt.xlabel('t, sec')
    plt.ylabel('Zreal, Ohm')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.savefig('Images/Z_convoluted.png', format='png', dpi=1000)
    np.savetxt('Field_solutions/Z_R_TimeDomain.csv', Signal_t_Zconv.real, delimiter=" ")
    np.savetxt('Field_solutions/Z_Im_TimeDomain.csv', Signal_t_Zconv.imag, delimiter=" ")

    return Signal_t_Zconv

#get electric potential in time over the axon 
def convolute_and_ifft(last_point,Ind_trunc1,trunc_method,post_truncation,i_axon,num_segments,N_freq,phi_shift,T,output):

    global solution_sort
    global Xs_signal_normalized
    global t_vect
    global FREQ_vector_signal
    
    for i_point in range(int(num_segments)):
        global_i_point=int(num_segments*i_axon+i_point)
        Xs_Tr=np.vectorize(complex)(solution_sort[(global_i_point)*N_freq:(global_i_point*N_freq+N_freq),3],solution_sort[global_i_point*N_freq:(global_i_point*N_freq+N_freq),4])         #real and im parts for the first point in VTA
            
        if trunc_method=='No Truncation':
            Xs_conv=Xs_signal_normalized*Xs_Tr
            
        if trunc_method=='Cutoff Method':
            cutoff = int(np.ceil((t_vect.shape[0]+1)/2.))             
            Xs_signal_full=np.complex(1.0,0.0)*np.zeros(cutoff,float)
            Xs_Tr_full=np.complex(1.0,0.0)*np.ones(cutoff,float)

            if post_truncation==0:
                np.put(Xs_Tr_full,np.arange(int(Ind_trunc1[1])),Xs_Tr)            
                np.put(Xs_signal_full,np.arange(int(Ind_trunc1[1])),Xs_signal_normalized)
                
            if post_truncation==1:
                Xs_signal_normalized_trunc=Xs_signal_normalized[np.arange(int(Ind_trunc1[1]))]
                Xs_Tr_trunc=Xs_Tr[np.arange(int(Ind_trunc1[1]))]
                
                np.put(Xs_Tr_full,np.arange(int(Ind_trunc1[1])),Xs_Tr_trunc)            
                np.put(Xs_signal_full,np.arange(int(Ind_trunc1[1])),Xs_signal_normalized_trunc)

            Xs_conv=Xs_signal_full*Xs_Tr_full
            
        if trunc_method=='High Amplitude Method':            

            cutoff = int(np.ceil((t_vect.shape[0]+1)/2.))             
            Xs_signal_full=np.complex(1.0,0.0)*np.zeros(cutoff,float)            
            Xs_Tr_full=np.complex(1.0,0.0)*np.ones(cutoff,float)
            
            if post_truncation==0:            
                np.put(Xs_Tr_full,Ind_trunc1.astype(int),Xs_Tr)            
                np.put(Xs_signal_full,Ind_trunc1.astype(int),Xs_signal_normalized)
            
            if post_truncation==1:
                Xs_signal_normalized_trunc=Xs_signal_normalized[Ind_trunc1.astype(int)]
                Xs_Tr_trunc=Xs_Tr[Ind_trunc1.astype(int)]
                
                np.put(Xs_Tr_full,Ind_trunc1.astype(int),Xs_Tr_trunc)            
                np.put(Xs_signal_full,Ind_trunc1.astype(int),Xs_signal_normalized_trunc)
                
        
            Xs_conv=Xs_signal_full*Xs_Tr_full
                      
        if np.mod(t_vect.shape[0], 2):  # if the FT vector is odd
            fv_conj = np.conjugate(Xs_conv[-1:0:-1])
        else:  # if the FT vector is even
            fv_conj = np.conjugate(Xs_conv[-2:0:-1])
          
        Y = np.concatenate((Xs_conv, fv_conj), axis=0)
        Signal_t_conv=np.fft.ifft(Y).real
        
        if global_i_point+last_point==0:
            plt.figure(11111234)
            plt.plot(t_vect,Signal_t_conv)
            plt.xlim(0.000,T*5)
            plt.grid(True)
            plt.xlabel('t, sec')
            plt.ylabel('Potential, V')
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            plt.savefig('Images/Signal_convoluted_1st_point.png', format='png', dpi=1000)
    
        np.save('Points_in_time/Signal_t_conv'+str(global_i_point+last_point), Signal_t_conv.real)
    
    output.put(i_axon)

#get electric potential in time over the axon when using octave band method for the frequency spectrum approximation
def convolute_and_ifft_octaves(last_point,i_axon,num_segments,N_freq,N_freq_octv,phi_shift,Fr_corresp_ar,Fr_octave_vec,i_start_octv,T,output):

    global solution_sort_octv
    global Xs_signal_normalized
    global t_vect
    global FREQ_vector_signal
    global FR_vec_sign_octv
    
    for i_point in range(int(num_segments)):
        global_i_point=int(num_segments*i_axon+i_point)
        Xs_Tr=np.vectorize(complex)(solution_sort_octv[(global_i_point)*N_freq_octv:(global_i_point*N_freq_octv+N_freq_octv),3],solution_sort_octv[global_i_point*N_freq_octv:(global_i_point*N_freq_octv+N_freq_octv),4])         #real and im parts for the first point in VTA
        Xs_Tr_full_real=np.zeros(FREQ_vector_signal.shape[0],float)
        Xs_Tr_full_imag=np.zeros(FREQ_vector_signal.shape[0],float)
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
        if global_i_point==0:

            np.savetxt('Field_solutions/Xs_Tr_full_real.csv', Xs_Tr_full_real, delimiter=" ")
            np.savetxt('Field_solutions/Xs_Tr_full_imag.csv', Xs_Tr_full_imag, delimiter=" ")
            
        Xs_Tr_full_complex=np.vectorize(complex)(Xs_Tr_full_real,Xs_Tr_full_imag)
    
        Xs_conv=Xs_signal_normalized*Xs_Tr_full_complex

        if np.mod(t_vect.shape[0], 2):  # if the FT vector is odd
            fv_conj = np.conjugate(Xs_conv[-1:0:-1])
        else:  # if the FT vector is even
            fv_conj = np.conjugate(Xs_conv[-2:0:-1])
        
        Y = np.concatenate((Xs_conv, fv_conj), axis=0)
        
        Signal_t_conv=np.fft.ifft(Y).real
        
        if global_i_point+last_point==0:
            plt.figure(11111234)
            plt.plot(t_vect,Signal_t_conv)
            plt.xlim(0.000,T*5)
            plt.grid(True)
            plt.xlabel('t, sec')
            plt.ylabel('Potential, V')
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            plt.savefig('Images/Signal_convoluted_1st_point.png', format='png', dpi=500)
       
        np.save('Points_in_time/Signal_t_conv'+str(global_i_point+last_point), Signal_t_conv.real)
    
    output.put(i_axon)

def convolute_signal_with_field_and_compute_ifft(d,XS_signal,models_in_population,number_of_segments,FR_vec_sign,t_vec,A,name_sol,inx_st_oct=0,dif_axons=False,last_point=0):
        
    # here I we don't want to pass large arrays so we use global variables    
    start_ifft=time_lib.time()
        
    if d["spectrum_trunc_method"]=='Octave Band Method':
        Fr_corresp_arr = np.genfromtxt('Stim_Signal/Fr_corresp_array'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')
        Fr_octave_vect = np.genfromtxt('Stim_Signal/Fr_octave_vector_'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')
        FR_vec_oct = np.genfromtxt('Stim_Signal/FR_vector_signal_octaves'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')
        
        Fr_corresp_arr=np.round(Fr_corresp_arr,6)
        Fr_octave_vect=np.round(Fr_octave_vect,6)
        
        global FR_vec_sign_octv
        FR_vec_sign_octv=FR_vec_oct
        global solution_sort_octv    
        
        if dif_axons==True:
            hf = h5py.File(name_sol[:-4]+'.h5', 'r')
            solution_sort_octv_full = hf.get('dataset_1')
            solution_sort_octv_full = np.array(solution_sort_octv_full)
            hf.close()                
            solution_sort_octv=solution_sort_octv_full[last_point*FR_vec_sign_octv.shape[0]:(last_point+models_in_population*number_of_segments)*FR_vec_sign_octv.shape[0],:]
        else:
            hf = h5py.File(name_sol[:-4]+'.h5', 'r')
            solution_sort_octv = hf.get('dataset_1')
            solution_sort_octv = np.array(solution_sort_octv)
            hf.close()   
    else:        
        global solution_sort

        if dif_axons==True:
            hf = h5py.File(name_sol[:-4]+'.h5', 'r')
            solution_sort_full = hf.get('dataset_1')
            solution_sort_full = np.array(solution_sort_full)
            hf.close()
            solution_sort=solution_sort_full[last_point*FR_vec_sign.shape[0]:(last_point+models_in_population*number_of_segments)*FR_vec_sign.shape[0],:]
        else:
            hf = h5py.File(name_sol[:-4]+'.h5', 'r')
            solution_sort = hf.get('dataset_1')
            solution_sort = np.array(solution_sort)
            hf.close()
            
        if d["spectrum_trunc_method"]=='Cutoff Method':
            Ind_trunc = np.genfromtxt('Stim_Signal/Indices_cutoff.csv', delimiter=' ')
            
        elif d["spectrum_trunc_method"]=='High Amplitude Method':
            Ind_trunc = np.genfromtxt('Stim_Signal/Indices_high_ampl.csv', delimiter=' ')
        else:
            Ind_trunc=0

    global FREQ_vector_signal   
    global t_vect
    global Xs_signal_normalized
    Xs_signal_normalized=XS_signal
    FREQ_vector_signal=FR_vec_sign
    t_vect=t_vec    
    N_freq=(FREQ_vector_signal.shape[0])

    j=0 #index of the process
    point_step=1
    i=0 #index of the population
    
    #just to check the progress
    axons_quart=[int(models_in_population/4.0),int(2*models_in_population/4.0),int(3*models_in_population/4.0)]
    
    while i<(models_in_population):
        output = mp.Queue()         #defines an output queue    
        proc=[] 
        while j<d["number_of_processors"] and i<(models_in_population):

            if d["spectrum_trunc_method"]=='Octave Band Method':
                N_freq_octv=(FR_vec_sign_octv.shape[0])
                processes=mp.Process(target=convolute_and_ifft_octaves,args=(last_point,i,number_of_segments,N_freq,N_freq_octv,d["phi"],Fr_corresp_arr,Fr_octave_vect,inx_st_oct,d["T"],output))
            else:
                processes=mp.Process(target=convolute_and_ifft,args=(last_point,Ind_trunc,d["spectrum_trunc_method"],d["Truncate_the_obtained_full_solution"],i,number_of_segments,N_freq,d["phi"],d["T"],output))
            proc.append(processes)
            j=j+1
            i=i+point_step
            if i in axons_quart:
                print(int(i*100/models_in_population)+1,"% of neuron models were processed")
        j=0
        for p in proc:       
            p.start()
        for p in proc:
            p.join()
        convoluted_lists=[output.get() for p in proc]#the files are already created
    
    if d["spectrum_trunc_method"]=='Octave Band Method':
        del solution_sort_octv
        if dif_axons==True:
            del solution_sort_octv_full
    else:
        del solution_sort
        if dif_axons==True:
            del solution_sort_full
    
    last_point=last_point+i*number_of_segments

    minutes=int((time_lib.time() - start_ifft)/60)
    secnds=int(time_lib.time() - start_ifft)-minutes*60
    print("----- Signal scaling and IFFT and took: ",minutes," min ",secnds," s -----\n")  
    
    return last_point




