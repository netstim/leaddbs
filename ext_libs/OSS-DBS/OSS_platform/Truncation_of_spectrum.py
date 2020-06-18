# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 15:53:14 2019

@author: butenko
"""


import numpy as np
#from pandas import read_csv
import matplotlib.pyplot as plt
import os


#get_freqs_for_calc is the manager function (called in Launcher)

def get_specific_freq(mode,freq,freq_param,FR_vector_signal,Xs_signal_norm,t_vect,signal_width):
    print("----- Truncating the frequency spectrum -----")
    #solution_sort_get=read_csv('Field_solutions/sorted_solution.csv', delimiter=' ', header=None)
    #solution_sort=solution_sort_get.values

    if mode=='High Amplitude Method':     #choosing freq. where FT has the highest amplitudes, here freq_param is the number of frequencies in the pack

        ind_max=np.sort(np.argpartition(np.absolute(Xs_signal_norm),-1*freq_param)[-1*freq_param:])
        #print ind_max
                 
        FR_vector_signal_new=FR_vector_signal[ind_max]
        Xs_signal_norm_new=Xs_signal_norm[ind_max]
        
#        with open('Stim_Signal/Indices_high_ampl.csv','w') as f_handle:
#            np.savetxt(f_handle,ind_max) 
        np.savetxt('Stim_Signal/Indices_high_ampl.csv', ind_max, delimiter=" ")
        
        
        '''To create the truncated signal'''
        Xs_signal_full=np.complex(1.0,0.0)*np.zeros(Xs_signal_norm.shape[0],float)
        np.put(Xs_signal_full,ind_max.astype(int),Xs_signal_norm)
        
        Xs_ones=np.complex(1.0,0.0)*np.ones(Xs_signal_norm.shape[0],float) 
        Xs_signal_unit=Xs_signal_full*Xs_ones
        
        if np.mod(t_vect.shape[0], 2):  # if the FT vector is odd
            fv_conj = np.conjugate(Xs_signal_unit[-1:0:-1])
        else:  # if the FT vector is even
            fv_conj = np.conjugate(Xs_signal_unit[-2:0:-1])
          
        Y = np.concatenate((Xs_signal_unit, fv_conj), axis=0)
        Signal_t_conv=np.fft.ifft(Y).real
        
        plt.figure(111112341)
        plt.plot(t_vect,Signal_t_conv)
        plt.xlim(0.000,signal_width*5)
        plt.grid(True)
        plt.xlabel('t, sec')
        plt.ylabel('Signal amplitude (A or V)')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.savefig('Images/Signal_convoluted_truncated.png', format='png', dpi=1000)
        
        
        return FR_vector_signal_new,Xs_signal_norm_new
        

            
    if mode=='Cutoff Method':         #setting up a new cutoff freq. Here freq_param is the cutoff frequency
    
#        inx_cutoff=-1
#        for fr in FR_vector_signal:
#            inx_cutoff=inx_cutoff+1
#            if fr>freq_param:
#                inx_cutoff=inx_cutoff     
#                break
    
        inx_cutoff=freq_param
            
        FR_vector_signal_new=FR_vector_signal[:inx_cutoff]
        Xs_signal_norm_new=Xs_signal_norm[:inx_cutoff]
        
#        with open('Stim_Signal/Indices_cutoff.csv','w') as f_handle:
#            np.savetxt(f_handle,np.array([Xs_signal_norm.shape[0],inx_cutoff]))
        np.savetxt('Stim_Signal/Indices_cutoff.csv', np.array([Xs_signal_norm.shape[0],inx_cutoff]), delimiter=" ")
        
        
        Xs_signal_full=np.complex(1.0,0.0)*np.zeros(Xs_signal_norm.shape[0],float)
        np.put(Xs_signal_full,np.arange(int(inx_cutoff)),Xs_signal_norm)
        
        Xs_ones=np.complex(1.0,0.0)*np.ones(Xs_signal_norm.shape[0],float) 
        Xs_signal_unit=Xs_signal_full*Xs_ones
        
        if np.mod(t_vect.shape[0], 2):  # if the FT vector is odd
            fv_conj = np.conjugate(Xs_signal_unit[-1:0:-1])
        else:  # if the FT vector is even
            fv_conj = np.conjugate(Xs_signal_unit[-2:0:-1])
          
        Y = np.concatenate((Xs_signal_unit, fv_conj), axis=0)
        Signal_t_conv=np.fft.ifft(Y).real
        
        plt.figure(111112341)
        plt.plot(t_vect,Signal_t_conv)
        plt.xlim(0.000,signal_width*5)
        plt.grid(True)
        plt.xlabel('t, sec')
        plt.ylabel('Signal amplitude (A or V)')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.savefig('Images/Signal_convoluted_truncated.png', format='png', dpi=1000)
        
        return FR_vector_signal_new,Xs_signal_norm_new
        
    if mode=='Octave Band Method':         # the octave band will be formed after the frequency set up in freq_param
        FR_vector_signal_new=-1*np.ones((FR_vector_signal.shape[0]),float)
        fr_indices=[]
        Fr_octave_vector=[]

        
        inx_cutoff=-1
        for fr in FR_vector_signal:
            inx_cutoff=inx_cutoff+1
            if fr<=freq_param:
                fr_indices.append(inx_cutoff)
                
            if fr>freq_param:
                octave_scale=1
                FR_vector_signal_new[:inx_cutoff]=FR_vector_signal[:inx_cutoff]
                base_octaves=FR_vector_signal_new[inx_cutoff-1]
                
                while (base_octaves+(freq*octave_scale)/np.sqrt(2.0))<FR_vector_signal[-1]:
                    Octave_freq=base_octaves+(freq*octave_scale)/np.sqrt(2.0)     #we will take a base not from freq_param, but from the closest freq from the left in FR_vector_signal
                    Fr_octave_vector.append(Octave_freq)
                    FR_vector_signal_new[inx_cutoff]=Octave_freq
                    inx_cutoff=inx_cutoff+1
                    #print ("fr for octave: ", Octave_freq)
                    next_Octave=base_octaves+(freq*2*octave_scale)/np.sqrt(2.0)
                    #print ("next octave: ",next_Octave)
                    Fr_in_octave=[]
                    for fr_init in FR_vector_signal:
                        if fr_init>Octave_freq and fr_init<next_Octave:
                            Fr_in_octave.append(fr_init)
                            #print "The frequency added from the octave: ",fr_init
                        
                    Fr_corresp_array=np.zeros((len(Fr_in_octave),2),float)
                    for i in range(Fr_corresp_array.shape[0]):
                        Fr_corresp_array[i,0]=Octave_freq
                        Fr_corresp_array[i,1]=Fr_in_octave[i]
                        
                    octave_scale=octave_scale*2
                        
#                    with open('Stim_Signal/Fr_corresp_array.csv','a') as f_handle:
#                        np.savetxt(f_handle,Fr_corresp_array)
                    #np.savetxt('Stim_Signal/Fr_corresp_array.csv','a', np.array([Xs_signal_norm.shape[0],inx_cutoff]), delimiter=" ")
                    f=open('Stim_Signal/Fr_corresp_array'+str(freq_param*1.0)+'.csv','ab')
                    np.savetxt(f, Fr_corresp_array, delimiter=" ")
                    f.close()
                    
                break
                
        #FR_vector_signal_new=FR_vector_signal_new.flatten()
        FR_vector_signal_new=FR_vector_signal_new[FR_vector_signal_new!=-1.0]
        print("Number of frequencies after truncation with the octave method: ", FR_vector_signal_new.shape[0])
                  
        Fr_octave_vector=np.asarray(Fr_octave_vector)    
        np.savetxt('Stim_Signal/Fr_octave_vector_'+str(freq_param*1.0)+'.csv', Fr_octave_vector, delimiter=" ")
        
        #print "number of frequencies: ", FR_vector_signal_new.shape[0]
        #print("Frequency vector: ", FR_vector_signal_new
        
        return FR_vector_signal_new
    
    return False
    
def get_freqs_for_calc(d,FR_vector_signal,Xs_signal_norm,t_vector):
            
    if d["spectrum_trunc_method"]=='Octave Band Method' and not(os.path.isfile('Stim_Signal/FR_vector_signal_octaves'+str(d["trunc_param"]*1.0)+'.csv')):
        
        from Truncation_of_spectrum import get_specific_freq
        
        #for mode=3 we do not get XS vector. For others we do.
        FR_vector_signal_new = get_specific_freq(d["spectrum_trunc_method"],d["freq"],d["trunc_param"],FR_vector_signal,Xs_signal_norm,t_vector,d["T"])
        
        Fr_corresp_array = np.genfromtxt('Stim_Signal/Fr_corresp_array'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')
        Fr_octave_vector = np.genfromtxt('Stim_Signal/Fr_octave_vector_'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')
        
        Fr_corresp_array=np.round(Fr_corresp_array,6)
        Fr_octave_vector=np.round(Fr_octave_vector,6)
        
        inx_start_octv_rslt=np.where(np.round(FR_vector_signal_new,6)==Fr_octave_vector[0])
        inx_start_octv=inx_start_octv_rslt[0][0]
    
        Xs_unit=3*(np.absolute(Xs_signal_norm)).max()*np.ones(FR_vector_signal_new.shape[0],float)
    
        plt.figure(161)
        plt.stem(FR_vector_signal, np.absolute(Xs_signal_norm), markerfmt=" ",linefmt='C0',basefmt="C0-")
        plt.stem(FR_vector_signal_new, Xs_unit, markerfmt=" ",linefmt="C1--",basefmt="C0-")      #we need to scale appropriately. Think about this image
        plt.xscale("log")
        plt.xlabel('Frequency, Hz')
        #plt.ylabel('Magnitude')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.savefig('Images/FT_octave'+str(d["trunc_param"]*1.0)+'.eps', format='eps', dpi=1000)  
     
        np.savetxt('Stim_Signal/FR_vector_signal_octaves'+str(d["trunc_param"]*1.0)+'.csv', FR_vector_signal_new, delimiter=" ")
        
        print('New frequency vector can be found in Stim_Signal/FR_vector_signal_octaves'+str(d["trunc_param"]*1.0)+'.csv')
        print(" ")
        
        return FR_vector_signal_new,inx_start_octv
    
    if d["spectrum_trunc_method"]=='Octave Band Method' and (os.path.isfile('Stim_Signal/FR_vector_signal_octaves'+str(d["trunc_param"]*1.0)+'.csv')):
        FR_vector_signal_new = np.genfromtxt('Stim_Signal/FR_vector_signal_octaves'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')
        Fr_corresp_array = np.genfromtxt('Stim_Signal/Fr_corresp_array'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')
        Fr_octave_vector = np.genfromtxt('Stim_Signal/Fr_octave_vector_'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')
        
        Fr_corresp_array=np.round(Fr_corresp_array,6)
        Fr_octave_vector=np.round(Fr_octave_vector,6)
        
        inx_start_octv_rslt=np.where(np.round(FR_vector_signal_new,6)==Fr_octave_vector[0])
        inx_start_octv=inx_start_octv_rslt[0][0]
        
        return FR_vector_signal_new,inx_start_octv
        
    if d["spectrum_trunc_method"]=='High Amplitude Method' and not(os.path.isfile('Stim_Signal/FR_vector_signal_high_ampl'+str(d["trunc_param"])+'.csv')):
        
        from Truncation_of_spectrum import get_specific_freq        
        FR_vector_signal_new,Xs_signal_norm_new = get_specific_freq(d["spectrum_trunc_method"],d["freq"],d["trunc_param"],FR_vector_signal,Xs_signal_norm,t_vector,d["T"])
        #with open('Stim_Signal/FR_vector_signal_high_ampl'+str(d["trunc_param"])+'.csv','w') as f_handle:
        #    np.savetxt(f_handle,FR_vector_signal_new)
        np.savetxt('Stim_Signal/FR_vector_signal_high_ampl'+str(d["trunc_param"])+'.csv', FR_vector_signal_new, delimiter=" ")
        
        print('New frequency vector can be found in Stim_Signal/FR_vector_signal_high_ampl'+str(d["trunc_param"])+'.csv')
        
        Xs_storage_trunc=np.zeros((np.real(Xs_signal_norm_new).shape[0],2),float)
        
        Xs_storage_trunc[:,0]=np.real(Xs_signal_norm_new)
        Xs_storage_trunc[:,1]=np.imag(Xs_signal_norm_new)
        #with open('Stim_Signal/Xs_storage_high_ampl'+str(d["trunc_param"])+'.csv','w') as f_handle:
        #    np.savetxt(f_handle,Xs_storage_trunc)
        np.savetxt('Stim_Signal/Xs_storage_high_ampl'+str(d["trunc_param"])+'.csv', Xs_storage_trunc, delimiter=" ")
        print(" ")
        
        return FR_vector_signal_new,Xs_signal_norm_new
            
    if d["spectrum_trunc_method"]=='High Amplitude Method' and (os.path.isfile('Stim_Signal/FR_vector_signal_high_ampl'+str(d["trunc_param"])+'.csv')):
        Xs_recovered_new = np.genfromtxt('Stim_Signal/Xs_storage_high_ampl'+str(d["trunc_param"])+'.csv', delimiter=' ')
        Xs_signal_norm_new=np.vectorize(complex)(Xs_recovered_new[:,0],Xs_recovered_new[:,1])
        FR_vector_signal_new = np.genfromtxt('Stim_Signal/FR_vector_signal_high_ampl'+str(d["trunc_param"])+'.csv', delimiter=' ')
        
        return FR_vector_signal_new,Xs_signal_norm_new
    
    if d["spectrum_trunc_method"]=='Cutoff Method' and not(os.path.isfile('Stim_Signal/FR_vector_signal_cutoff'+str(d["trunc_param"])+'.csv')):
        from Truncation_of_spectrum import get_specific_freq        
        FR_vector_signal_new,Xs_signal_norm_new = get_specific_freq(d["spectrum_trunc_method"],d["freq"],d["trunc_param"],FR_vector_signal,Xs_signal_norm,t_vector,d["T"])
        #with open('Stim_Signal/FR_vector_signal_cutoff'+str(d["trunc_param"])+'.csv','w') as f_handle:
        #    np.savetxt(f_handle,FR_vector_signal_new)
        np.savetxt('Stim_Signal/FR_vector_signal_cutoff'+str(d["trunc_param"])+'.csv', FR_vector_signal_new, delimiter=" ")
        
        print('New frequency vector can be found in Stim_Signal/FR_vector_signal_cutoff'+str(d["trunc_param"])+'.csv')
        
        Xs_storage_trunc=np.zeros((np.real(Xs_signal_norm_new).shape[0],2),float)
        
        Xs_storage_trunc[:,0]=np.real(Xs_signal_norm_new)
        Xs_storage_trunc[:,1]=np.imag(Xs_signal_norm_new)
        #with open('Stim_Signal/Xs_storage_cutoff'+str(d["trunc_param"])+'.csv','w') as f_handle:
        #    np.savetxt(f_handle,Xs_storage_trunc)
        np.savetxt('Stim_Signal/Xs_storage_cutoff'+str(d["trunc_param"])+'.csv', Xs_storage_trunc, delimiter=" ")
        print(" ")
        
        return FR_vector_signal_new,Xs_signal_norm_new
            
    if d["spectrum_trunc_method"]=='Cutoff Method' and (os.path.isfile('Stim_Signal/FR_vector_signal_cutoff'+str(d["trunc_param"])+'.csv')):
        Xs_recovered_new = np.genfromtxt('Stim_Signal/Xs_storage_cutoff'+str(d["trunc_param"])+'.csv', delimiter=' ')
        Xs_signal_norm_new=np.vectorize(complex)(Xs_recovered_new[:,0],Xs_recovered_new[:,1])
        FR_vector_signal_new = np.genfromtxt('Stim_Signal/FR_vector_signal_cutoff'+str(d["trunc_param"])+'.csv', delimiter=' ')
    
        return FR_vector_signal_new,Xs_signal_norm_new
    
    
    
    
    
    