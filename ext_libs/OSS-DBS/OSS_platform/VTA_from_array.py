#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 14:09:42 2020

@author: butenko
"""

import h5py
import numpy as np
from dolfin import *
import matplotlib.pyplot as plt
from pandas import read_csv
import time as time_lib

import nibabel as nib

#This script allows to use VTA arrays of points instead of axons in OSS-DBS

def create_VTA_array(Xt,Yt,Zt):    #in mm, in MRI space
    [__,__,__,__,__,__,__,__,__,MRI_voxel_size_x,MRI_voxel_size_y,MRI_voxel_size_z]=np.genfromtxt('/opt/Patient/MRI_DTI_derived_data/MRI_misc.csv', delimiter=' ')
    
    VTA_res=max(MRI_voxel_size_x,MRI_voxel_size_y,MRI_voxel_size_z)
    #VTA_res=0.5
    VTA_box_length=20.0       # can be adjusted

    x_vector=np.arange(Xt-VTA_box_length/2.0,Xt+VTA_box_length/2.0+VTA_res,VTA_res)
    y_vector=np.arange(Yt-VTA_box_length/2.0,Yt+VTA_box_length/2.0+VTA_res,VTA_res)
    z_vector=np.arange(Zt-VTA_box_length/2.0,Zt+VTA_box_length/2.0+VTA_res,VTA_res)
    
    VTA_array=np.zeros((x_vector.shape[0]*y_vector.shape[0]*z_vector.shape[0],3),float)
    
    total_counter=0
    for i in range(x_vector.shape[0]):
       for j in range(y_vector.shape[0]): 
           for k in range(z_vector.shape[0]):   
               VTA_array[total_counter,:]=(x_vector[i],y_vector[j],z_vector[k])
               total_counter+=1
           
    #hf = h5py.File('VTA_default_array.h5', 'a')
    #hf.create_dataset('VTA_default', data=VTA_array)
    #hf.close()      
    
    np.savetxt('/opt/Patient/VTA_default_array.csv', VTA_array, delimiter=" ")
           
    return(x_vector.shape[0],'VTA_default_array.csv',VTA_res)

def resave_as_verts(array_name):    #in mm, in MRI space

    mesh = Mesh("/opt/Patient/Meshes/Mesh_unref.xml")
    [__,__,__,x_min,y_min,z_min,__,__,__,MRI_voxel_size_x,MRI_voxel_size_y,MRI_voxel_size_z]=np.genfromtxt('/opt/Patient/MRI_DTI_derived_data/MRI_misc.csv', delimiter=' ')

    arrays_shapes=[]

    Array_coord_get=read_csv('/opt/Patient/'+array_name, delimiter=' ', header=None)   
    Array_coord=Array_coord_get.values  
    
    for j in range(Array_coord.shape[0]):
        pnt=Point(Array_coord[j,0]-x_min,Array_coord[j,1]-y_min,Array_coord[j,2]-z_min)
        if not(mesh.bounding_box_tree().compute_first_entity_collision(pnt)<mesh.num_cells()*10): 
            Array_coord[j,:]=-100000000.0
#                    
    Array_coord=Array_coord[~np.all(Array_coord==-100000000.0,axis=1)] 
    arrays_shapes.append(Array_coord.shape[0])        #save to use later to recognize the array. Make sure you know the order!

    # shift to the positive octant space
    Array_coord[:,0]=Array_coord[:,0]-x_min
    Array_coord[:,1]=Array_coord[:,1]-y_min
    Array_coord[:,2]=Array_coord[:,2]-z_min    

        
#    if array_name[-3:]=='.h5':
#        hf = h5py.File(array_name, 'r')
#        lst=list(hf.keys())
#        result_total=[]
#        arrays_shapes=[]
#        
#        for i in lst:
#            a=hf.get(i)
#            a=np.array(a)
#            
#            for j in range(a.shape[0]):
#                pnt=Point(a[j,0]-x_min,a[j,1]-y_min,a[j,2]-z_min)
#                if not(mesh.bounding_box_tree().compute_first_entity_collision(pnt)<mesh.num_cells()*10): 
#                    a[j,:]=-100000000.0
#                    
#            a=a[~np.all(a==-100000000.0,axis=1)] 
#            arrays_shapes.append(a.shape[0])        #save to use later to recognize the array. Make sure you know the order!
#            
#            result_total.append(a)  
#        
#        Array_coord=np.concatenate(result_total)
#        hf.close()
#      
#        # shift to the positive octant space
#        Array_coord[:,0]=Array_coord[:,0]-x_min
#        Array_coord[:,1]=Array_coord[:,1]-y_min
#        Array_coord[:,2]=Array_coord[:,2]-z_min    
        
    
    np.savetxt('/opt/Patient/Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', Array_coord, delimiter=" ")  

    return arrays_shapes

#here we need to add function that will check that these vertices do not intersect with CSF, encap and so on

# do not use parallelization here


# This is outdated
def ifft_on_VTA_array(name_sol,d,FREQ_vector_signal,Xs_signal_normalized,t_vect,T,i_start_octv,arrays_shape):    #in mm, in MRI space

    start_IFFT=time_lib.time()
    
    num_segments=sum(arrays_shape)

    Max_signal_for_point=np.zeros(num_segments,float)

    hf = h5py.File('/opt/Patient/'+name_sol[:-4]+'.h5', 'r')
    solution_sort_octv = hf.get('dataset_1')
    solution_sort_octv = np.array(solution_sort_octv)
    hf.close() 

    Fr_corresp_ar = np.genfromtxt('/opt/Patient/Stim_Signal/Fr_corresp_array'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')
    #Fr_octave_vect = np.genfromtxt('Stim_Signal/Fr_octave_vector_'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')
    FR_vec_sign_octv = np.genfromtxt('/opt/Patient/Stim_Signal/FR_vector_signal_octaves'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')
        
    Fr_corresp_ar=np.round(Fr_corresp_ar,6)
    #Fr_octave_vect=np.round(Fr_octave_vect,6)
    N_freq_octv=(FR_vec_sign_octv.shape[0])
        
    
    for i_point in range(num_segments):
        Xs_Tr=np.vectorize(complex)(solution_sort_octv[(i_point)*N_freq_octv:(i_point*N_freq_octv+N_freq_octv),3],solution_sort_octv[i_point*N_freq_octv:(i_point*N_freq_octv+N_freq_octv),4])         #real and im parts for the first point in VTA
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
        
        Max_signal_for_point[i_point]=abs(max(Signal_t_conv[:], key=abs))
        
        if i_point==1:
            plt.figure(11111234)
            plt.plot(t_vect,Signal_t_conv)
            plt.xlim(0.000,T*5)
            plt.grid(True)
            plt.xlabel('t, sec')
            plt.ylabel('Potential, V')
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            plt.savefig('/opt/Patient/Images/Signal_convoluted_1st_point.png', format='png', dpi=500)
       
        np.save('/opt/Patient/Points_in_time/Signal_t_conv'+str(i_point), Signal_t_conv.real)

    minutes=int((time_lib.time() - start_IFFT)/60)
    secnds=int(time_lib.time() - start_IFFT)-minutes*60
    print("----- IFFT took ",minutes," min ",secnds," s -----") 
        
    return(Max_signal_for_point)
   
def get_VTA(d,array_full_name,Max_signal_for_point,arrays_shape,vox_along_axis,VTA_res):

    VTA_Vertices_get=read_csv('/opt/Patient/Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', delimiter=' ', header=None)    # get only physiologically correct neuron models 
    VTA_Vertices=VTA_Vertices_get.values    

    Array_full_coord_get=read_csv('/opt/Patient/'+array_full_name, delimiter=' ', header=None)    # get only physiologically correct neuron models 
    Array_full_coord=Array_full_coord_get.values  

#    if array_full_name[-3:]=='.h5':
#        hf = h5py.File(array_full_name, 'r')
#        lst=list(hf.keys())
#        result_total=[]
#        arrays_shapes=[]
#        
#        for i in lst:
#            a=hf.get(i)
#            a=np.array(a)
#            
#            arrays_shapes.append(a.shape[0])        #save to use later to recognize the array. Make sure you know the order!            
#            result_total.append(a)  
#        
#        Array_full_coord=np.concatenate(result_total)
#        hf.close()
    
    
    VTA_affected=np.zeros((VTA_Vertices.shape[0],4),float)
    VTA_affected[:,:3]=VTA_Vertices
    
    VTA_size=0.0
    
    for i in range(VTA_Vertices.shape[0]):
        if abs(Max_signal_for_point[i])>=d["Activation_threshold_VTA"]:
            VTA_affected[i,3]=1.0
            VTA_size+=VTA_res**3
#            if i<arrays_shape[0]:
#                VTA_EPN += VTA_res**3
#            else:
#                VTA_outside += VTA_res**3
                
    print("VTA_size: ",VTA_size)

    [__,__,__,x_min,y_min,z_min,__,__,__,MRI_voxel_size_x,MRI_voxel_size_y,MRI_voxel_size_z]=np.genfromtxt('/opt/Patient/MRI_DTI_derived_data/MRI_misc.csv', delimiter=' ')            
    shift_to_MRI_space=np.array([x_min,y_min,z_min])

    VTA_affected_MRI_space=np.zeros((VTA_Vertices.shape[0],4),float)
    VTA_affected_MRI_space[:,:3]=VTA_affected[:,:3]+shift_to_MRI_space
    VTA_affected_MRI_space[:,3]=VTA_affected[:,3]
                
    np.savetxt('/opt/Patient/Field_solutions/VTA_affected.csv', VTA_affected, delimiter=" ")
    
    VTA_nifti=np.zeros((vox_along_axis,vox_along_axis,vox_along_axis),int)
    E_field_nifti=np.zeros((vox_along_axis,vox_along_axis,vox_along_axis),float)

    #abs_xv_yv_index=xv_index+yv_index*vox_along_axis
    #abs_xv_yv_zv_index=xv_index+yv_index*vox_along_axis+vox_along_axis*vox_along_axis*zv_index

    # will throw an error, because we need to have the same number of points (no extractions)
    counter_truncated=0
    
    print("vox_along_axis :",vox_along_axis)
    
    for k in range(vox_along_axis):  #go over all voxels
        for j in range(vox_along_axis):  #go over all voxels
            for i in range(vox_along_axis):  #go over all voxels
                total_counter=i+j*vox_along_axis+k*vox_along_axis*vox_along_axis
                
                if np.all(np.round(VTA_affected_MRI_space[counter_truncated,:3],6)==np.round(Array_full_coord[total_counter,:],6)):            # if coordinates match, then
                    VTA_nifti[i,j,k]=int(VTA_affected_MRI_space[counter_truncated,3])
                    E_field_nifti[i,j,k]=Max_signal_for_point[counter_truncated]
                    counter_truncated+=1
                else:
                    VTA_nifti[i,j,k]=0
                    E_field_nifti[i,j,k]=0.0

    if counter_truncated!=VTA_affected_MRI_space.shape[0]:
        print("Hasn't iterated over whole VTA_affected_MRI_space, check the algorithm")
        raise SystemExit
                
                
    affine_info=np.eye(4)    
    affine_info[0,0]=VTA_res   # always isotropic voxels for VTA array
    affine_info[1,1]=VTA_res
    affine_info[2,2]=VTA_res
    affine_info[:3,3]=VTA_affected_MRI_space[0,:3]


    import os
    example_filename = os.path.join('/opt/Patient/'+d['MRI_data_name'])
    img = nib.load(example_filename)
        
    img3 = nib.Nifti1Image(VTA_nifti, affine_info,img.header)
    nib.save(img3, '/opt/Patient/VTA_solution.nii.gz')
    
    img4 = nib.Nifti1Image(E_field_nifti, affine_info,img.header)
    nib.save(img4, '/opt/Patient/E_field_solution.nii.gz')
                
    return True
                
        
#def VTA_from_
