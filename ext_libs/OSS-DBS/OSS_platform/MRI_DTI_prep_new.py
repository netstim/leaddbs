# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 14:56:19 2018

@author: Butenko K.
"""

#obtain_MRI_class and obtain_DTI_class are manager functions (called in Launcher)
#The functions in the script import nifti voxel arrays or grid voxel arrays (COMSOL format) and resave them in the platform's format.


from dolfin import *
import numpy as np
import time
import pickle

class MRI_info:
    def __init__(self,MRI_n,Mx,My,Mz,voxel_size_x,voxel_size_y,voxel_size_z,x_min,y_min,z_min,x_max,y_max,z_max,x_shift,y_shift,z_shift):
        self.name=MRI_n               #name of the external MRI file
        self.M_x=Mx                   #number of voxels along x-axis
        self.M_y=My
        self.M_z=Mz
        self.x_vox_size=voxel_size_x        #size of voxel along x-axis
        self.y_vox_size=voxel_size_y
        self.z_vox_size=voxel_size_z
        self.x_min=x_min                    #x coordinate of the first voxel
        self.y_min=y_min
        self.z_min=z_min
        self.x_max=x_max
        self.y_max=y_max
        self.z_max=z_max
        self.x_shift=x_shift                #defined as -xmin. The shift is required to have MRI data in positive octant starting in (0,0,0)
        self.y_shift=y_shift
        self.z_shift=z_shift
        
class DTI_info:
    def __init__(self,MRI_n,Mx,My,Mz,voxel_size_x,voxel_size_y,voxel_size_z,x_st,y_st,z_st):
        self.name=MRI_n
        self.M_x=Mx
        self.M_y=My
        self.M_z=Mz
        self.x_vox_size=voxel_size_x
        self.y_vox_size=voxel_size_y
        self.z_vox_size=voxel_size_z
        self.x_start=x_st                 #relative x coordinate of the first voxel in DTI data (relative to MRI data)
        self.y_start=y_st
        self.z_start=z_st


def map_MRI(MRI_name,MRI_data_in_m,default_material,CSF_inx,WM_inx,GM_inx,from_grid_txt):      # exctracts MRI data from the grid txt file (COMSOL format)
    start_voxel=time.clock()

    if from_grid_txt==True:
        #kicks out strings with comments 
        infile = open('/opt/Patient/'+MRI_name,'r').readlines()
        with open('/opt/Patient/MRI_DTI_derived_data/Filtered_'+MRI_name,'w') as outfile:   
            for index,line in enumerate(infile):
                #if index != 0 and index != 4 and index != 5:      #Waxholm atlas space
                if index != 0 and index != 4:           #here, 1st and 5th line are comments (tissue_full.txt)
                    outfile.write(line)
    
        x_vector=[]
        y_vector=[]
        z_vector=[]
        spc=' '
            
        #extracts vectors and values in MRI slices
        with open('/opt/Patient/MRI_DTI_derived_data/Filtered_'+MRI_name, 'r') as f:
            for index,line in enumerate(f):
                if index==0:
                    line=line.rstrip()
                    x_vector=line
                    x_arr = np.fromstring(x_vector, dtype=np.float, sep=' ')
                if index==1:
                    line=line.rstrip()
                    y_vector=line
                    y_arr = np.fromstring(y_vector, dtype=np.float, sep=' ')
                if index==2:
                    line=line.rstrip()
                    z_vector=line
                    z_arr = np.fromstring(z_vector, dtype=np.float, sep=' ')
                if index>=3:
                    line=line.rstrip()
                    if index==3:
                        voxel_val=line+spc
                    else:
                        voxel_val=voxel_val_old+line+spc     #if voxel data is devided into slices, this will collect all of them
                    voxel_val_old=voxel_val
      
            voxel_array_temp = np.fromstring(voxel_val, dtype=np.float, sep=' ')
    
        Mx,My,Mz=(x_arr.shape[0],y_arr.shape[0],z_arr.shape[0])       #number of voxels along axes    
        
        if Mx*My*Mz!=voxel_array_temp.shape[0]:
            print("Error while processing MRI, maybe not all slices were extracted")
            raise SystemExit
    else:
        import nibabel as nib      #nibabel should be installed
        import os
        
        example_filename = os.path.join('/opt/Patient/'+MRI_name)
        img = nib.load(example_filename)
        img.shape
        tissue_array = img.get_fdata()
        voxel_array_temp=tissue_array.flatten('F')

        Mx,My,Mz=(tissue_array.shape[0],tissue_array.shape[1],tissue_array.shape[2])       #number of voxels along axes    
        
        if MRI_data_in_m==1:        #switch to mm if the MRI data is in m
            voxel_size_x=img.header.get_zooms()[0]*1000
            voxel_size_y=img.header.get_zooms()[1]*1000 
            voxel_size_z=img.header.get_zooms()[2]*1000
            img_start_x,img_start_y,img_start_z=(img.affine[0,3]*1000,img.affine[1,3]*1000,img.affine[2,3]*1000)
        else:
            voxel_size_x=img.header.get_zooms()[0]
            voxel_size_y=img.header.get_zooms()[1]
            voxel_size_z=img.header.get_zooms()[2]   
            img_start_x,img_start_y,img_start_z=(img.affine[0,3],img.affine[1,3],img.affine[2,3])
        
        x_arr=np.arange(img_start_x,img_start_x+voxel_size_x*Mx,voxel_size_x)
        y_arr=np.arange(img_start_y,img_start_y+voxel_size_y*My,voxel_size_y)
        z_arr=np.arange(img_start_z,img_start_z+voxel_size_z*Mz,voxel_size_z)

    voxel_arr=np.zeros(voxel_array_temp.shape[0],int)
    voxel_arr[voxel_array_temp==CSF_inx]=1        #changes indices to the internal notation
    voxel_arr[voxel_array_temp==WM_inx]=2
    voxel_arr[voxel_array_temp==GM_inx]=3
    
    i=0         #counter for voxels

    ## For Lead-DBS we assume that the MRI data is well defined    
    # for z_i in z_arr:
    #     for y_i in y_arr:
    #         for x_i in x_arr:
    
    #             '''If a voxel has an unassigned material or background (0), algorithm will take average from the surrounding cube (at least 14 out of 26 filled voxels are required)'''      
            
    #             if voxel_arr[i]==0:
    #                 arr_around=np.array([i-Mx*My-Mx-1,i-Mx*My-Mx,i-Mx*My-Mx+1,i-Mx*My-1,i-Mx*My,i-Mx*My+1,i-Mx*My+Mx-1,i-Mx*My+Mx,i-Mx*My+Mx+1, #lower 9
    #                                      i-Mx-1,i-Mx,i-Mx+1,i-1,i+1,i+Mx-1,i+Mx,i+Mx+1,                                                         #middle 8
    #                                      i+Mx*My-Mx-1,i+Mx*My-Mx,i+Mx*My-Mx+1,i+Mx*My-1,i+Mx*My,i+Mx*My+1,i+Mx*My+Mx-1,i+Mx*My+Mx,i+Mx*My+Mx+1])#upper 9
    #                 counter=0.0
    #                 summa=0.0

    #                 for indx in arr_around:
    #                     if indx<Mx*My*Mz and indx>=0:                        
    #                         if voxel_arr[indx]!=0:
    #                             summa=summa+voxel_arr[indx]
    #                             counter=counter+1
    #                 if counter>14:      #how many voxels are needed for interpolation/approximation
    #                     voxel_arr[i]=int(round(summa/counter))
    #                 else:
    #                     voxel_arr[i]=int(default_material)                #if no material was averaged, put to default
       
    #             i=i+1
    
    # i=0
    
    #np.savetxt('MRI_DTI_derived_data/Tissue_array_MRI.csv', voxel_arr.astype(int), fmt='%i', delimiter=" ")
    np.save('/opt/Patient/MRI_DTI_derived_data/Tissue_array_MRI', voxel_arr.astype('b'), allow_pickle=False, fix_imports=False)
    del voxel_arr,voxel_array_temp
    
    if MRI_data_in_m==1:        #switch to mm if the MRI data is in m
        x_arr[:] = [x_i * 1000.0 for x_i in x_arr]
        y_arr[:] = [y_i * 1000.0 for y_i in y_arr]
        z_arr[:] = [z_i * 1000.0 for z_i in z_arr]
        
    voxel_size_x=abs(round(x_arr[1]-x_arr[0],6))    #size of voxels along x-axis
    voxel_size_y=abs(round(y_arr[1]-y_arr[0],6))
    voxel_size_z=abs(round(z_arr[1]-z_arr[0],6))
    
    x_vector_MRI_Box=np.zeros(x_arr.shape[0],float)
    y_vector_MRI_Box=np.zeros(y_arr.shape[0],float)
    z_vector_MRI_Box=np.zeros(z_arr.shape[0],float)
    
    #shift vectors to the positive octant starting in (0,0,0)
    x_vector_MRI_Box[:] = [round(x_i-min(x_arr)+voxel_size_x,6) for x_i in x_arr]   
    y_vector_MRI_Box[:] = [round(y_i-min(y_arr)+voxel_size_y,6) for y_i in y_arr]
    z_vector_MRI_Box[:] = [round(z_i-min(z_arr)+voxel_size_z,6) for z_i in z_arr]
    
    np.savetxt('/opt/Patient/MRI_DTI_derived_data/x_vector_MRI_Box.csv', x_vector_MRI_Box, delimiter=" ")
    np.savetxt('/opt/Patient/MRI_DTI_derived_data/y_vector_MRI_Box.csv', y_vector_MRI_Box, delimiter=" ")
    np.savetxt('/opt/Patient/MRI_DTI_derived_data/z_vector_MRI_Box.csv', z_vector_MRI_Box, delimiter=" ")

    print("----- Preparation of MRI data took %s seconds -----" % (time.clock() - start_voxel))

    return (Mx,My,Mz,round(min(x_arr),6),round(min(y_arr),6),round(min(z_arr),6),round(max(x_arr),6),round(max(y_arr),6),round(max(z_arr),6),voxel_size_x,voxel_size_y,voxel_size_z)


def map_DTI(DTI_name,DTI_data_in_m,from_grid_txt):        # exctracts Tensor data from the grid txt file (COMSOL format)
    
    start_voxel=time.clock()
    
    if from_grid_txt==True:
    
        infile = open('/opt/Patient/'+DTI_name,'r').readlines()
        with open('/opt/Patient/MRI_DTI_derived_data/Filtered_'+DTI_name,'w') as outfile:
    
            #first, we check the length of z and y vectors (one line of DTI data corresponds to the same y and z coordinate, directions are separated by a comment)
            for index,line in enumerate(infile):
                if index==2:
                    line=line.rstrip()
                    y_vector=line
                    y_stps = np.fromstring(y_vector, dtype=np.float, sep=' ')
                    y_length=y_stps.shape[0]
                if index==3:
                    line=line.rstrip()
                    z_vector=line
                    z_stps = np.fromstring(z_vector, dtype=np.float, sep=' ')
                    z_length=z_stps.shape[0]
    
            #jumps over line with comments and writes only vectors
            for index,line in enumerate(infile):
                if index != 0 and index != 4 and index != 4+(y_length*z_length+1) and index != 4+2*(y_length*z_length+1) and (index < 4+3*(y_length*z_length+1) or index > 4+4*(y_length*z_length+1)) and index!=4+5*(y_length*z_length+1) and (index < 4+6*(y_length*z_length+1) or index > 4+8*(y_length*z_length+1)):
                    outfile.write(line)
    
        x_vector=[]
        y_vector=[]
        z_vector=[]
        spc=' '
    
        with open('/opt/Patient/MRI_DTI_derived_data/Filtered_'+DTI_name, 'r') as f:
            for index,line in enumerate(f):
                if index==0:
                    line=line.rstrip()
                    x_vector=line
                    x_arr = np.fromstring(x_vector, dtype=np.float, sep=' ')
                if index==1:
                    line=line.rstrip()
                    y_vector=line
                    y_arr = np.fromstring(y_vector, dtype=np.float, sep=' ')
                if index==2:
                    line=line.rstrip()
                    z_vector=line
                    z_arr = np.fromstring(z_vector, dtype=np.float, sep=' ')
                    
                '''Here we assumed the order c11,c21,c31,c22,c32,c33'''    
                if index>=3 and index<3+y_length*z_length:      #corresponds to collection of slices (separated or combined) for c11
                    line=line.rstrip()
                    if index==3:
                        voxel_val_c11=line+spc
                    else:
                        voxel_val_c11=voxel_val_old_c11+line+spc
                    voxel_val_old_c11=voxel_val_c11
                
                if index>=3+y_length*z_length and index<3+2*y_length*z_length:
                    line=line.rstrip()
                    if index==3+y_length*z_length:
                        voxel_val_c21=line+spc
                    else:
                        voxel_val_c21=voxel_val_old_c21+line+spc
                    voxel_val_old_c21=voxel_val_c21
                    
                if index>=3+2*y_length*z_length and index<3+3*y_length*z_length:
                    line=line.rstrip()
                    if index==3+2*y_length*z_length:
                        voxel_val_c31=line+spc
                    else:
                        voxel_val_c31=voxel_val_old_c31+line+spc
                    voxel_val_old_c31=voxel_val_c31
                    
                if index>=3+3*y_length*z_length and index<3+4*y_length*z_length:
                    line=line.rstrip()
                    if index==3+3*y_length*z_length:
                        voxel_val_c22=line+spc
                    else:
                        voxel_val_c22=voxel_val_old_c22+line+spc
                    voxel_val_old_c22=voxel_val_c22
                    
                if index>=3+4*y_length*z_length and index<3+5*y_length*z_length:
                    line=line.rstrip()
                    if index==3+4*y_length*z_length:
                        voxel_val_c32=line+spc
                    else:
                        voxel_val_c32=voxel_val_old_c32+line+spc
                    voxel_val_old_c32=voxel_val_c32
                    
                if index>=3+5*y_length*z_length and index<3+6*y_length*z_length:
                    line=line.rstrip()
                    if index==3+5*y_length*z_length:
                        voxel_val_c33=line+spc
                    else:
                        voxel_val_c33=voxel_val_old_c33+line+spc
                    voxel_val_old_c33=voxel_val_c33
                    
             
            voxel_arr_c11 = np.fromstring(voxel_val_c11, dtype=np.float, sep=' ')
            voxel_arr_c21 = np.fromstring(voxel_val_c21, dtype=np.float, sep=' ')
            voxel_arr_c31 = np.fromstring(voxel_val_c31, dtype=np.float, sep=' ')
            voxel_arr_c22 = np.fromstring(voxel_val_c22, dtype=np.float, sep=' ')
            voxel_arr_c32 = np.fromstring(voxel_val_c32, dtype=np.float, sep=' ')
            voxel_arr_c33 = np.fromstring(voxel_val_c33, dtype=np.float, sep=' ')
                
        
        Tensor_array=np.zeros((voxel_arr_c11.shape[0],6),float)
        
        Mx,My,Mz=(x_arr.shape[0],y_arr.shape[0],z_arr.shape[0])       #number of voxels along axes    
        
        if Mx*My*Mz*6!=voxel_arr_c11.shape[0]+voxel_arr_c21.shape[0]+voxel_arr_c31.shape[0]+voxel_arr_c22.shape[0]+voxel_arr_c32.shape[0]+voxel_arr_c33.shape[0]:
            print("Error while processing DTI data, maybe not all slices were extracted")
            raise SystemExit
            
        if DTI_data_in_m==1:                  #if DTI data is in m
            x_arr[:] = [x_i * 1000.0 for x_i in x_arr]
            y_arr[:] = [y_i * 1000.0 for y_i in y_arr]
            z_arr[:] = [z_i * 1000.0 for z_i in z_arr]
            
        voxel_size_x=abs(round(x_arr[1]-x_arr[0],6))    #size of voxel along x-axis
        voxel_size_y=abs(round(y_arr[1]-y_arr[0],6))
        voxel_size_z=abs(round(z_arr[1]-z_arr[0],6))
    else:
        import nibabel as nib      #nibabel should be installed
        import os
        
        example_filename = os.path.join('/opt/Patient/'+DTI_name)
        img = nib.load(example_filename)
        img.shape
        tissue_array = img.get_fdata()
    
        voxel_arr_c11=tissue_array[:,:,:,0].flatten('F')
        voxel_arr_c21=tissue_array[:,:,:,1].flatten('F')
        voxel_arr_c31=tissue_array[:,:,:,2].flatten('F')
        voxel_arr_c22=tissue_array[:,:,:,3].flatten('F')
        voxel_arr_c32=tissue_array[:,:,:,4].flatten('F')
        voxel_arr_c33=tissue_array[:,:,:,5].flatten('F')
        
        Mx,My,Mz=(tissue_array.shape[0],tissue_array.shape[1],tissue_array.shape[2])       #number of voxels along axes    
        
        if Mx*My*Mz*6!=voxel_arr_c11.shape[0]+voxel_arr_c21.shape[0]+voxel_arr_c31.shape[0]+voxel_arr_c22.shape[0]+voxel_arr_c32.shape[0]+voxel_arr_c33.shape[0]:
            print("Error while processing DTI data, maybe not all slices were extracted")
            raise SystemExit
            
        Tensor_array=np.zeros((voxel_arr_c11.shape[0],6),float)
            
        if DTI_data_in_m==1:        #switch to mm if the MRI data is in m
            voxel_size_x=img.header.get_zooms()[0]*1000
            voxel_size_y=img.header.get_zooms()[1]*1000 
            voxel_size_z=img.header.get_zooms()[2]*1000
            img_start_x,img_start_y,img_start_z=(img.affine[0,3]*1000,img.affine[1,3]*1000,img.affine[2,3]*1000)
        else:
            voxel_size_x=img.header.get_zooms()[0]
            voxel_size_y=img.header.get_zooms()[1]
            voxel_size_z=img.header.get_zooms()[2]   
            img_start_x,img_start_y,img_start_z=(img.affine[0,3],img.affine[1,3],img.affine[2,3])
        
        x_arr=np.arange(img_start_x,img_start_x+voxel_size_x*Mx,voxel_size_x)
        y_arr=np.arange(img_start_y,img_start_y+voxel_size_y*My,voxel_size_y)
        z_arr=np.arange(img_start_z,img_start_z+voxel_size_z*Mz,voxel_size_z)        
    
    x_vector_DTI_Box=np.zeros(x_arr.shape[0],float)
    y_vector_DTI_Box=np.zeros(y_arr.shape[0],float)
    z_vector_DTI_Box=np.zeros(z_arr.shape[0],float)
    
    x_vector_DTI_Box[:] = [round(x_i-min(x_arr)+voxel_size_x,6) for x_i in x_arr] #shift vectors to positive octant staring in 0,0,0 (DTI will be repositioned relatively to MRI in Tissue_marking.py)
    y_vector_DTI_Box[:] = [round(y_i-min(y_arr)+voxel_size_y,6) for y_i in y_arr]
    z_vector_DTI_Box[:] = [round(z_i-min(z_arr)+voxel_size_z,6) for z_i in z_arr]
    
    np.savetxt('/opt/Patient/MRI_DTI_derived_data/x_vector_DTI_Box.csv', x_vector_DTI_Box, delimiter=" ")
    np.savetxt('/opt/Patient/MRI_DTI_derived_data/y_vector_DTI_Box.csv', y_vector_DTI_Box, delimiter=" ")
    np.savetxt('/opt/Patient/MRI_DTI_derived_data/z_vector_DTI_Box.csv', z_vector_DTI_Box, delimiter=" ")
    
    
    i=0
    for z_i in z_arr:
        for y_i in y_arr:
            for x_i in x_arr:
                
                Tensor_array[i,0]=round(voxel_arr_c11[i],8)
                Tensor_array[i,1]=round(voxel_arr_c21[i],8)
                Tensor_array[i,2]=round(voxel_arr_c31[i],8)
                Tensor_array[i,3]=round(voxel_arr_c22[i],8)
                Tensor_array[i,4]=round(voxel_arr_c32[i],8)
                Tensor_array[i,5]=round(voxel_arr_c33[i],8)
          
                i=i+1
                
        
    #'''Initialyly the data should be in ascending order'''
    #'''To ensure that coordinate vectors go in ascending order'''       
    #voxel_array = voxel_array[voxel_array[:,0].argsort()]
    #voxel_array = voxel_array[voxel_array[:,1].argsort(kind='mergesort')]
    #voxel_array = voxel_array[voxel_array[:,2].argsort(kind='mergesort')]
    
    #np.savetxt('MRI_DTI_derived_data/Tensor_array_DTI.csv', Tensor_array, delimiter=" ")
    np.save('/opt/Patient/MRI_DTI_derived_data/Tensor_array_DTI', Tensor_array)

    del Tensor_array,voxel_arr_c11,voxel_arr_c21,voxel_arr_c31,voxel_arr_c22,voxel_arr_c32,voxel_arr_c33

    print("----- Preparation of DTI data took %s seconds -----" % (time.clock() - start_voxel))
    #print("File MRI_DTI_derived_data/voxel_array_DTI.csv was created and can be visualized in Paraview to check the processed DTI data")
  
    return (Mx,My,Mz,round(min(x_arr),6),round(min(y_arr),6),round(min(z_arr),6),voxel_size_x,voxel_size_y,voxel_size_z)


def obtain_MRI_class(inp_dict):
    
    if inp_dict["MRI_data_name"]==0:
        print("MRI_data_name was not provided, exiting")
        raise Exception('exit')
    
    if inp_dict["voxel_arr_MRI"]==0:		# 1 if MRI data were already processed by the platform and corresp. meta data were created
        print("--- processing provided MRI data")
        if inp_dict["MRI_data_name"][-3:]=='nii' or inp_dict["MRI_data_name"][-6:]=='nii.gz':
            [Mx_mri,My_mri,Mz_mri,x_min,y_min,z_min,x_max,y_max,z_max,MRI_voxel_size_x,MRI_voxel_size_y,MRI_voxel_size_z]=map_MRI(inp_dict["MRI_data_name"],inp_dict["MRI_in_m"],inp_dict["default_material"],inp_dict["CSF_index"],inp_dict["WM_index"],inp_dict["GM_index"],False)       #will also prepare voxel_array!        
        else:
            [Mx_mri,My_mri,Mz_mri,x_min,y_min,z_min,x_max,y_max,z_max,MRI_voxel_size_x,MRI_voxel_size_y,MRI_voxel_size_z]=map_MRI(inp_dict["MRI_data_name"],inp_dict["MRI_in_m"],inp_dict["default_material"],inp_dict["CSF_index"],inp_dict["WM_index"],inp_dict["GM_index"],True)       #will also prepare voxel_array!

        '''Save meta data for the future simulations with the current MRI data set'''
        MRI_misc=np.array([Mx_mri,My_mri,Mz_mri,x_min,y_min,z_min,x_max,y_max,z_max,MRI_voxel_size_x,MRI_voxel_size_y,MRI_voxel_size_z])
        np.savetxt('/opt/Patient/MRI_DTI_derived_data/MRI_misc.csv', MRI_misc, delimiter=" ")
        print("--- MRI meta data were created\n")
    else:
        [Mx_mri,My_mri,Mz_mri,x_min,y_min,z_min,x_max,y_max,z_max,MRI_voxel_size_x,MRI_voxel_size_y,MRI_voxel_size_z]=np.genfromtxt('/opt/Patient/MRI_DTI_derived_data/MRI_misc.csv', delimiter=' ')
        print("--- MRI meta data were loaded\n")
    
    x_shift,y_shift,z_shift=(-1*(x_min),-1*(y_min),-1*(z_min))  #shift of MRI to have it in the positive octant and start in (0,0,0)    
        
    MRI_param=MRI_info(inp_dict["MRI_data_name"],Mx_mri,My_mri,Mz_mri,MRI_voxel_size_x,MRI_voxel_size_y,MRI_voxel_size_z,x_min,y_min,z_min,x_max,y_max,z_max,x_shift,y_shift,z_shift)    
    with open('/opt/Patient/MRI_DTI_derived_data/MRI_class.file', "wb") as f:
        pickle.dump(MRI_param, f, pickle.HIGHEST_PROTOCOL)
        
    return MRI_param

def obtain_DTI_class(inp_dict,MRI_param):

    if inp_dict["voxel_arr_DTI"]==0:       #1 if DTI data were already processed by the platform and corresp. meta data were created          
        if inp_dict["DTI_data_name"][-3:]=='nii' or inp_dict["DTI_data_name"][-6:]=='nii.gz':
            [Mx_dti,My_dti,Mz_dti,x_min_dti,y_min_dti,z_min_dti,DTI_voxel_size_x,DTI_voxel_size_y,DTI_voxel_size_z]=map_DTI(inp_dict["DTI_data_name"],inp_dict["MRI_in_m"],False)        
        else:
            [Mx_dti,My_dti,Mz_dti,x_min_dti,y_min_dti,z_min_dti,DTI_voxel_size_x,DTI_voxel_size_y,DTI_voxel_size_z]=map_DTI(inp_dict["DTI_data_name"],inp_dict["MRI_in_m"],True)
                    
        x_start_dti=x_min_dti-MRI_param.x_min               #DTI can be shifted from the MRI origin (0,0,0) (but only to the positive direction).
        y_start_dti=y_min_dti-MRI_param.y_min
        z_start_dti=z_min_dti-MRI_param.z_min
    
        if  x_start_dti<0 or y_start_dti<0 or z_start_dti<0:
            print("DTI data is not in positive octant, exiting.")
            raise Exception('exit')
            
        DTI_misc=np.array([Mx_dti,My_dti,Mz_dti,x_min_dti,y_min_dti,z_min_dti,DTI_voxel_size_x,DTI_voxel_size_y,DTI_voxel_size_z,x_start_dti,y_start_dti,z_start_dti])
        '''Save meta data for the future simulations with the current MRI data set'''
        np.savetxt('/opt/Patient/MRI_DTI_derived_data/DTI_misc.csv', DTI_misc, delimiter=" ")
        print("--- DTI meta data were created\n")

    if inp_dict["voxel_arr_DTI"]==1:
        [Mx_dti,My_dti,Mz_dti,x_min_dti,y_min_dti,z_min_dti,DTI_voxel_size_x,DTI_voxel_size_y,DTI_voxel_size_z,x_start_dti,y_start_dti,z_start_dti]=np.genfromtxt('/opt/Patient/MRI_DTI_derived_data/DTI_misc.csv', delimiter=' ')
        print("--- DTI meta data were loaded\n")

    DTI_param=DTI_info(inp_dict["DTI_data_name"],Mx_dti,My_dti,Mz_dti,DTI_voxel_size_x,DTI_voxel_size_y,DTI_voxel_size_z,x_start_dti,y_start_dti,z_start_dti)        
    with open('/opt/Patient/MRI_DTI_derived_data/DTI_class.file', "wb") as f:
        pickle.dump(DTI_param, f, pickle.HIGHEST_PROTOCOL)
        
    return DTI_param
