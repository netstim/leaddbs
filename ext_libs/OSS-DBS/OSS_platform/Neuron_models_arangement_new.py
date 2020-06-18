# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 14:06:45 2018

@author: konstantin
"""
import h5py

import numpy as np
from dolfin import *
from pandas import read_csv
from Axon_files.axon import Axon
import time
import os
import pickle

# Butenko: the functions in this script work, but it is difficult to follow and too much of redunancy. Will be reritten in the next version.


class Neuron_info:      #contains spatial data on full neuron array
    def __init__(self,pattern_model_name,X_coord,Y_coord,Z_coord,YZ_angles=[0.0],ZX_angles=[0.0],XY_angles=[0.0],OX_angles=[0.0],OY_angles=[0.0],OZ_angles=[0.0],N_models_plane=0):
        self.name=pattern_model_name
        self.x_loc=X_coord
        self.y_loc=Y_coord
        self.z_loc=Z_coord
        self.alpha_array=YZ_angles
        self.beta_array=ZX_angles
        self.gamma_array=XY_angles
        self.alpha_array_global=OX_angles
        self.beta_array_global=OY_angles
        self.gamma_array_global=OZ_angles
        self.N_models_in_plane=N_models_plane



def rotate_globally(x_pos,y_pos,z_pos,alpha,beta,gamma):
           
    alpha_matrix=np.array([[1,0,0],
                           [0,cos(alpha),-1*sin(alpha)],
                           [0,sin(alpha),cos(alpha)]])
                           
    beta_matrix=np.array([[cos(beta),0,sin(beta)],
                           [0,1,0],
                           [-1*sin(beta),0,cos(beta)]])
    
                           
    gamma_matrix=np.array([[cos(gamma),-1*sin(gamma),0],
                           [sin(gamma),cos(gamma),0],
                           [0,0,1]])
    
    #for inx in xrange(Neuron_coord.shape[0]):
        
    xyz_alpha=np.array([x_pos,y_pos,z_pos])
    [x_pos_new,y_pos_new,z_pos_new]=alpha_matrix.dot(xyz_alpha)
        
    xyz_beta=np.array([x_pos_new,y_pos_new,z_pos_new])
    [x_pos_new,y_pos_new,z_pos_new]=beta_matrix.dot(xyz_beta)
                
    xyz_gamma=np.array([x_pos_new,y_pos_new,z_pos_new])
    [x_pos_new,y_pos_new,z_pos_new]=gamma_matrix.dot(xyz_gamma)

    Point_coords=np.zeros((1,3),float)
    Point_coords[0,:]=(x_pos_new,y_pos_new,z_pos_new)
    
    return Point_coords

def create_structured_array(d):
    
    x_coord=[]
    y_coord=[]
    z_coord=[]
       
    x_start_point=0-(d["x_step"]*(d["x_steps"])/2)
    y_start_point=0-(d["y_step"]*(d["y_steps"])/2)
    z_start_point=0-(d["z_step"]*(d["z_steps"])/2)
    
    for x_ind in range(int(d["x_steps"])+1):
        for y_ind in range(int(d["y_steps"])+1):
            for z_ind in range(int(d["z_steps"])+1):
                x_coord.append(x_start_point+x_ind*d["x_step"])
                y_coord.append(y_start_point+y_ind*d["y_step"])
                z_coord.append(z_start_point+z_ind*d["z_step"]) 
     
    for inx_angle in range(len(d["alpha_array_glob"])):           
        for inx in range(len(x_coord)):
           
            Rotated_point=rotate_globally(x_coord[inx],y_coord[inx],z_coord[inx],d["alpha_array_glob"][inx_angle],d["beta_array_glob"][inx_angle],d["gamma_array_glob"][inx_angle])
            
            if not('Array_coord_str' in locals()):
                Array_coord_str=Rotated_point
            else:
                Array_coord_str=np.append(Array_coord_str,Rotated_point,axis=0)
 
    x_coord=Array_coord_str[:,0]
    y_coord=Array_coord_str[:,1]
    z_coord=Array_coord_str[:,2]
         
    return (x_coord,y_coord,z_coord)
    
# takes a pattern model, rotates it around X,Y,Z and shifts its center to (x_loc,y_loc,z_loc)
def place_neuron(Arr_coord,alpha,beta,gamma,x_loc,y_loc,z_loc):
    
    Neuron_coord=np.zeros((Arr_coord.shape[0],3),float)
       
    alpha_matrix=np.array([[1,0,0],
                           [0,cos(alpha),-1*sin(alpha)],
                           [0,sin(alpha),cos(alpha)]])
                           
    beta_matrix=np.array([[cos(beta),0,sin(beta)],
                           [0,1,0],
                           [-1*sin(beta),0,cos(beta)]])
                               
    gamma_matrix=np.array([[cos(gamma),-1*sin(gamma),0],
                           [sin(gamma),cos(gamma),0],
                           [0,0,1]])
    
    for inx in range(Neuron_coord.shape[0]):
        
        xyz_alpha=np.array([Arr_coord[inx,0],Arr_coord[inx,1],Arr_coord[inx,2]])
        [Neuron_coord[inx,0],Neuron_coord[inx,1],Neuron_coord[inx,2]]=alpha_matrix.dot(xyz_alpha)
        
        xyz_beta=np.array([Neuron_coord[inx,0],Neuron_coord[inx,1],Neuron_coord[inx,2]])
        [Neuron_coord[inx,0],Neuron_coord[inx,1],Neuron_coord[inx,2]]=beta_matrix.dot(xyz_beta)
                
        xyz_gamma=np.array([Neuron_coord[inx,0],Neuron_coord[inx,1],Neuron_coord[inx,2]])
        [Neuron_coord[inx,0],Neuron_coord[inx,1],Neuron_coord[inx,2]]=gamma_matrix.dot(xyz_gamma)
    
    Neuron_coord[:,0]=np.round(Neuron_coord[:,0]+x_loc,6)
    Neuron_coord[:,1]=np.round(Neuron_coord[:,1]+y_loc,6)
    Neuron_coord[:,2]=np.round(Neuron_coord[:,2]+z_loc,6)
    
    return Neuron_coord
    

def generate_pattern_model(name,N_Ranv,axon_param,Axon_Model_Type):
    
    x_center=0.0
    y_center=0.0
    z_center=0.0    
    ZY_angle=0.0
    
    [ranvier_nodes, para1_nodes, para2_nodes, inter_nodes, ranvier_length, para1_length, para2_length, inter_length, deltax, fiberD]=(axon_param[:])  
    if Axon_Model_Type == 'McIntyre2002': 
        n_comp=((ranvier_nodes-1)+inter_nodes+para1_nodes+para2_nodes)/(ranvier_nodes-1)
        n_segments=int((N_Ranv-1)*n_comp+1)        #overall number of points on Axon incl. internodal
    elif Axon_Model_Type == 'Reilly2016': 
        n_comp=2
        n_segments=int((N_Ranv-1)*n_comp+1)        #overall number of points on Axon incl. internodal

    Array_coord_pattern=np.zeros((n_segments,3),float)
    z_start=z_center-(0.001*deltax*(N_Ranv-1)/2.0)*sin(ZY_angle)  #sin. theorem
    y_start=y_center-(0.001*deltax*(N_Ranv-1)/2.0)*sin(np.pi/2-ZY_angle)
    
    Array_coord_pattern[0,0]=x_center
    Array_coord_pattern[0,1]=y_start
    Array_coord_pattern[0,2]=z_start
    
    loc_inx=1       #because first node (with index 0) was already seeded
    if Axon_Model_Type == 'McIntyre2002': 
        if fiberD>3.0:
    
            for inx in range(1,n_segments):
                if loc_inx==0:
                    l_step=(para1_length+ranvier_length)/2000   #switch to mm from µm
                if loc_inx==1 or loc_inx==11:
                    l_step=(ranvier_length+para1_length)/2000   #switch to mm from µm
                if loc_inx==2 or loc_inx==10:
                    l_step=(para1_length+para2_length)/2000   #switch to mm from µm
                if loc_inx==3 or loc_inx==9:
                    l_step=(para2_length+inter_length)/2000   #switch to mm from µm
                if loc_inx==4 or loc_inx==5 or loc_inx==6 or loc_inx==7 or loc_inx==8:
                    l_step=(inter_length)/1000   #switch to mm from µm
                
                loc_inx=loc_inx+1
                Array_coord_pattern[inx,0]=x_center
                Array_coord_pattern[inx,1]=Array_coord_pattern[inx-1,1]+(l_step)*sin(np.pi/2-ZY_angle)
                Array_coord_pattern[inx,2]=Array_coord_pattern[inx-1,2]+(l_step)*sin(ZY_angle)
                if inx%n_comp==0:
                    loc_inx=0
                    
        else:
    
            for inx in range(1,n_segments):
                if loc_inx==0:
                    l_step=(para1_length+ranvier_length)/2000   #switch to mm from µm
                if loc_inx==1 or loc_inx==8:
                    l_step=(ranvier_length+para1_length)/2000   #switch to mm from µm
                if loc_inx==2 or loc_inx==7:
                    l_step=(para1_length+para2_length)/2000   #switch to mm from µm
                if loc_inx==3 or loc_inx==6:
                    l_step=(para2_length+inter_length)/2000   #switch to mm from µm
                if loc_inx==4 or loc_inx==5:
                    l_step=(inter_length)/1000   #switch to mm from µm
                
                loc_inx=loc_inx+1
                Array_coord_pattern[inx,0]=x_center
                Array_coord_pattern[inx,1]=Array_coord_pattern[inx-1,1]+(l_step)*sin(np.pi/2-ZY_angle)
                Array_coord_pattern[inx,2]=Array_coord_pattern[inx-1,2]+(l_step)*sin(ZY_angle)
                if inx%n_comp==0:
                    loc_inx=0
    elif Axon_Model_Type == 'Reilly2016':
        l_step=deltax/2000.0   # node - internodal - node only
        for inx in range(1,n_segments):
            Array_coord_pattern[inx,0]=x_center
            Array_coord_pattern[inx,1]=Array_coord_pattern[inx-1,1]+(l_step)*sin(np.pi/2-ZY_angle)
            Array_coord_pattern[inx,2]=Array_coord_pattern[inx-1,2]+(l_step)*sin(ZY_angle)
        
    
    Array_coord_pattern=np.round(Array_coord_pattern,8)
    
    np.savetxt('Neuron_model_arrays/'+name, Array_coord_pattern, delimiter=" ") 
        
    return n_segments

# builds or resaves the full neuron array in the positive octant, computes its extent (distance between the implantation point and the furtherst compartment).
def get_neuron_models_dims(Full_model_ready,X_imp,Y_imp,Z_imp,Neuron_param,axon_param,plane_rot):
    
    start_neuron_models=time.clock()
      
    if Full_model_ready==0:
        [ranvier_nodes, para1_nodes, para2_nodes, inter_nodes]=(axon_param[:4]) 
        if Neuron_param.name != 'pattern.csv':      #get a pattern from the externally provided file
            Array_coord_load_get=read_csv(Neuron_param.name, delimiter=' ', header=None)
            Array_coord_load=Array_coord_load_get.values    #
            
        else:   #or get a pattern from the internally created file
            Array_coord_load_get=read_csv('Neuron_model_arrays/'+Neuron_param.name, delimiter=' ', header=None)
            Array_coord_load=Array_coord_load_get.values
    
        n_segments=Array_coord_load.shape[0]        #number of compartments per neuron (axon)

        Array_coord_temp=Array_coord_load         #template model
        del Array_coord_load

        if plane_rot==0:    # if manual placement
            for inx in range(len(Neuron_param.x_loc)):  #goes through all seeding (central) nodes of neurons (axons) 
                for alpha_angle in Neuron_param.alpha_array:
                    for beta_angle in Neuron_param.beta_array:
                        for gamma_angle in Neuron_param.gamma_array:                
                            New_model_coord=place_neuron(Array_coord_temp,alpha_angle,beta_angle,gamma_angle,Neuron_param.x_loc[inx],Neuron_param.y_loc[inx],Neuron_param.z_loc[inx])
                            if not('Array_coord' in locals()):
                                Array_coord=New_model_coord
                            else:
                                Array_coord=np.append(Array_coord,New_model_coord,axis=0)
            
        if plane_rot==1:    # if placement in the ordered array (as for VTA)
            for inx in range(len(Neuron_param.x_loc)):     #goes through all seeding (central) nodes of neurons (axons)                                
                New_model_coord=place_neuron(Array_coord_temp,Neuron_param.alpha_array_global[int(inx/Neuron_param.N_models_in_plane)],Neuron_param.beta_array_global[int(inx/Neuron_param.N_models_in_plane)],Neuron_param.gamma_array_global[int(inx/Neuron_param.N_models_in_plane)],Neuron_param.x_loc[inx],Neuron_param.y_loc[inx],Neuron_param.z_loc[inx])               
                if not('Array_coord' in locals()):
                    Array_coord=New_model_coord
                else:
                    Array_coord=np.append(Array_coord,New_model_coord,axis=0)
                    
    if Full_model_ready==1:         # load the provided neuron array (already preprocessed by cut_models_by_domain)
    
        # load MRI meta data to check the shift of the coordinate system
        [Mx_mri,My_mri,Mz_mri,x_min,y_min,z_min,x_max,y_max,z_max,MRI_voxel_size_x,MRI_voxel_size_y,MRI_voxel_size_z]=np.genfromtxt('MRI_DTI_derived_data/MRI_misc.csv', delimiter=' ')
        
        if os.path.isfile('Neuron_model_arrays/Prepared_models_full_filtered.csv'):              # if the the neuron array was provided in .csv   
            Array_coord_get=read_csv('Neuron_model_arrays/Prepared_models_full_filtered.csv', delimiter=' ', header=None)        #
            Array_coord=Array_coord_get.values
                
        if os.path.isfile('Neuron_model_arrays/Prepared_models_full_filtered.h5'):              # if the the neuron array was provided in .h5
            hf = h5py.File('Neuron_model_arrays/Prepared_models_full_filtered.h5', 'r')
            lst=list(hf.keys())
            result_total=[]
            for i in lst:
                a=hf.get(i)
                a=np.array(a)
                result_total.append(a)               
                b=np.copy(a)
                
                hf3 = h5py.File('Neuron_model_arrays/All_neuron_models_by_populations.h5', 'a')     # resave to .h5
                b[:,0]=b[:,0]-x_min           # shift to the positive octant space
                b[:,1]=b[:,1]-y_min
                b[:,2]=b[:,2]-z_min
                hf3.create_dataset(i, data=b)
                hf3.close()
                
            Array_coord=np.concatenate(result_total)
            hf.close()
      
        # shift to the positive octant space
        Array_coord[:,0]=Array_coord[:,0]-x_min
        Array_coord[:,1]=Array_coord[:,1]-y_min
        Array_coord[:,2]=Array_coord[:,2]-z_min
                
    np.savetxt('Neuron_model_arrays/All_neuron_models.csv', Array_coord, delimiter=" ") 
        
    print("Initial neuron models can be visualized from Neuron_model_arrays/All_neuron_models.csv in Paraview\n")

    # find the maximum distance from the meuron compartment to the implantation site
    max_length=0.0
    i_dist_comp=0
    for i in range(Array_coord.shape[0]):        
        if sqrt((Array_coord[i,0]-X_imp)**2+(Array_coord[i,1]-Y_imp)**2+(Array_coord[i,2]-Z_imp)**2)>max_length:
            max_length=sqrt((Array_coord[i,0]-X_imp)**2+(Array_coord[i,1]-Y_imp)**2+(Array_coord[i,2]-Z_imp)**2)
            i_dist_comp=i
    print("Max distance from a compartment in the neuron array to the implantation site: ", max_length)
    #print("The compartment index: ", i_dist_comp)
    return max_length
       

def generate_neuron_models(N_Ranv,Full_model_ready,Domains,MRI_param,Neuron_param,axon_param,plane_rot):
    
    import os
    start_neuron_models=time.clock()
        
    if Full_model_ready==0:             # if the neuron array was created internally
        same_axon_type=True             # for now, can't create different neuron morphologies internally
        if Neuron_param.name!='pattern.csv':
            Array_coord_load_get=read_csv(Neuron_param.name, delimiter=' ', header=None)
            Array_coord_load=Array_coord_load_get.values        
        else:
            Array_coord_load_get=read_csv('Neuron_model_arrays/'+Neuron_param.name, delimiter=' ', header=None)
            Array_coord_load=Array_coord_load_get.values
    
        n_segments=Array_coord_load.shape[0]

        Array_coord_get=read_csv('Neuron_model_arrays/All_neuron_models.csv', delimiter=' ', header=None)        #
        Array_coord=Array_coord_get.values
                    
    if Full_model_ready == 1:             # if the neuron array was provided
                         
        if os.path.isfile('Neuron_model_arrays/All_neuron_models.csv'):
            same_axon_type=True            
            Array_coord_get=read_csv('Neuron_model_arrays/All_neuron_models.csv', delimiter=' ', header=None)        #
            Array_coord=Array_coord_get.values
            [Mx_mri,My_mri,Mz_mri,x_min,y_min,z_min,x_max,y_max,z_max,MRI_voxel_size_x,MRI_voxel_size_y,MRI_voxel_size_z]=np.genfromtxt('MRI_DTI_derived_data/MRI_misc.csv', delimiter=' ')
            [__, __, __, __, __, __, __, __, __, __,__,__,n_segments]=np.genfromtxt('Neuron_model_arrays/Neuron_model_misc.csv', delimiter=' ')
                           
        if os.path.isfile('Neuron_model_arrays/All_neuron_models_by_populations.h5'):
            
            same_axon_type=False            
            
            hf = h5py.File('Neuron_model_arrays/All_neuron_models_by_populations.h5', 'r')
            lst=list(hf.keys())
            List_of_arrays=[]
            List_of_empty=[]
            List_of_placed=[]
                       
            [Mx_mri,My_mri,Mz_mri,x_min,y_min,z_min,x_max,y_max,z_max,MRI_voxel_size_x,MRI_voxel_size_y,MRI_voxel_size_z]=np.genfromtxt('MRI_DTI_derived_data/MRI_misc.csv', delimiter=' ')
            
            for i in lst:
                Array_coord=hf.get(i)
                Array_coord=np.array(Array_coord)      
                
                if Array_coord.size == 0:
                    List_of_empty.append(i)
                    List_of_arrays.append(0)                               
                else:
                    List_of_arrays.append(Array_coord)                    
                    List_of_placed.append(Array_coord) 
                                          
                del Array_coord
            hf.close()
            
            n_segments_fib_diam_array=np.load('Neuron_model_arrays/Neuron_populations_misc.npy')
            N_models=np.zeros((len(List_of_arrays)),int)    # number of models per population  

    n_segments=int(n_segments)  # if loaded from .csv

    # now we will folter out unphysiological neurins    
    mesh = Mesh("Meshes/Mesh_unref.xml")
    subdomains_assigned=MeshFunction('size_t',mesh,"Meshes/Mesh_unref_physical_region.xml")

    #first, the neuron compartments should not pass through CSF
    import os      
    if os.path.isfile('MRI_DTI_derived_data/'+MRI_param.name[:-4]+'_voxel_array_CSF.npy') or os.path.isfile('MRI_DTI_derived_data/'+MRI_param.name[:-7]+'_voxel_array_CSF.npy'):        #if array was already prepared
        if MRI_param.name[-2:]=='gz':
            voxel_array_CSF=np.load('MRI_DTI_derived_data/'+MRI_param.name[:-7]+'_voxel_array_CSF.npy')
        else:   
            voxel_array_CSF=np.load('MRI_DTI_derived_data/'+MRI_param.name[:-4]+'_voxel_array_CSF.npy') 
        print("voxel_array_CSF is loaded")    
    else:   #otherwise prepare an array that stores coordinated of all voxels with CSF in the vicinity of the neurons    
        Tissue_array=np.load('MRI_DTI_derived_data/Tissue_array_MRI.npy')
        
        x_vect=np.genfromtxt('MRI_DTI_derived_data/x_vector_MRI_Box.csv', delimiter=' ')
        y_vect=np.genfromtxt('MRI_DTI_derived_data/y_vector_MRI_Box.csv', delimiter=' ')
        z_vect=np.genfromtxt('MRI_DTI_derived_data/z_vector_MRI_Box.csv', delimiter=' ')
                
        voxel_size_x=MRI_param.x_vox_size
        voxel_size_y=MRI_param.y_vox_size
        voxel_size_z=MRI_param.z_vox_size
        
        voxel_array_CSF=np.zeros((Tissue_array.shape[0],3),float)      #array to store all CSF voxels in the specified ROI
        
        Mx_mri=MRI_param.M_x
        My_mri=MRI_param.M_y
        Mz_mri=MRI_param.M_z
        
        if same_axon_type==True:
            x_neuron_max,y_neuron_max,z_neuron_max=(max(Array_coord[:,0]),max(Array_coord[:,1]),max(Array_coord[:,2]))
            x_neuron_min,y_neuron_min,z_neuron_min=(min(Array_coord[:,0]),min(Array_coord[:,1]),min(Array_coord[:,2]))
        else:
            max_values=np.zeros((len(List_of_placed),3),float)
            min_values=np.zeros((len(List_of_placed),3),float)
            for i in range(len(List_of_placed)):
                max_values[i,:]=(max(List_of_placed[i][:,0]),max(List_of_placed[i][:,1]),max(List_of_placed[i][:,2]))
                min_values[i,:]=(min(List_of_placed[i][:,0]),min(List_of_placed[i][:,1]),min(List_of_placed[i][:,2]))

            x_neuron_max,y_neuron_max,z_neuron_max=(max(max_values[:,0]),max(max_values[:,1]),max(max_values[:,2]))
            x_neuron_min,y_neuron_min,z_neuron_min=(min(min_values[:,0]),min(min_values[:,0]),min(min_values[:,0]))
        
        space_from_neurons=1.0          #here we do not need to check further away
        for z_coord in z_vect:
            for y_coord in y_vect:
                for x_coord in x_vect:
                    
                    x_pos=x_coord-voxel_size_x/2.0
                    y_pos=y_coord-voxel_size_y/2.0
                    z_pos=z_coord-voxel_size_z/2.0
                    
                    if (x_pos<=x_neuron_max+space_from_neurons and x_pos>=x_neuron_min-space_from_neurons and y_pos<=y_neuron_max+space_from_neurons and y_pos>=y_neuron_min-space_from_neurons and z_pos<=z_neuron_max+space_from_neurons and z_pos>=z_neuron_min-space_from_neurons):
                    
                        xv_mri=int((x_coord)/voxel_size_x-0.000000001)                                  #defines number of steps to get to the voxels containing x[0] coordinate
                        yv_mri=(int((y_coord)/voxel_size_y-0.000000001))*Mx_mri                  #defines number of steps to get to the voxels containing x[0] and x[1] coordinates
                        zv_mri=(int((z_coord)/voxel_size_z-0.000000001))*Mx_mri*My_mri           #defines number of steps to get to the voxels containing x[0], x[1] and x[2] coordinates
                        
                        glob_index=xv_mri+yv_mri+zv_mri
                        glob_index=int(glob_index)
                                        
                        if Tissue_array[glob_index]==1:
                            voxel_array_CSF[glob_index,0]=x_pos
                            voxel_array_CSF[glob_index,1]=y_pos
                            voxel_array_CSF[glob_index,2]=z_pos
                        
        voxel_array_CSF=voxel_array_CSF[~np.all(voxel_array_CSF==0.0,axis=1)]  #deletes all zero enteries

        if MRI_param.name[-2:]=='gz':
            np.save('MRI_DTI_derived_data/'+MRI_param.name[:-7]+'_voxel_array_CSF', voxel_array_CSF) 
        else:
            np.save('MRI_DTI_derived_data/'+MRI_param.name[:-4]+'_voxel_array_CSF', voxel_array_CSF)
                             
        del Tissue_array,x_vect,y_vect,z_vect
        print("voxel_array_CSF (contains CSF voxels close to the neuron array) is prepared")

    voxel_array_CSF_shifted=np.zeros((voxel_array_CSF.shape[0],3),float)
    voxel_array_CSF_shifted[:,0]=voxel_array_CSF[:,0]+MRI_param.x_vox_size/2
    voxel_array_CSF_shifted[:,1]=voxel_array_CSF[:,1]+MRI_param.y_vox_size/2
    voxel_array_CSF_shifted[:,2]=voxel_array_CSF[:,2]+MRI_param.z_vox_size/2
    points_csf=0
    points_encap=0
    points_outside=0
    del voxel_array_CSF

    #also neurons should not pass through the encapsulation layer and floating contacts
    subdomains_enc = MeshFunction("size_t", mesh, mesh.topology().dim())
    subdomains_enc.set_all(0)
    
    for i in range(len(Domains.Encup_index)):
        subdomains_enc.array()[subdomains_assigned.array()==Domains.Encup_index[i]]=1
        
    subdomains_enc.array()[subdomains_assigned.array()==Domains.Float_contacts]=1       #neuron models cannot be located in floating conductors

    submesh_encup=SubMesh(mesh,subdomains_enc,1)

    inx=0
    if same_axon_type==True:
       
        while inx < Array_coord.shape[0]:
            #print(inx)
            pnt=Point(Array_coord[inx,0],Array_coord[inx,1],Array_coord[inx,2])
            
            if (submesh_encup.bounding_box_tree().compute_first_entity_collision(pnt)<submesh_encup.num_cells()):        #this is a condition to check whether the point is inside encap. layer or floating conductor
                points_encap=points_encap+1
                inx_start=int(inx/n_segments)*n_segments
                Array_coord[inx_start:inx_start+n_segments,:]=-100000000.0
                inx=inx_start+n_segments
            else:    
                if not(mesh.bounding_box_tree().compute_first_entity_collision(pnt)<mesh.num_cells()*100):        #if one point of the neural model is absent, the whole model is disabled
                   points_outside=points_outside+1
                   inx_start=int(inx/n_segments)*n_segments
                   Array_coord[inx_start:inx_start+n_segments,:]=-100000000.0
                   inx=inx_start+n_segments
                else:                                                                                               # finally checks whether the neuron compartment is inside the CSF voxel
                    check1_1=(voxel_array_CSF_shifted[:,0]-Array_coord[inx,0]<=MRI_param.x_vox_size)
                    check1_2=(voxel_array_CSF_shifted[:,1]-Array_coord[inx,1]<=MRI_param.y_vox_size)
                    check1_3=(voxel_array_CSF_shifted[:,2]-Array_coord[inx,2]<=MRI_param.z_vox_size)
                    check2_1=(voxel_array_CSF_shifted[:,0]>=Array_coord[inx,0])
                    check2_2=(voxel_array_CSF_shifted[:,1]>=Array_coord[inx,1])
                    check2_3=(voxel_array_CSF_shifted[:,2]>=Array_coord[inx,2])

                    check3=np.logical_and(np.logical_and(check1_1,check2_1),np.logical_and(np.logical_and(check1_2,check2_2),np.logical_and(check1_3,check2_3)))
                    a=np.where((check3 == (True)))                    
                    if str(a)!='(array([], dtype=int64),)':     
                        points_csf=points_csf+1
                        inx_start=int(inx/n_segments)*n_segments
                        Array_coord[inx_start:inx_start+n_segments,:]=-100000000.0
                        inx=inx_start+n_segments
                    else:
                        inx=inx+1

        print("Points in CSF, encapsulation layer (and floating conductors) and outside (and intersecting with the electrode): ",points_csf,points_encap,points_outside)       
        
        inx=0
        
        del voxel_array_CSF_shifted       
        
        Array_coord=Array_coord[~np.all(Array_coord==-100000000.0,axis=1)]  #deletes all -100000000.0 enteries          
        Array_coord=np.round(Array_coord,8)
        
        np.savetxt('Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', Array_coord, delimiter=" ")   
        print("Adjusted neuron models can be visualized from Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv in Paraview")
 
        import os
        import subprocess
        with open(os.devnull, 'w') as FNULL: subprocess.call('python Visualization_files/Paraview_csv_neurons.py', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
           
        N_models=int(Array_coord.shape[0]/n_segments)
        
        np.savetxt('Neuron_model_arrays/Adjusted_neuron_array_info.csv', np.array([N_models,points_csf,points_encap,points_outside]), delimiter=" ") 
        print("Number of placed neuron models: ",N_models)

                        
    if same_axon_type==False:
        
        list_of_connections=[]
        
        if os.path.isfile('Neuron_model_arrays/Vert_of_Neural_model_NEURON_by_populations.h5'):
            os.remove('Neuron_model_arrays/Vert_of_Neural_model_NEURON_by_populations.h5')
            
        number_of_points_filtered=0
        N_models=np.zeros((len(List_of_arrays)),int)
        
        for i in range(len(List_of_arrays)):
            
            if not(type(List_of_arrays[i]) is np.ndarray):
                hf3 = h5py.File('Neuron_model_arrays/Vert_of_Neural_model_NEURON_by_populations.h5', 'a')
                hf3.create_dataset(lst[i], data=0)
                hf3.close()                
            else:            
                Array_coord=List_of_arrays[i]
                
                while inx < Array_coord.shape[0]:

                    pnt=Point(Array_coord[inx,0],Array_coord[inx,1],Array_coord[inx,2])
                    
                    if (submesh_encup.bounding_box_tree().compute_first_entity_collision(pnt)<mesh.num_cells()):        #this is a condition to check whether the point is inside encap. layer or floating conductor
                        points_encap=points_encap+1
                        inx_start=int(inx/int(n_segments_fib_diam_array[i,0]))*int(n_segments_fib_diam_array[i,0])
                        Array_coord[inx_start:inx_start+int(n_segments_fib_diam_array[i,0]),:]=-100000000.0  
                        inx=inx_start+int(n_segments_fib_diam_array[i,0])
                    else:    
                        if not(mesh.bounding_box_tree().compute_first_entity_collision(pnt)<mesh.num_cells()*100):        #if one point of the neural model is absent, the whole model is disabled
                           points_outside=points_outside+1
                           inx_start=int(inx/int(n_segments_fib_diam_array[i,0]))*int(n_segments_fib_diam_array[i,0])
                           Array_coord[inx_start:inx_start+int(n_segments_fib_diam_array[i,0]),:]=-100000000.0
                           inx=inx_start+int(n_segments_fib_diam_array[i,0])
                        else:                                                                                               # finally checks whether the neuron compartment is inside the CSF voxel
                            check1_1=(voxel_array_CSF_shifted[:,0]-Array_coord[inx,0]<=MRI_param.x_vox_size)
                            check1_2=(voxel_array_CSF_shifted[:,1]-Array_coord[inx,1]<=MRI_param.y_vox_size)
                            check1_3=(voxel_array_CSF_shifted[:,2]-Array_coord[inx,2]<=MRI_param.z_vox_size)
                            check2_1=(voxel_array_CSF_shifted[:,0]>=Array_coord[inx,0])
                            check2_2=(voxel_array_CSF_shifted[:,1]>=Array_coord[inx,1])
                            check2_3=(voxel_array_CSF_shifted[:,2]>=Array_coord[inx,2])
                            #print(check1_1)
                            
                            check3=np.logical_and(np.logical_and(check1_1,check2_1),np.logical_and(np.logical_and(check1_2,check2_2),np.logical_and(check1_3,check2_3)))
                            a=np.where((check3 == (True)))                            
                            if str(a)!='(array([], dtype=int64),)':
                                points_csf=points_csf+1
                                inx_start=int(inx/int(n_segments_fib_diam_array[i,0]))*int(n_segments_fib_diam_array[i,0])
                                Array_coord[inx_start:inx_start+int(n_segments_fib_diam_array[i,0]),:]=-100000000.0
                                inx=inx_start+int(n_segments_fib_diam_array[i,0])
                            else:
                                inx=inx+1
    
                Array_coord=Array_coord[~np.all(Array_coord==-100000000.0,axis=1)]  #deletes all zero enteries          

                inx=0
                
                Array_coord=np.round(Array_coord,8)
                N_models[i]=int(Array_coord.shape[0]/int(n_segments_fib_diam_array[i,0]))
                
                hf3 = h5py.File('Neuron_model_arrays/Vert_of_Neural_model_NEURON_by_populations.h5', 'a')
                hf3.create_dataset(lst[i], data=Array_coord)
                hf3.close()
                        
                np.savetxt('Neuron_model_arrays/'+lst[i]+'.csv', Array_coord, delimiter=" ") 
                list_of_connections.append(lst[i])

                number_of_points_filtered=number_of_points_filtered+Array_coord.shape[0]
        
        hf = h5py.File('Neuron_model_arrays/Vert_of_Neural_model_NEURON_by_populations.h5', 'r')
        lst=list(hf.keys())
        result_total=[]        
        for i in lst:
            if not(i in List_of_empty):                         
                #print(i)
                a=hf.get(i)
                a=np.array(a)
                result_total.append(a)

        Array_coord_total=np.concatenate(result_total)
        hf.close()

        np.savetxt('Neuron_model_arrays/Adjusted_neuron_array_info.csv', N_models, delimiter=" ")
        np.savetxt('Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', Array_coord_total, delimiter=" ")   
        print("Points in CSF, encapsulation layer (and floating conductors) and outside (and intersecting with the electrode): ",points_csf,points_encap,points_outside)
        print("Adjusted neuron models can be visualized from Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv in Paraview")
        
        #might not work like this
        #from Visualization_files.Paraview_connections_processed import show_connections
        #show_connections(list_of_connections)
        
        from Parameter_insertion import paste_paraview_connections_vis
        paste_paraview_connections_vis(list_of_connections)
        import os
        import subprocess
        with open(os.devnull, 'w') as FNULL: subprocess.call('python Visualization_files/Paraview_connections_processed.py', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

 
        print("Number of placed neuron models per population: ",N_models)
        del Array_coord_total,voxel_array_CSF_shifted         

    if same_axon_type==False:
        del List_of_arrays
    del Array_coord

    minutes=int((time.clock() - start_neuron_models)/60)
    secnds=int(time.clock() - start_neuron_models)-minutes*60
    print("----- Adjustment of the neuron models took ",minutes," min ",secnds," s -----")    
    
    return N_models

# creates or resaves the full neuron array in the positive octant space, estimates neurons' extent for ROI size, stores meta data.
def build_neuron_models(d,MRI_param):
    
    if d["Global_rot"]==1:          # for VTA
        
        #center node of the center neuron is located at the seed (the seeding coordinate is in the positive octant space)  
        x_seeding=MRI_param.x_shift+d["x_seed"]
        y_seeding=MRI_param.y_shift+d["y_seed"]
        z_seeding=MRI_param.z_shift+d["z_seed"]
        
        #number of models, seeded in one plane. This parameter is 0 by default in Neuron_info class
        N_models_in_plane=(d["x_steps"]+1)*(d["y_steps"]+1)*(d["z_steps"]+1)    
                   
        [X_coord_old,Y_coord_old,Z_coord_old]=create_structured_array(d)    
        X_coord=[xs+x_seeding for xs in X_coord_old]
        Y_coord=[ys+y_seeding for ys in Y_coord_old]
        Z_coord=[zs+z_seeding for zs in Z_coord_old]  
    
    else:                           # manual allocation
                
        X_coord=[xs+MRI_param.x_shift for xs in d["X_coord_old"]]
        Y_coord=[ys+MRI_param.y_shift for ys in d["Y_coord_old"]]
        Z_coord=[zs+MRI_param.z_shift for zs in d["Z_coord_old"]]
        
        N_models_in_plane=0     #not need for this method

    #get morphology depending on the neuron model
    if d["Axon_Model_Type"] == 'McIntyre2002':        
        param_ax={
        'centered':True,
        'diameter':d["diam_fib"]
        }
        a=Axon(param_ax)
        nr=Axon.get_axonparams(a)    
        param_axon=[nr["ranvier_nodes"], nr["para1_nodes"], nr["para2_nodes"], nr["inter_nodes"], nr["ranvier_length"], nr["para1_length"], nr["para2_length"], nr["inter_length"], nr["deltax"], d["diam_fib"]]
    elif d["Axon_Model_Type"] == 'Reilly2016':      # lumped internodal compartments
        if d["diam_fib"]>=5.0 and d["diam_fib"]<=10.0:      # if I understand correctly, this is the only parameter that is influenced by the change of the diameter in Reilly2016 by Carnevale
            internodel_length=d["diam_fib"]*200.0     # from 1 to 2 mm
            param_axon=[d["n_Ranvier"], 0, 0, d["n_Ranvier"]-1, 0.0, 0.0, 0.0, 0.0, internodel_length, d["diam_fib"]]
        else:
            print("Wrong fiber diameter for Reilly2016, exiting")
            raise SystemExit
    else:
        print("The neuron model "+str(d["Axon_Model_Type"])+" is not implemented, exiting")
        raise SystemExit        
    
    if d["pattern_model_name"]==0:       #we can build axon model pattern
        pattern_model_name=('pattern.csv')      # name with .csv, angle in rad                    
        n_segments=generate_pattern_model(pattern_model_name,d["n_Ranvier"],param_axon,d["Axon_Model_Type"])
    else:
        pattern_model_name=d["pattern_model_name"]      # or get a pattern from the external file
        Array_coord_load_get=read_csv(pattern_model_name, delimiter=' ', header=None)
        Array_coord_load=Array_coord_load_get.values        
        n_segments=Array_coord_load.shape[0]            #all segments should be in the pattern model
        del Array_coord_load
       
    Neuron_param=Neuron_info(pattern_model_name,X_coord,Y_coord,Z_coord,d["YZ_angles"],d["ZX_angles"],d["XY_angles"],d["alpha_array_glob"],d["beta_array_glob"],d["gamma_array_glob"],N_models_in_plane)                
    with open('Neuron_model_arrays/Neuron_param_class.file', "wb") as f:
        pickle.dump(Neuron_param, f, pickle.HIGHEST_PROTOCOL)
    
    # implantation coordinated in the positive octant space
    Xt_new=MRI_param.x_shift+d["Implantation_coordinate_X"]
    Yt_new=MRI_param.y_shift+d["Implantation_coordinate_Y"]
    Zt_new=MRI_param.z_shift+d["Implantation_coordinate_Z"]
    
    "definition of ROI is required for adaptive mesh refinement. ROI is a sphere encompassing all neuron models"
    ROI_radius=get_neuron_models_dims(d["Neuron_model_array_prepared"],Xt_new,Yt_new,Zt_new,Neuron_param,param_axon,d["Global_rot"])
    
    if d["Axon_Model_Type"] == 'McIntyre2002': 
        Neuron_model_misc=np.array([nr["ranvier_nodes"], nr["para1_nodes"], nr["para2_nodes"], nr["inter_nodes"], nr["ranvier_length"], nr["para1_length"], nr["para2_length"], nr["inter_length"], nr["deltax"],d["diam_fib"],d["n_Ranvier"],ROI_radius,n_segments])
    elif d["Axon_Model_Type"] == 'Reilly2016':
        Neuron_model_misc=np.array([d["n_Ranvier"], 0, 0, d["n_Ranvier"]-1, 0.0, 0.0, 0.0, 0.0, internodel_length,d["diam_fib"],d["n_Ranvier"],ROI_radius,n_segments])
    
    '''Save meta data for the future simulations with the current Neuron model data set'''
    np.savetxt('Neuron_model_arrays/Neuron_model_misc.csv', Neuron_model_misc, delimiter=" ")
    
    print("Initial neuron models and corresponding meta data were created\n")

    return Neuron_param,param_axon,ROI_radius,n_segments

def create_meta_data_for_predefined_models(d,MRI_param):
    
    # shift to the positive octant space
    Xt_new=MRI_param.x_shift+d["Implantation_coordinate_X"]
    Yt_new=MRI_param.y_shift+d["Implantation_coordinate_Y"]
    Zt_new=MRI_param.z_shift+d["Implantation_coordinate_Z"]
    
    if d["Name_prepared_neuron_array"][-4:]=='.csv':
        if d["Axon_Model_Type"] == 'McIntyre2002': 
            from Axon_files.axon import Axon
            param_ax={
            'centered':True,
            'diameter':d["diam_fib"]
            }
            a=Axon(param_ax)
            nr=Axon.get_axonparams(a)
            
            param_axon=[nr["ranvier_nodes"], nr["para1_nodes"], nr["para2_nodes"], nr["inter_nodes"], nr["ranvier_length"], nr["para1_length"], nr["para2_length"], nr["inter_length"], nr["deltax"], d["diam_fib"]]
            n_comp=((nr["ranvier_nodes"]-1)+nr["inter_nodes"]+nr["para1_nodes"]+nr["para2_nodes"])/(nr["ranvier_nodes"]-1)
            n_segments=int((d["n_Ranvier"]-1)*n_comp+1)        #overall number of points on Axon incl. internodal
        elif d["Axon_Model_Type"] == 'Reilly2016':
            if d["diam_fib"]>=5.0 and d["diam_fib"]<=10.0:      # if I understand correctly, this is the only parameter that is influenced by the change of the diameter in Reilly2016 by Carnevale
                internodel_length=d["diam_fib"]*200.0     # from 1 to 2 mm
                param_axon=[d["n_Ranvier"], 0, 0, d["n_Ranvier"]-1, 0.0, 0.0, 0.0, 0.0, internodel_length, d["diam_fib"]]        
            else:
                print("Wrong fiber diameter for Reilly2016, exiting")
                raise SystemExit
        else:
            print("The neuron model "+str(d["Axon_Model_Type"])+" is not implemented, exiting")
            raise SystemExit    
        
        #from Neuron_models_arangement import get_neuron_models_dims
        ROI_radius=get_neuron_models_dims(d["Neuron_model_array_prepared"],Xt_new,Yt_new,Zt_new,0,param_axon,0)   #here we do not need Neuron_param class
        if d["Axon_Model_Type"] == 'McIntyre2002':
            Neuron_model_misc=np.array([nr["ranvier_nodes"], nr["para1_nodes"], nr["para2_nodes"], nr["inter_nodes"], nr["ranvier_length"], nr["para1_length"], nr["para2_length"], nr["inter_length"], nr["deltax"],d["diam_fib"],d["n_Ranvier"],ROI_radius,n_segments])
        elif d["Axon_Model_Type"] == 'Reilly2016':
            Neuron_model_misc=np.array([d["n_Ranvier"], 0, 0, d["n_Ranvier"]-1, 0.0, 0.0, 0.0, 0.0, internodel_length,d["diam_fib"],d["n_Ranvier"],ROI_radius,n_segments])
        '''Save meta data for the future simulations with the current Neuron model data set'''
        np.savetxt('Neuron_model_arrays/Neuron_model_misc.csv', Neuron_model_misc, delimiter=" ")
    
        print("Meta data for predefined neuron models was created\n")
        
    if d["Name_prepared_neuron_array"][-3:]=='.h5':
        
        param_axon=[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1]          #not needed here

        n_segments_fib_diam_array=np.load('Neuron_model_arrays/Neuron_populations_misc.npy')
        ROI_radius=get_neuron_models_dims(d["Neuron_model_array_prepared"],Xt_new,Yt_new,Zt_new,0,0,0)
        
        Neuron_model_misc=np.array([-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ROI_radius,-1])          #not needed here
        np.savetxt('Neuron_model_arrays/Neuron_model_misc.csv', Neuron_model_misc, delimiter=" ")
        
        n_segments=n_segments_fib_diam_array[:,0]
        n_segments=n_segments.astype(int)
        
    return ROI_radius,param_axon,n_segments

def adjust_neuron_models(d,MRI_param,Domains,Neuron_param,param_axon):
    
    #from Neuron_models_arangement import generate_neuron_models
    print("----- Adjusting neuron models to the geometry and MRI data -----")
    N_models=generate_neuron_models(d["n_Ranvier"],d["Neuron_model_array_prepared"],Domains,MRI_param,Neuron_param,param_axon,d["Global_rot"])         #creates file with vertices and also returns numbver of nodes in one neural model      
    
    print("Neuron models were adjusted and corresp. meta data were created\n")
    
    return N_models

# if brain approximation is used, the function will drop the models outside of the computational domain
# otherwise will just store the meta data for .h5 and resave the neuron array as Prepared_models_full_filtered
# the reason to cut for the brain approximation is possible wrong choice of the dimensions that leads to problems with ROI
# if an imported brain geometry is use, we assume the users made sure that it encompases their neuron array (non-physiologically placed neurons will be removed later!)
def cut_models_by_domain(d,Brain_shape_name,name_prepared_neuron_array):
    
    if Brain_shape_name=='Brain_substitute.brep':       # load the brain approximation
        mesh_brain_sub = Mesh("Meshes/Mesh_brain_substitute_max_ROI.xml")

    if name_prepared_neuron_array[-4:]=='.csv':      #if we have an array where all axons have the same length
        if d["Axon_Model_Type"] == 'McIntyre2002':
            param_ax={
            'centered':True,
            'diameter':d["diam_fib"]
            }
            a=Axon(param_ax)
            nr=Axon.get_axonparams(a)
            n_comp=((nr["ranvier_nodes"]-1)+nr["inter_nodes"]+nr["para1_nodes"]+nr["para2_nodes"])/(nr["ranvier_nodes"]-1)
            n_segments=int((d["n_Ranvier"]-1)*n_comp+1)        #overall number of points on Axon incl. internodal
        elif d["Axon_Model_Type"] == 'Reilly2016':
            n_comp=2        #only nodes and one internodal per segment
            n_segments=int((d["n_Ranvier"]-1)*n_comp+1)
        else:
            print("The neuron model "+str(d["Axon_Model_Type"])+" is not implemented, exiting")
            raise SystemExit  
      
        Array_coord_get=read_csv(name_prepared_neuron_array, delimiter=' ', header=None)        #
        Array_coord_in_MRI=Array_coord_get.values
        
        if Brain_shape_name=='Brain_substitute.brep':
            n_models_before=int(Array_coord_in_MRI.shape[0]/n_segments)
            print("Initial number of neuron models: ",n_models_before)                        
            points_outside=0
            
            for inx in range(Array_coord_in_MRI.shape[0]):                                        
                pnt=Point(Array_coord_in_MRI[inx,0],Array_coord_in_MRI[inx,1],Array_coord_in_MRI[inx,2])            
                if not(mesh_brain_sub.bounding_box_tree().compute_first_entity_collision(pnt)<mesh_brain_sub.num_cells()*100):        #if one point of the neural model is absent, the whole model is subtracted
                    points_outside=points_outside+1
                    inx_start=int(inx/n_segments)*n_segments
                    Array_coord_in_MRI[inx_start:inx_start+n_segments,:]=-100000000.0
                   
            Array_coord_in_MRI=Array_coord_in_MRI[~np.all(Array_coord_in_MRI==-100000000.0,axis=1)]  #deletes all zero enteries          
            n_models_after=int(Array_coord_in_MRI.shape[0]/n_segments)
            print(n_models_before-n_models_after, " models were outside of the approximating geometrical domain")            
            Array_coord_in_MRI=np.round(Array_coord_in_MRI,8)
        
        np.savetxt('Neuron_model_arrays/Prepared_models_full_filtered.csv', Array_coord_in_MRI, delimiter=" ")
        del Array_coord_in_MRI
        return True
    
    if name_prepared_neuron_array[-3:]=='.h5':  

        if isinstance(d["diam_fib"],list):
            diam_of_populations=d["diam_fib"]
        else:
            diam_of_populations=[d["diam_fib"]]
            d["n_Ranvier"]=[d["n_Ranvier"]]

        hf = h5py.File(name_prepared_neuron_array, 'r')
        lst=list(hf.keys())

        population_index=0
        n_segments_fib_diam_array=np.zeros((len(lst),2),float)        #colums are n_segments and fib_diam
        
        for i in lst:
            Array_coord_in_MRI=hf.get(i)
            Array_coord_in_MRI=np.array(Array_coord_in_MRI)
            
            if d["Axon_Model_Type"] == 'McIntyre2002':
                 param_ax={
                 'centered':True,
                 'diameter':diam_of_populations[population_index]
                 }
                 a=Axon(param_ax)
                 nr=Axon.get_axonparams(a)
                 n_comp=((nr["ranvier_nodes"]-1)+nr["inter_nodes"]+nr["para1_nodes"]+nr["para2_nodes"])/(nr["ranvier_nodes"]-1)
                 n_segments=int((d["n_Ranvier"][population_index]-1)*n_comp+1)        #overall number of points on Axon incl. internodal
            elif d["Axon_Model_Type"] == 'Reilly2016':
                n_comp=2        #only nodes and one internodal per segment
                n_segments=int((d["n_Ranvier"][population_index]-1)*n_comp+1)
            else:
                print("The neuron model "+str(d["Axon_Model_Type"])+" is not implemented, exiting")
                raise SystemExit
                       
            n_segments_fib_diam_array[population_index,0]=n_segments
            n_segments_fib_diam_array[population_index,1]=diam_of_populations[population_index]
            
            if Brain_shape_name=='Brain_substitute.brep':   # cut models if use brain approximation                
                n_models_before=int(Array_coord_in_MRI.shape[0]/n_segments)
                print("Initial number of neuron models in",str(i),": ",n_models_before)                                
                points_outside=0
                
                for inx in range(Array_coord_in_MRI.shape[0]):                                            
                    pnt=Point(Array_coord_in_MRI[inx,0],Array_coord_in_MRI[inx,1],Array_coord_in_MRI[inx,2])                
                    if not(mesh_brain_sub.bounding_box_tree().compute_first_entity_collision(pnt)<mesh_brain_sub.num_cells()*100):        #if one point of the neural model is absent, the whole model is disabled
                        points_outside=points_outside+1
                        inx_start=int(inx/n_segments)*n_segments
                        Array_coord_in_MRI[inx_start:inx_start+n_segments,:]=-100000000.0
                       
                Array_coord_in_MRI=Array_coord_in_MRI[~np.all(Array_coord_in_MRI==-100000000.0,axis=1)]  #deletes all -10000000 enteries          
                n_models_after=int(Array_coord_in_MRI.shape[0]/n_segments)
                if int(n_models_before-n_models_after)!=0:
                    print("In ",str(i), n_models_before-n_models_after, " models were outside of the approximating geometrical domain")
                
                Array_coord_in_MRI=np.round(Array_coord_in_MRI,8)

            hf2 = h5py.File('Neuron_model_arrays/Prepared_models_full_filtered.h5', 'a')
            hf2.create_dataset(i, data=Array_coord_in_MRI)
            hf2.close()

            population_index=population_index+1

        hf.close()
        
        '''Here the meta data is different as we have different structure of the input array'''
        np.save('Neuron_model_arrays/Neuron_populations_misc', n_segments_fib_diam_array)

        return True
    
