from dolfin import *
from pandas import read_csv
#from MRI_DTI_prep import map_MRI
#from MRI_DTI_prep import map_DTI
from tissue_dielectrics import DielectricProperties
import numpy as np
import time
import os
set_log_active(False)   #turns off debugging info
'''this script will create CellFunction on the provided mesh and mark tissues with different indices'''
'''this is not a script just to plug into any simulation and run'''
'''it rather gives an example of how MRI (DTI) can be prepared for FEniCS simulations (the implementation is adapted from "Tensor-weighted Poisson")'''



def get_cellmap(mesh,subdomains_assigned,Domains,MRI_param,default_material):
    
    #Tissue_array_get=read_csv('MRI_DTI_derived_data/Tissue_array_MRI.csv', delimiter=' ', header=None)
    #Tissue_array=Tissue_array_get.values
    Tissue_array=np.load('MRI_DTI_derived_data/Tissue_array_MRI.npy')
    
    x_vect=np.genfromtxt('MRI_DTI_derived_data/x_vector_MRI_Box.csv', delimiter=' ')
    y_vect=np.genfromtxt('MRI_DTI_derived_data/y_vector_MRI_Box.csv', delimiter=' ')
    z_vect=np.genfromtxt('MRI_DTI_derived_data/z_vector_MRI_Box.csv', delimiter=' ')
    
    #[Mx_mri,My_mri,Mz_mri]=map_MRI()         #don't need to run know, the voxel_array already formed
    
    #voxel_array_get=read_csv('voxel_array_ROI.csv', delimiter=' ', header=None)
    #voxel_array_get=read_csv('MRI_DTI_files/voxel_array_'+MRI_param.name[:-4]+'.csv', delimiter=' ', header=None)
    #voxel_array=voxel_array_get.values
    #print voxel_array[0:3,:]
    '''structure of voxel array'''
    '''rows: x,y,z,material index'''
    '''first voxel has to have coordinates: 0+voxelsize_x,0+voxelsize_y,0+voxelsize_z'''
    '''Indeed, such a structure of the array is not optimal, since x,y,z coordinates can be simply stored in coordinate vectors'''
    '''Any suggestions (incl. minimal working example) will be appreciated'''
    
    voxel_size_x=MRI_param.x_vox_size
    voxel_size_y=MRI_param.y_vox_size
    voxel_size_z=MRI_param.z_vox_size

    #print MRI_param.M_x,MRI_param.M_y,MRI_param.M_z

    Mx_mri=MRI_param.M_x
    My_mri=MRI_param.M_y
    Mz_mri=MRI_param.M_z

      
    subdomains = MeshFunction('size_t',mesh,3)
    subdomains.set_all(0)
    #subdomains = MeshFunctionSizet(mesh, 3, 0)
    
    for cell in cells(mesh):
        
        x_coord=cell.midpoint().x()
    
        y_coord=cell.midpoint().y()
        z_coord=cell.midpoint().z()
        
        
        xv_mri=int((x_coord)/voxel_size_x-0.000000001)                                  #defines number of steps to get to the voxels containing x[0] coordinate
        yv_mri=(int((y_coord)/voxel_size_y-0.000000001))*Mx_mri                  #defines number of steps to get to the voxels containing x[0] and x[1] coordinates
        zv_mri=(int((z_coord)/voxel_size_z-0.000000001))*Mx_mri*My_mri           #defines number of steps to get to the voxels containing x[0], x[1] and x[2] coordinates
        k_mri=xv_mri+yv_mri+zv_mri
	#print k_mri
        k_mri=int(k_mri)
        
        x_vect_index=(int((x_coord)/voxel_size_x-0.000000001))
        y_vect_index=(int((y_coord)/voxel_size_y-0.000000001))
        z_vect_index=(int((z_coord)/voxel_size_z-0.000000001))
        
        if x_vect_index>=Mx_mri or y_vect_index>=My_mri or z_vect_index>=Mz_mri or x_coord<0.0 or y_coord<0.0 or z_coord<0.0:
            x_vect_index,y_vect_index,z_vect_index=0,0,0        #cell is outside MRI data, assign first voxel of MRI, which will not fullfil the following condition, and the cell will be remain marked with 0 (default)
        
        if k_mri>=Mx_mri*My_mri*Mz_mri or k_mri<0:
            k_mri=0         #cell is outside MRI data, assign first voxel of MRI, which will not fullfil the following condition, and the cell will be remain marked with 0 (default)
        

        '''First, find the corresponding voxel in MRI, then check the value of the voxel and assign a corresponding conductivity'''    
        if (round(y_vect[y_vect_index]-y_coord,8)<=voxel_size_y and round(y_vect[y_vect_index]-y_coord,8)>=0 and (round(x_vect[x_vect_index]-x_coord,8)<=voxel_size_x and round(x_vect[x_vect_index]-x_coord,8)>=0) and round(z_vect[z_vect_index]-z_coord)<=voxel_size_z and round(z_vect[z_vect_index]-z_coord,8)>=0):

            if int(Tissue_array[k_mri])==3:
                #cell.mark(subdomains, 1)
                subdomains[cell]=3
                
            if int(Tissue_array[k_mri])==2:
                subdomains[cell]=2
                
            if int(Tissue_array[k_mri])==1:
                subdomains[cell]=1
                
            if int(Tissue_array[k_mri])==0:
                subdomains[cell]=default_material
#old way 
#    if Domains.Float_contacts!=-1:
#        subdomains.array()[subdomains_assigned.array()==Domains.Float_contacts]=5          #5 is index for float contact (very high cond and perm)            

    if Domains.Float_contacts!=-1:
        if isinstance(Domains.Float_contacts,int):
            subdomains.array()[subdomains_assigned.array()==Domains.Float_contacts]=5
        else:
            for i in range(len(Domains.Float_contacts)):
                subdomains.array()[subdomains_assigned.array()==Domains.Float_contacts[i]]=5          #5 is index for float contact (very high cond and perm)            
  
         
    for i in range(len(Domains.Encup_index)):
        subdomains.array()[subdomains_assigned.array()==Domains.Encup_index[i]]=4          #4 is index of encap


    
#    boundaries = MeshFunction('size_t',mesh,'Meshes/Mesh_unref_facet_region.xml')
#
#    for cell in cells(mesh):
#        for f in facets(cell):
#            #print(boundaries[f]!=0)
#            if boundaries[f]!=0:
#                subdomains[cell]=4          #cells on the interface should be marked with encap, otherwise no current on the floating contact       
#                for v in vertices(f):       #also all cells which vertices are on the active boundaries should be marked
#                    for c in cells(v):
#                        subdomains[c]=4
#                #cell_markers[cell] = True


    
    
    
    #file=File('subdomains.pvd')
    #file<<subdomains
            
    del Tissue_array
    
    return (subdomains)          #cond, perm are cell functions and can be directly used 


def get_cellmap_tensors(mesh,subdomains_assigned,Domains,MRI_param,DTI_param,default_material):
    
    #Tissue_array_get=read_csv('MRI_DTI_derived_data/Tissue_array_MRI.csv', delimiter=' ', header=None)
    #Tissue_array=Tissue_array_get.values
    Tissue_array=np.load('MRI_DTI_derived_data/Tissue_array_MRI.npy')
    
    x_vect=np.genfromtxt('MRI_DTI_derived_data/x_vector_MRI_Box.csv', delimiter=' ')
    y_vect=np.genfromtxt('MRI_DTI_derived_data/y_vector_MRI_Box.csv', delimiter=' ')
    z_vect=np.genfromtxt('MRI_DTI_derived_data/z_vector_MRI_Box.csv', delimiter=' ')
    
    #DTI_array_get=read_csv('MRI_DTI_derived_data/Tensor_array_DTI.csv', delimiter=' ', header=None)
    #DTI_array=DTI_array_get.values
    DTI_array=np.load('MRI_DTI_derived_data/Tensor_array_DTI.npy')
    
    x_vect_DTI=np.genfromtxt('MRI_DTI_derived_data/x_vector_DTI_Box.csv', delimiter=' ')
    y_vect_DTI=np.genfromtxt('MRI_DTI_derived_data/y_vector_DTI_Box.csv', delimiter=' ')
    z_vect_DTI=np.genfromtxt('MRI_DTI_derived_data/z_vector_DTI_Box.csv', delimiter=' ')

    #[Mx_dti,My_dti,Mz_dti]=map_DTI()
    
    #voxel_array_get=read_csv('voxel_array_ROI.csv', delimiter=' ', header=None)

    #print voxel_array[0:3,:]
    '''structure of voxel array'''
    '''rows: x,y,z,material index'''
    '''first voxel has to have coordinates: 0+voxelsize_x,0+voxelsize_y,0+voxelsize_z'''
    '''Indeed, such a structure of the array is not optimal, since x,y,z coordinates can be simply stored in coordinate vectors'''
    '''Any suggestions (incl. minimal working example) will be appreciated'''
    
    '''ROUND DTI DATA!'''
        
    #[aa,ab,ac,ad,ae,af,ag,x_start_DTI,y_start_DTI,z_start_DTI]=np.genfromtxt('MRI_DTI_files/DTI_misc.csv', delimiter=' ')    
    

    voxel_size_x=MRI_param.x_vox_size
    voxel_size_y=MRI_param.y_vox_size
    voxel_size_z=MRI_param.z_vox_size
    
    Mx_mri=MRI_param.M_x
    My_mri=MRI_param.M_y
    Mz_mri=MRI_param.M_z

    voxel_size_x_DTI=DTI_param.x_vox_size
    voxel_size_y_DTI=DTI_param.y_vox_size
    voxel_size_z_DTI=DTI_param.z_vox_size

    Mx_dti=DTI_param.M_x
    My_dti=DTI_param.M_y
    Mz_dti=DTI_param.M_z
      
    
    x_vect_DTI[:]=x_vect_DTI[:]+DTI_param.x_start
    y_vect_DTI[:]=y_vect_DTI[:]+DTI_param.y_start
    z_vect_DTI[:]=z_vect_DTI[:]+DTI_param.z_start
        
   # subdomains =MeshFunctionSizet(mesh, 3, 0)
    subdomains = MeshFunction('size_t',mesh,3)
    subdomains.set_all(0)
    
    #Create mesh functions for c00, c01, c11
    c00 = MeshFunction("double", mesh, 3, 1.0)
    c01 = MeshFunction("double", mesh, 3, 0.0)
    c02 = MeshFunction("double", mesh, 3, 0.0)
    c11 = MeshFunction("double", mesh, 3, 1.0)
    c12 = MeshFunction("double", mesh, 3, 0.0)
    c22 = MeshFunction("double", mesh, 3, 1.0)    
    
    
    for cell in cells(mesh):
        
        x_coord=cell.midpoint().x()
    
        y_coord=cell.midpoint().y()
        z_coord=cell.midpoint().z()
        
        
        xv_mri=int(x_coord/voxel_size_x-0.000000001)                                  #defines number of steps to get to the voxels containing x[0] coordinate
        yv_mri=xv_mri+(int(y_coord/voxel_size_y-0.000000001))*Mx_mri                  #defines number of steps to get to the voxels containing x[0] and x[1] coordinates
        zv_mri=yv_mri+(int(z_coord/voxel_size_z-0.000000001))*Mx_mri*My_mri           #defines number of steps to get to the voxels containing x[0], x[1] and x[2] coordinates
        k_mri=zv_mri
        
        k_mri=int(k_mri)
        
        if k_mri>=Mx_mri*My_mri*Mz_mri or k_mri<0:
            k_mri=0         #cell is outside MRI data, assign first voxel of MRI, which will not fullfil the following condition, and the cell will be remain marked with 0 (default)
  
        x_vect_index=(int((x_coord)/voxel_size_x-0.000000001))
        y_vect_index=(int((y_coord)/voxel_size_y-0.000000001))
        z_vect_index=(int((z_coord)/voxel_size_z-0.000000001))
        
        if x_vect_index>=Mx_mri or y_vect_index>=My_mri or z_vect_index>=Mz_mri or x_coord<0.0 or y_coord<0.0 or z_coord<0.0:
            x_vect_index,y_vect_index,z_vect_index=0,0,0        #cell is outside MRI data, assign first voxel of MRI, which will not fullfil the following condition, and the cell will be remain marked with 0 (default)


               
        xv_dti=int((x_coord-DTI_param.x_start)/voxel_size_x_DTI-0.00000001)
        yv_dti=xv_dti+(int((y_coord-DTI_param.y_start)/voxel_size_y_DTI-0.00000001))*Mx_dti
        zv_dti=yv_dti+(int((z_coord-DTI_param.z_start)/voxel_size_z_DTI-0.00000001))*Mx_dti*My_dti
        k_dti=zv_dti
        
        k_dti=int(k_dti)        
        
        if k_dti>=Mx_dti*My_dti*Mz_dti or k_dti<0:
            k_dti=0         #cell is outside DTI data, assign first voxel of DTI, which will not fullfil the following condition, and the cell will be considered as isotropic
        ''' otherwise comment it out'''
        
        x_vect_DTI_index=(int((x_coord-DTI_param.x_start)/voxel_size_x_DTI-0.000000001))
        y_vect_DTI_index=(int((y_coord-DTI_param.y_start)/voxel_size_y_DTI-0.000000001))
        z_vect_DTI_index=(int((z_coord-DTI_param.z_start)/voxel_size_z_DTI-0.000000001))
        
        if x_vect_DTI_index>=Mx_dti or y_vect_DTI_index>=My_dti or z_vect_DTI_index>=Mz_dti or x_coord<0.0 or y_coord<0.0 or z_coord<0.0 or z_vect_DTI_index<0 or y_vect_DTI_index<0 or x_vect_DTI_index<0:
            x_vect_DTI_index,y_vect_DTI_index,z_vect_DTI_index=0,0,0        #cell is outside MRI data, assign first voxel of MRI, which will not fullfil the following condition, and the cell will be remain marked with 0 (default)

        '''First, find it in DTI'''    
        if round(y_vect_DTI[y_vect_DTI_index]-y_coord,8)<=voxel_size_y_DTI and round(y_vect_DTI[y_vect_DTI_index]-y_coord,8)>=0 and round(x_vect_DTI[x_vect_DTI_index]-x_coord,8)<=voxel_size_x_DTI and round(x_vect_DTI[x_vect_DTI_index]-x_coord,8)>=0 and (round(z_vect_DTI[z_vect_DTI_index]-z_coord,8)<=voxel_size_z_DTI and round(z_vect_DTI[z_vect_DTI_index]-z_coord,8)>=0):
            if int(Tissue_array[k_mri])==3.0 and (round(y_vect[y_vect_index]-y_coord,8)<=voxel_size_y and round(y_vect[y_vect_index]-y_coord,8)>=0 and (round(x_vect[x_vect_index]-x_coord,8)<=voxel_size_x and round(x_vect[x_vect_index]-x_coord,8)>=0) and round(z_vect[z_vect_index]-z_coord)<=voxel_size_z and round(z_vect[z_vect_index]-z_coord,8)>=0):
                subdomains[cell]=3
            if int(Tissue_array[k_mri])==2.0 and (round(y_vect[y_vect_index]-y_coord,8)<=voxel_size_y and round(y_vect[y_vect_index]-y_coord,8)>=0 and (round(x_vect[x_vect_index]-x_coord,8)<=voxel_size_x and round(x_vect[x_vect_index]-x_coord,8)>=0) and round(z_vect[z_vect_index]-z_coord)<=voxel_size_z and round(z_vect[z_vect_index]-z_coord,8)>=0):
                subdomains[cell]=2
            if int(Tissue_array[k_mri])==1.0 and (round(y_vect[y_vect_index]-y_coord,8)<=voxel_size_y and round(y_vect[y_vect_index]-y_coord,8)>=0 and (round(x_vect[x_vect_index]-x_coord,8)<=voxel_size_x and round(x_vect[x_vect_index]-x_coord,8)>=0) and round(z_vect[z_vect_index]-z_coord)<=voxel_size_z and round(z_vect[z_vect_index]-z_coord,8)>=0):
                subdomains[cell]=1
            if int(Tissue_array[k_mri])==0.0 and (round(y_vect[y_vect_index]-y_coord,8)<=voxel_size_y and round(y_vect[y_vect_index]-y_coord,8)>=0 and (round(x_vect[x_vect_index]-x_coord,8)<=voxel_size_x and round(x_vect[x_vect_index]-x_coord,8)>=0) and round(z_vect[z_vect_index]-z_coord)<=voxel_size_z and round(z_vect[z_vect_index]-z_coord,8)>=0):
                #print "Unknown material, assigning GM marker"
                subdomains[cell]=int(default_material)
                
            
            c00[cell] = DTI_array[k_dti,0]      #we don't scale tensor with cond. values, because they depend on frequency, will do it directly in EQS_function2
            c01[cell] = DTI_array[k_dti,1]
            c02[cell] = DTI_array[k_dti,2]
            c11[cell] = DTI_array[k_dti,3]
            c12[cell] = DTI_array[k_dti,4]
            c22[cell] = DTI_array[k_dti,5]            
            
            
            '''if MRI data set is not for the full computational domain, the following condition will be triggered (cell does not correspond to MRI voxel)'''
            if not (round(y_vect[y_vect_index]-y_coord,8)<=voxel_size_y and round(y_vect[y_vect_index]-y_coord,8)>=0 and (round(x_vect[x_vect_index]-x_coord,8)<=voxel_size_x and round(x_vect[x_vect_index]-x_coord,8)>=0) and round(z_vect[z_vect_index]-z_coord)<=voxel_size_z and round(z_vect[z_vect_index]-z_coord,8)>=0):
                subdomains[cell]=int(default_material)

              
        
        '''If DTI voxel wasn't found (we might have DTI data only for a part of MRI data)'''
        #else:
        if not(round(y_vect_DTI[y_vect_DTI_index]-y_coord,8)<=voxel_size_y_DTI and round(y_vect_DTI[y_vect_DTI_index]-y_coord,8)>=0 and round(x_vect_DTI[x_vect_DTI_index]-x_coord,8)<=voxel_size_x_DTI and round(x_vect_DTI[x_vect_DTI_index]-x_coord,8)>=0 and (round(z_vect_DTI[z_vect_DTI_index]-z_coord,8)<=voxel_size_z_DTI and round(z_vect_DTI[z_vect_DTI_index]-z_coord,8)>=0)):
            if int(Tissue_array[k_mri])==3.0 and (round(y_vect[y_vect_index]-y_coord,8)<=voxel_size_y and round(y_vect[y_vect_index]-y_coord,8)>=0 and (round(x_vect[x_vect_index]-x_coord,8)<=voxel_size_x and round(x_vect[x_vect_index]-x_coord,8)>=0) and round(z_vect[z_vect_index]-z_coord)<=voxel_size_z and round(z_vect[z_vect_index]-z_coord,8)>=0):
                subdomains[cell]=3
            if int(Tissue_array[k_mri])==2.0 and (round(y_vect[y_vect_index]-y_coord,8)<=voxel_size_y and round(y_vect[y_vect_index]-y_coord,8)>=0 and (round(x_vect[x_vect_index]-x_coord,8)<=voxel_size_x and round(x_vect[x_vect_index]-x_coord,8)>=0) and round(z_vect[z_vect_index]-z_coord)<=voxel_size_z and round(z_vect[z_vect_index]-z_coord,8)>=0):
                subdomains[cell]=2
            if int(Tissue_array[k_mri])==1.0 and (round(y_vect[y_vect_index]-y_coord,8)<=voxel_size_y and round(y_vect[y_vect_index]-y_coord,8)>=0 and (round(x_vect[x_vect_index]-x_coord,8)<=voxel_size_x and round(x_vect[x_vect_index]-x_coord,8)>=0) and round(z_vect[z_vect_index]-z_coord)<=voxel_size_z and round(z_vect[z_vect_index]-z_coord,8)>=0):
                subdomains[cell]=1
            if int(Tissue_array[k_mri])==0.0 and (round(y_vect[y_vect_index]-y_coord,8)<=voxel_size_y and round(y_vect[y_vect_index]-y_coord,8)>=0 and (round(x_vect[x_vect_index]-x_coord,8)<=voxel_size_x and round(x_vect[x_vect_index]-x_coord,8)>=0) and round(z_vect[z_vect_index]-z_coord)<=voxel_size_z and round(z_vect[z_vect_index]-z_coord,8)>=0):
                subdomains[cell]=int(default_material)
                
            c00[cell] = 1.0      #no anisotropy defined for this voxel
            c01[cell] = 0.0
            c02[cell] = 0.0
            c11[cell] = 1.0
            c12[cell] = 0.0
            c22[cell] = 1.0
            
            
            
            #'''if MRI data set is not for the full computational domain, the following condition will be triggered (cell does not correspond to MRI voxel)'''
            #if not (round(voxel_array[k_mri,1]-y_coord,8)<=voxelsize and round(voxel_array[k_mri,1]-y_coord,8)>=0) and (round(voxel_array[k_mri,0]-x_coord,8)<=voxelsize and round(voxel_array[k_mri,0]-x_coord,8)>=0) and (round(voxel_array[k_mri,2]-z_coord,8)<=voxelsize and round(voxel_array[k_mri,2]-z_coord,8)>=0):
            #    subdomains[cell]=int(default_material)
                
    
    for i in range(len(Domains.Encup_index)):
        subdomains.array()[subdomains_assigned.array()==Domains.Encup_index[i]]=4          #mark every cell with value 3 from xml to index

        c00.array()[subdomains_assigned.array()==Domains.Encup_index[i]]=1.0          #no anisotropy in encap, index 3 is from xml
        c01.array()[subdomains_assigned.array()==Domains.Encup_index[i]]=0.0          #no anisotropy in encap 
        c02.array()[subdomains_assigned.array()==Domains.Encup_index[i]]=0.0          #no anisotropy in encap 
        c11.array()[subdomains_assigned.array()==Domains.Encup_index[i]]=1.0          #no anisotropy in encap 
        c12.array()[subdomains_assigned.array()==Domains.Encup_index[i]]=0.0          #no anisotropy in encap 
        c22.array()[subdomains_assigned.array()==Domains.Encup_index[i]]=1.0          #no anisotropy in encap   

    if Domains.Float_contacts!=-1:
        if isinstance(Domains.Float_contacts,int):
            subdomains.array()[subdomains_assigned.array()==Domains.Float_contacts]=5          #4 is index of encap
            c00.array()[subdomains_assigned.array()==Domains.Float_contacts]=1.0          #no anisotropy in encap, index 3 is from xml
            c01.array()[subdomains_assigned.array()==Domains.Float_contacts]=0.0          #no anisotropy in encap 
            c02.array()[subdomains_assigned.array()==Domains.Float_contacts]=0.0          #no anisotropy in encap 
            c11.array()[subdomains_assigned.array()==Domains.Float_contacts]=1.0          #no anisotropy in encap 
            c12.array()[subdomains_assigned.array()==Domains.Float_contacts]=0.0          #no anisotropy in encap 
            c22.array()[subdomains_assigned.array()==Domains.Float_contacts]=1.0          #no anisotropy in encap 
        else:
            for i in range(len(Domains.Float_contacts)):
                subdomains.array()[subdomains_assigned.array()==Domains.Float_contacts[i]]=5
                c00.array()[subdomains_assigned.array()==Domains.Float_contacts[i]]=1.0          #no anisotropy in encap, index 3 is from xml
                c01.array()[subdomains_assigned.array()==Domains.Float_contacts[i]]=0.0          #no anisotropy in encap 
                c02.array()[subdomains_assigned.array()==Domains.Float_contacts[i]]=0.0          #no anisotropy in encap 
                c11.array()[subdomains_assigned.array()==Domains.Float_contacts[i]]=1.0          #no anisotropy in encap 
                c12.array()[subdomains_assigned.array()==Domains.Float_contacts[i]]=0.0          #no anisotropy in encap 
                c22.array()[subdomains_assigned.array()==Domains.Float_contacts[i]]=1.0          #no anisotropy in encap 


    #subdomains.array()[subdomains_assigned.array()==3]=4          #mark every cell with value 3 from xml to index encap as 4

    #here we also keep a reference to the number of elements as the mesh identifier
    hdf = HDF5File(mesh.mpi_comm(), 'Results_adaptive/Tensors_to_solve_num_el_'+str(mesh.num_cells())+'.h5', 'w')
    hdf.write(c00, "/c00")
    hdf.write(c01, "/c01")
    hdf.write(c02, "/c02")
    hdf.write(c11, "/c11")
    hdf.write(c12, "/c12")
    hdf.write(c22, "/c22")
    hdf.close()
    
    file=File('Tensors/c00_unscaled.pvd')
    file<<c00,mesh      
    file=File('Tensors/c01_unscaled.pvd')
    file<<c01,mesh         
    file=File('Tensors/c02_unscaled.pvd')
    file<<c02,mesh
    file=File('Tensors/c11_unscaled.pvd')
    file<<c11,mesh
    file=File('Tensors/c12_unscaled.pvd')
    file<<c12,mesh
    file=File('Tensors/c22_unscaled.pvd')
    file<<c22,mesh
        
    #Store to file
    #mesh_file = File("mesh.xml.gz")
    #c00_file = File("Tensors/c00.xml.gz")
    #c01_file = File("Tensors/c01.xml.gz")
    #c02_file = File("Tensors/c02.xml.gz")
    #c11_file = File("Tensors/c11.xml.gz")
    #c12_file = File("Tensors/c12.xml.gz")
    #c22_file = File("Tensors/c22.xml.gz")
    
    ##mesh_file << mesh
    #c00_file << c00
    #c01_file << c01
    #c02_file << c02
    #c11_file << c11
    #c12_file << c12
    #c22_file << c22
  

            
    return (subdomains)          #cond, perm are cell functions and can be directly used 
#interactive()
