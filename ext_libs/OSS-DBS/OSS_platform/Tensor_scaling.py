
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 12:11:41 2020
@author: scaling algorithms by A.Andree, parallelization by K.Butenko
"""

import os
import nibabel as nib
import matplotlib.pyplot as plt

from multiprocessing import sharedctypes,cpu_count,Pool
from functools import partial

import numpy as np
import itertools

import sys

eps = 1e-12
def theta_star(w12, w13):
   """
   Nonlinear analytic expression to approximate the mapping of eigenvalues
   based on Howell_McIntyre_2016 for load preservation method.
   """
   v1 = 2.15
   u1 = 1.21
   m = 8e-1
   v2 = 1.85
   u2 = 1.12
   n = 8e-1
   
   theta = (v1/(np.power(u1/(w12+eps),m)+1))*\
           (v2/(np.power(u2/(w13+eps),n)+1))
   
   return np.round(theta, -int(np.log10(eps))) # to supress round-off
   
def fill_out_in_parallel(z_ind_vector,tensor_order,scaling_method,args):
    i,j=args
    tmp = np.ctypeslib.as_array(shared_array)
    #tmp_DTITK = np.ctypeslib.as_array(shared_array_DTITK)

    data_reshape=np.zeros((DTI_data.shape[0],DTI_data.shape[1],DTI_data.shape[2],6),float)




    for k in z_ind_vector:

        if len(DTI_data.shape)>4:
            data_reshape[i,j,k,:]=DTI_data[i,j,k,0,:]
        elif len(DTI_data.shape)==4:
            data_reshape[i,j,k,:]=DTI_data[i,j,k,:]



        if np.all(data_reshape[i,j,k,:]==0.0):
            tmp[i,j,k,:]=np.array([1.0,0.0,0.0,1.0,0.0,1.0])
        else:
            if tensor_order=="NIFTI":
                matrix_from_array=np.array([[data_reshape[i,j,k,0],data_reshape[i,j,k,1],data_reshape[i,j,k,3]],   #if tensor is ordered xx,yx,yy,zx,zy,zz (NIFTI standard)
                                            [data_reshape[i,j,k,1],data_reshape[i,j,k,2],data_reshape[i,j,k,4]],
                                            [data_reshape[i,j,k,3],data_reshape[i,j,k,4],data_reshape[i,j,k,5]]])
            elif tensor_order=="DSI_studio":
                 matrix_from_array=np.array([[data_reshape[i,j,k,0],data_reshape[i,j,k,3],data_reshape[i,j,k,4]],   #if tensor is ordered xx,yy,zz,yx,zx,zy (DSI Studio)
                                             [data_reshape[i,j,k,3],data_reshape[i,j,k,1],data_reshape[i,j,k,5]],
                                             [data_reshape[i,j,k,4],data_reshape[i,j,k,5],data_reshape[i,j,k,2]]])
            elif tensor_order=="FSL":
                 matrix_from_array=np.array([[data_reshape[i,j,k,0],data_reshape[i,j,k,1],data_reshape[i,j,k,2]],   #if tensor is ordered xx,yx,zx,yy,zy,zz (FSL)
                                             [data_reshape[i,j,k,1],data_reshape[i,j,k,3],data_reshape[i,j,k,4]],
                                             [data_reshape[i,j,k,2],data_reshape[i,j,k,4],data_reshape[i,j,k,5]]])
            elif tensor_order=="Johnson_Wistar":
                 matrix_from_array=np.array([[data_reshape[i,j,k,3],data_reshape[i,j,k,1],data_reshape[i,j,k,2]],   #if tensor is ordered yy,yx,zy,xx,zx,zz (Johnson Wistar)
                                             [data_reshape[i,j,k,1],data_reshape[i,j,k,0],data_reshape[i,j,k,4]],
                                             [data_reshape[i,j,k,2],data_reshape[i,j,k,4],data_reshape[i,j,k,5]]])

            #compute eigenvalues and eigenvectors
            eigVals,eigVecs = np.linalg.eig(matrix_from_array)
            if np.any(eigVals<0):
                 print("Warning, no negative eigenvalues should be present for DTI voxels by definition, taking an absolute value. But check your DTI data!")
                 #raise SystemExit
                 eigVals=abs(eigVals)

            # ein einfacher Weg um zu überprüfen, ob die Umrechnugen stimmen, ist sich den max Wert von data (Diffusion=0.003 mm²/s)und data_reshape (Conductivity=2.5 S/m) für CSF anzusehen

            if scaling_method=='Tuch':
                #CRP approach as in Tuch el al. 2001 linear fit
                k_DTI=0.844e3 #S*s/mm³ Umrechnung in (S/m)*(s/mm²) daher e3
                d_epsilon=0.124e-3 #micrometer^2/ms extracellular diffusivity Umrechnung in mm²/s daher e-3
                eigVals_scaled=k_DTI*(eigVals-d_epsilon)
#
#            elif scaling_method=='Nordin':
#                #CRP Nordin approach as in Nordin el al. 2019 BUT Do Not Use EigenValues that is wrong use Dxx Dyy Dzz
#                eigVals_scaled=eigVals/((data_reshape[i,j,k,0]+data_reshape[i,j,k,2]+data_reshape[i,j,k,5])/3)
#                

            elif scaling_method=='NormMapping':
            ##Normalized MAPPING approach as in Güllmar et al./Schmidt el al.
                eigVals_scaled=eigVals/(eigVals[0]*eigVals[1]*eigVals[2])**(1/3.0)

            ##Load preservation method as in Howell, B., McIntyre, C.C., 2016.
            elif scaling_method=='LoadPreservation':
                w12 = eigVals[0]/(eigVals[1]+eps)
                w13 = eigVals[0]/(eigVals[1]+eps)
                theta = theta_star(w12, w13)
                eigVals_scaled = np.array([1., 1./(w12+eps), 1./(w13+eps)])*theta

            if np.any(eigVals_scaled<=0): # if there are still negative eigenvalues put eigVals<=0.0000001
                print("Error, no negative eigenvalues are allowed by definition!")
                raise SystemExit
            #Using eigendecomposition of an SPD tensor (A=Q*lambda*Qtransposed)
            #A is the SPD tensor, which here is the diffusion tensor (which is scaled to become the conductivity tensor)
            #Q is an orthogonal matrix whose columns are eigenvectors of A
            #lambda is the diagonal matrix with the eigenvalues as entries

            elif scaling_method=='Nordin':
                #CRP Nordin approach as in Nordin el al. 2019
                tensor=matrix_from_array/((data_reshape[i,j,k,0]+data_reshape[i,j,k,2]+data_reshape[i,j,k,5])/3)
            else:
                eigVals_matrix=np.diag(eigVals_scaled) #HAS TO BE COMMENTED FOR CRP AS IN ASTROEM!
                tensor=eigVecs.dot(eigVals_matrix).dot(eigVecs.T) #HAS TO BE COMMENTED FOR CRP AS IN ASTROEM!



            eigVals_tensor,eigVecs_tensor = np.linalg.eig(tensor)
            if np.any(eigVals_tensor<=0): # if there are still negative eigenvalues put eigVals<=0.0000001
                 print("Error, no negative eigenvalues are allowed by definition!")
                 raise SystemExit

            # define a lower boundary for the used conductivity (mean value for all tissues in the brain):
#                Sigma_iso_low=1.28e-1 # 1.28e-1 [S/m] IT'IS data base for low frequency conductivities; WM across fibers mean value
#                Sigma_iso_low=1.1e-1 # 1.1e-1 [S/m] IT'IS data base for low frequency conductivities; Brain minimum value
            Sigma_iso_low=0.027512 # 0.027512 [S/m] Gabriel et al. for 10 Hz, WM
            Sigma_iso_lowerBoundary=Sigma_iso_low/20 # 10% of sigma_iso
#                Sigma_iso_lowerBoundary=0.000025

            #define a lower boundary for the eigenvalues:
            if np.any(eigVals_tensor*Sigma_iso_low<Sigma_iso_lowerBoundary):
#                if np.any((eigVals*Sigma_iso_low)<Sigma_iso_lowerBoundary):
                print("Error, lower boundary detected at Voxel:")
                print(i,j,k)
                raise SystemExit

            tmp[i,j,k,:]=np.array([tensor[0][0],tensor[1][0],tensor[2][0],tensor[1][1],tensor[2][1],tensor[2][2]]) #we need to have it as xx,yx,zx,yy,zy,zz (which is FSL standard saving procedure of DTI)
            #for visualization with DTI TK
            #tmp_DTITK[i_ind,j_ind,z_ind,:]=np.array([tensor[0][0],tensor[1][0],tensor[1][1],tensor[2][0],tensor[2][1],tensor[2][2]]) #we need to have it as xx,yx,yy,zx,zy,zz (NIFTI standard)

    #%%

def main_part(tensor_order,scaling_method):
    global shared_array
    #global shared_array_DTITK #visulazation with DTI TK

    Mx,My,Mz=(DTI_data.shape[0],DTI_data.shape[1],DTI_data.shape[2])

    normalized_DTI=np.ctypeslib.as_ctypes(np.zeros((Mx,My,Mz,6),float))
    #normalized_DTITK=np.ctypeslib.as_ctypes(np.zeros((Mx,My,Mz,6),float)) #visulazation with DTI TK
    #np.ctypeslib.as_ctypes = Create and return a ctypes object from a numpy array. Actually anything that exposes the __array_interface__ is accepted.
    shared_array = sharedctypes.RawArray(normalized_DTI._type_, normalized_DTI)
    #shared_array_DTITK = sharedctypes.RawArray(normalized_DTITK._type_, normalized_DTITK) #visulazation with DTI TK
    # multiprocessing.sharedctypes.RawArray(typecode_or_type, size_or_initializer) = Returns a ctypes array allocated from shared memory.

    #i_vector=np.arange(Mx)
    #j_vector=np.arange(My)
    k_vector=np.arange(Mz)

    # we will iterate over z inside the parallelized function
    # this will create all combinations of x,y indices
    window_idxs = [(i, j) for i, j in
               itertools.product(range(0, Mx),    #Combinatoric iterators: product() = cartesian product, equivalent to a nested for-loop [(0,0),(0,1),..,(0,len(My)),..(len(Mx),len(My))]
                                 range(0, My))]
    p = Pool(cpu_count()-1)
    p.map(partial(fill_out_in_parallel,k_vector,tensor_order,scaling_method),window_idxs)
    #p.map =  map(func, iterable[, chunksize]) This method chops the iterable into a number of chunks which it submits to the process pool as separate tasks.
    #The (approximate) size of these chunks can be specified by setting chunksize to a positive integer.
    #partial(func,/,*args,**keywords)
    p.terminate()
    normalized_DTI = np.ctypeslib.as_array(shared_array)


def scale_tensor_data(tensor_data_name,scaling_method='NormMapping',tensor_order='NIFTI'):

    global DTI_data
    #DTI_data=np.zeros((18,21,18,6),float)
    #load DTI data
    filepath = os.path.realpath(tensor_data_name)
    img = nib.load(filepath)
    # img.shape
    DTI_data = img.get_fdata()
    if np.any(np.isnan(DTI_data))
        print("NaN detected in the DTI, please remove them!")
        raise SystemExit

#    #plot DTI data as test
#    fig = plt.figure()
#    a = fig.add_subplot(1, 3, 1)
#    img_ax = np.rot90(DTI_data[..., 90,0,0])
#    imgplot = plt.imshow(img_ax)
#    a.axis('off')
#    a.set_title('Axial_org')
#    a = fig.add_subplot(1, 3, 2)
#    img_cor = np.rot90(DTI_data[:, 100, :,0,0])
#    imgplot = plt.imshow(img_cor)
#    a.axis('off')
#    a.set_title('Coronal_org')
#    a = fig.add_subplot(1, 3, 3)
#    img_cor = np.rot90(DTI_data[90, :, :,0,0])
#    imgplot = plt.imshow(img_cor)
#    a.axis('off')
#    a.set_title('Sagittal_org')

    main_part(tensor_order,scaling_method)
    normalized_DTI=np.ctypeslib.as_array(shared_array)
    #normalized_DTITK=np.ctypeslib.as_array(shared_array_DTITK)


    #save data
    img3 = nib.Nifti1Image(normalized_DTI, img.affine)
    nib.save(img3, nib.filename_parser.splitext_addext(filepath)[0]+'_'+scaling_method+'.nii.gz')          #Has to be changed either _CRPTuch2.nii.gz or _CRPTuch1.nii.gz or _CRPAstroem.nii.gz or _CRPNordin.nii.gz or _normalizedMapping.nii.gz _loadPreservation.nii.gz
    #nib.save(img3, tensor_data_name+'_'+scaling_method+'.nii.gz')

    # save for visualization with DTI TK
    #img4 = nib.Nifti1Image(normalized_DTITK, img.affine)
    #nib.save(img4, 'DTI_TK_IIT3mean_tensor_normalizedMapping.nii.gz')          #Has to be changed either _CRP.nii.gz or _CRPAstroem.nii.gz or _CRPNordin.nii.gz or _normalizedMapping.nii.gz _loadPreservation.nii.gz


    #test plot
#    fig = plt.figure()
#    a = fig.add_subplot(1, 3, 1)
#    img_ax = np.rot90(normalized_DTI[..., 90,0])
#    imgplot = plt.imshow(img_ax)
#    a.axis('off')
#    a.set_title('Axial_norm')
#    a = fig.add_subplot(1, 3, 2)
#    img_cor = np.rot90(normalized_DTI[:, 100, :,0])
#    imgplot = plt.imshow(img_cor)
#    a.axis('off')
#    a.set_title('Coronal_norm')
#    a = fig.add_subplot(1, 3, 3)
#    img_cor = np.rot90(normalized_DTI[90, :, :,0])
#    imgplot = plt.imshow(img_cor)
#    a.axis('off')
#    a.set_title('Sagittal_norm')

if __name__ == '__main__':
    scale_tensor_data(*sys.argv[1:])
