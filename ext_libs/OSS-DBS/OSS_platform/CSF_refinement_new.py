'''Written by K.Butenko'''
'''The script refines provided mesh in the regions of specified material (CSF)'''


import logging
logging.getLogger('UFL').setLevel(logging.WARNING)
logging.getLogger('FFC').setLevel(logging.WARNING)

from dolfin import *

from pandas import read_csv
from tissue_dielectrics import DielectricProperties
import numpy as np
import os

import time as tim

#launch_CSF_refinement is the manager function (called in Launcher)

parameters['linear_algebra_backend']='PETSc'
parameters["refinement_algorithm"] = "plaza_with_parent_facets"
parameters["allow_extrapolation"] = True;


class Field_calc_parameters:
    def __init__(self,default_material,element_order,anisotropy,c_c,CPE,refinement_frequency,Laplace_formulation,external_grounding,Solver):
        self.default_material=default_material
        self.element_order=element_order
        self.anisotropy=anisotropy
        self.c_c=c_c
        self.CPE=CPE
        self.frequenc=refinement_frequency      #list
        self.EQS_mode=Laplace_formulation
        self.external_grounding=external_grounding
        self.Solver_type=Solver

def save_mesh_and_subdomains_to_h5(mesh_to_h5,subdomains_to_h5,subdomains_assigned_to_h5,boundaries_to_h5,Scaling):

    hdf = HDF5File(mesh_to_h5.mpi_comm(), os.environ['PATIENTDIR']+'/CSF_ref/Mesh_to_solve_scaling_'+str(Scaling)+'.h5', 'w')
    hdf.write(mesh_to_h5, "/mesh")
    hdf.write(subdomains_to_h5, "/subdomains")
    hdf.write(subdomains_assigned_to_h5, "/subdomains_assigned")
    hdf.write(boundaries_to_h5, "/boundaries")
    hdf.close()

    return True

def load_mesh_and_subdomains_from_h5(Scaling):

    mesh_from_h5 = Mesh()
    hdf = HDF5File(mesh_from_h5.mpi_comm(), os.environ['PATIENTDIR']+'/CSF_ref/Mesh_to_solve_scaling_'+str(Scaling)+'.h5', "r")
    hdf.read(mesh_from_h5, "/mesh", False)
    subdomains_from_h5 = MeshFunction("size_t", mesh_from_h5, 3)
    hdf.read(subdomains_from_h5, "/subdomains")
    subdomains_assigned_from_h5 = MeshFunction("size_t", mesh_from_h5, 3)
    hdf.read(subdomains_assigned_from_h5, "/subdomains_assigned")
    boundaries_from_h5 = MeshFunction("size_t", mesh_from_h5, 2)
    hdf.read(boundaries_from_h5, "/boundaries")
    hdf.close()

    return mesh_from_h5,subdomains_from_h5,subdomains_assigned_from_h5,boundaries_from_h5

def Dummy_CSF():        #if we want to skip adaptive mesh refinement

    mesh = Mesh(os.environ['PATIENTDIR']+"/Meshes/Mesh_unref.xml")
    boundaries = MeshFunction('size_t',mesh,os.environ['PATIENTDIR']+'/Meshes/Mesh_unref_facet_region.xml')
    subdomains_assigned=MeshFunction('size_t',mesh,os.environ['PATIENTDIR']+"/Meshes/Mesh_unref_physical_region.xml")


    # print("Before ground ref: ",mesh.num_cells())
    # from Ext_ground_preref import prerefine_ground

    # mesh,subdomains_assigned,boundaries=prerefine_ground(mesh,subdomains_assigned,boundaries,Domains)
    # mesh,subdomains_assigned,boundaries=prerefine_ground(mesh,subdomains_assigned,boundaries,Domains)
    # mesh,subdomains_assigned,boundaries=prerefine_ground(mesh,subdomains_assigned,boundaries,Domains)
    # mesh,subdomains_assigned,boundaries=prerefine_ground(mesh,subdomains_assigned,boundaries,Domains)
    # mesh,subdomains_assigned,boundaries=prerefine_ground(mesh,subdomains_assigned,boundaries,Domains)
    # #mesh,subdomains_assigned,boundaries=prerefine_ground(mesh,subdomains_assigned,boundaries,Domains)
    # #mesh,subdomains_assigned,boundaries=prerefine_ground(mesh,subdomains_assigned,boundaries,Domains)
    # print("After ground ref: ",mesh.num_cells())

    mesh_file=File(os.environ['PATIENTDIR']+'/Results_adaptive/mesh_adapt.xml.gz')
    boundaries_file = File(os.environ['PATIENTDIR']+'/Results_adaptive/boundaries_adapt.xml')
    subdomains_assigned_file=File(os.environ['PATIENTDIR']+'/Results_adaptive/subdomains_assigned_adapt.xml')

    mesh_file<<mesh
    boundaries_file<<boundaries
    subdomains_assigned_file<<subdomains_assigned

    return True

def mesh_refiner(mesh_old,boundaries,subdomains_assigned,cell_markers,Domains,cc_multicontact):
    parameters['linear_algebra_backend']='PETSc'
    parameters["refinement_algorithm"] = "plaza_with_parent_facets"
    parameters["allow_extrapolation"] = True;

    facets_old = MeshFunction('size_t',mesh_old,2)
    facets_old.set_all(0)

    # take the first contact to check whether the mesh refinement was done correctly
    facets_old.array()[boundaries.array() == Domains.Active_contacts[0]]=1

    if cc_multicontact==True and Domains.Amp_vector[0]!=0.0:            #because ground contact is always subtracted from the mesh
        dsSSS=Measure("dS",domain=mesh_old,subdomain_data=facets_old)
        An_surface_size_old=assemble(1.0*dsSSS(1))
    else:
        dss=Measure("ds",domain=mesh_old,subdomain_data=facets_old)
        An_surface_size_old=assemble(1.0*dss(1))

    # refine marked cells on the old mesh and adapt other entities
    mesh_new = refine(mesh_old, cell_markers)
    subdomains_assigned_new=adapt(subdomains_assigned,mesh_new)
    boundaries_new = adapt(boundaries,mesh_new) # put function space

    facets = MeshFunction('size_t',mesh_new,2)
    facets.set_all(0)
    facets.array()[boundaries_new.array()==Domains.Active_contacts[0]]=1

    if cc_multicontact==True and Domains.Amp_vector[0]!=0.0:
        dsS_new=Measure("dS",domain=mesh_new,subdomain_data=facets)
        An_surface_size_new=assemble(1.0*dsS_new(1))
    else:
        dss_new=Measure("ds",domain=mesh_new,subdomain_data=facets)
        An_surface_size_new=assemble(1.0*dss_new(1))

    # if the mesh refinement failed, it will be detected here
    if (An_surface_size_new-An_surface_size_old)/An_surface_size_new>0.01:
        #(An_surface_size_new-An_surface_size_old)/An_surface_size_new
        logging.critical("Refinement broke the imposed B.C.!")
        exit()

    return (mesh_new,boundaries_new,subdomains_assigned_new)

#refine cells which contain CSF voxels (listed in index_array)
#scaler is the maximum allowed element edge/voxel size ratio
def index_cell_marker(mesh, index_array ,MRI_param, Scaler):

    cell_ref = MeshFunction('bool',mesh,3)
    cell_ref.set_all(False)
    cell_to_ref=0
    cell_processed=0
    c00 = MeshFunction("double", mesh, 3)       #to check, which cells will be refined
    for cell in cells(mesh):
        cell_processed=cell_processed+1
        smallest_edge=min([MRI_param.voxel_dims[0],MRI_param.voxel_dims[1],MRI_param.voxel_dims[2]])   # of the voxel
        if np.any(np.isin(index_array,cell.index())) and cell.h()>Scaler*smallest_edge:
            cell_ref[cell] = True
            cell_to_ref=cell_to_ref+1
            c00[cell]=1.0
        else:
            c00[cell]=0.0

    return cell_ref


def get_CSF_voxels(MRI_param, Array_coord, Array_Type):
    if Array_Type == 'Imported':
        List_of_placed = Array_coord

    # first, the neuron compartments should not pass through CSF, maybe put it as a separate function, maybe even in CSF_refinement
    import os
    if os.path.isfile(os.environ['PATIENTDIR']+'/MRI_DTI_derived_data/' + MRI_param.name[:-4] + '_voxel_array_CSF.npy') or os.path.isfile(
            os.environ['PATIENTDIR']+'/MRI_DTI_derived_data/' + MRI_param.name[:-7] + '_voxel_array_CSF.npy'):  # if array was already prepared
        if MRI_param.name[-2:] == 'gz':
            voxel_array_CSF = np.load(os.environ['PATIENTDIR']+'/MRI_DTI_derived_data/' + MRI_param.name[:-7] + '_voxel_array_CSF.npy')
        else:
            voxel_array_CSF = np.load(os.environ['PATIENTDIR']+'/MRI_DTI_derived_data/' + MRI_param.name[:-4] + '_voxel_array_CSF.npy')
        logging.critical("voxel_array_CSF is loaded")
    else:  # otherwise prepare an array that stores coordinated of all voxels with CSF in the vicinity of the neurons
        Tissue_array = np.load(os.environ['PATIENTDIR']+'/MRI_DTI_derived_data/Tissue_array_MRI.npy')

        voxel_array_CSF = np.zeros((Tissue_array.shape[0], 3),
                                   float)  # array to store all CSF voxels in the specified ROI

        if Array_Type != 'Imported':
            x_neuron_max, y_neuron_max, z_neuron_max = (
            max(Array_coord[:, 0]), max(Array_coord[:, 1]), max(Array_coord[:, 2]))
            x_neuron_min, y_neuron_min, z_neuron_min = (
            min(Array_coord[:, 0]), min(Array_coord[:, 1]), min(Array_coord[:, 2]))
        else:
            max_values = np.zeros((len(List_of_placed), 3), float)
            min_values = np.zeros((len(List_of_placed), 3), float)
            for i in range(len(List_of_placed)):
                max_values[i, :] = (
                max(List_of_placed[i][:, 0]), max(List_of_placed[i][:, 1]), max(List_of_placed[i][:, 2]))
                min_values[i, :] = (
                min(List_of_placed[i][:, 0]), min(List_of_placed[i][:, 1]), min(List_of_placed[i][:, 2]))

            x_neuron_max, y_neuron_max, z_neuron_max = (
            max(max_values[:, 0]), max(max_values[:, 1]), max(max_values[:, 2]))
            x_neuron_min, y_neuron_min, z_neuron_min = (
            min(min_values[:, 0]), min(min_values[:, 1]), min(min_values[:, 2]))

        space_from_neurons = 1.0  # here we do not need to check further away

        # to create a nifti of CSF (for a visual check of intersections)
        CSF_nifti = np.zeros((MRI_param.N_voxels[0], MRI_param.N_voxels[1], MRI_param.N_voxels[2]), int)

        for z_i in range(int(MRI_param.N_voxels[2])):
            for y_i in range(int(MRI_param.N_voxels[1])):
                for x_i in range(int(MRI_param.N_voxels[0])):
                    # this is defined in O(0,0,0) of the voxel
                    coords_reals = np.dot(MRI_param.affine_MRI, np.array([x_i, y_i, z_i, 1.0]))

                    # this is defined in the center of the voxel
                    x_pos = coords_reals[0] + MRI_param.MRI_shift[0] + MRI_param.voxel_dims[0]/ 2.0
                    y_pos = coords_reals[1] + MRI_param.MRI_shift[1] + MRI_param.voxel_dims[1] / 2.0
                    z_pos = coords_reals[2] + MRI_param.MRI_shift[2] + MRI_param.voxel_dims[2] / 2.0

                    # this is defined in the opposite corner from O(0,0,0)
                    x_coord, y_coord, z_coord = (x_pos + MRI_param.voxel_dims[0] / 2.0, y_pos + MRI_param.voxel_dims[1] / 2.0, z_pos + MRI_param.voxel_dims[2] / 2.0)

                    if x_neuron_max + space_from_neurons >= x_pos >= x_neuron_min - space_from_neurons and y_neuron_max + space_from_neurons >= y_pos >= y_neuron_min - space_from_neurons and z_neuron_max + space_from_neurons >= z_pos >= z_neuron_min - space_from_neurons:

                        glob_index = x_i + y_i*int(MRI_param.N_voxels[0]) + z_i*int(MRI_param.N_voxels[0])*int(MRI_param.N_voxels[1])
                        glob_index = int(glob_index)

                        if Tissue_array[glob_index] == 1:
                            voxel_array_CSF[glob_index, 0] = x_pos
                            voxel_array_CSF[glob_index, 1] = y_pos
                            voxel_array_CSF[glob_index, 2] = z_pos

                            CSF_nifti[x_i,y_i,z_i] = 1





        # import nibabel as nib
        # import os
        # example_filename = os.path.join(os.environ['PATIENTDIR'] + '/' + MRI_param.name)
        # img = nib.load(example_filename)
        # img3 = nib.Nifti1Image(CSF_nifti, MRI_param.affine_MRI, img.header)
        # nib.save(img3, os.environ['PATIENTDIR'] + '/Neuron_model_arrays/CSF_vox_nearby.nii')

        voxel_array_CSF = voxel_array_CSF[~np.all(voxel_array_CSF == 0.0, axis=1)]  # deletes all zero enteries

        if MRI_param.name[-2:] == 'gz':
            np.save(os.environ['PATIENTDIR']+'/MRI_DTI_derived_data/' + MRI_param.name[:-7] + '_voxel_array_CSF', voxel_array_CSF)
        else:
            np.save(os.environ['PATIENTDIR']+'/MRI_DTI_derived_data/' + MRI_param.name[:-4] + '_voxel_array_CSF', voxel_array_CSF)

        del Tissue_array
        logging.critical("voxel_array_CSF (contains CSF voxels close to the neuron array) is prepared")

    voxel_array_CSF_shifted = np.zeros((voxel_array_CSF.shape[0], 3), float)
    voxel_array_CSF_shifted[:, 0] = voxel_array_CSF[:, 0] + MRI_param.voxel_dims[0] / 2
    voxel_array_CSF_shifted[:, 1] = voxel_array_CSF[:, 1] + MRI_param.voxel_dims[1] / 2
    voxel_array_CSF_shifted[:, 2] = voxel_array_CSF[:, 2] + MRI_param.voxel_dims[2] / 2
    del voxel_array_CSF

    return voxel_array_CSF_shifted


def Refine_CSF(MRI_param,DTI_param,Scaling,Domains,Field_calc_param,rel_div,CSF_frac_div,CSF_ref_add,EQS_mode,cc_multicontact,ref_freq,Best_scaling=0,scaling_old=0):

    start_CSF_refinement=tim.time()

    if cc_multicontact==True:
        from Math_module_hybrid_floating import compute_field_with_superposition,get_field_on_points
    else:
        from Math_module_hybrid import get_field,get_field_on_points

    scaling_old=float(scaling_old)
    Scaling=float(Scaling)
    Best_scaling=float(Best_scaling)

    if Field_calc_param.anisotropy==1:
        from Tissue_marking_new import get_cellmap_tensors
    else:
        from Tissue_marking_new import get_cellmap

    '''load results and mesh from the most refined iteration if available'''
    if Best_scaling!=0:
        Phi_amp_on_neuron_old_get=read_csv(os.environ['PATIENTDIR']+'/CSF_ref/Field_on_points'+str(Best_scaling)+'.csv', delimiter=' ', header=None)
        Phi_amp_on_neuron_old=Phi_amp_on_neuron_old_get.values

        mesh = Mesh(os.environ['PATIENTDIR']+'/CSF_ref/mesh_adapt_CSF'+str(scaling_old)+'.xml.gz')
        boundaries = MeshFunction('size_t',mesh,os.environ['PATIENTDIR']+'/CSF_ref/boundaries_adapt_CSF'+str(scaling_old)+'.xml')
        subdomains_assigned=MeshFunction('size_t',mesh,os.environ['PATIENTDIR']+'/CSF_ref/subdomains_assigned_adapt_CSF'+str(scaling_old)+'.xml')
    else:
        if Field_calc_param.external_grounding==False:
            mesh = Mesh(os.environ['PATIENTDIR']+"/Meshes/Mesh_unref.xml")
            boundaries = MeshFunction('size_t',mesh,os.environ['PATIENTDIR']+'/Meshes/Mesh_unref_facet_region.xml')
            subdomains_assigned=MeshFunction('size_t',mesh,os.environ['PATIENTDIR']+"/Meshes/Mesh_unref_physical_region.xml")
        else:
            mesh = Mesh(os.environ['PATIENTDIR']+'/Results_adaptive/mesh_adapt.xml.gz')
            boundaries = MeshFunction('size_t',mesh,os.environ['PATIENTDIR']+'/Results_adaptive/boundaries_adapt.xml')
            subdomains_assigned=MeshFunction('size_t',mesh,os.environ['PATIENTDIR']+'/Results_adaptive/subdomains_assigned_adapt.xml')

    # load neuron compartments
    Vertices_neur_get=read_csv(os.environ['PATIENTDIR']+'/Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', delimiter=' ', header=None)
    Vertices_neur=Vertices_neur_get.values

    if Best_scaling==0:
        mesh_file=File(os.environ['PATIENTDIR']+'/CSF_ref/mesh_adapt_CSF'+str(Best_scaling)+'.xml.gz')
        boundaries_file = File(os.environ['PATIENTDIR']+'/CSF_ref/boundaries_adapt_CSF'+str(Best_scaling)+'.xml')
        subdomains_assigned_file=File(os.environ['PATIENTDIR']+'/CSF_ref/subdomains_assigned_adapt_CSF'+str(Best_scaling)+'.xml')
        mesh_file<<mesh
        boundaries_file<<boundaries
        subdomains_assigned_file<<subdomains_assigned

        logging.critical("Field calculation on the initial mesh")
        if Field_calc_param.anisotropy==1:
            subdomains=get_cellmap_tensors(mesh,subdomains_assigned,Domains,MRI_param,DTI_param,Field_calc_param.default_material)
        else:
            subdomains=get_cellmap(mesh,subdomains_assigned,Domains,MRI_param,Field_calc_param.default_material)

        save_mesh_and_subdomains_to_h5(mesh,subdomains,subdomains_assigned,boundaries,Best_scaling)

        logging.critical("CSF_Subdomains_unref file was created")
        file=File(os.environ['PATIENTDIR']+'/CSF_ref/CSF_Subdomains_unref.pvd')
        file<<subdomains,mesh

        if cc_multicontact==True:
            Phi_r_old,Phi_im_old,Field_r_old,Field_im_old,max_E_old,J_r_old,J_im_old,j_dens_real,j_dens_im=compute_field_with_superposition(mesh,Domains,subdomains_assigned,subdomains,boundaries,Field_calc_param)
        else:
            Phi_r_old,Phi_im_old,Field_r_old,Field_im_old,max_E_old,J_r_old,J_im_old,j_dens_real,j_dens_im=get_field(mesh,Domains,subdomains,boundaries,Field_calc_param)

        file=File(os.environ['PATIENTDIR']+'/CSF_ref/Field_r_'+str(Best_scaling)+'.pvd')
        file<<Field_r_old

        Phi_amp_on_neuron_old = get_field_on_points(Phi_r_old,Phi_im_old,Field_calc_param.c_c,J_r_old,J_im_old)
        np.savetxt(os.environ['PATIENTDIR']+'/CSF_ref/Field_on_points'+str(Best_scaling)+'.csv', Phi_amp_on_neuron_old, delimiter=" ")

    # load array with the CSF voxels if available, otherwise extract it from the MRI data
    if os.path.isfile(os.environ['PATIENTDIR']+'/MRI_DTI_derived_data/'+MRI_param.name[:-4]+'_voxel_array_CSF_'+str(CSF_ref_add)+'.npy') or os.path.isfile(os.environ['PATIENTDIR']+'/MRI_DTI_derived_data/'+MRI_param.name[:-7]+'_voxel_array_CSF_'+str(CSF_ref_add)+'.npy'):
        if MRI_param.name[-2:]=='gz':
            voxel_array_CSF=np.load(os.environ['PATIENTDIR']+'/MRI_DTI_derived_data/'+MRI_param.name[:-7]+'_voxel_array_CSF_'+str(CSF_ref_add)+'.npy')
        else:
            voxel_array_CSF=np.load(os.environ['PATIENTDIR']+'/MRI_DTI_derived_data/'+MRI_param.name[:-4]+'_voxel_array_CSF_'+str(CSF_ref_add)+'.npy')

        logging.critical("voxel_array_CSF in {} mm vicinity is loaded".format(str(CSF_ref_add)))
    else:
        start_voxel_array_CSF=tim.time()

        Tissue_array=np.load(os.environ['PATIENTDIR']+'/MRI_DTI_derived_data/Tissue_array_MRI.npy')
        voxel_array_CSF=np.zeros((Tissue_array.shape[0],3),float)      #array to store all CSF voxels in the specified ROI

        bb = mesh.bounding_box_tree()
        #check the extent of the neuron array (loof for CSF only there + vicinity defined by CSF_ref_add)
        x_neuron_max=max(Vertices_neur[:,0])
        y_neuron_max=max(Vertices_neur[:,1])
        z_neuron_max=max(Vertices_neur[:,2])
        x_neuron_min=min(Vertices_neur[:,0])
        y_neuron_min=min(Vertices_neur[:,1])
        z_neuron_min=min(Vertices_neur[:,2])

        #go over all voxels and check whether it contains CSF and intersect with the mesh
        
        
        affine=MRI_param.affine_MRI

        glob_index=0

        for z_i in range(int(MRI_param.N_voxels[2])):
            for y_i in range(int(MRI_param.N_voxels[1])):
                for x_i in range(int(MRI_param.N_voxels[0])):

                    coords_reals = np.dot(affine,np.array([x_i,y_i,z_i,1.0]))

                    x_pos=coords_reals[0]+MRI_param.MRI_shift[0]+MRI_param.voxel_dims[0]/2.0
                    y_pos=coords_reals[1]+MRI_param.MRI_shift[1]+MRI_param.voxel_dims[1]/2.0
                    z_pos=coords_reals[2]+MRI_param.MRI_shift[2]+MRI_param.voxel_dims[2]/2.0

                    if (x_pos<=x_neuron_max+CSF_ref_add and x_pos>=x_neuron_min-CSF_ref_add and y_pos<=y_neuron_max+CSF_ref_add and y_pos>=y_neuron_min-CSF_ref_add and z_pos<=z_neuron_max+CSF_ref_add and z_pos>=z_neuron_min-CSF_ref_add):
                        pnt = Point(x_pos,y_pos,z_pos)

                        glob_index = x_i + y_i*int(MRI_param.N_voxels[0]) + z_i*int(MRI_param.N_voxels[0])*int(MRI_param.N_voxels[1])
                        glob_index = int(glob_index)

                        if Tissue_array[glob_index] == 1 and bb.compute_first_entity_collision(pnt)<mesh.num_cells()*100:
                            voxel_array_CSF[glob_index,0] = x_pos
                            voxel_array_CSF[glob_index,1] = y_pos
                            voxel_array_CSF[glob_index,2] = z_pos

        voxel_array_CSF=voxel_array_CSF[~np.all(voxel_array_CSF==0.0,axis=1)]  #deletes all zero enteries

        if MRI_param.name[-2:]=='gz':
            np.save(os.environ['PATIENTDIR']+'/MRI_DTI_derived_data/'+MRI_param.name[:-7]+'_voxel_array_CSF_'+str(CSF_ref_add), voxel_array_CSF)
        else:
            np.save(os.environ['PATIENTDIR']+'/MRI_DTI_derived_data/'+MRI_param.name[:-4]+'_voxel_array_CSF_'+str(CSF_ref_add), voxel_array_CSF)

        del Tissue_array
        logging.critical("----- voxel_array_CSF for {} mm vicinity was prepared in {} seconds -----".format(str(CSF_ref_add),tim.time() - start_voxel_array_CSF))


    '''Here we pre-refine mesh on elements with CSF voxels'''
    csf_ref=0
    loaded_from_h5=0
    logging.critical("refining CSF voxels with scaling: {}".format(int(Scaling)))
    num_cell_old=mesh.num_cells()

    if os.path.isfile(os.environ['PATIENTDIR']+'/CSF_ref/Mesh_to_solve_scaling_'+str(Scaling)+'.h5'):
        logging.critical("Mesh and tissue mapping was loaded from the refinement at the previous frequency")
        mesh,subdomains,subdomains_assigned,boundaries=load_mesh_and_subdomains_from_h5(Scaling)
        loaded_from_h5=1
    elif os.path.isfile(os.environ['PATIENTDIR']+'/CSF_ref/mesh_adapt_CSF'+str(Scaling)+'.xml.gz'):   #this and not the case above could be only triggered if no computations took place at this scalling at the last frequency (because of the same mesh size)
        if Best_scaling==0:
            logging.critical("No elements were refined during the full CSF refinement, skipping to adaptive mesh refinement (consider decreasing Minimum Element to Voxel Ratio)")
            return 1
        logging.critical("skipping scaling {}".format(Scaling))
        mesh_file=File(os.environ['PATIENTDIR']+'/CSF_ref/mesh_adapt_CSF'+str(Scaling)+'.xml.gz')
        mesh_file<<mesh
        boundaries_file = File(os.environ['PATIENTDIR']+'/CSF_ref/boundaries_adapt_CSF'+str(Scaling)+'.xml')
        subdomains_assigned_file=File(os.environ['PATIENTDIR']+'/CSF_ref/subdomains_assigned_adapt_CSF'+str(Scaling)+'.xml')
        boundaries_file<<boundaries
        subdomains_assigned_file<<subdomains_assigned
        return 0
    else:
        while csf_ref==0:       #refine mesh until get the required edge size to voxel ratio
            inx_pnt=0       #to check, how much voxels were processed
            cell_index_list=[]
            bb = mesh.bounding_box_tree()

            for i_csf in range(voxel_array_CSF.shape[0]):       #find cells which contain the CSF voxels
                pnt=Point(voxel_array_CSF[i_csf,0],voxel_array_CSF[i_csf,1],voxel_array_CSF[i_csf,2])
                inx_pnt=inx_pnt+1
                cell_index_list.append(bb.compute_first_entity_collision(pnt))
            cell_index_array=np.asarray(cell_index_list)
            cell_ref=index_cell_marker(mesh, cell_index_array, MRI_param, Scaling)

            if not(cell_ref.where_equal(True)):     #if any cell was marked for refinement, will return True
                csf_ref=1

                mesh_file=File(os.environ['PATIENTDIR']+'/CSF_ref/mesh_adapt_CSF'+str(Scaling)+'.xml.gz')
                mesh_file<<mesh
                boundaries_file = File(os.environ['PATIENTDIR']+'/CSF_ref/boundaries_adapt_CSF'+str(Scaling)+'.xml')
                subdomains_assigned_file=File(os.environ['PATIENTDIR']+'/CSF_ref/subdomains_assigned_adapt_CSF'+str(Scaling)+'.xml')
                boundaries_file<<boundaries
                subdomains_assigned_file<<subdomains_assigned

                logging.critical("Number of cells after CSF refinement iteration: {}".format(mesh.num_cells()))

                if num_cell_old==mesh.num_cells():
                    if Best_scaling==0:
                        logging.critical("No elements were refined during the full CSF refinement, skipping to adaptive mesh refinement (consider decreasing Minimum Element to Voxel Ratio)")
                        return 1
                    logging.critical("skipping scaling ".format(Scaling))
                    return 0
            else:

                mesh_file=File(os.environ['PATIENTDIR']+'/CSF_ref/mesh_adapt_CSF'+str(Scaling)+'_old.xml.gz')
                mesh_file<<mesh
                boundaries_file = File(os.environ['PATIENTDIR']+'/CSF_ref/boundaries_adapt_CSF'+str(Scaling)+'_old.xml')
                subdomains_assigned_file=File(os.environ['PATIENTDIR']+'/CSF_ref/subdomains_assigned_adapt_CSF'+str(Scaling)+'_old.xml')
                boundaries_file<<boundaries
                subdomains_assigned_file<<subdomains_assigned

                mesh,boundaries,subdomains_assigned=mesh_refiner(mesh,boundaries,subdomains_assigned,cell_ref,Domains,cc_multicontact)

                if mesh.num_cells()>10000000:       #users can adjust for their hardware
                    logging.critical("Mesh is too large, will have to check with bigger scaling")
                    csf_refined=-1
                    return csf_refined

    if loaded_from_h5==1:
        logging.critical("CSF_Subdomains_refinement file with scaling was loaded".format(int(Scaling)))
    else:
        if Field_calc_param.anisotropy==1:
            subdomains=get_cellmap_tensors(mesh,subdomains_assigned,Domains,MRI_param,DTI_param,Field_calc_param.default_material)
        else:
            subdomains=get_cellmap(mesh,subdomains_assigned,Domains,MRI_param,Field_calc_param.default_material)

        logging.critical("CSF_Subdomains_refinement file with scaling was created".format(int(Scaling)))
        file=File(os.environ['PATIENTDIR']+'/CSF_ref/CSF_Subdomains_refinement_'+str(int(Scaling))+'.pvd')
        file<<subdomains,mesh

    save_mesh_and_subdomains_to_h5(mesh,subdomains,subdomains_assigned,boundaries,Scaling)

    if cc_multicontact==True:
        Phi_r,Phi_im,Field_r,Field_im,max_E,J_r,J_im,j_dens_real,j_dens_im=compute_field_with_superposition(mesh,Domains,subdomains_assigned,subdomains,boundaries,Field_calc_param)
    else:
        Phi_r,Phi_im,Field_r,Field_im,max_E,J_r,J_im,j_dens_real,j_dens_im=get_field(mesh,Domains,subdomains,boundaries,Field_calc_param)

    if Scaling==1:
        file=File(os.environ['PATIENTDIR']+'/CSF_ref/Field_r_'+str(Scaling)+'.pvd')
        file<<Field_r
        logging.critical("CSF_Subdomains full refinement was created")
        file=File(os.environ['PATIENTDIR']+'/CSF_ref/CSF_Subdomains_full_ref.pvd')
        file<<subdomains,mesh

        #import subprocess
        #subprocess.call('python Visualization_files/Paraview_CSFref.py', shell=True)

    Phi_amp_on_neuron =get_field_on_points(Phi_r,Phi_im,Field_calc_param.c_c,J_r,J_im)
    np.savetxt(os.environ['PATIENTDIR']+'/CSF_ref/Field_on_points'+str(Scaling)+'.csv', Phi_amp_on_neuron, delimiter=" ")

    csf_refined=1
    max_div=0.0

    # define the absolute error threshold
    if Field_calc_param.c_c == True:
        if EQS_mode=='EQS':
            #not the best approach, but both should be on the Dirichlet BCs
            max_phi_r=max(Phi_r.vector()[:])
            max_phi_im=max(Phi_im.vector()[:])

            min_phi_r=min(Phi_r.vector()[:])
            min_phi_im=min(Phi_im.vector()[:])

            phi_error=abs((np.sqrt((max_phi_r-min_phi_r)**2+(max_phi_im-min_phi_im)**2))*CSF_frac_div)
        else:
            phi_error=abs((max(Phi_r.vector()[:])-min(Phi_r.vector()[:]))*CSF_frac_div)   #should be scaled
    else:
        Phi_vector=[x for x in Domains.Amp_vector if x is not None]
        if Field_calc_param.external_grounding==True:
            Phi_vector.append(0.0)
        phi_error=abs((max(Phi_vector)-min(Phi_vector))*CSF_frac_div)      #Absolute potential error defined as a 1% of the maximum potential difference, VC case

    # compare solutions on the neuron compartments
    for inx in range(Phi_amp_on_neuron_old.shape[0]):
        if Best_scaling==0:             #first iteration
            delim=abs(Phi_amp_on_neuron[inx,3])
        else:
            delim=abs(Phi_amp_on_neuron_old[inx,3])

        if max_div<abs(Phi_amp_on_neuron_old[inx,3]-Phi_amp_on_neuron[inx,3]):
            max_div=abs(Phi_amp_on_neuron_old[inx,3]-Phi_amp_on_neuron[inx,3])

        if max_div> phi_error:
            logging.critical("Deviation threshold: {} V".format(phi_error))
            logging.critical("Deviation at least: {} V".format(max_div))
            logging.critical("At point: {} {} {}".format(Phi_amp_on_neuron_old[inx,0],Phi_amp_on_neuron_old[inx,1],Phi_amp_on_neuron_old[inx,2]))
            logging.critical("Need further refinement of CSF")
            csf_refined=0
            break

    if csf_refined==1:
        mesh_file=File(os.environ['PATIENTDIR']+'/CSF_ref/mesh_adapt_CSF'+str(Scaling)+'.xml.gz')
        mesh_file<<mesh
        boundaries_file = File(os.environ['PATIENTDIR']+'/CSF_ref/boundaries_adapt_CSF'+str(Scaling)+'.xml')
        subdomains_assigned_file=File(os.environ['PATIENTDIR']+'/CSF_ref/subdomains_assigned_adapt_CSF'+str(Scaling)+'.xml')
        boundaries_file<<boundaries
        subdomains_assigned_file<<subdomains_assigned

        logging.critical("Deviation at least: V".format(max_div))
        logging.critical("Deviation threshold: V".format(phi_error))
        logging.critical("CSF is refined enough")

    del voxel_array_CSF

    minutes=int((tim.time() - start_CSF_refinement)/60)
    secnds=int(tim.time() - start_CSF_refinement)-minutes*60
    logging.critical("----- CSF refinement iteration took {} min {} sec -----\n".format(minutes, secnds))

    return csf_refined


def launch_CSF_refinement(d,MRI_param,DTI_param,Domains,anisotrop,cc_multicontact,ref_freqs):

    # the threshold for deviation due to CSF is twice as for the adaptive refinement and limitied to 1%
    d["CSF_frac_div"]=2*d["Adaptive_frac_div"]
    if d["CSF_frac_div"]<0.01:
        d["CSF_frac_div"]=0.01

    el_order_for_CSF=1      #always first order due to the size of the mesh after CSF refinement. Users can increase if hardware allows.

    if cc_multicontact==True and el_order_for_CSF==1:   # for the multicontact current-controlled case we require accurate projection of the electric field
        el_order_for_CSF=2
        "element_order is 1, increasing to 2 for multicontact current-controlled stimulation"

    logging.critical("----- Conducting evaluation of CSF refinement -----")

    Min_Scaling=d["Min_Scaling"]    # the maxixum CSF element ref. criterion (element edge size/voxel site). 1 by default.
    Scaling_results=[]

    for freq in ref_freqs:          # conduct refinement at different frequencies
        logging.critical("At frequency: {}".format(freq))
        Field_calc_param=Field_calc_parameters(d["default_material"],el_order_for_CSF,anisotrop,d["current_control"],d["CPE_activ"],freq,d["EQS_core"],d["external_grounding"],d["Solver_Type"])

        csf_ref=-1
        '''csf_ref is 1, when further refinement of elements with CSF voxels does not significantly change the result'''
        '''in this loop csf_ref will be 0 if the calcualtions were performed, 1 if deviation below the thresold and -1 if calculations failed (too much elements)'''
        while csf_ref==-1:
            csf_ref=Refine_CSF(MRI_param,DTI_param,Min_Scaling,Domains,Field_calc_param,d["rel_div_CSF"],d["CSF_frac_div"],d["CSF_ref_reg"],d["EQS_core"],cc_multicontact,ref_freqs)          #return -1, when solution was not obtained (mesh is too large)
            if csf_ref==-1:
                Min_Scaling=Min_Scaling*2           #if the amount of elements is too large, we will have to refine less.

        if csf_ref==1:
            '''resave unref as CSF with scaling 0'''
            Scaling=0
        else:
            iteration_csf=0
            scaling_factors=[32*Min_Scaling,16*Min_Scaling,8*Min_Scaling,4*Min_Scaling,2*Min_Scaling,1*Min_Scaling]
            scaling_old=0
            while csf_ref==0:       #subsequently refine to smaller edge sizes until the solution is below the deviation threshold
                logging.critical('scaling factor: {}'.format(scaling_factors[iteration_csf]))
                Scaling=float(scaling_factors[iteration_csf])
                csf_ref=Refine_CSF(MRI_param,DTI_param,Scaling,Domains,Field_calc_param,d["rel_div_CSF"],d["CSF_frac_div"],d["CSF_ref_reg"],d["EQS_core"],cc_multicontact,ref_freqs,Min_Scaling,scaling_old)
                iteration_csf=iteration_csf+1
                scaling_old=Scaling

        Scaling_results.append(Scaling)

    # store to load the correct mesh for the next stage of the refinement
    Scaling_arr=np.array([min(Scaling_results)])
    np.savetxt(os.environ['PATIENTDIR']+'/CSF_ref/Scaling_CSF.csv', Scaling_arr, delimiter=" ")
    Scaling=float(Scaling)

    return Scaling

