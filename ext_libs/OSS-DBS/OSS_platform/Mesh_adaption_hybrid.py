# -*- coding: utf-8 -*-
"""
Created on Sun Jun 24 20:53:29 2018

@author: butenko
"""

#adapt_mesh is the manager function (called in Launcher)
import h5py     #testing


import time as time_lib

import logging
logging.getLogger('UFL').setLevel(logging.WARNING)
logging.getLogger('FFC').setLevel(logging.WARNING)

from dolfin import *
from pandas import read_csv
from tissue_dielectrics import DielectricProperties
import numpy as np
import os
import subprocess
import pickle

#from Math_module import get_field,get_field_on_points
#from Math_module_new import get_field,get_field_QS,get_field_dual,get_field_on_points

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

def load_mesh(mesh_designation):
    mesh = Mesh(os.environ['PATIENTDIR']+'/Results_adaptive/mesh_'+mesh_designation+'.xml.gz')
    boundaries = MeshFunction('size_t',mesh,os.environ['PATIENTDIR']+'/Results_adaptive/boundaries_'+mesh_designation+'.xml')
    subdomains = MeshFunction('size_t',mesh,os.environ['PATIENTDIR']+'/Results_adaptive/subdomains_assigned_'+mesh_designation+'.xml')

    return mesh, boundaries, subdomains

def save_mesh(mesh_designation,mesh,boundaries,subdomains):
    mesh_file=File(os.environ['PATIENTDIR']+"/Results_adaptive/mesh_"+mesh_designation+".xml.gz")
    boundaries_file = File(os.environ['PATIENTDIR']+'/Results_adaptive/boundaries_'+mesh_designation+'.xml')
    subdomains_assigned_file=File(os.environ['PATIENTDIR']+'/Results_adaptive/subdomains_assigned_'+mesh_designation+'.xml')
    mesh_file<<mesh
    boundaries_file<<boundaries
    subdomains_assigned_file<<subdomains

    return True

def save_mesh_and_kappa_to_h5(mesh_to_h5,subdomains_to_h5,boundaries_to_h5,Field_calc_param):
    logging.critical("Number of mesh elements: {}".format(mesh_to_h5.num_cells()))
    #due to the glitch with the ghost model and subdomains. Used only for c-c multicontact, so always with floating
    V0_r=FunctionSpace(mesh_to_h5,'DG',0)
    kappa_r=Function(V0_r)
    default_material=Field_calc_param.default_material

    [cond_GM, perm_GM]=DielectricProperties(3).get_dielectrics(Field_calc_param.frequenc)        #3 for grey matter and so on (numeration as in voxel_data)
    [cond_WM, perm_WM]=DielectricProperties(2).get_dielectrics(Field_calc_param.frequenc)
    [cond_CSF, perm_CSF]=DielectricProperties(1).get_dielectrics(Field_calc_param.frequenc)
    [cond_default,perm_default]=DielectricProperties(default_material).get_dielectrics(Field_calc_param.frequenc)

    from GUI_inp_dict import d as d_encap
    [cond_encap, perm_encap]=DielectricProperties(d_encap['encap_tissue_type']).get_dielectrics(Field_calc_param.frequenc)
    cond_encap=cond_encap*d_encap['encap_scaling_cond']
    perm_encap=perm_encap*d_encap['encap_scaling_perm']

    k_val_r=[cond_default*0.001,cond_CSF*0.001,cond_WM*0.001,cond_GM*0.001,cond_encap*0.001,1000.0]
    help = np.asarray(subdomains_to_h5.array(), dtype=np.int32)
    kappa_r.vector()[:] = np.choose(help, k_val_r)

    hdf = HDF5File(mesh_to_h5.mpi_comm(), os.environ['PATIENTDIR']+'/Results_adaptive/Mesh_to_solve.h5', 'w')
    hdf.write(mesh_to_h5, "/mesh")
    hdf.write(subdomains_to_h5, "/subdomains")
    hdf.write(boundaries_to_h5, "/boundaries")
    hdf.write(kappa_r, "/kappa_r")

    file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_subdomains_map.pvd')
    file<<subdomains_to_h5
    file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_conductivity_map.pvd')
    file<<kappa_r

    if Field_calc_param.EQS_mode == 'EQS':
        V0_i=FunctionSpace(mesh_to_h5,'DG',0)
        kappa_i=Function(V0_i)
        omega_eps0=2*np.pi*130.0*8.854e-12           #2*pi*f*eps0
        k_val_i=[omega_eps0*perm_default*0.001,omega_eps0*perm_CSF*0.001,omega_eps0*perm_WM*0.001,omega_eps0*perm_GM*0.001,1*omega_eps0*perm_encap*0.001,1000000000*omega_eps0]
        help = np.asarray(subdomains_to_h5.array(), dtype=np.int32)
        kappa_i.vector()[:] = np.choose(help, k_val_i)

        hdf.write(kappa_i, "/kappa_i")
        file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_permittivity_map.pvd')
        file<<kappa_i

    if Field_calc_param.anisotropy == 1:

        #we get unscaled tensors and scale them with conductivity here. This approach is needed only for MPI
        c00 = MeshFunction("double", mesh_to_h5, 3, 0.0)
        c01 = MeshFunction("double", mesh_to_h5, 3, 0.0)
        c02 = MeshFunction("double", mesh_to_h5, 3, 0.0)
        c11 = MeshFunction("double", mesh_to_h5, 3, 0.0)
        c12 = MeshFunction("double", mesh_to_h5, 3, 0.0)
        c22 = MeshFunction("double", mesh_to_h5, 3, 0.0)

        hdf2 = HDF5File(mesh_to_h5.mpi_comm(), os.environ['PATIENTDIR']+"/Results_adaptive/Tensors_to_solve_num_el_"+str(mesh_to_h5.num_cells())+".h5", "r")
        hdf2.read(c00, "/c00")
        hdf2.read(c01, "/c01")
        hdf2.read(c02, "/c02")
        hdf2.read(c11, "/c11")
        hdf2.read(c12, "/c12")
        hdf2.read(c22, "/c22")
        hdf2.close()

        cell_Anis = MeshFunction('bool',mesh_to_h5,3)
        cell_Anis.set_all(False)

        for cell in cells(mesh_to_h5):

            scale_cond=k_val_r[subdomains_to_h5[cell]]

            if c00[cell]!=1.0:
                cell_Anis[cell]=True

            c00[cell]=c00[cell]*scale_cond
            c01[cell]=c01[cell]*scale_cond
            c02[cell]=c02[cell]*scale_cond
            c11[cell]=c11[cell]*scale_cond
            c12[cell]=c12[cell]*scale_cond
            c22[cell]=c22[cell]*scale_cond


        file=File(os.environ['PATIENTDIR']+'/Tensors/c00_mapped.pvd')
        file<<c00,mesh_to_h5
        file=File(os.environ['PATIENTDIR']+'/Tensors/c01_mapped.pvd')
        file<<c01,mesh_to_h5
        file=File(os.environ['PATIENTDIR']+'/Tensors/c02_mapped.pvd')
        file<<c02,mesh_to_h5
        file=File(os.environ['PATIENTDIR']+'/Tensors/c11_mapped.pvd')
        file<<c11,mesh_to_h5
        file=File(os.environ['PATIENTDIR']+'/Tensors/c12_mapped.pvd')
        file<<c12,mesh_to_h5
        file=File(os.environ['PATIENTDIR']+'/Tensors/c22_mapped.pvd')
        file<<c22,mesh_to_h5

        file=File(os.environ['PATIENTDIR']+'/Tensors/Anis_cells.pvd')
        file<<cell_Anis,mesh_to_h5

        hdf.write(c00, "/c00")
        hdf.write(c01, "/c01")
        hdf.write(c02, "/c02")
        hdf.write(c11, "/c11")
        hdf.write(c12, "/c12")
        hdf.write(c22, "/c22")

    hdf.close()

    return True

def mark_cells(ext_ground,mesh_for_ref,ref_mode,Field_r,Field_im,Field_r_new,Field_im_new,Ampl_p_new,Ampl_p,max_E_best,rel_error,phi_error,subdomains_imp,Domains):

    if ext_ground==True and ref_mode==0:
        E_max_scaler=0.0004
    else:
        E_max_scaler=0.004

    cell_ref = MeshFunction('bool',mesh_for_ref,3)
    cell_ref.set_all(False)

    ROI_affected=0      #indicator, that refinement outside of ROI, affected results in ROI
    max_div=0.0
    Div_sum=0.0
    N_points_affected=0

    affected_points=np.zeros((Ampl_p.shape[0],4),float)
    aff_point_inx=0


    for point_inx in range(Ampl_p.shape[0]):
        #if (Ampl_p_new[point_inx,3]-Ampl_p[point_inx,3])/Ampl_p_new[point_inx,3]>0.2 and Ampl_p_new[point_inx,3]>=ampl_threshold:
        if abs(Ampl_p_new[point_inx,3]-Ampl_p[point_inx,3])>phi_error:
            ROI_affected=1
            affected_points[aff_point_inx,:3]=Ampl_p_new[point_inx,:3]
            affected_points[aff_point_inx,3]=abs(Ampl_p_new[point_inx,3]-Ampl_p[point_inx,3])
            aff_point_inx=aff_point_inx+1

        Div_sum=Div_sum+abs(Ampl_p_new[point_inx,3]-Ampl_p[point_inx,3])
        N_points_affected=N_points_affected+1

        if abs(Ampl_p_new[point_inx,3]-Ampl_p[point_inx,3])>max_div:
            max_div=abs(Ampl_p_new[point_inx,3]-Ampl_p[point_inx,3])

    affected_points=affected_points[~np.all(affected_points==0.0,axis=1)]  #deletes all zero enteries


    logging.critical("Av. dev. on the neuron segments: {} V".format(Div_sum/Ampl_p.shape[0]))
    with open(os.environ['PATIENTDIR']+'/Results_adaptive/Av_dev.txt','a') as f:
        f.write(str(Div_sum/float(Ampl_p.shape[0])))
        f.write('\n')

    logging.critical("Max dev. on the neuron segments: {} V".format(max_div))
    with open(os.environ['PATIENTDIR']+'/Results_adaptive/Max_dev.txt','a') as f:
        f.write(str(max_div))
        f.write('\n')

    logging.critical("Absolute error threshold on neuron compartments: {} V".format(phi_error))

    if np.any(affected_points!=0.0):
        np.savetxt(os.environ['PATIENTDIR']+'/Results_adaptive/affected_points.csv', affected_points, delimiter=" ")

    cell_counter=0
    cell_counter_total=mesh_for_ref.num_cells()

    cell_quart=[int(cell_counter_total/4.0),int(2*cell_counter_total/4.0),int(3*cell_counter_total/4.0),int(4*cell_counter_total/4.0)]

    for cell in cells(mesh_for_ref):

        if cell_counter in cell_quart:
            logging.critical("{}% of cells were processed".format(int(cell_counter*100/cell_counter_total)+1))

        cell_counter=cell_counter+1

        if ROI_affected==0:     #no need to refine if solution on neuron compartments was not affected
            break

        if (ref_mode==0 and not(subdomains_imp[cell] in Domains.Tis_index)) or (ref_mode==1 and not(subdomains_imp[cell]==Domains.Contact_index or subdomains_imp[cell]==Domains.Float_contacts)) or (ref_mode==2 and not(subdomains_imp[cell]==Domains.ROI_index)):
            wrong_domain=1
        else:
            if ref_mode==2:     #also mark cells where affected neuron compartments
                for aff_point_inx in range(affected_points.shape[0]):
                    pnt=Point(affected_points[aff_point_inx,0],affected_points[aff_point_inx,1],affected_points[aff_point_inx,2])
                    if (cell.collides(pnt)):
                        cell_ref[cell]=True
                        break

            x_coord=cell.midpoint().x()
            y_coord=cell.midpoint().y()
            z_coord=cell.midpoint().z()
            probe=Point(x_coord,y_coord,z_coord)

            '''Compute E_ampl as sqrt(ex_r*ex_r+ey_r*ey_r*ez_r*ez_r+ex_i*ex_i+ey_i*ey_i*ez_i*ez_i)'''
            #Ampl=sqrt(Field_r(probe)[0]*Field_r(probe)[0]+Field_r(probe)[1]*Field_r(probe)[1]+Field_r(probe)[2]*Field_r(probe)[2]+Field_im(probe)[0]*Field_im(probe)[0]+Field_im(probe)[1]*Field_im(probe)[1]+Field_im(probe)[2]*Field_im(probe)[2])
            #Ampl_new=sqrt(Field_r_new(probe)[0]*Field_r_new(probe)[0]+Field_r_new(probe)[1]*Field_r_new(probe)[1]+Field_r_new(probe)[2]*Field_r_new(probe)[2]+Field_im_new(probe)[0]*Field_im_new(probe)[0]+Field_im_new(probe)[1]*Field_im_new(probe)[1]+Field_im_new(probe)[2]*Field_im_new(probe)[2])

            '''Compute Phi_ampl as sqrt(phi_r*phi_r+phi_im*phi_im)'''
            #Ampl=sqrt(Field_r(probe)*Field_r(probe)+Field_im(probe)*Field_im(probe))

            Field_old=np.array([Field_r(probe)[0],Field_r(probe)[1],Field_r(probe)[2],Field_im(probe)[0],Field_im(probe)[1],Field_im(probe)[2]])
            Field_new=np.array([Field_r_new(probe)[0],Field_r_new(probe)[1],Field_r_new(probe)[2],Field_im_new(probe)[0],Field_im_new(probe)[1],Field_im_new(probe)[2]])
            Field_diff_L2 = np.linalg.norm(Field_new[:]-Field_old[:])

            #if Field_diff_L2/np.linalg.norm(Field_new[:])>rel_error and np.linalg.norm(Field_new[:])>max_E_best*0.002:
            if Field_diff_L2/np.linalg.norm(Field_new[:])>rel_error and np.linalg.norm(Field_new[:])>max_E_best*E_max_scaler:
                cell_ref[cell]=True

    file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_Marked_cells.pvd')
    file<<cell_ref

    return cell_ref

def mark_cells_start(mesh_for_ref,ref_mode,subdomains_imp,Domains):

    cell_ref = MeshFunction('bool',mesh_for_ref,3)
    cell_ref.set_all(False)

    for cell in cells(mesh_for_ref):

        if ref_mode==0: #outside of ROI
            if subdomains_imp[cell] in Domains.Tis_index:
                cell_ref[cell]=True
        elif ref_mode==1:
            if subdomains_imp[cell]==Domains.Contact_index or subdomains_imp[cell]==Domains.Float_contacts:
                cell_ref[cell]=True   #around contacts in encap. layer
        elif ref_mode==2:
            if subdomains_imp[cell]==Domains.ROI_index:
                if ref_mode==2: #inside ROI
                    cell_ref[cell]=True

    file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Marked_cells_start_ref_mode'+str(ref_mode)+'.pvd')
    file<<cell_ref

    return cell_ref

def check_current_conv(J_real,J_imag,J_real_new,J_imag_new,rel_error_current,info_only=0):
    J_ampl=sqrt(J_real**2+J_imag**2)
    J_ampl_new=sqrt(J_real_new**2+J_imag_new**2)

    #in some cases this current is not scaled
    #print("Sum(J_ampl): ", J_ampl)
    #print("Sum(J_ampl_new): ", J_ampl_new)


    logging.critical("Current deviation: {}%".format(abs((J_ampl_new-J_ampl)/J_ampl_new)*100))
    logging.critical("Relative error threshold for current convergence: {}%".format(rel_error_current*100.0))
    if info_only==0:
        if abs(J_ampl_new-J_ampl)/J_ampl_new>rel_error_current:
            logging.critical("Current deviation is high, refining...")
            return True
        else:
            return False
    else:
        return True

def mark_cell_loc_J(ext_ground,subdomains_imp,j_real,j_im,j_real_new,j_im_new,mesh_old,mesh_new,ref_mode,ref_iteration,Domains,rel_error_current):
    #this function will find mapping between cells of the old and the new meshes and compare currents in the corresponding elements

    a=list(mesh_new.data().array("parent_cell", mesh_new.topology().dim()))
    b=list(map(int,a))

    count = MeshFunction('size_t',mesh_old, 3)
    count.set_all(0)

    for idx in b:
        count[idx] += 1		#count[cell] contains a number of cells in the new mesh

    dx = Measure("dx",domain=mesh_old)

    cl_number=-1
    J_el_old=np.zeros((mesh_old.num_cells()),dtype=float)

    for cell in cells(mesh_old):
        cl_number=cl_number+1
        J_el_old[cl_number]=sqrt(abs(assemble_local((inner(j_real,j_real)+inner(j_im,j_im))*dx,cell)))/cell.volume()

    max_normalized_current=max(J_el_old[:])

    if ext_ground==True and ref_mode==0:
        threshold_current_in_element=0.000000000001
        rel_error_current=0.65
        #threshold_current_in_element=0.0000000001
        #rel_error_current=0.5
        #threshold_current_in_element=0.00000000001
        #rel_error_current=0.25
        #threshold_current_in_element=0.0000000001    #was
        #rel_error_current=0.001
    else:
        threshold_current_in_element=0.25

    #logging.critical("threshold_current_in_element: {}".format(threshold_current_in_element))

    int_cell_ref= MeshFunction('bool',mesh_old,3)       #cell marker for mesh_old
    int_cell_ref.set_all(False)

    int_cell_ref_new = MeshFunction('bool',mesh_new,3)  #cell marker for mesh_new
    int_cell_ref_new.set_all(False)
    #print("ref_iteration: ", ref_iteration)
    cl_number=-1

    if ref_mode==0:
        Domain_index=Domains.Tis_index    #a list because rst+encap_out
    elif ref_mode==1:
        Domain_index=[Domains.Contact_index,-1]             #just to make it iterable
    elif ref_mode==2:
        Domain_index=[Domains.ROI_index,-1]

    dx_new = Measure("dx",domain=mesh_new)
    cell_counter_total=mesh_old.num_cells()
    cell_quart=[int(cell_counter_total/4.0),int(2*cell_counter_total/4.0),int(3*cell_counter_total/4.0),int(4*cell_counter_total/4.0)]

    for cell in cells(mesh_old):
        if cl_number in cell_quart:
            logging.critical("{}% of cells were processed".format(int(cl_number*100/cell_counter_total)+1))

        cl_number=cl_number+1

        if (subdomains_imp[cell] in Domain_index) and J_el_old[cl_number]>=threshold_current_in_element*max_normalized_current: #or subdomains_imp[cell]==Domain_add:
            J_el_sum=0.0
            new_cells_indices=[i for i, e in enumerate(b) if e == cl_number]
            for new_cells_index in new_cells_indices:
                J_el_sum=J_el_sum+sqrt(abs(assemble_local((inner(j_real_new,j_real_new)+inner(j_im_new,j_im_new))*dx_new,Cell(mesh_new,new_cells_index))))/cell.volume()

            if J_el_sum!=0.0:
                if abs((J_el_sum-J_el_old[cl_number])/J_el_sum)>rel_error_current:
                    if ref_iteration==1:
                        int_cell_ref[cell]=True
                    else:
                        for new_cells_index in new_cells_indices:
                            int_cell_ref_new[new_cells_index]=True

    if ref_iteration==1:
        #mesh_old_new = refine(mesh_old, int_cell_ref)
        #print("N cells in new mesh after refinement for current: ", mesh_old_new.num_cells())
        file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_marked_cell_current_ref.pvd')
        file<<int_cell_ref
        return int_cell_ref
    else:
        file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_marked_cell_current_ref_on_new_mesh.pvd')
        file<<int_cell_ref_new
        return int_cell_ref_new


def adapt_mesh(region,mesh_initial,boundaries_initial,subdomains_assigned_initial,MRI_param,DTI_param,Domains,d,cc_multicontact,num_proc,anisotropy,Field_calc_param,previous_results):

    logging.critical("----- Conducting refinement in {} -----".format(region))
    start_adapt_region=time_lib.time()

    save_mesh('initial_for_the_step',mesh_initial,boundaries_initial,subdomains_assigned_initial)

    if d["FEniCS_MPI"]==True:
        if cc_multicontact==True:
            from Math_module_floating_MPI import get_solutions,get_field_on_points
        else:
            from Math_module_MPI_only_FEniCS import get_solutions,get_field_on_points
    else:
        if cc_multicontact==True:
            from Math_module_hybrid_floating import compute_field_with_superposition,get_field_on_points
        else:
            from Math_module_hybrid import get_field,get_field_on_points


    from CSF_refinement_new import mesh_refiner

    if anisotropy==1:
        from Tissue_marking_new import get_cellmap_tensors
        subdomains=get_cellmap_tensors(mesh_initial,subdomains_assigned_initial,Domains,MRI_param,DTI_param,d["default_material"])
    else:
        from Tissue_marking_new import get_cellmap
        subdomains=get_cellmap(mesh_initial,subdomains_assigned_initial,Domains,MRI_param,d["default_material"])

    if d["current_control"]==1:     #in this case we need to check current convergence
        current_checked=0       #phi_error will be defined after phi evaluation on the initial mesh
    else:
        current_checked=1
        Phi_vector=[x for x in d["Phi_vector"] if x is not None]
        if Field_calc_param.external_grounding==True:
            Phi_vector.append(0.0)
        phi_error=abs((max(Phi_vector)-min(Phi_vector))*d["Adaptive_frac_div"])


    if region == 'it_outside_ROI':
        ref_mode = 0
    elif region == 'it_on_contact':
        ref_mode = 1
        current_checked=0           #always check current when refining around contacts
        if d["current_control"] == 0:     #in this case we need to check current convergence
            d["rel_div_current"] = 0.012
            logging.critical("Although VC mode is used, current convergence will be checked during refinement around contacts with 1.2% rel. deviation. You can change the threshold in Mesh_adaption_hybrid.py")

        if d["rel_div"]>=0.01:
            d["rel_div"] = d["rel_div"]
        else:
            d["rel_div"] = 0.01          #always fixed to 1% to avoid flickering
            logging.critical("Rel. error threshold during refinement around contacts is set to 1% (to discard flickering effect)")
    elif region == 'it_in_ROI':
        ref_mode = 2

    ref_it=1    #first iteration (yes, here we start from 1)

    if previous_results[0] == -1:    #no previous results available
        #Calculate on the initial mesh
        if d["FEniCS_MPI"]==True:
            save_mesh_and_kappa_to_h5(mesh_initial,subdomains,boundaries_initial,Field_calc_param)
            if cc_multicontact==True:
                subprocess.call(["mpirun", "-np", "{}".format(num_proc), "python3", "Math_module_floating_MPI.py"])
            else:
                subprocess.call(["mpirun", "-np", "{}".format(num_proc), "python3", "Math_module_MPI_only_FEniCS.py"])

            Phi_r,Phi_im,Field_real,Field_imag,max_E,J_r,J_im,j_dens_real,j_dens_im=get_solutions(d["EQS_core"],Field_calc_param.frequenc,Field_calc_param.element_order)
        else:
            if cc_multicontact==True:
                [Phi_r,Phi_im,Field_real,Field_imag,max_E,J_r,J_im,j_dens_real,j_dens_im]=compute_field_with_superposition(mesh_initial,Domains,subdomains_assigned_initial,subdomains,boundaries_initial,Field_calc_param)
            else:
                [Phi_r,Phi_im,Field_real,Field_imag,max_E,J_r,J_im,j_dens_real,j_dens_im]=get_field(mesh_initial,Domains,subdomains,boundaries_initial,Field_calc_param)

        Phi_amp_on_neuron =get_field_on_points(Phi_r,Phi_im,d["current_control"],J_r,J_im)   #To calculate the second norm of the Field
        previous_results = [Phi_r,Phi_im,Field_real,Field_imag,J_r,J_im,j_dens_real,j_dens_im,Phi_amp_on_neuron,max_E]
    else:
        Phi_r,Phi_im,Field_real,Field_imag,J_r,J_im,j_dens_real,j_dens_im,Phi_amp_on_neuron,max_E=previous_results[:]

    if d["current_control"] == True:
        if d["EQS_core"]=='EQS':
            #not the best approach, but both should be on the Dirichlet BCs
            max_phi_r=max(Phi_r.vector()[:])
            max_phi_im=max(Phi_im.vector()[:])

            min_phi_r=min(Phi_r.vector()[:])
            min_phi_im=min(Phi_im.vector()[:])

            phi_error=abs((np.sqrt((max_phi_r-min_phi_r)**2+(max_phi_im-min_phi_im)**2))*d["Adaptive_frac_div"])   #should be scaled
        else:
            phi_error=abs((max(Phi_r.vector()[:])-min(Phi_r.vector()[:]))*d["Adaptive_frac_div"])   #should be scaled

    logging.critical("Absolute error threshold on neuron compartments: {} V".format(phi_error))

    # mesh refined uniformly in the specified region by ref_mode
    logging.critical("\n--- Initial uniform refinement step for {}".format(region))
    cells_ref=mark_cells_start(mesh_initial,ref_mode,subdomains_assigned_initial, Domains)
    [mesh_new,boundaries_new,subdomains_assigned_new]=mesh_refiner(mesh_initial,boundaries_initial,subdomains_assigned_initial,cells_ref,Domains,cc_multicontact)

    #here we start to refine
    while (cells_ref.where_equal(True)):        # if True, then there are some cells marked for refinement

        if anisotropy==1:
            subdomains_new=get_cellmap_tensors(mesh_new,subdomains_assigned_new,Domains,MRI_param,DTI_param,d["default_material"])
        else:
            subdomains_new=get_cellmap(mesh_new,subdomains_assigned_new,Domains,MRI_param,d["default_material"])

        if d["FEniCS_MPI"]==True:
            save_mesh_and_kappa_to_h5(mesh_new,subdomains_new,boundaries_new,Field_calc_param)
            if cc_multicontact == True:
                subprocess.call(["mpirun", "-np", "{}".format(num_proc), "python3", "Math_module_floating_MPI.py"])
            else:
                subprocess.call(["mpirun", "-np", "{}".format(num_proc), "python3", "Math_module_MPI_only_FEniCS.py"])
            Phi_r_new,Phi_im_new,Field_real_new,Field_imag_new,max_E_new,J_r_new,J_im_new,j_dens_real_new,j_dens_im_new=get_solutions(d["EQS_core"],Field_calc_param.frequenc,Field_calc_param.element_order)
        else:
            if cc_multicontact == True:
                [Phi_r_new,Phi_im_new,Field_real_new,Field_imag_new,max_E_new,J_r_new,J_im_new,j_dens_real_new,j_dens_im_new]=compute_field_with_superposition(mesh_new,Domains,subdomains_assigned_new,subdomains_new,boundaries_new,Field_calc_param)
            else:
                [Phi_r_new,Phi_im_new,Field_real_new,Field_imag_new,max_E_new,J_r_new,J_im_new,j_dens_real_new,j_dens_im_new]=get_field(mesh_new,Domains,subdomains_new,boundaries_new,Field_calc_param)



        Phi_amp_on_neuron_new=get_field_on_points(Phi_r_new,Phi_im_new,d["current_control"],J_r_new,J_im_new)


        #check E-field convergence first (always, even if only current did not converge)
        if ref_it==1:       #on the first iteration we will mark cells on the initial mesh
            cells_ref=mark_cells(Field_calc_param.external_grounding,mesh_initial,ref_mode,Field_real,Field_imag,Field_real_new,Field_imag_new,Phi_amp_on_neuron_new,Phi_amp_on_neuron,max_E_new,d["rel_div"],phi_error,subdomains_assigned_initial,Domains)
        else:
            cells_ref=mark_cells(Field_calc_param.external_grounding,mesh_new,ref_mode,Field_real,Field_imag,Field_real_new,Field_imag_new,Phi_amp_on_neuron_new,Phi_amp_on_neuron,max_E_new,d["rel_div"],phi_error,subdomains_assigned_new,Domains)

        if (cells_ref.where_equal(True)):
            ref_due_to_phi_dev=1

        ###### check current convergence if necessary #####
        if current_checked==0 and not(cells_ref.where_equal(True)):
            logging.critical("Checking current convergence for {}".format(region))
            if check_current_conv(J_r,J_im,J_r_new,J_im_new,d["rel_div_current"]):           # True if current deviation above rel_div_current
                if ref_it==1:       #marking will be on the initial mesh
                    cells_ref=mark_cell_loc_J(Field_calc_param.external_grounding,subdomains_assigned_initial,j_dens_real,j_dens_im,j_dens_real_new,j_dens_im_new,mesh_initial,mesh_new,ref_mode,1,Domains,d["rel_div_current"])
                elif ref_it==2:
                    cells_ref=mark_cell_loc_J(Field_calc_param.external_grounding,subdomains_assigned_initial,j_dens_real,j_dens_im,j_dens_real_new,j_dens_im_new,mesh_initial,mesh_new,ref_mode,2,Domains,d["rel_div_current"])
                else:
                    mesh_old = Mesh(os.environ['PATIENTDIR']+"/Results_adaptive/mesh_adapt.xml.gz")
                    subdomains_assigned_old = MeshFunction('size_t',mesh_old,os.environ['PATIENTDIR']+'/Results_adaptive/subdomains_assigned_adapt.xml')
                    if ref_due_to_phi_dev==0:   #if previous refinement due to current, refine on mesh_new
                        cells_ref=mark_cell_loc_J(Field_calc_param.external_grounding,subdomains_assigned_old,j_dens_real,j_dens_im,j_dens_real_new,j_dens_im_new,mesh_old,mesh_new,ref_mode,2,Domains,d["rel_div_current"])
                    else:                       #else, refine on mesh
                        cells_ref=mark_cell_loc_J(Field_calc_param.external_grounding,subdomains_assigned_old,j_dens_real,j_dens_im,j_dens_real_new,j_dens_im_new,mesh_old,mesh_new,ref_mode,1,Domains,d["rel_div_current"])

                if not (cells_ref.where_equal(True)):
                    current_checked=1
                    logging.critical("Current did not converge, but nothing marked to refine! Consider decreasing threshold_current_in_element in Mesh_adaptation.py")
                else:
                    if ref_it>2:
                        if ref_due_to_phi_dev==1:   #if this is a refinement due to current right after phi_dev check, refine on mesh (resaved as mesh_new). Otherwise, on mesh_new
                            Phi_r_new,Phi_im_new,J_r_new,J_im_new,j_dens_real_new,j_dens_im_new,Field_real_new,Field_imag_new,Phi_amp_on_neuron_new,max_E_new=(Phi_r,Phi_im,J_r,J_im,j_dens_real,j_dens_im,Field_real,Field_imag,Phi_amp_on_neuron,max_E)
                            del mesh_new,boundaries_new,subdomains_assigned_new
                            mesh_new,boundaries_new,subdomains_assigned_new=load_mesh('adapt')

                    ref_due_to_phi_dev=0   #the refinement will be conducted due to current deviation

            else:
                current_checked=1
                logging.critical("Current converged\n")


        if not(cells_ref.where_equal(True)) and ref_it==2:       # written on the envelope.   if we refine adaptively only once, we want to evaluate mesh that we just refined. Otherwise, we will evaluate already saved mesh
            save_mesh('adapt_it2',mesh_new,boundaries_new,subdomains_assigned_new)
            Phi_r,Phi_im,J_r,J_im,j_dens_real,j_dens_im,Field_real,Field_imag,Phi_amp_on_neuron,max_E=(Phi_r_new,Phi_im_new,J_r_new,J_im_new,j_dens_real_new,j_dens_im_new,Field_real_new,Field_imag_new,Phi_amp_on_neuron_new,max_E_new)
            previous_results = [Phi_r_new,Phi_im_new,Field_real_new,Field_imag_new,J_r_new,J_im_new,j_dens_real_new,j_dens_im_new,Phi_amp_on_neuron_new,max_E_new]

        # if we got local convergence, let us check with uniform refinement again and mark cells, where high deviation occurs
        if not (cells_ref.where_equal(True)) and ref_it>1:
            logging.critical("\n--- Uniform refinement step for {}".format(region))
            if ref_it==2:
                mesh,boundaries,subdomains_assigned=load_mesh('adapt_it2')
            else:
                mesh,boundaries,subdomains_assigned=load_mesh('adapt')

            #uniform refinement
            cells_ref=mark_cells_start(mesh,ref_mode,subdomains_assigned, Domains)
            [mesh_uni,boundaries_uni,subdomains_assigned_uni]=mesh_refiner(mesh,boundaries,subdomains_assigned,cells_ref,Domains,cc_multicontact)

            if not(os.path.isfile(os.environ['PATIENTDIR']+'/Results_adaptive/cells_to_ref_after_uni_'+str(mesh_uni.num_cells())+'.h5')):     #check whether calculations for this uniform refinement were already conducted
                if anisotropy==1:
                    subdomains_uni=get_cellmap_tensors(mesh_uni,subdomains_assigned_uni,Domains,MRI_param,DTI_param,d["default_material"])
                else:
                    subdomains_uni=get_cellmap(mesh_uni,subdomains_assigned_uni,Domains,MRI_param,d["default_material"])

                if d["FEniCS_MPI"]==True:
                    save_mesh_and_kappa_to_h5(mesh_uni,subdomains_uni,boundaries_uni,Field_calc_param)
                    if cc_multicontact == True:
                        subprocess.call(["mpirun", "-np", "{}".format(num_proc), "python3", "Math_module_floating_MPI.py"])
                    else:
                        subprocess.call(["mpirun", "-np", "{}".format(num_proc), "python3", "Math_module_MPI_only_FEniCS.py"])
                    [Phi_r_uni,Phi_im_uni,Field_real_uni,Field_imag_uni,max_E_uni,J_r_uni,J_im_uni,j_dens_real_uni,j_dens_im_uni]=get_solutions(d["EQS_core"],Field_calc_param.frequenc,Field_calc_param.element_order)
                else:
                    if cc_multicontact == True:
                        [Phi_r_uni,Phi_im_uni,Field_real_uni,Field_imag_uni,max_E_uni,J_r_uni,J_im_uni,j_dens_real_uni,j_dens_im_uni]=compute_field_with_superposition(mesh_uni,Domains,subdomains_assigned_uni,subdomains_uni,boundaries_uni,Field_calc_param)
                    else:
                        [Phi_r_uni,Phi_im_uni,Field_real_uni,Field_imag_uni,max_E_uni,J_r_uni,J_im_uni,j_dens_real_uni,j_dens_im_uni]=get_field(mesh_uni,Domains,subdomains_uni,boundaries_uni,Field_calc_param)

                Phi_amp_on_neuron_uni=get_field_on_points(Phi_r_uni,Phi_im_uni,d["current_control"],J_r_uni,J_im_uni)
                cells_ref=mark_cells(Field_calc_param.external_grounding,mesh,ref_mode,Field_real,Field_imag,Field_real_uni,Field_imag_uni,Phi_amp_on_neuron_uni,Phi_amp_on_neuron,max_E_uni,d["rel_div"],phi_error,subdomains_assigned,Domains)

                hdf = HDF5File(mesh.mpi_comm(), os.environ['PATIENTDIR']+'/Results_adaptive/cells_to_ref_after_uni_'+str(mesh_uni.num_cells())+'.h5', 'w')
                hdf.write(cells_ref, "/cells_ref_after_uni")
                hdf.close()
                np.savetxt(os.environ['PATIENTDIR']+'/Results_adaptive/current_in_uni_ref_'+str(mesh_uni.num_cells())+'.csv', np.array([J_r_uni,J_im_uni]), delimiter=" ")
            else:       #the field was already computed for this uniform refinement
                logging.critical("Cells for refinement after uniform check were loaded from the previous iteration\n")
                hdf = HDF5File(mesh.mpi_comm(), os.environ['PATIENTDIR']+'/Results_adaptive/cells_to_ref_after_uni_'+str(mesh_uni.num_cells())+'.h5', 'r')
                cells_ref = MeshFunction('bool', mesh,3)
                hdf.read(cells_ref, "/cells_ref_after_uni")
                hdf.close()

                [J_r_uni,J_im_uni]=np.genfromtxt(os.environ['PATIENTDIR']+'/Results_adaptive/current_in_uni_ref_'+str(mesh_uni.num_cells())+'.csv', delimiter=' ')

            mesh_check=refine(mesh, cells_ref)   #we need to check whether it leads to mesh_new
            num_cells_mesh_check=mesh_check.num_cells()

            file=File(os.environ['PATIENTDIR']+'/Results_adaptive/cells_ref_after_uni.pvd')
            file<<cells_ref

            del mesh_check

            if mesh_new.num_cells()<num_cells_mesh_check or (ref_due_to_phi_dev==0 and mesh_new.num_cells()!=num_cells_mesh_check):   #use mesh to refine further (for this you will have to resave it as mesh_new)
                Phi_r_new,Phi_im_new,J_r_new,J_im_new,j_dens_real_new,j_dens_im_new,Field_real_new,Field_imag_new,Phi_amp_on_neuron_new,max_E_new=(Phi_r,Phi_im,J_r,J_im,j_dens_real,j_dens_im,Field_real,Field_imag,Phi_amp_on_neuron,max_E)
                del mesh_new,boundaries_new,subdomains_assigned_new
                hdf = HDF5File(mesh.mpi_comm(), os.environ['PATIENTDIR']+'/Results_adaptive/cells_to_ref_after_uni_'+str(mesh_uni.num_cells())+'.h5', 'r')
                cells_ref = MeshFunction('bool', mesh,3)
                hdf.read(cells_ref, "/cells_ref_after_uni")
                hdf.close()
                if ref_it==2:
                    mesh_new,boundaries_new,subdomains_assigned_new=load_mesh('adapt_it2')
                else:
                    mesh_new,boundaries_new,subdomains_assigned_new=load_mesh('adapt')
            elif mesh_new.num_cells()==num_cells_mesh_check and cells_ref.where_equal(True):   #if mesh_uni and mesh_new are the same after current-based ref, switch status to phi-based
                ref_due_to_phi_dev=1

            if check_current_conv(J_r,J_im,J_r_uni,J_im_uni,d["rel_div_current"]):
                logging.critical("Current ref. will be conducted only if deviation on neuron compartments was high")
                current_checked=0
            else:
                if ref_mode==1:   #
                    if cells_ref.where_equal(True):
                        logging.critical("After uniform check, the current converged, the iteration is complete (only when refining around the contacts)")
                        cells_ref.set_all(False)

            ###
            #if uniform refinement on mesh revealed deviation, but number of marked cells is <= than in mesh_new, then uniformly refine mesh_new. Do this step only if the last local refinement was due to phi ref. If current ref, refine on mesh (no danger of repetition as we have a different criterion)
            if (cells_ref.where_equal(True) and ref_it!=2 and (mesh_new.num_cells()>=num_cells_mesh_check) and ref_due_to_phi_dev==1):       #if mesh and mesh_uni have a high deviation in the solution, then create new mesh_uni from mesh_new using uniform refinement (otherwise it might stuck in the loop mesh-mesh_uni-mesh_new-mesh)
                cells_ref=mark_cells_start(mesh_new,ref_mode,subdomains_assigned_new, Domains)
                logging.critical("Deviation is high, now uniformly refining mesh_new")
                save_mesh('adapt',mesh_new,boundaries_new,subdomains_assigned_new)

                [mesh_uni,boundaries_uni,subdomains_assigned_uni]=mesh_refiner(mesh_new,boundaries_new,subdomains_assigned_new,cells_ref,Domains,cc_multicontact)

                if anisotropy==1:
                    subdomains_uni=get_cellmap_tensors(mesh_uni,subdomains_assigned_uni,Domains,MRI_param,DTI_param,d["default_material"])
                else:
                    subdomains_uni=get_cellmap(mesh_uni,subdomains_assigned_uni,Domains,MRI_param,d["default_material"])

                if d["FEniCS_MPI"]==True:
                    save_mesh_and_kappa_to_h5(mesh_uni,subdomains_uni,boundaries_uni,Field_calc_param)
                    if cc_multicontact==True:
                        subprocess.call(["mpirun", "-np", "{}".format(num_proc), "python3", "Math_module_floating_MPI.py"])
                    else:
                        subprocess.call(["mpirun", "-np", "{}".format(num_proc), "python3", "Math_module_MPI_only_FEniCS.py"])
                    [Phi_r_uni,Phi_im_uni,Field_real_uni,Field_imag_uni,max_E_uni,J_r_uni,J_im_uni,j_dens_real_uni,j_dens_im_uni]=get_solutions(d["EQS_core"],Field_calc_param.frequenc,Field_calc_param.element_order)
                else:
                    if cc_multicontact==True:
                        [Phi_r_uni,Phi_im_uni,Field_real_uni,Field_imag_uni,max_E_uni,J_r_uni,J_im_uni,j_dens_real_uni,j_dens_im_uni]=compute_field_with_superposition(mesh_uni,Domains,subdomains_assigned_uni,subdomains_uni,boundaries_uni,Field_calc_param)
                    else:
                        [Phi_r_uni,Phi_im_uni,Field_real_uni,Field_imag_uni,max_E_uni,J_r_uni,J_im_uni,j_dens_real_uni,j_dens_im_uni]=get_field(mesh_uni,Domains,subdomains_uni,boundaries_uni,Field_calc_param)

                Phi_amp_on_neuron_uni=get_field_on_points(Phi_r_uni,Phi_im_uni,d["current_control"],J_r_uni,J_im_uni)

                #for stability
                del mesh_new,boundaries_new,subdomains_assigned_new
                mesh_new,boundaries_new,subdomains_assigned_new=load_mesh('adapt')

                cells_ref=mark_cells(Field_calc_param.external_grounding,mesh_new,ref_mode,Field_real_new,Field_imag_new,Field_real_uni,Field_imag_uni,Phi_amp_on_neuron_uni,Phi_amp_on_neuron_new,max_E_uni,d["rel_div"],phi_error,subdomains_assigned_new,Domains)

                if d["current_control"]==1 or ref_mode==1:
                    if check_current_conv(J_r_new,J_im_new,J_r_uni,J_im_uni,d["rel_div_current"]):
                        logging.critical("Current ref. will be conducted only if deviation on neuron compartments was high")
                        current_checked=0

                #save field solution
                hdf = HDF5File(mesh.mpi_comm(), os.environ['PATIENTDIR']+'/Results_adaptive/cells_to_ref_after_uni_'+str(mesh_uni.num_cells())+'.h5', 'w')
                hdf.write(cells_ref, "/cells_ref_after_uni")
                hdf.close()
                np.savetxt(os.environ['PATIENTDIR']+'/Results_adaptive/current_in_uni_ref_'+str(mesh_uni.num_cells())+'.csv', np.array([J_r_uni,J_im_uni]), delimiter=" ")

                if not (cells_ref.where_equal(True)):
                    previous_results = [Phi_r_new,Phi_im_new,Field_real_new,Field_imag_new,J_r_new,J_im_new,j_dens_real_new,j_dens_im_new,Phi_amp_on_neuron_new,max_E_new]



        #exiting from here if no cells marked, otherwise a new iteration
        if not (cells_ref.where_equal(True)):
            if ref_it==1:
                logging.critical("mesh meets the requirements (or no cells were marked for refinement) in {} was completed\n".format(region))
                save_mesh(region,mesh_initial,boundaries_initial,subdomains_assigned_initial)
                if ref_mode == 2:
                    save_mesh('adapt',mesh_initial,boundaries_initial,subdomains_assigned_initial)
            elif ref_it==2:
                logging.critical("mesh_new meets the requirements (or no cells were marked for refinement) in {} was completed\n".format(region))
                mesh_2it,boundaries_2it,subdomains_assigned_2it=load_mesh('adapt_it2')
                save_mesh(region,mesh_2it,boundaries_2it,subdomains_assigned_2it)
                if ref_mode == 2:
                    save_mesh('adapt',mesh_initial,boundaries_initial,subdomains_assigned_initial)
            else:
                logging.critical("mesh_new meets the requirements (or no cells were marked for refinement) in {} was completed\n".format(region))
                mesh_adapted,boundaries_adapted,subdomains_assigned_adapted=load_mesh('adapt')
                save_mesh(region,mesh_adapted,boundaries_adapted,subdomains_assigned_adapted)

        else:
            logging.critical("--- Local refinement step for {}".format(region))
            if ref_it==1:
                mesh_initial,boundaries_initial,subdomains_assigned_initial=load_mesh('initial_for_the_step')   #for stability
                [mesh_new,boundaries_new,subdomains_assigned_new]=mesh_refiner(mesh_initial,boundaries_initial,subdomains_assigned_initial,cells_ref,Domains,cc_multicontact)
                logging.critical("mesh size after adaptive refinement: {}".format(mesh_new.num_cells()))
                ref_it=ref_it+1
            else:
                Phi_r,Phi_im,J_r,J_im,j_dens_real,j_dens_im,Field_real,Field_imag,Phi_amp_on_neuron,max_E=(Phi_r_new,Phi_im_new,J_r_new,J_im_new,j_dens_real_new,j_dens_im_new,Field_real_new,Field_imag_new,Phi_amp_on_neuron_new,max_E_new)
                previous_results = [Phi_r_new,Phi_im_new,Field_real_new,Field_imag_new,J_r_new,J_im_new,j_dens_real_new,j_dens_im_new,Phi_amp_on_neuron_new,max_E_new]
                ref_it=ref_it+1
                save_mesh('adapt',mesh_new,boundaries_new,subdomains_assigned_new)
                [mesh_new,boundaries_new,subdomains_assigned_new]=mesh_refiner(mesh_new,boundaries_new,subdomains_assigned_new,cells_ref,Domains,cc_multicontact)

    minutes=int((time_lib.time() - start_adapt_region)/60)
    secnds=int(time_lib.time() - start_adapt_region)-minutes*60
    logging.critical("----- Mesh in {} was adapted in {} min {} sec -----\n".format(region,minutes, secnds))

    status=1        #for now always 1

    return previous_results,status


def mesh_refiner(mesh_old,boundaries,subdomains_assigned,cell_markers,Domains,cc_multicontact):
    parameters['linear_algebra_backend']='PETSc'
    parameters["refinement_algorithm"] = "plaza_with_parent_facets"
    parameters["allow_extrapolation"] = True;

    facets_old = MeshFunction('size_t',mesh_old,2)
    facets_old.set_all(0)

    facets_old.array()[boundaries.array() == Domains.Active_contacts[0]]=1

    if cc_multicontact==True and Domains.Amp_vector[0]!=0.0:            #because ground contact is always subtracted from the mesh
        dsSSS=Measure("dS",domain=mesh_old,subdomain_data=facets_old)
        An_surface_size_old=assemble(1.0*dsSSS(1))
    else:
        dss=Measure("ds",domain=mesh_old,subdomain_data=facets_old)
        An_surface_size_old=assemble(1.0*dss(1))

    mesh_new = refine(mesh_old, cell_markers)
    subdomains_assigned_new=adapt(subdomains_assigned,mesh_new)
    boundaries_new = adapt(boundaries,mesh_new) # put function space

    facets = MeshFunction('size_t',mesh_new,2)
    facets.set_all(0)
    facets.array()[boundaries_new.array()== Domains.Active_contacts[0]]=1
    #facets.array()[boundaries_new.array()==Domains.Active_contacts[1]]=2

    if cc_multicontact==True and Domains.Amp_vector[0]!=0.0:
        dsS_new=Measure("dS",domain=mesh_new,subdomain_data=facets)
        An_surface_size_new=assemble(1.0*dsS_new(1))
    else:
        dss_new=Measure("ds",domain=mesh_new,subdomain_data=facets)
        An_surface_size_new=assemble(1.0*dss_new(1))


    if (An_surface_size_new-An_surface_size_old)/An_surface_size_new>0.02:
        #print((An_surface_size_new-An_surface_size_old)/An_surface_size_new)
        logging.critical("Refinement broke the imposed B.C.!")
        exit()

    return (mesh_new,boundaries_new,subdomains_assigned_new)


def stupid_prerefiner(mesh_old,boundaries,subdomains_assigned):
    parameters['linear_algebra_backend']='PETSc'
    parameters["refinement_algorithm"] = "plaza_with_parent_facets"
    parameters["allow_extrapolation"] = True;

    sub_domains = MeshFunction("bool", mesh_old, mesh_old.topology().dim())
    sub_domains.set_all(0)
    inflow = Ref_box()
    inflow.mark(sub_domains, True)

    file=File(os.environ['PATIENTDIR']+'/Results_adaptive/sub_domains.pvd')
    file<<sub_domains

    mesh_new = refine(mesh_old, sub_domains)
    subdomains_assigned_new=adapt(subdomains_assigned,mesh_new)
    boundaries_new = adapt(boundaries,mesh_new) # put function s

    return mesh_new,boundaries_new,subdomains_assigned_new


def mesh_adapter(MRI_param,DTI_param,Scaling,Domains,d,anisotropy,cc_multicontact,ref_freqs):
    logging.critical("----- Conducting mesh convergence study -----")
    start_adapt=time_lib.time()

    for i in range(len(ref_freqs)):     # go over the refinement frequencies
        logging.critical("At frequency: {} Hz".format(ref_freqs[i]))
        if i == 0:    # load mesh after the CSF refinement
            mesh = Mesh(os.environ['PATIENTDIR']+'/CSF_ref/mesh_adapt_CSF'+str(Scaling)+'.xml.gz')
            boundaries = MeshFunction('size_t',mesh,os.environ['PATIENTDIR']+'/CSF_ref/boundaries_adapt_CSF'+str(Scaling)+'.xml')
            subdomains_assigned = MeshFunction('size_t',mesh,os.environ['PATIENTDIR']+'/CSF_ref/subdomains_assigned_adapt_CSF'+str(Scaling)+'.xml')
            save_mesh('for_ref',mesh,boundaries,subdomains_assigned)
        else:       # load mesh after the adaptive ref. at the previous frequency
            mesh,boundaries,subdomains_assigned=load_mesh('adapt')          #the output accepted mesh is always saved with 'adapt'

        Field_calc_param=Field_calc_parameters(d["default_material"],d["el_order"],anisotropy,d["current_control"],d["CPE_activ"],ref_freqs[i],d["EQS_core"],d["external_grounding"],d["Solver_Type"])
        if d['FEniCS_MPI']==True:
            with open(os.environ['PATIENTDIR']+'/Results_adaptive/Field_calc_param.file', "wb") as f:
                pickle.dump(Field_calc_param, f, pickle.HIGHEST_PROTOCOL)

        ref_regions=['it_outside_ROI','it_on_contact','it_in_ROI']
        regions_to_skip=[None]          # important: 'it_in_ROI' should not be skipped
        convergence_status=[]   # 1 (converged), -1 (did not converge, but no cells were marked for refinement), 0 (region was skipped)
        previous_results=[-1]     # here we will store Field and Current functions from previous iteration so that we do not need to recalculate. -1 - no previous results available

        for region in ref_regions:
            if region in regions_to_skip:
                logging.critical("{} was skipped. Make sure you are working with the desired mesh (i.e. initial or after some ref. iterations)".format(region))
                #mesh,boundaries,subdomain=load_mesh('adapt_CSF'+str(Scaling))     #loading the initial mesh here
                convergence_status.append(0)
            else:       #also saves the mesh with the name of the refined region. The last iteration will also save it as mesh_adapt.xml.gz
                previous_results,convergence = adapt_mesh(region,mesh,boundaries,subdomains_assigned,MRI_param,DTI_param,Domains,d,cc_multicontact,d['number_of_processors'],anisotropy,Field_calc_param,previous_results)
                convergence_status.append(convergence)
                #convergence here can be only 1 (converged) or -1 (did not converge, but no cells were marked for refinement)

                #reloading mesh to ensure stability using the name of the last refined region
                mesh,boundaries,subdomains_assigned=load_mesh(region)

    logging.critical("Mesh was adapted, stored as mesh_adapt.xml.gz in Results_adaptive/")
    mesh_fin = Mesh(os.environ['PATIENTDIR']+"/Results_adaptive/mesh_it_in_ROI.xml.gz")
    logging.critical("number of elements: {}".format(mesh_fin.num_cells()))
    Field_real = previous_results[0]
    file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Adapted_Field_real.pvd')
    file<<Field_real

    minutes=int((time_lib.time() - start_adapt)/60)
    secnds=int(time_lib.time() - start_adapt)-minutes*60
    logging.critical("----- Mesh adaptation took {} min {} sec -----\n".format(minutes, secnds))

    return previous_results[-2]     #Ampl of the potential on the neuron compartments
