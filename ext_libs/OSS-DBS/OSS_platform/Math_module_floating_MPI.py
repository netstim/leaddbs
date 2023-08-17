# -*- coding: utf-8 -*-
"""
Created on Sun Aug 19 13:23:54 2018

@author: butenko
"""
import warnings
warnings.filterwarnings("ignore")

from dolfin import *
from pandas import read_csv
import numpy as np
import os
import subprocess
import pickle
import time as tm

from tissue_dielectrics import DielectricProperties

parameters['linear_algebra_backend']='PETSc'
set_log_active(False)   #turns off debugging info

parameters['ghost_mode'] = 'shared_vertex'

def get_CPE_impedance(beta,K_Area,freq):
    Z_CPE=K_Area/((1j*2*np.pi*freq)**(beta))
    return Z_CPE


def get_tissue_impedance(V_drop,J_r_contact,J_im_contact,C=0):         #adapted from FanPy

    z_tis = V_drop/(J_r_contact+1j*J_im_contact)
    return z_tis

def get_solutions(EQS_form,frequency,el_order):

    start_reassamble=tm.time()

    mesh_sol = Mesh()
    f = HDF5File(mesh_sol.mpi_comm(),os.environ['PATIENTDIR']+"/Results_adaptive/Solution_"+str(np.round(frequency,6))+".h5",'r')
    f.read(mesh_sol,"mesh_sol", False)

    if EQS_form == 'EQS':
        Er = FiniteElement("Lagrange", mesh_sol.ufl_cell(),el_order)
        Ei = FiniteElement("Lagrange", mesh_sol.ufl_cell(),el_order)
        Ec = Er * Ei
        V = FunctionSpace(mesh_sol, Ec)
        phi_sol=Function(V)
        f.read(phi_sol,'solution_phi_full')
        phi_r_sol,phi_i_sol=phi_sol.split(deepcopy=True)

        if el_order>1:
            W =VectorFunctionSpace(mesh_sol,'DG',el_order-1)
            W_i =VectorFunctionSpace(mesh_sol,'DG',el_order-1)
            V_normE=FunctionSpace(mesh_sol,"CG",el_order-1)
        else:
            W =VectorFunctionSpace(mesh_sol,'DG',el_order)
            W_i =VectorFunctionSpace(mesh_sol,'DG',el_order)
            V_normE=FunctionSpace(mesh_sol,"CG",el_order)

        E_field = Function(W)
        E_field_im = Function(W_i)
        f.read(E_field,'solution_E_field')
        f.read(E_field_im,'solution_E_field_im')

        j_dens_real = Function(W)
        j_dens_im = Function(W_i)
        f.read(j_dens_real, "solution_j_real")
        f.read(j_dens_im, "solution_j_im")

        J_Vector=PETScVector(MPI.comm_world,2)
        f.read(J_Vector, "J_Vector",False)
        J_real,J_im=J_Vector[:]

        E_norm=project(sqrt(inner(E_field,E_field)+inner(E_field_im,E_field_im)),V_normE,solver_type="cg", preconditioner_type="amg")
        max_E=E_norm.vector().max()

        file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Phi_r_field_EQS.pvd')
        file<<phi_r_sol,mesh_sol
        file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Phi_im_field_EQS.pvd')
        file<<phi_i_sol,mesh_sol
        file=File(os.environ['PATIENTDIR']+'/Results_adaptive/E_norm_EQS.pvd')
        file<<E_norm,mesh_sol
    elif EQS_form == 'QS':
        V = FunctionSpace(mesh_sol, "Lagrange",el_order)
        phi_r_sol=Function(V)
        phi_i_sol=Function(V)
        f.read(phi_r_sol,'solution_phi_full')
        phi_i_sol.vector()[:]=0.0

        if el_order>1:
            W =VectorFunctionSpace(mesh_sol,'DG',el_order-1)
            V_normE=FunctionSpace(mesh_sol,"CG",el_order-1)
        else:
            W =VectorFunctionSpace(mesh_sol,'DG',el_order)
            V_normE=FunctionSpace(mesh_sol,"CG",el_order)

        E_field = Function(W)
        f.read(E_field,'solution_E_field')
        E_field_im=Function(W)
        E_field_im.vector()[:] = 0.0        #fake

        j_dens_real = Function(W)
        j_dens_im = Function(W)
        f.read(j_dens_real, "solution_j_real")
        j_dens_im.vector()[:] = 0.0        #fake

        J_Vector=PETScVector(MPI.comm_world,2)
        f.read(J_Vector, "J_Vector",False)
        J_real,J_im=J_Vector[:]

        E_norm=project(sqrt(inner(E_field,E_field)),V_normE,solver_type="cg", preconditioner_type="amg")
        max_E=E_norm.vector().max()
        file=File(os.environ['PATIENTDIR']+'/Results_adaptive/E_norm_QS.pvd')
        file<<E_norm,mesh_sol

        file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Phi_r_field_QS.pvd')
        file<<phi_r_sol,mesh_sol

    f.close()

    #if we want to get the potential magnitude on the neuron compartments
    Vertices_get=read_csv(os.environ['PATIENTDIR']+'/Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', delimiter=' ', header=None)
    Vertices_array=Vertices_get.values
    Phi_ROI=np.zeros((Vertices_array.shape[0],4),float)

    for inx in range(Vertices_array.shape[0]):
        pnt=Point(Vertices_array[inx,0],Vertices_array[inx,1],Vertices_array[inx,2])

        Phi_ROI[inx,0]=Vertices_array[inx,0]
        Phi_ROI[inx,1]=Vertices_array[inx,1]
        Phi_ROI[inx,2]=Vertices_array[inx,2]

        Phi_ROI[inx,3]=np.sqrt(phi_r_sol(pnt)*phi_r_sol(pnt)+phi_i_sol(pnt)*phi_i_sol(pnt))

    np.savetxt(os.environ['PATIENTDIR']+'/Results_adaptive/Phi_'+str(frequency)+'.csv',  Phi_ROI, delimiter=" ")


    print("Quasi impedance (to check current convergence): ",J_real,J_im)

    minutes=int((tm.time() - start_reassamble)/60)
    secnds=int(tm.time() - start_reassamble)-minutes*60
    print("--- solution reassambled in ",minutes," min ",secnds," s ---")

    return phi_r_sol,phi_i_sol,E_field,E_field_im,max_E,J_real,J_im,j_dens_real,j_dens_im

def get_field_on_points(phi_r,phi_i,c_c,J_r,J_i):

    Vertices_neur_get=read_csv(os.environ['PATIENTDIR']+'/Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', delimiter=' ', header=None)
    Vertices_neur=Vertices_neur_get.values

    Ampl_ROI=np.zeros((Vertices_neur.shape[0],4),float)


    for inx in range(Vertices_neur.shape[0]):
        pnt=Point(Vertices_neur[inx,0],Vertices_neur[inx,1],Vertices_neur[inx,2])

        Ampl_ROI[inx,3]=sqrt(phi_r(pnt)*phi_r(pnt)+phi_i(pnt)*phi_i(pnt))
        Ampl_ROI[inx,0]=Vertices_neur[inx,0]
        Ampl_ROI[inx,1]=Vertices_neur[inx,1]
        Ampl_ROI[inx,2]=Vertices_neur[inx,2]

    Ampl_ROI=Ampl_ROI[~np.all(Ampl_ROI==0.0,axis=1)]  #deletes all zero enteries

    return Ampl_ROI



def compute_field_with_superposition(mesh_sol,Domains,subdomains,boundaries_sol,kappa,Field_calc_param):

    start_math=tm.time()

    from Math_module_hybrid import choose_solver_for_me
    if Field_calc_param.Solver_type=='Default':
        Solver_type=choose_solver_for_me(Field_calc_param.EQS_mode,Domains.Float_contacts)    #choses solver basing on the Laplace formulation and whether the floating conductors are used
    else:
        Solver_type=Field_calc_param.Solver_type      # just get the solver directly
    #IMPORTANT: for get_field_with_floats when solving EQS we always use direct solver MUMPS for stability issues (multiple floating conductors)

    if Field_calc_param.element_order==1 and MPI.comm_world.rank==1:
        print("Selected element_order (1st) is too low for current-controlled stimulation, increasing to 2nd")
        Field_calc_param.element_order=2

    if MPI.comm_world.rank==1:
        #print("Computing field with superposition on mesh with ",mesh_sol.num_cells(), " elements")
        print(len(Domains.Amp_vector)," computations are required for the iteration")

    contacts_with_current=[x for x in Domains.Amp_vector if x != 0.0]       #0.0 are grounded contacts
    phi_r_floating=np.zeros((len(contacts_with_current),len(contacts_with_current)-1),float)       #stores real potential field in the virtual floating contacts (where current is actually assigned)
    J_real_current_contacts=np.zeros(len(contacts_with_current),float)                  #currents computed on the contacts when we solve "one active contact vs ground" system (other contacts are floating)
    contact_amplitude=np.zeros(len(contacts_with_current),float)                        #stores assigned amplitudes of the currents
    fl_ind=np.zeros((len(contacts_with_current),len(contacts_with_current)-1),float)    ##strores relative ind of floats (if three current contacts, it will store [[1,2][0,2],[0,1]])
    fl_contacts_rel_ind=np.arange(len(contacts_with_current))

    if Field_calc_param.EQS_mode == 'EQS':
        phi_i_floating=np.zeros((len(contacts_with_current),len(contacts_with_current)-1),float)
        J_im_current_contacts=np.zeros(len(contacts_with_current),float)

    glob_counter=0
    for i in range(len(Domains.Amp_vector)):
        for j in range(len(Domains.Float_on_lead)):
            if Domains.Active_on_lead[i] == Domains.Float_on_lead[j]:      # find the index of the floating conductor (in .med/.msh file) for the active contact (i)

                # to solve "one active contact (i) vs ground" system, get potentials on the rest of the active contacts (which are put to floating condcutors), get current on the active contact
                from Math_module_hybrid_floating import get_field_with_floats
                if Field_calc_param.EQS_mode == 'EQS':
                    phi_r_floating[glob_counter,:],phi_i_floating[glob_counter,:],J_real_current_contacts[glob_counter],J_im_current_contacts[glob_counter]=get_field_with_floats(Field_calc_param.external_grounding,mesh_sol,i,Domains,subdomains,boundaries_sol,Field_calc_param.default_material,Field_calc_param.element_order,Field_calc_param.anisotropy,Field_calc_param.frequenc,Field_calc_param.EQS_mode,Solver_type,calc_with_MPI=True,kappa=kappa)
                else:
                    phi_r_floating[glob_counter,:],__,J_real_current_contacts[glob_counter],__=get_field_with_floats(Field_calc_param.external_grounding,mesh_sol,i,Domains,subdomains,boundaries_sol,Field_calc_param.default_material,Field_calc_param.element_order,Field_calc_param.anisotropy,Field_calc_param.frequenc,Field_calc_param.EQS_mode,Solver_type,calc_with_MPI=True,kappa=kappa)

                fl_ind[glob_counter,:]=fl_contacts_rel_ind[np.arange(len(fl_contacts_rel_ind))!=glob_counter]   # if three current contacts, it will store [[1,2][0,2],[0,1]]
                contact_amplitude[glob_counter]=Domains.Amp_vector[i]

                glob_counter=glob_counter+1

    glob_counter=0

    V_r_BC_for_current=np.zeros(len(contacts_with_current),float)       #real potential for the contact to match the given current
    V_im_BC_for_current=np.zeros(len(contacts_with_current),float)      # only for the EQS formulation

    for i in range(V_r_BC_for_current.shape[0]):
        floating_ind=np.argwhere(fl_ind==i)
        contact_amplitude_others=contact_amplitude[np.arange(len(contact_amplitude))!=i]
        if Field_calc_param.EQS_mode == 'EQS':
            phi_float_vector=(phi_r_floating[floating_ind[:,0],floating_ind[:,1]]+1j*phi_i_floating[floating_ind[:,0],floating_ind[:,1]])   #vector of values of the floating potentials at the contact
            J_others_vector=(J_real_current_contacts[np.arange(len(J_real_current_contacts))!=i]+1j*J_im_current_contacts[np.arange(len(J_im_current_contacts))!=i])
            V_r_BC_for_current[i]=np.real(contact_amplitude[i]*contact_amplitude[i]/(J_real_current_contacts[i]+1j*J_im_current_contacts[i])+np.sum(phi_float_vector*contact_amplitude_others/J_others_vector))
            V_im_BC_for_current[i]=np.imag(contact_amplitude[i]*contact_amplitude[i]/(J_real_current_contacts[i]+1j*J_im_current_contacts[i])+np.sum(phi_float_vector*contact_amplitude_others/J_others_vector))
        else:
            phi_float_vector=(phi_r_floating[floating_ind[:,0],floating_ind[:,1]])   #vector of values of the floating potentials at the contact
            J_others_vector=(J_real_current_contacts[np.arange(len(J_real_current_contacts))!=i])
            V_r_BC_for_current[i]=np.real(contact_amplitude[i]*contact_amplitude[i]/(J_real_current_contacts[i])+np.sum(phi_float_vector*contact_amplitude_others/J_others_vector))

    # not an elegant way but just for the maximum transparency
    if Field_calc_param.EQS_mode == 'EQS':
        scaled_phi=np.complex(1.0,0.0)*np.zeros(len(Domains.Amp_vector),float)
        for i in range(len(Domains.Amp_vector)):
            if Domains.Amp_vector[i]==0.0:
                scaled_phi[i]=0.0+1j*0.0
            else:
                scaled_phi[i]=V_r_BC_for_current[glob_counter]+1j*V_im_BC_for_current[glob_counter]
                glob_counter=glob_counter+1
    else:
        scaled_phi=np.zeros(len(Domains.Amp_vector),float)
        for i in range(len(Domains.Amp_vector)):
            if Domains.Amp_vector[i]==0.0:
                scaled_phi[i]=0.0
            else:
                scaled_phi[i]=V_r_BC_for_current[glob_counter]
                glob_counter=glob_counter+1

    glob_counter=0

    #the results are stored in h5 file, check get_solutions above
    from Math_module_hybrid_floating import get_field_with_scaled_BC
    get_field_with_scaled_BC(Field_calc_param.external_grounding,mesh_sol,Domains,scaled_phi,subdomains,boundaries_sol,Field_calc_param.default_material,Field_calc_param.element_order,Field_calc_param.EQS_mode,Field_calc_param.anisotropy,Field_calc_param.frequenc,Solver_type,calc_with_MPI=True,kappa=kappa)

    minutes=int((tm.time() - start_math)/60)
    secnds=int(tm.time() - start_math)-minutes*60
    if MPI.comm_world.rank==1:
        print("--- Field with superposition was calculated in ",minutes," min ",secnds," s ")
        print("__________________________________")

    return True

if __name__ == '__main__':

    with open(os.environ['PATIENTDIR']+'/Meshes/Mesh_ind.file', "rb") as f:
        Domains = pickle.load(f)

    with open(os.environ['PATIENTDIR']+'/Results_adaptive/Field_calc_param.file', "rb") as f:
        Field_calc_param = pickle.load(f)

    mesh = Mesh()
    hdf = HDF5File(mesh.mpi_comm(), os.environ['PATIENTDIR']+"/Results_adaptive/Mesh_to_solve.h5", "r")
    hdf.read(mesh, "/mesh", False)
    subdomains = MeshFunction("size_t", mesh, 3)
    hdf.read(subdomains, "/subdomains")
    boundaries = MeshFunction("size_t", mesh, 2)
    hdf.read(boundaries, "/boundaries")

    V0_r=FunctionSpace(mesh,'DG',0)
    V0_i=FunctionSpace(mesh,'DG',0)
    kappa_r=Function(V0_r)
    kappa_i=Function(V0_i)
    hdf.read(kappa_r, "/kappa_r")
    kappa=[kappa_r]

    if Field_calc_param.EQS_mode=='EQS':
        hdf.read(kappa_i, "/kappa_i")
        kappa=[kappa_r,kappa_i]

    #anisotropy will be read in at the site

    hdf.close()

    compute_field_with_superposition(mesh,Domains,subdomains,boundaries,kappa,Field_calc_param)
