# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 10:39:20 2017

@author: konstantin
"""
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=FutureWarning)
    import h5py

import logging
logging.getLogger('UFL').setLevel(logging.WARNING)
logging.getLogger('FFC').setLevel(logging.WARNING)

from dolfin import *
import numpy as np
import time
import pickle
import os

from tissue_dielectrics import DielectricProperties

parameters['linear_algebra_backend']='PETSc'
set_log_active(False)   #turns off debugging info

def get_E_field(mesh,element_order,Laplace_eq,phi_real,phi_imag):

    if element_order>1:
        W = VectorFunctionSpace(mesh,'DG',element_order-1)
        W_i = VectorFunctionSpace(mesh,'DG',element_order-1)
    else:
        W = VectorFunctionSpace(mesh,'DG',element_order)
        W_i = VectorFunctionSpace(mesh,'DG',element_order)

    #Explicit E-field projection
    w = TestFunction(W)
    Pv = TrialFunction(W)
    E_field = Function(W)
    a_local = inner(w, Pv) * dx
    L_local = inner(w, -grad(phi_real)) * dx
    A_local, b_local = assemble_system(a_local, L_local, bcs=[])
    local_solver = PETScKrylovSolver('bicgstab')
    local_solver.solve(A_local,E_field.vector(),b_local)

    E_field_im = Function(W)
    E_field_im.vector()[:] = 0.0        #fake
    if Laplace_eq == 'EQS':
        w_i = TestFunction(W_i)
        Pv_i = TrialFunction(W_i)
        E_field_im = Function(W_i)
        a_local = inner(w_i, Pv_i) * dx
        L_local = inner(w_i, -grad(phi_imag)) * dx
        A_local, b_local = assemble_system(a_local, L_local, bcs=[])
        local_solver = PETScKrylovSolver('bicgstab')
        local_solver.solve(A_local,E_field_im.vector(),b_local)

    return E_field,E_field_im

def get_current_on_multiple_contacts(ext_grounding,facets_ground,mesh,boundaries,Laplace_eq,Contacts_indices,Phi_vector,E_field,E_field_im,kappa,Cond_tensor):

    facets=MeshFunction('size_t',mesh, 2)
    facets.set_all(0)

    # mark all contacts where we have Diriclet BC
    for bc_i in range(len(Contacts_indices)):
        facets.array()[boundaries.array()==Contacts_indices[bc_i]]=bc_i+1

    ds=Measure("ds",domain=mesh,subdomain_data=facets)
    dS_int=Measure("dS",domain=mesh,subdomain_data=facets)

    if ext_grounding==False:
        J_currents_real=np.zeros(len(Contacts_indices),float)
        J_currents_imag=np.zeros(len(Contacts_indices),float)
    else:
        J_currents_real=np.zeros(len(Contacts_indices)+1,float)
        J_currents_imag=np.zeros(len(Contacts_indices)+1,float)

    n = FacetNormal(mesh)
    # getting currets where we have Dirichlet BC. For all contacts but the ground we have to integrate inside the computational domain (thus different expressions)
    for bc_i in range(len(Contacts_indices)):
        if Phi_vector[bc_i]==0.0:
            if Cond_tensor!=False:
                if Laplace_eq == 'EQS':
                    j_dens_real = dot(Cond_tensor*E_field,-1*n)*ds(bc_i+1)-dot(kappa[1]*E_field_im,-1*n)*ds(bc_i+1)
                    j_dens_im= dot(Cond_tensor*E_field_im,-1*n)*ds(bc_i+1)+dot(kappa[1]*E_field,-1*n)*ds(bc_i+1)
                else:
                    j_dens_real = dot(Cond_tensor*E_field,-1*n)*ds(bc_i+1)
            else:
                if Laplace_eq == 'EQS':
                    j_dens_real = dot(kappa[0]*E_field,-1*n)*ds(bc_i+1)-dot(kappa[1]*E_field_im,-1*n)*ds(bc_i+1)
                    j_dens_im= dot(kappa[0]*E_field_im,-1*n)*ds(bc_i+1)+dot(kappa[1]*E_field,-1*n)*ds(bc_i+1)
                else:
                    j_dens_real = dot(kappa[0]*E_field,-1*n)*ds(bc_i+1)
            J_currents_real[bc_i]=assemble(j_dens_real)
            if Laplace_eq == 'EQS':
                J_currents_imag[bc_i]=assemble(j_dens_im)
        else:
            if Cond_tensor!=False:
                if Laplace_eq == 'EQS':
                    j_dens_real = dot(Cond_tensor*E_field,-1*n)('-')*dS_int(bc_i+1)-dot(kappa[1]*E_field_im,-1*n)('-')*dS_int(bc_i+1)
                    j_dens_im= dot(Cond_tensor*E_field_im,-1*n)('-')*dS_int(bc_i+1)+dot(kappa[1]*E_field,-1*n)('-')*dS_int(bc_i+1)
                else:
                    j_dens_real = dot(Cond_tensor*E_field,-1*n)('-')*dS_int(bc_i+1)

            else:
                if Laplace_eq == 'EQS':
                    j_dens_real = dot(kappa[0]*E_field,-1*n)('-')*dS_int(bc_i+1)-dot(kappa[1]*E_field_im,-1*n)('-')*dS_int(bc_i+1)
                    j_dens_im= dot(kappa[0]*E_field_im,-1*n)('-')*dS_int(bc_i+1)+dot(kappa[1]*E_field,-1*n)('-')*dS_int(bc_i+1)
                else:
                    j_dens_real = dot(kappa[0]*E_field,-1*n)('-')*dS_int(bc_i+1)
            J_currents_real[bc_i]=assemble(j_dens_real)
            if Laplace_eq == 'EQS':
                J_currents_imag[bc_i]=assemble(j_dens_im)

    if ext_grounding == True:
        ds_ext=Measure("ds",domain=mesh,subdomain_data=facets_ground)
        if Cond_tensor!=False:
            if Laplace_eq == 'EQS':
                j_dens_real = dot(Cond_tensor*E_field,-1*n)*ds_ext(1)-dot(kappa[1]*E_field_im,-1*n)*ds_ext(1)
                j_dens_im= dot(Cond_tensor*E_field_im,-1*n)*ds_ext(1)+dot(kappa[1]*E_field,-1*n)*ds_ext(1)
            else:
                j_dens_real = dot(Cond_tensor*E_field,-1*n)*ds_ext(1)
        else:
            if Laplace_eq == 'EQS':
                j_dens_real = dot(kappa[0]*E_field,-1*n)*ds_ext(1)-dot(kappa[1]*E_field_im,-1*n)*ds_ext(1)
                j_dens_im= dot(kappa[0]*E_field_im,-1*n)*ds_ext(1)+dot(kappa[1]*E_field,-1*n)*ds_ext(1)
            else:
                j_dens_real = dot(kappa[0]*E_field,-1*n)*ds_ext(1)

        J_currents_real[-1]=assemble(j_dens_real)
        if Laplace_eq == 'EQS':
            J_currents_imag[-1]=assemble(j_dens_im)

    return J_currents_real,J_currents_imag

def get_field_with_floats(Sim_setup,active_index,Domains,Solver_type):

    set_log_active(False)   #turns off debugging info
    parameters['linear_algebra_backend']='PETSc'

    # to get conductivity (and permittivity if EQS formulation) mapped accrodingly to the subdomains. k_val_r is just a list of conductivities (S/mm!) in a specific order to scale the cond. tensor
    from FEM_in_spectrum import get_dielectric_properties_from_subdomains
    kappa,k_val_r=get_dielectric_properties_from_subdomains(Sim_setup.mesh,Sim_setup.subdomains,Sim_setup.Laplace_eq,Domains.Float_contacts,Sim_setup.conductivities,Sim_setup.rel_permittivities,Sim_setup.sine_freq)

    # to get tensor scaled by the conductivity map
    if Sim_setup.anisotropy==1:
        from FEM_in_spectrum import get_scaled_cond_tensor
        Cond_tensor=get_scaled_cond_tensor(Sim_setup.mesh,Sim_setup.subdomains,Sim_setup.sine_freq,Sim_setup.signal_freq,Sim_setup.unscaled_tensor,k_val_r)
    else:
        Cond_tensor=False  #just to initialize

    from FEM_in_spectrum import get_solution_space_and_Dirichlet_BC
    V_space,facets=get_solution_space_and_Dirichlet_BC(Sim_setup.external_grounding,1,Sim_setup.mesh,Sim_setup.subdomains,Sim_setup.boundaries,Sim_setup.element_order,Sim_setup.Laplace_eq,Domains.Active_contacts,Domains.Amp_vector,only_space=True)

    facets_active = MeshFunction('size_t',Sim_setup.mesh,2)
    facets_active.set_all(0)

    # here we have a custom way to assign Dirichlet BC
    dirichlet_bc=[]
    float_surface=2
    active_floats=0

    # assign only the chosen active contact and the ground, other should be just marked on the mesh (starting from 2)
    for bc_i in range(len(Domains.Active_contacts)):
        if bc_i==active_index or Domains.Amp_vector[bc_i]==0.0:
            if Sim_setup.Laplace_eq == 'EQS':
                dirichlet_bc.append(DirichletBC(V_space.sub(0), Domains.Amp_vector[bc_i], Sim_setup.boundaries,Domains.Active_contacts[bc_i]))
                dirichlet_bc.append(DirichletBC(V_space.sub(1), Constant(0.0), Sim_setup.boundaries,Domains.Active_contacts[bc_i]))
            else:
                dirichlet_bc.append(DirichletBC(V_space, Domains.Amp_vector[bc_i], Sim_setup.boundaries,Domains.Active_contacts[bc_i]))

            if bc_i==active_index:
                facets_active.array()[Sim_setup.boundaries.array()==Domains.Active_contacts[bc_i]]=1
        else:
            facets_active.array()[Sim_setup.boundaries.array()==Domains.Active_contacts[bc_i]]=float_surface    #it will not be assigned to always floating contacts
            float_surface+=1
            active_floats+=1

    if Sim_setup.external_grounding==True:
        if Sim_setup.Laplace_eq == 'EQS':
            dirichlet_bc.append(DirichletBC(V_space.sub(0),0.0,facets,1))
            dirichlet_bc.append(DirichletBC(V_space.sub(1),0.0,facets,1))
        else:
            dirichlet_bc.append(DirichletBC(V_space,0.0,facets,1))

    #definitions for integrators
    dx = Measure("dx",domain=Sim_setup.mesh)
    dsS=Measure("ds",domain=Sim_setup.mesh,subdomain_data=facets_active)
    dsS_int=Measure("dS",domain=Sim_setup.mesh,subdomain_data=facets_active)

    # to solve the Laplace equation div(kappa*grad(phi))=0   (variational form: a(u,v)=L(v))
    from FEM_in_spectrum import define_variational_form_and_solve
    if Sim_setup.Laplace_eq == 'EQS':
        phi_sol=define_variational_form_and_solve(V_space,dirichlet_bc,kappa,Sim_setup.Laplace_eq,Cond_tensor,'MUMPS')      #fixed to MUMPS here
    else:
        phi_sol=define_variational_form_and_solve(V_space,dirichlet_bc,kappa,Sim_setup.Laplace_eq,Cond_tensor,Solver_type)

    if Sim_setup.Laplace_eq=='EQS':
        (phi_r,phi_i)=phi_sol.split(deepcopy=True)
    else:
        phi_r=phi_sol
        phi_i=Function(V_space)
        phi_i.vector()[:] = 0.0

    Float_potentials_real=np.zeros(active_floats,float)
    Float_potentials_imag=np.zeros(active_floats,float)

    float_surface=2         #float indicies start from 2

    # assess the floating potential by integrating over the contact's surface
    for fl_in in range(active_floats):
        Float_surface_size=assemble(1.0*dsS_int(float_surface))
        Float_potentials_real[float_surface-2]=assemble(phi_r*dsS_int(float_surface))/Float_surface_size
        if Sim_setup.Laplace_eq=='EQS':
            Float_potentials_imag[float_surface-2]=assemble(phi_i*dsS_int(float_surface))/Float_surface_size

        float_surface+=1

    # to get manual projections for the E-field
    E_field,E_field_im=get_E_field(Sim_setup.mesh,Sim_setup.element_order,Sim_setup.Laplace_eq,phi_r,phi_i)
    #if QS, E_field_im is a null function

    n = FacetNormal(Sim_setup.mesh)
    if Sim_setup.Laplace_eq == 'EQS':
        if Cond_tensor!=False:
            j_dens_real_contact = dot(Cond_tensor*E_field,-1*n)('-')*dsS_int(1)-dot(kappa[1]*E_field_im,-1*n)('-')*dsS_int(1)
            j_dens_im_contact= dot(Cond_tensor*E_field_im,-1*n)('-')*dsS_int(1)+dot(kappa[1]*E_field,-1*n)('-')*dsS_int(1)
        else:
            j_dens_real_contact = dot(kappa[0]*E_field,-1*n)('-')*dsS_int(1)-dot(kappa[1]*E_field_im,-1*n)('-')*dsS_int(1)
            j_dens_im_contact= dot(kappa[0]*E_field_im,-1*n)('-')*dsS_int(1)+dot(kappa[1]*E_field,-1*n)('-')*dsS_int(1)

        J_real=assemble(j_dens_real_contact)
        J_im=assemble(j_dens_im_contact)

        return Float_potentials_real,Float_potentials_imag,J_real,J_im
    else:
        if Cond_tensor!=False:
            j_dens_real_contact = dot(Cond_tensor*E_field,-1*n)('-')*dsS_int(1)
        else:
            j_dens_real_contact = dot(kappa[0]*E_field,-1*n)('-')*dsS_int(1)

        J_real=assemble(j_dens_real_contact)

        return Float_potentials_real,0,J_real,0


#see the manuscript (Design and Implementation/Physics of Volume Conductor Model)
def scale_bc_potentials_with_superposition(Sim_setup,Domains,Solver_type):

    set_log_active(False)   #turns off debugging info
    if Sim_setup.element_order==1:
        logging.critical("element_order is 1, increasing to 2 for current-controlled stimulation")
        Sim_setup.element_order=2

    contacts_with_current=[x for x in Domains.Amp_vector if x != 0.0]       #0.0 are grounded contacts
    phi_r_floating=np.zeros((len(contacts_with_current),len(contacts_with_current)-1),float)       #stores real potential field in the virtual floating contacts (where current is actually assigned)
    J_real_current_contacts=np.zeros(len(contacts_with_current),float)                              #currents computed on the contacts when we solve "one active contact vs ground" system (other contacts are floating)
    contact_amplitude=np.zeros(len(contacts_with_current),float)                                   #stores assigned amplitudes of the currents
    fl_ind=np.zeros((len(contacts_with_current),len(contacts_with_current)-1),float)               #strores relative ind of floats (if three current contacts, it will store [[1,2][0,2],[0,1]])
    fl_contacts_rel_ind=np.arange(len(contacts_with_current))

    if Sim_setup.Laplace_eq=='EQS':
        phi_i_floating=np.zeros((len(contacts_with_current),len(contacts_with_current)-1),float)
        J_im_current_contacts=np.zeros(len(contacts_with_current),float)

    glob_counter=0
    for i in range(len(Domains.Amp_vector)):
        for j in range(len(Domains.Float_on_lead)):
            if Domains.Active_on_lead[i] == Domains.Float_on_lead[j]:      #find the index of the floating conductor (in .med/.msh file) for the active contact (i)

                # to solve "one active contact (i) vs ground" system, get potentials on the rest of the active contacts (which are put to floating condcutors), get current on the active contact
                if Sim_setup.Laplace_eq=='EQS':
                    phi_r_floating[glob_counter,:],phi_i_floating[glob_counter,:],J_real_current_contacts[glob_counter],J_im_current_contacts[glob_counter]=get_field_with_floats(Sim_setup,i,Domains,Solver_type)
                else:
                    phi_r_floating[glob_counter,:],__,J_real_current_contacts[glob_counter],__=get_field_with_floats(Sim_setup,i,Domains,Solver_type)

                fl_ind[glob_counter,:]=fl_contacts_rel_ind[np.arange(len(fl_contacts_rel_ind))!=glob_counter]   # if three current contacts, it will store [[1,2][0,2],[0,1]]
                contact_amplitude[glob_counter]=Domains.Amp_vector[i]

                glob_counter=glob_counter+1

    glob_counter=0

    V_r_BC_for_current=np.zeros(len(contacts_with_current),float)       #real potential for the contact to match the given current
    V_im_BC_for_current=np.zeros(len(contacts_with_current),float)          # only for the EQS formulation

    for i in range(V_r_BC_for_current.shape[0]):
        floating_ind=np.argwhere(fl_ind==i)
        contact_amplitude_others=contact_amplitude[np.arange(len(contact_amplitude))!=i]
        if Sim_setup.Laplace_eq=='EQS':
            phi_float_vector=(phi_r_floating[floating_ind[:,0],floating_ind[:,1]]+1j*phi_i_floating[floating_ind[:,0],floating_ind[:,1]])   #vector of values of the floating potentials at the contact
            J_others_vector=(J_real_current_contacts[np.arange(len(J_real_current_contacts))!=i]+1j*J_im_current_contacts[np.arange(len(J_im_current_contacts))!=i])
            V_r_BC_for_current[i]=np.real(contact_amplitude[i]*contact_amplitude[i]/(J_real_current_contacts[i]+1j*J_im_current_contacts[i])+np.sum(phi_float_vector*contact_amplitude_others/J_others_vector))
            V_im_BC_for_current[i]=np.imag(contact_amplitude[i]*contact_amplitude[i]/(J_real_current_contacts[i]+1j*J_im_current_contacts[i])+np.sum(phi_float_vector*contact_amplitude_others/J_others_vector))
        else:
            phi_float_vector=(phi_r_floating[floating_ind[:,0],floating_ind[:,1]])   #vector of values of the floating potentials at the contact
            J_others_vector=(J_real_current_contacts[np.arange(len(J_real_current_contacts))!=i])
            V_r_BC_for_current[i]=np.real(contact_amplitude[i]*contact_amplitude[i]/(J_real_current_contacts[i])+np.sum(phi_float_vector*contact_amplitude_others/J_others_vector))

    # not an elegant way but just for the maximum transparency
    if Sim_setup.Laplace_eq=='EQS':
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

    return scaled_phi


#def solve_Laplace_multicontact(Full_IFFT,frequenc,freq_signal,mesh,boundaries,subdomains,Domains,core,Vertices_array,cond_vector,element_order,anisotropy,c_c,c00_unscaled,c01_unscaled,c02_unscaled,c11_unscaled,c12_unscaled,c22_unscaled,CPE,CPE_param,output):
def solve_Laplace_multicontact(Sim_setup,Solver_type,Vertices_array,Domains,core,VTA_IFFT,output):

    Sim_setup.external_grounding=True

    set_log_active(False)   #turns off debugging info
    tol=1e-14           #tolerance

    Sim_setup.mesh.coordinates()
    Sim_setup.mesh.init()

    #get scaled potentials on the contacts that provide the desired currents through the contacts
    Phi_scaled=scale_bc_potentials_with_superposition(Sim_setup,Domains,Solver_type)

    #Dirichlet_bc was scaled to match the desired current (using system's linearity).
    from FEM_in_spectrum import get_solution_space_and_Dirichlet_BC
    V_space,facets=get_solution_space_and_Dirichlet_BC(Sim_setup.external_grounding,Sim_setup.c_c,Sim_setup.mesh,Sim_setup.subdomains,Sim_setup.boundaries,Sim_setup.element_order,Sim_setup.Laplace_eq,Domains.Active_contacts,Phi_scaled,only_space=True)


    Dirichlet_bc_scaled=[]
    for bc_i in range(len(Domains.Active_contacts)):
        if Sim_setup.Laplace_eq == 'EQS':
            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(0), np.real(Phi_scaled[bc_i]), Sim_setup.boundaries,Domains.Active_contacts[bc_i]))
            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(1), np.imag(Phi_scaled[bc_i]), Sim_setup.boundaries,Domains.Active_contacts[bc_i]))
        else:
            Dirichlet_bc_scaled.append(DirichletBC(V_space, Phi_scaled[bc_i], Sim_setup.boundaries,Domains.Active_contacts[bc_i]))

    if Sim_setup.external_grounding==True:
        if Sim_setup.Laplace_eq == 'EQS':
            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(0),0.0,facets,1))
            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(1),0.0,facets,1))
        else:
            Dirichlet_bc_scaled.append(DirichletBC(V_space,0.0,facets,1))

    # to get conductivity (and permittivity if EQS formulation) mapped accrodingly to the subdomains. k_val_r is just a list of conductivities (S/mm!) in a specific order to scale the cond. tensor
    from FEM_in_spectrum import get_dielectric_properties_from_subdomains
    kappa,k_val_r=get_dielectric_properties_from_subdomains(Sim_setup.mesh,Sim_setup.subdomains,Sim_setup.Laplace_eq,Domains.Float_contacts,Sim_setup.conductivities,Sim_setup.rel_permittivities,Sim_setup.sine_freq)
    if int(Sim_setup.sine_freq)==int(Sim_setup.signal_freq):
        file=File(os.environ['PATIENTDIR']+'/Field_solutions/Conductivity_map_'+str(Sim_setup.signal_freq)+'Hz.pvd')
        file<<kappa[0]
        if Sim_setup.Laplace_eq == 'EQS':
            file=File(os.environ['PATIENTDIR']+'/Field_solutions/Permittivity_map_'+str(Sim_setup.signal_freq)+'Hz.pvd')
            file<<kappa[1]

    # to get tensor scaled by the conductivity map
    if Sim_setup.anisotropy==1:
        from FEM_in_spectrum import get_scaled_cond_tensor
        Cond_tensor=get_scaled_cond_tensor(Sim_setup.mesh,Sim_setup.subdomains,Sim_setup.sine_freq,Sim_setup.signal_freq,Sim_setup.unscaled_tensor,k_val_r)
    else:
        Cond_tensor=False  #just to initialize

    # to solve the Laplace equation div(kappa*grad(phi))=0   (variational form: a(u,v)=L(v))
    from FEM_in_spectrum import define_variational_form_and_solve
    phi_sol=define_variational_form_and_solve(V_space,Dirichlet_bc_scaled,kappa,Sim_setup.Laplace_eq,Cond_tensor,Solver_type)

    if Sim_setup.Laplace_eq=='EQS':
        (phi_r,phi_i)=phi_sol.split(deepcopy=True)
    else:
        phi_r=phi_sol
        phi_i=Function(V_space)
        phi_i.vector()[:] = 0.0

    # the system was solved for already scaled potential, but to check the scaling accuracy, we will assess the currents through the contacts
    if int(Sim_setup.sine_freq)==int(Sim_setup.signal_freq):
        file=File(os.environ['PATIENTDIR']+'/Field_solutions/Phi_real_scaled_'+str(Sim_setup.signal_freq)+'Hz.pvd')
        file<<phi_r,Sim_setup.mesh
        logging.critical("DoFs on the mesh for {} : {}".format(Sim_setup.Laplace_eq, max(V_space.dofmap().dofs()) + 1))

        # to get function space and manual projections for the E-field
        E_field,E_field_im=get_E_field(Sim_setup.mesh,Sim_setup.element_order,Sim_setup.Laplace_eq,phi_r,phi_i)
        #if QS, E_field_im is 0

        # to get current on the active contacts (including the ground)
        J_currents_real,J_currents_imag = get_current_on_multiple_contacts(Sim_setup.external_grounding,facets,Sim_setup.mesh,Sim_setup.boundaries,Sim_setup.Laplace_eq,Domains.Active_contacts,Phi_scaled,E_field,E_field_im,kappa,Cond_tensor)
        # J_currents_imag is a zero array if 'QS' mode
        logging.critical("Complex currents on contacts at the base frequency (signal repetition rate):")
        for j in range(J_currents_real.shape[0]):
            logging.critical("{} {}".format(J_currents_real[j],J_currents_imag[j]))
        file=File(os.environ['PATIENTDIR']+'/Field_solutions/E_real_scaled_'+str(Sim_setup.sine_freq)+'Hz.pvd')
        file<<E_field,Sim_setup.mesh

    #if Full_IFFT==1:
    #    Hdf=HDF5File(Sim_setup.mesh.mpi_comm(), os.environ['PATIENTDIR']+"/Field_solutions_functions/solution"+str(np.round(Sim_setup.sine_freq,6))+".h5", "w")
    #    Hdf.write(Sim_setup.mesh, "mesh")
    #    Hdf.write(phi_sol, "solution_full")
    #    Hdf.close()

    if VTA_IFFT==1:
        Sim_type='Astrom' #   fixed for now
        E_field_real,E_field_im=get_E_field(Sim_setup.mesh,Sim_setup.element_order,Sim_setup.Laplace_eq,phi_r,phi_i)
        if Sim_type=='Astrom':
            W_amp=FunctionSpace(Sim_setup.mesh,'DG',Sim_setup.element_order-1)
            w_amp = TestFunction(W_amp)
            Pv_amp = TrialFunction(W_amp)
            E_amp_real = Function(W_amp)
            a_local = inner(w_amp, Pv_amp) * dx
            L_local = inner(w_amp, sqrt(dot(E_field_real,E_field_real))) * dx
            A_local, b_local = assemble_system(a_local, L_local, bcs=[])

            local_solver = PETScKrylovSolver('bicgstab')
            local_solver.solve(A_local,E_amp_real.vector(),b_local)

            #E_amp_real.vector()[:]=E_amp_real.vector()

            E_amp_imag = Function(W_amp)
            a_local = inner(w_amp, Pv_amp) * dx
            L_local = inner(w_amp, sqrt(dot(E_field_im,E_field_im))) * dx
            A_local, b_local = assemble_system(a_local, L_local, bcs=[])

            local_solver = PETScKrylovSolver('bicgstab')
            local_solver.solve(A_local,E_amp_imag.vector(),b_local)

        Phi_ROI=np.zeros((Vertices_array.shape[0],5),float)

        #VTA=0.0

        for inx in range(Vertices_array.shape[0]):
            pnt=Point(Vertices_array[inx,0],Vertices_array[inx,1],Vertices_array[inx,2])
            if Sim_setup.mesh.bounding_box_tree().compute_first_entity_collision(pnt)<Sim_setup.mesh.num_cells()*100:

                Phi_ROI[inx,0]=Vertices_array[inx,0]
                Phi_ROI[inx,1]=Vertices_array[inx,1]
                Phi_ROI[inx,2]=Vertices_array[inx,2]

                #if Sim_setup.c_c==1:
                if Sim_type=='Butson':
                    Phi_ROI[inx,3]=Second_deriv(pnt)   # already scaled
                    Phi_ROI[inx,4]=Second_deriv_imag(pnt)    # already scaled
                elif Sim_type=='Astrom':
                    Phi_ROI[inx,3]=E_amp_real(pnt)   # already scaled
                    Phi_ROI[inx,4]=E_amp_imag(pnt)    # already scaled

                #if Sim_setup.sine_freq==Sim_setup.signal_freq and abs(Phi_ROI[inx,3])>=0.3:
                #    VTA+=0.1**3

            else:       # we assign 0.0 here
                Phi_ROI[inx,3]=0.0
                Phi_ROI[inx,4]=0.0

                #print("Couldn't probe the potential at the point ",Vertices_array[inx,0],Vertices_array[inx,1],Vertices_array[inx,2])
                #print("check the neuron array, exiting....")
                #raise SystemExit

        fre_vector=[Sim_setup.sine_freq]*Phi_ROI.shape[0]
        comb=np.vstack((Phi_ROI[:,0],Phi_ROI[:,1],Phi_ROI[:,2],Phi_ROI[:,3],Phi_ROI[:,4],fre_vector)).T

        f = h5py.File(os.environ['PATIENTDIR']+'/Field_solutions/sol_cor'+str(core)+'.h5','a')
        f.create_dataset(str(Sim_setup.sine_freq), data=comb)
        f.close()

    else:
        Phi_ROI=np.zeros((Vertices_array.shape[0],5),float)

        for inx in range(Vertices_array.shape[0]):
            pnt=Point(Vertices_array[inx,0],Vertices_array[inx,1],Vertices_array[inx,2])
            if Sim_setup.mesh.bounding_box_tree().compute_first_entity_collision(pnt)<Sim_setup.mesh.num_cells()*100:

                Phi_ROI[inx,0]=Vertices_array[inx,0]
                Phi_ROI[inx,1]=Vertices_array[inx,1]
                Phi_ROI[inx,2]=Vertices_array[inx,2]
                Phi_ROI[inx,3]=phi_r(pnt)
                Phi_ROI[inx,4]=phi_i(pnt)

            else:
                logging.critical("Couldn't probe the potential at the point {}, {}, {}".format(Vertices_array[inx,0],Vertices_array[inx,1],Vertices_array[inx,2]))
                logging.critical("check the neuron array, exitting....")
                raise SystemExit

        fre_vector=[Sim_setup.sine_freq]*Phi_ROI.shape[0]
    ###    freq_vector=[frequenc]*coordinates.shape[0]
    #####    com=np.vstack((coordinat[0:10,0],coordinat[0:10,1],coordinat[0:10,2],real_par[0:10],image_par[0:10],fre_vector[0:10])).T
        comb=np.vstack((Phi_ROI[:,0],Phi_ROI[:,1],Phi_ROI[:,2],Phi_ROI[:,3],Phi_ROI[:,4],fre_vector)).T

        if int(Sim_setup.sine_freq) == int(Sim_setup.signal_freq):
            np.savetxt(os.environ['PATIENTDIR']+'/Phi_real_on_Axons_' + str(Sim_setup.signal_freq) + '_PO.csv', Phi_ROI[:, :4], delimiter=' ')

        f = h5py.File(os.environ['PATIENTDIR']+'/Field_solutions/sol_cor'+str(core)+'.h5','a')
        f.create_dataset(str(Sim_setup.sine_freq), data=comb)
        f.close()

        output.put(1)
#=============================================================================#




