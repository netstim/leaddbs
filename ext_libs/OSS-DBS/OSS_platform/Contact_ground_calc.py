'''
    By K. Butenko
    This script computes Contact-Ground solution for each contact and stores V on the neuron points as well as the floating potentials at active contacts
'''
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=FutureWarning)
    import h5py

import logging
logging.getLogger('UFL').setLevel(logging.WARNING)
logging.getLogger('FFC').setLevel(logging.WARNING)

from dolfin import *
from pandas import read_csv
import numpy as np
import os.path
import time as tm

import warnings
warnings.filterwarnings("ignore")

#from tissue_dielectrics import DielectricProperties

parameters['linear_algebra_backend']='PETSc'
set_log_active(False)   #turns off debugging info

# to compute Contact-Ground solution we use the same routine as for floats, but with minor adjustments
# in the context of this function Domains.Amp_vector is all 1.0 (grounding is external), thus Float_contacts should be empty



#if calculating with MPI, the dielectic properties (kappa) and the scaled tensor were already prepared
def get_field_with_floats(external_grounding,mesh_sol,active_index,Domains,subdomains,boundaries_sol,element_order,anisotropy,frequenc,Laplace_mode,Solver_type,Vertices_array,conductivities,rel_permittivities,unscaled_tensor,VTA_IFFT,calc_with_MPI=False,kappa=False):

    set_log_active(False)   #turns off debugging info
    parameters['linear_algebra_backend']='PETSc'
    mesh_sol.coordinates()
    mesh_sol.init()

    from FEM_in_spectrum import get_dielectric_properties_from_subdomains
    kappa,k_val_r=get_dielectric_properties_from_subdomains(mesh_sol,subdomains,Laplace_mode,Domains.Float_contacts,conductivities,rel_permittivities,frequenc)

    if anisotropy==1:
        from FEM_in_spectrum import get_scaled_cond_tensor
        Cond_tensor=get_scaled_cond_tensor(mesh_sol,subdomains,frequenc,130.0,unscaled_tensor,k_val_r)
    else:
        Cond_tensor=False  #just to initialize

    from FEM_in_spectrum import get_solution_space_and_Dirichlet_BC
    V_space,facets=get_solution_space_and_Dirichlet_BC(external_grounding,1,mesh_sol,subdomains,boundaries_sol,element_order,Laplace_mode,Domains.Active_contacts,Domains.Amp_vector,only_space=True)

    facets_active = MeshFunction('size_t',mesh_sol,2)
    facets_active.set_all(0)

    # here we have a custom way to assign Dirichlet BC
    dirichlet_bc=[]
    float_surface=2
    active_floats=0

    # assign only the chosen active contact and the ground, other should be just marked on the mesh (starting from 2)
    #print("Contacts: ",Domains.Active_contacts)
    #print("fi: ",Domains.Amp_vector)
    #print("active_index: ",active_index)
    #print("ponteial: ", Domains.Amp_vector[active_index])

    for bc_i in range(len(Domains.Active_contacts)):
        if bc_i==active_index or Domains.Amp_vector[bc_i]==0.0:
            if Laplace_mode == 'EQS':
                dirichlet_bc.append(DirichletBC(V_space.sub(0), Domains.Amp_vector[bc_i], boundaries_sol,Domains.Active_contacts[bc_i]))
                dirichlet_bc.append(DirichletBC(V_space.sub(1), Constant(0.0), boundaries_sol,Domains.Active_contacts[bc_i]))
            else:
                dirichlet_bc.append(DirichletBC(V_space, Domains.Amp_vector[bc_i], boundaries_sol,Domains.Active_contacts[bc_i]))

            if bc_i==active_index:
                facets_active.array()[boundaries_sol.array()==Domains.Active_contacts[bc_i]]=1
                #print("phi to solve for: ",Domains.Amp_vector[bc_i])
                #print(active_index)
        else:
            facets_active.array()[boundaries_sol.array()==Domains.Active_contacts[bc_i]]=float_surface    #it will not be assigned to always floating contacts
            float_surface=float_surface+1
            active_floats=active_floats+1


    if external_grounding==True:
        if Laplace_mode == 'EQS':
            dirichlet_bc.append(DirichletBC(V_space.sub(0),0.0,facets,1))
            dirichlet_bc.append(DirichletBC(V_space.sub(1),0.0,facets,1))
        else:
            dirichlet_bc.append(DirichletBC(V_space,0.0,facets,1))

    #definitions for integrators
    dx = Measure("dx",domain=mesh_sol)
    dsS=Measure("ds",domain=mesh_sol,subdomain_data=facets_active)
    dsS_int=Measure("dS",domain=mesh_sol,subdomain_data=facets_active)

    #An_surface_size=assemble(1.0*dsS_int(1))
    #Cat_surface_size=assemble(1.0*dsS_int(2))
    #print "Anode size: "
    #print An_surface_size

   # to solve the Laplace equation div(kappa*grad(phi))=0   (variational form: a(u,v)=L(v))

    from FEM_in_spectrum import define_variational_form_and_solve
    phi_sol=define_variational_form_and_solve(V_space,dirichlet_bc,kappa,Laplace_mode,Cond_tensor,Solver_type)      # with multiple floats MUMPS is the most stable though slow

    if Laplace_mode == 'EQS':
        (phi_r_sol,phi_i_sol)=phi_sol.split(deepcopy=True)
    else:
        phi_r_sol=phi_sol
        phi_i_sol=Function(V_space)
        phi_i_sol.vector()[:] = 0.0


    # to get manual projections for the E-field
    from FEM_in_spectrum_multicontact import get_E_field
    E_field,E_field_im=get_E_field(mesh_sol,element_order,Laplace_mode,phi_r_sol,phi_i_sol)
    #if QS, E_field_im is a null function

    n = FacetNormal(mesh_sol)
    if Laplace_mode == 'EQS':
        if Cond_tensor!=False:
            j_dens_real_contact = dot(Cond_tensor*E_field,-1*n)('-')*dsS_int(1)-dot(kappa[1]*E_field_im,-1*n)('-')*dsS_int(1)
            j_dens_im_contact= dot(Cond_tensor*E_field_im,-1*n)('-')*dsS_int(1)+dot(kappa[1]*E_field,-1*n)('-')*dsS_int(1)
        else:
            j_dens_real_contact = dot(kappa[0]*E_field,-1*n)('-')*dsS_int(1)-dot(kappa[1]*E_field_im,-1*n)('-')*dsS_int(1)
            j_dens_im_contact= dot(kappa[0]*E_field_im,-1*n)('-')*dsS_int(1)+dot(kappa[1]*E_field,-1*n)('-')*dsS_int(1)

        J_real=assemble(j_dens_real_contact)
        J_im=assemble(j_dens_im_contact)

        #return Float_potentials_real,Float_potentials_imag,J_real,J_im
    else:
        if Cond_tensor!=False:
            j_dens_real_contact = dot(Cond_tensor*E_field,-1*n)('-')*dsS_int(1)
        else:
            j_dens_real_contact = dot(kappa[0]*E_field,-1*n)('-')*dsS_int(1)

        J_real=assemble(j_dens_real_contact)





    # we will save for a unit current
    if VTA_IFFT==1:
        Sim_type='Astrom' #   fixed for now
        if Sim_type=='Astrom' or Sim_type=='Butson':
            if Laplace_mode == 'QS':
                phi_i_check=Function(V_space)
                phi_i_check.vector()[:] = 0.0
                phi_r_check=Function(V_space)
                phi_r_check.vector()[:]=1.0*phi_r_sol.vector()[:]/J_real
#            else:
#            # Solve for rescaled
#                Dirichlet_bc_scaled=[]
#                for bc_i in range(len(Domains.Active_contacts)):          #CPE estimation is valid only for one activa and one ground contact configuration
#                    if Laplace_mode == 'EQS':
#                        if Domains.Amp_vector[bc_i]==0.0:
#                            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(0), 0.0, Sim_setup.boundaries,Domains.Active_contacts[bc_i]))
#                            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(1), 0.0, Sim_setup.boundaries,Domains.Active_contacts[bc_i]))
#                        else:
#                            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(0), np.real((Domains.Amp_vector[bc_i])/J_ground),Sim_setup.boundaries,Domains.Active_contacts[bc_i]))
#                            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(1), np.imag((Domains.Amp_vector[bc_i])/J_ground),Sim_setup.boundaries,Domains.Active_contacts[bc_i]))
#                    else:
#                        if Domains.Amp_vector[bc_i]==0.0:
#                            Dirichlet_bc_scaled.append(DirichletBC(V_space, 0.0, Sim_setup.boundaries,Domains.Active_contacts[bc_i]))
#                        else:
#                            Dirichlet_bc_scaled.append(DirichletBC(V_space, np.real((Domains.Amp_vector[bc_i])/J_ground),Sim_setup.boundaries,Domains.Active_contacts[bc_i]))
#
#                if Sim_setup.external_grounding==True:
#                    if Laplace_mode == 'EQS':
#                        Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(0),0.0,facets,1))
#                        Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(1),0.0,facets,1))
#                    else:
#                        Dirichlet_bc_scaled.append(DirichletBC(V_space,0.0,facets,1))
#
#                phi_sol_check=define_variational_form_and_solve(V_space,Dirichlet_bc_scaled,kappa,Sim_setup.Laplace_eq,Cond_tensor,Solver_type)
#
#                if Laplace_mode=='EQS':
#                    (phi_r_check,phi_i_check)=phi_sol_check.split(deepcopy=True)
#                else:
#                    phi_r_check=phi_sol_check
#                    phi_i_check=Function(V_space)
#                    phi_i_check.vector()[:] = 0.0

            #from FEM_in_spectrum import get_current
            #J_ground,E_field_r,E_field_im=get_current(mesh_sol,facets,boundaries_sol,element_order,Laplace_mode,Domains.Active_contacts,kappa,Cond_tensor,phi_r_check,phi_i_check,ground_index,get_E_field=True)

            # to get manual projections for the E-field
            from FEM_in_spectrum_multicontact import get_E_field
            E_field,E_field_im=get_E_field(mesh_sol,element_order,Laplace_mode,phi_r_check,phi_i_check)
            #if QS, E_field_im is a null function

            if Sim_type=='Astrom':
                W_amp=FunctionSpace(mesh_sol,'DG',element_order-1)
                w_amp = TestFunction(W_amp)
                Pv_amp = TrialFunction(W_amp)
                E_amp_real = Function(W_amp)
                a_local = inner(w_amp, Pv_amp) * dx
                L_local = inner(w_amp, sqrt(dot(E_field,E_field))) * dx
                A_local, b_local = assemble_system(a_local, L_local, bcs=[])

                local_solver = PETScKrylovSolver('bicgstab')
                local_solver.solve(A_local,E_amp_real.vector(),b_local)

                #E_amp_real.vector()[:]=E_amp_real.vector()

#                E_amp_imag = Function(W_amp)
#                a_local = inner(w_amp, Pv_amp) * dx
#                L_local = inner(w_amp, sqrt(dot(E_field_im,E_field_im))) * dx
#                A_local, b_local = assemble_system(a_local, L_local, bcs=[])
#
#                local_solver = PETScKrylovSolver('bicgstab')
#                local_solver.solve(A_local,E_amp_imag.vector(),b_local)
            elif Sim_type=='Butson':
                from ufl import nabla_div

                W_amp=FunctionSpace(mesh_sol,'DG',element_order-1)
                w_amp = TestFunction(W_amp)
                Pv_amp = TrialFunction(W_amp)
                Second_deriv= Function(W_amp)
                a_local = inner(w_amp, Pv_amp) * dx
                L_local = inner(w_amp, nabla_div(E_field)) * dx
                A_local, b_local = assemble_system(a_local, L_local, bcs=[])

                local_solver = PETScKrylovSolver('bicgstab')
                local_solver.solve(A_local,Second_deriv.vector(),b_local)

#                W_amp=FunctionSpace(mesh_sol,'DG',element_order-1)
#                w_amp = TestFunction(W_amp)
#                Pv_amp = TrialFunction(W_amp)
#                Second_deriv_imag= Function(W_amp)
#                a_local = inner(w_amp, Pv_amp) * dx
#                L_local = inner(w_amp, nabla_div(E_field_im)) * dx
#                A_local, b_local = assemble_system(a_local, L_local, bcs=[])
#
#                local_solver = PETScKrylovSolver('bicgstab')
#                local_solver.solve(A_local,Second_deriv_imag.vector(),b_local)


            Phi_ROI=np.zeros((Vertices_array.shape[0]),float)

            for inx in range(Vertices_array.shape[0]):
                pnt=Point(Vertices_array[inx,0],Vertices_array[inx,1],Vertices_array[inx,2])
                if mesh_sol.bounding_box_tree().compute_first_entity_collision(pnt)<mesh_sol.num_cells()*100:

                    #Phi_ROI[inx,0]=Vertices_array[inx,0]
                    #Phi_ROI[inx,1]=Vertices_array[inx,1]
                    #Phi_ROI[inx,2]=Vertices_array[inx,2]

                    #if Sim_setup.c_c==1:
                    if Sim_type=='Butson':
                        Phi_ROI[inx]=Second_deriv(pnt)
                        #Phi_ROI[inx,4]=Second_deriv_imag(pnt)
                    elif Sim_type=='Astrom':
                        Phi_ROI[inx]=E_amp_real(pnt)  # if VC, they are already scaled here and the signal will be unit
                        #Phi_ROI[inx,4]=E_amp_imag(pnt)  # if CC, they will be scaled as the signal (only one contact and ground here, so ok)

#            fre_vector=[Sim_setup.sine_freq]*Phi_ROI.shape[0]
#            comb=np.vstack((Phi_ROI[:,0],Phi_ROI[:,1],Phi_ROI[:,2],Phi_ROI[:,3],Phi_ROI[:,4],fre_vector)).T
#
#            f = h5py.File(os.environ['PATIENTDIR']+'/Field_solutions/sol_cor'+str(core)+'.h5','a')
#            f.create_dataset(str(Sim_setup.sine_freq), data=comb)
#            f.close()


    else:
        Phi_ROI=np.zeros((Vertices_array.shape[0]),float)
        for inx in range(Vertices_array.shape[0]):
            pnt=Point(Vertices_array[inx,0],Vertices_array[inx,1],Vertices_array[inx,2])
            if mesh_sol.bounding_box_tree().compute_first_entity_collision(pnt)<mesh_sol.num_cells()*100:

                Phi_ROI[inx]=phi_r_sol(pnt)/J_real

            else:
                logging.critical("Couldn't probe the potential at the point {},{},{}".format (Vertices_array[inx,0],Vertices_array[inx,1],Vertices_array[inx,2]))
                logging.critical("check the neuron array, exiting....")
                raise SystemExit


    Float_potentials_real=np.zeros(active_floats,float)
    Float_potentials_imag=np.zeros(active_floats,float)

    float_surface=2      #float indicies start from 2

    # assess the floating potential by integrating over the contact's surface
    for fl_in in range(active_floats):
        Float_surface_size=assemble(1.0*dsS_int(float_surface))
        Float_potentials_real[float_surface-2]=assemble(phi_r_sol*dsS_int(float_surface))/Float_surface_size
        if Laplace_mode == 'EQS':
            Float_potentials_imag[float_surface-2]=assemble(phi_i_sol*dsS_int(float_surface))/Float_surface_size

        float_surface=float_surface+1

    return Float_potentials_real,Float_potentials_imag,Phi_ROI



def compute_fields_from_unit_currents(Field_calc_param,Solver_type,Vertices,Domains,core,VTA_IFFT,output):

    #start_math=tm.time()

    ##load the neuron array
    #Vertices_get=read_csv('Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', delimiter=' ', header=None)
    #Vertices=Vertices_get.values


    if Field_calc_param.element_order==1:
        if Field_calc_param.sine_freq==Field_calc_param.signal_freq:
            logging.critical("Selected element_order (1st) is too low for current-controlled stimulation, increasing to 2nd")
        Field_calc_param.element_order=2

    if Field_calc_param.sine_freq==Field_calc_param.signal_freq:

        logging.critical("Computing field with superposition on mesh with {} elements".format(Field_calc_param.mesh.num_cells()))
        logging.critical("{} computations are required for the iteration".format(len(Domains.Amp_vector)))

    #contacts_with_current=[x for x in Domains.Amp_vector if x != 0.0]       #0.0 are grounded contacts
    contacts_with_current=[x for x in Domains.Amp_vector if x != 0.0]       #0.0 are grounded contacts

    phi_r_floating=np.zeros((len(contacts_with_current),len(contacts_with_current)-1),float)       #stores real potential field in the virtual floating contacts (where current is actually assigned)
    J_real_current_contacts=np.zeros(len(contacts_with_current),float)                  #currents computed on the contacts when we solve "one active contact vs ground" system (other contacts are floating)
    contact_amplitude=np.zeros(len(contacts_with_current),float)                        #stores assigned amplitudes of the currents
    fl_ind=np.zeros((len(contacts_with_current),len(contacts_with_current)-1),float)    ##strores relative ind of floats (if three current contacts, it will store [[1,2][0,2],[0,1]])
    fl_contacts_rel_ind=np.arange(len(contacts_with_current))

    if Field_calc_param.Laplace_eq == 'EQS':
        phi_i_floating=np.zeros((len(contacts_with_current),len(contacts_with_current)-1),float)
        J_im_current_contacts=np.zeros(len(contacts_with_current),float)

    #print("Active_on_lead: ",Domains.Active_on_lead)
    #print("Float_on_lead: ",Domains.Float_on_lead)
    #print(len(contacts_with_current))

    if VTA_IFFT==1:
        logging.critical("VTA from E-field metrics is not yet supported for current superposition mode")
        raise SystemExit
    else:
        Solutions_on_points=np.zeros((Vertices.shape[0],len(contacts_with_current)),float)

    glob_counter=0
    for i in range(len(Domains.Amp_vector)):
        for j in range(len(Domains.Float_on_lead)):     # this iterator does not exist here
            if Domains.Active_on_lead[i] == Domains.Float_on_lead[j]:       #find the index of the floating conductor (in .med/.msh file) for the active contact (i)

                # to solve "one active contact (i) vs ground" system, get potentials on the rest of the active contacts (which are put to floating condcutors), get current on the active contact
                if Field_calc_param.Laplace_eq == 'EQS':
                    phi_r_floating[glob_counter,:],phi_i_floating[glob_counter,:],Solutions_on_points[:,glob_counter]=get_field_with_floats(Field_calc_param.external_grounding,Field_calc_param.mesh,i,Domains,Field_calc_param.subdomains,Field_calc_param.boundaries,Field_calc_param.element_order,Field_calc_param.anisotropy,Field_calc_param.sine_freq,Field_calc_param.Laplace_eq,Solver_type,Vertices,Field_calc_param.conductivities,Field_calc_param.rel_permittivities,Field_calc_param.unscaled_tensor,VTA_IFFT)
                else:
                    phi_r_floating[glob_counter,:],__,Solutions_on_points[:,glob_counter]=get_field_with_floats(Field_calc_param.external_grounding,Field_calc_param.mesh,i,Domains,Field_calc_param.subdomains,Field_calc_param.boundaries,Field_calc_param.element_order,Field_calc_param.anisotropy,Field_calc_param.sine_freq,Field_calc_param.Laplace_eq,Solver_type,Vertices,Field_calc_param.conductivities,Field_calc_param.rel_permittivities,Field_calc_param.unscaled_tensor,VTA_IFFT)

                #fl_ind[glob_counter,:]=fl_contacts_rel_ind[np.arange(len(fl_contacts_rel_ind))!=glob_counter]   # if three current contacts, it will store [[1,2][0,2],[0,1]]
                #contact_amplitude[glob_counter]=Domains.Amp_vector[i]

                glob_counter=glob_counter+1

    glob_counter=0


    # here load in loop all the solutions on points into one array
    fre_vector=[Field_calc_param.sine_freq]*Vertices.shape[0]
    # 8 contact case for now

    N_contacts=len(contacts_with_current)
    if N_contacts==8:
        comb=np.vstack((Solutions_on_points[:,0],Solutions_on_points[:,1],Solutions_on_points[:,2],Solutions_on_points[:,3],Solutions_on_points[:,4],Solutions_on_points[:,5],Solutions_on_points[:,6],Solutions_on_points[:,7],fre_vector)).T
    elif N_contacts==4:
        comb=np.vstack((Solutions_on_points[:,0],Solutions_on_points[:,1],Solutions_on_points[:,2],Solutions_on_points[:,3],fre_vector)).T
    else:
        logging.critical("Electrode type is not supported or contacts were not counted correctly")    #comb=np.vstack((Solutions_on_points[:,:],fre_vector)).T




    f = h5py.File(os.environ['PATIENTDIR']+'/Field_solutions/sol_per_contact_cor'+str(core)+'.h5','a')
    f.create_dataset(str(Field_calc_param.sine_freq), data=comb)
    f.close()

    #only store the real part for now
    fre_vector_2=[Field_calc_param.sine_freq]*(len(contacts_with_current))
    #comb_fl=np.vstack((phi_r_floating[0,:],phi_r_floating[1,:],phi_r_floating[2,:],phi_r_floating[3,:],phi_r_floating[4,:],phi_r_floating[5,:],phi_r_floating[6,:],phi_r_floating[7,:],fre_vector_2)).T

    if N_contacts==8:
        comb_fl=np.vstack((phi_r_floating[:,0],phi_r_floating[:,1],phi_r_floating[:,2],phi_r_floating[:,3],phi_r_floating[:,4],phi_r_floating[:,5],phi_r_floating[:,6],fre_vector_2)).T
    elif N_contacts==4:
        comb_fl=np.vstack((phi_r_floating[:,0],phi_r_floating[:,1],phi_r_floating[:,2],fre_vector_2)).T

    # we will need this file to recompute BCs on contacts
    f = h5py.File(os.environ['PATIENTDIR']+'/Field_solutions/sol_fl_contacts'+str(core)+'.h5','a')
    f.create_dataset(str(Field_calc_param.sine_freq), data=comb_fl)
    f.close()

    output.put(1)


    # this is only a solution for one frequency. It should be stored i


    # V_r_BC_for_current=np.zeros(len(contacts_with_current),float)       #real potential for the contact to match the given current
    # V_im_BC_for_current=np.zeros(len(contacts_with_current),float)      # only for the EQS formulation

    # for i in range(V_r_BC_for_current.shape[0]):
    #     floating_ind=np.argwhere(fl_ind==i)
    #     contact_amplitude_others=contact_amplitude[np.arange(len(contact_amplitude))!=i]
    #     if Field_calc_param.EQS_mode == 'EQS':
    #         phi_float_vector=(phi_r_floating[floating_ind[:,0],floating_ind[:,1]]+1j*phi_i_floating[floating_ind[:,0],floating_ind[:,1]])   #vector of values of the floating potentials at the contact
    #         J_others_vector=(J_real_current_contacts[np.arange(len(J_real_current_contacts))!=i]+1j*J_im_current_contacts[np.arange(len(J_im_current_contacts))!=i])
    #         V_r_BC_for_current[i]=np.real(contact_amplitude[i]*contact_amplitude[i]/(J_real_current_contacts[i]+1j*J_im_current_contacts[i])+np.sum(phi_float_vector*contact_amplitude_others/J_others_vector))
    #         V_im_BC_for_current[i]=np.imag(contact_amplitude[i]*contact_amplitude[i]/(J_real_current_contacts[i]+1j*J_im_current_contacts[i])+np.sum(phi_float_vector*contact_amplitude_others/J_others_vector))
    #     else:
    #         phi_float_vector=(phi_r_floating[floating_ind[:,0],floating_ind[:,1]])   #vector of values of the floating potentials at the contact
    #         J_others_vector=(J_real_current_contacts[np.arange(len(J_real_current_contacts))!=i])
    #         V_r_BC_for_current[i]=np.real(contact_amplitude[i]*contact_amplitude[i]/(J_real_current_contacts[i])+np.sum(phi_float_vector*contact_amplitude_others/J_others_vector))

    # # not an elegant way but just for the maximum transparency
    # if Field_calc_param.EQS_mode == 'EQS':
    #     scaled_phi=np.complex(1.0,0.0)*np.zeros(len(Domains.Amp_vector),float)
    #     for i in range(len(Domains.Amp_vector)):
    #         if Domains.Amp_vector[i]==0.0:
    #             scaled_phi[i]=0.0+1j*0.0
    #         else:
    #             scaled_phi[i]=V_r_BC_for_current[glob_counter]+1j*V_im_BC_for_current[glob_counter]
    #             glob_counter=glob_counter+1
    # else:
    #     scaled_phi=np.zeros(len(Domains.Amp_vector),float)
    #     for i in range(len(Domains.Amp_vector)):
    #         if Domains.Amp_vector[i]==0.0:
    #             scaled_phi[i]=0.0
    #         else:
    #             scaled_phi[i]=V_r_BC_for_current[glob_counter]
    #             glob_counter=glob_counter+1

    # # quasi_imp is a metric we use to assess the current convergence (an analog to unscaled current).
    # phi_r_sol,phi_i_sol,Field_real,Field_imag,max_E,quasi_imp_real,quasi_imp_im,j_dens_real,j_dens_im=get_field_with_scaled_BC(mesh_sol,Domains,scaled_phi,subdomains,boundaries_sol,Field_calc_param.default_material,Field_calc_param.element_order,Field_calc_param.EQS_mode,Field_calc_param.anisotropy,Field_calc_param.frequenc,Solver_type)

    # minutes=int((tm.time() - start_math)/60)
    # secnds=int(tm.time() - start_math)-minutes*60
    # print("--- Field with superposition was calculated in ",minutes," min ",secnds," s ")
    # print("__________________________________")

    # return phi_r_sol,phi_i_sol,Field_real,Field_imag,max_E,quasi_imp_real,quasi_imp_im,j_dens_real,j_dens_im
