# -*- coding: utf-8 -*-
"""
Created on Sun Aug 19 13:23:54 2018

@author: butenko
"""

import logging
logging.getLogger('UFL').setLevel(logging.WARNING)
logging.getLogger('FFC').setLevel(logging.WARNING)

from dolfin import *
from pandas import read_csv
import numpy as np
import os
import time as tm

from tissue_dielectrics import DielectricProperties

parameters["allow_extrapolation"]=True
parameters['linear_algebra_backend']='PETSc'
set_log_level(LogLevel.CRITICAL)
if MPI.comm_world.rank == 0:
  set_log_level(LogLevel.CRITICAL)


def choose_solver_for_me(EQS_mode,float_conductors):
    if float_conductors != -1:   #that means we have floating conductors
        if EQS_mode=='EQS':
            return('MUMPS')    # maybe for QS with only one floating conductor we could use GMRES
        else:
            return('GMRES')
    else:
        if EQS_mode=='EQS':
            return('BiCGSTAB')
        else:
            return('GMRES')

def get_current_density(mesh,element_order,EQS_mode,kappa,Cond_tensor,E_field_real,E_field_imag):

    if element_order>1:
        W =VectorFunctionSpace(mesh,'DG',element_order-1)
        W_i =VectorFunctionSpace(mesh,'DG',element_order-1)
    else:
        W =VectorFunctionSpace(mesh,'DG',element_order)
        W_i =VectorFunctionSpace(mesh,'DG',element_order)

    w = TestFunction(W)
    Pv = TrialFunction(W)
    j_dens_real = Function(W)
    a_local = inner(w, Pv) * dx

    if EQS_mode == 'EQS':
        w_i = TestFunction(W_i)
        Pv_i = TrialFunction(W_i)
        j_dens_im = Function(W_i)
        a_local_imag = inner(w_i, Pv_i) * dx

        if Cond_tensor!=False:
            L_local = inner(w, (Cond_tensor*E_field_real-kappa[1]*E_field_imag)) * dx
            L_local_imag = inner(w_i, (Cond_tensor*E_field_imag+kappa[1]*E_field_real)) * dx
        else:
            L_local = inner(w, (kappa[0]*E_field_real-kappa[1]*E_field_imag)) * dx
            L_local_imag = inner(w_i, (kappa[0]*E_field_imag+kappa[1]*E_field_real)) * dx

        A_local_imag, b_local_imag = assemble_system(a_local_imag, L_local_imag, bcs=[])
        local_solver = PETScKrylovSolver('bicgstab')
        local_solver.solve(A_local_imag,j_dens_im.vector(),b_local_imag)
    else:
        j_dens_im=Function(W)
        j_dens_im.vector()[:] = 0.0

        if Cond_tensor!=False:
            L_local = inner(w, (Cond_tensor*E_field_real)) * dx
        else:
            L_local = inner(w, (kappa[0]*E_field_real)) * dx

    A_local, b_local = assemble_system(a_local, L_local, bcs=[])
    local_solver = PETScKrylovSolver('bicgstab')
    local_solver.solve(A_local,j_dens_real.vector(),b_local)

    return j_dens_real,j_dens_im


def get_field(mesh_sol,Domains,subdomains,boundaries_sol,Field_calc_param):

    set_log_active(False)   #turns off debugging info
    logging.critical("_________________________")
    parameters['linear_algebra_backend']='PETSc'

    [cond_GM, perm_GM]=DielectricProperties(3).get_dielectrics(Field_calc_param.frequenc)        #3 for grey matter and so on (numeration as in voxel_data)
    [cond_WM, perm_WM]=DielectricProperties(2).get_dielectrics(Field_calc_param.frequenc)
    [cond_CSF, perm_CSF]=DielectricProperties(1).get_dielectrics(Field_calc_param.frequenc)

    [cond_default,perm_default]=DielectricProperties(Field_calc_param.default_material).get_dielectrics(Field_calc_param.frequenc)

    from GUI_inp_dict import d as d_encap
    [cond_encap, perm_encap]=DielectricProperties(d_encap['encap_tissue_type']).get_dielectrics(Field_calc_param.frequenc)
    cond_encap=cond_encap*d_encap['encap_scaling_cond']
    perm_encap=perm_encap*d_encap['encap_scaling_perm']

    if Field_calc_param.Solver_type=='Default':
        Solver_type=choose_solver_for_me(Field_calc_param.EQS_mode,Domains.Float_contacts)    #choses solver basing on the Laplace formulation and whether the floating conductors are used
    else:
        Solver_type=Field_calc_param.Solver_type      # just get the solver directly

    conductivities=[cond_default,cond_GM,cond_WM,cond_CSF,cond_encap]
    rel_permittivities=[perm_default,perm_GM,perm_WM,perm_CSF,perm_encap]

    # to get conductivity (and permittivity if EQS formulation) mapped accrodingly to the subdomains. k_val_r is just a list of conductivities (S/mm!) in a specific order to scale the cond. tensor
    from FEM_in_spectrum import get_dielectric_properties_from_subdomains
    kappa,k_val_r=get_dielectric_properties_from_subdomains(mesh_sol,subdomains,Field_calc_param.EQS_mode,Domains.Float_contacts,conductivities,rel_permittivities,Field_calc_param.frequenc)
    file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_subdomains_map.pvd')
    file<<subdomains
    file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_conductivity_map.pvd')
    file<<kappa[0]
    if Field_calc_param.EQS_mode == 'EQS':
        file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_permittivity_map.pvd')
        file<<kappa[1]


    if Field_calc_param.anisotropy==1:
        # order xx,xy,xz,yy,yz,zz
        c00 = MeshFunction("double", mesh_sol, 3, 0.0)
        c01 = MeshFunction("double", mesh_sol, 3, 0.0)
        c02 = MeshFunction("double", mesh_sol, 3, 0.0)
        c11 = MeshFunction("double", mesh_sol, 3, 0.0)
        c12 = MeshFunction("double", mesh_sol, 3, 0.0)
        c22 = MeshFunction("double", mesh_sol, 3, 0.0)

        # load the diffusion tensor (should be normalized beforehand)
        hdf = HDF5File(mesh_sol.mpi_comm(), os.environ['PATIENTDIR']+"/Results_adaptive/Tensors_to_solve_num_el_"+str(mesh_sol.num_cells())+".h5", "r")
        hdf.read(c00, "/c00")
        hdf.read(c01, "/c01")
        hdf.read(c02, "/c02")
        hdf.read(c11, "/c11")
        hdf.read(c12, "/c12")
        hdf.read(c22, "/c22")
        hdf.close()

        unscaled_tensor=[c00,c01,c02,c11,c12,c22]

        # to get tensor scaled by the conductivity map (twice send Field_calc_param.frequenc to always get unscaled ellipsoid tensor for visualization)
        from FEM_in_spectrum import get_scaled_cond_tensor
        Cond_tensor=get_scaled_cond_tensor(mesh_sol,subdomains,Field_calc_param.frequenc,Field_calc_param.frequenc,unscaled_tensor,k_val_r,plot_tensors=True)
    else:
        Cond_tensor=False  #just to initialize


    #In case of current-controlled stimulation, Dirichlet_bc or the whole potential distribution will be scaled afterwards (due to the system's linearity)
    from FEM_in_spectrum import get_solution_space_and_Dirichlet_BC
    V_space,Dirichlet_bc,ground_index,facets=get_solution_space_and_Dirichlet_BC(Field_calc_param.external_grounding,Field_calc_param.c_c,mesh_sol,subdomains,boundaries_sol,Field_calc_param.element_order,Field_calc_param.EQS_mode,Domains.Active_contacts,Domains.Amp_vector)
    #ground index refers to the ground in .med/.msh file

    logging.critical("dofs: {}".format(max(V_space.dofmap().dofs())+1))
    logging.critical("N of elements: {}".format(mesh_sol.num_cells()))

    #facets = MeshFunction('size_t',mesh_sol,2)
    #facets.set_all(0)
    if Field_calc_param.external_grounding==False:       # otherwise we have it already from get_solution_space_and_Dirichlet_BC()
        facets.array()[boundaries_sol.array()==Domains.Active_contacts[ground_index]]=1
    dsS=Measure("ds",domain=mesh_sol,subdomain_data=facets)
    Ground_surface_size=assemble(1.0*dsS(1))
    dx = Measure("dx",domain=mesh_sol)

    # to solve the Laplace equation div(kappa*grad(phi))=0   (variational form: a(u,v)=L(v))
    start_math=tm.time()
    from FEM_in_spectrum import define_variational_form_and_solve
    phi_sol=define_variational_form_and_solve(V_space,Dirichlet_bc,kappa,Field_calc_param.EQS_mode,Cond_tensor,Solver_type)
    minutes=int((tm.time() - start_math)/60)
    secnds=int(tm.time() - start_math)-minutes*60
    logging.critical("--- assembled and solved in {} min {} sec ---".format(minutes, secnds))

    if Field_calc_param.EQS_mode == 'EQS':
        (phi_r_sol,phi_i_sol)=phi_sol.split(deepcopy=True)
    else:
        phi_r_sol=phi_sol
        phi_i_sol=Function(V_space)
        phi_i_sol.vector()[:] = 0.0

    # get current flowing through the grounded contact and the electric field in the whole domain
    from FEM_in_spectrum import get_current
    J_ground,E_field,E_field_im = get_current(mesh_sol,facets,boundaries_sol,Field_calc_param.element_order,Field_calc_param.EQS_mode,Domains.Active_contacts,kappa,Cond_tensor,phi_r_sol,phi_i_sol,ground_index,get_E_field=True)

    #print("J_ground_unscaled: ",J_ground)
    # If EQS, J_ground is a complex number. If QS, E_field_im is a null function

    # to get current density function which is required for mesh refinement when checking current convergence
    j_dens_real,j_dens_im = get_current_density(mesh_sol,Field_calc_param.element_order,Field_calc_param.EQS_mode,kappa,Cond_tensor,E_field,E_field_im)
    # If QS, j_dens_im is null function

    # will be used for mesh refinement
    j_dens_real_unscaled=j_dens_real.copy(deepcopy=True)
    j_dens_im_unscaled=j_dens_im.copy(deepcopy=True)            #  null function if QS
    import copy
    J_real_unscaled=copy.deepcopy(np.real(J_ground))
    J_im_unscaled=copy.deepcopy(np.imag(J_ground))  # 0 if QS
    # to project the E-field magnitude
    if Field_calc_param.element_order>1:
        V_normE=FunctionSpace(mesh_sol,"CG",Field_calc_param.element_order-1)
    else:
        V_normE=FunctionSpace(mesh_sol,"CG",Field_calc_param.element_order)

    #V_across=max(Domains.Amp_vector[:], key=abs)    #actually, not across, but against ground!!!

    if Field_calc_param.external_grounding==True and (Field_calc_param.c_c==1 or len(Domains.Amp_vector)==1):
        V_max=max(Domains.Amp_vector[:], key=abs)
        V_min=0.0
    elif -1*Domains.Amp_vector[0]==Domains.Amp_vector[1]:     # V_across is needed only for 2 active contact systems
        V_min=-1*abs(Domains.Amp_vector[0])
        V_max=abs(Domains.Amp_vector[0])
    else:
        V_min=min(Domains.Amp_vector[:], key=abs)
        V_max=max(Domains.Amp_vector[:], key=abs)
    V_across=V_max-V_min   # this can be negative


    Vertices_get=read_csv(os.environ['PATIENTDIR']+'/Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', delimiter=' ', header=None)
    Vertices_array=Vertices_get.values

    Phi_ROI=np.zeros((Vertices_array.shape[0],4),float)

    for inx in range(Vertices_array.shape[0]):
        pnt=Point(Vertices_array[inx,0],Vertices_array[inx,1],Vertices_array[inx,2])

        Phi_ROI[inx,0]=Vertices_array[inx,0]
        Phi_ROI[inx,1]=Vertices_array[inx,1]
        Phi_ROI[inx,2]=Vertices_array[inx,2]
        if Field_calc_param.c_c==1:
            phi_r_sol_scaled_on_point=V_across*np.real((phi_r_sol(pnt)+1j*phi_i_sol(pnt))/J_ground)
            phi_i_sol_scaled_on_point=V_across*np.imag((phi_r_sol(pnt)+1j*phi_i_sol(pnt))/J_ground)
            Phi_ROI[inx,3]=np.sqrt(phi_r_sol_scaled_on_point*phi_r_sol_scaled_on_point+phi_i_sol_scaled_on_point*phi_i_sol_scaled_on_point)
        else:
            Phi_ROI[inx,3]=np.sqrt(phi_r_sol(pnt)*phi_r_sol(pnt)+phi_i_sol(pnt)*phi_i_sol(pnt))

    np.savetxt(os.environ['PATIENTDIR']+'/Results_adaptive/Phi_'+str(Field_calc_param.frequenc)+'.csv',  Phi_ROI, delimiter=" ")      # this is amplitude, actually

#    #Probe_of_potential
#    probe_z=np.zeros((100,4),float)
#    for inx in range(100):
#        pnt=Point(75.5,78.5,27.865+inx/10.0)
#        probe_z[inx,0]=75.5
#        probe_z[inx,1]=78.5
#        probe_z[inx,2]=27.865+inx/10.0
#        if Field_calc_param.c_c==1:
#            phi_r_sol_scaled_on_point=V_across*np.real((phi_r_sol(pnt)+1j*phi_i_sol(pnt))/(J_real_unscaled+1j*J_im_unscaled))
#            phi_i_sol_scaled_on_point=V_across*np.imag((phi_r_sol(pnt)+1j*phi_i_sol(pnt))/(J_real_unscaled+1j*J_im_unscaled))
#            probe_z[inx,3]=np.sqrt(phi_r_sol_scaled_on_point*phi_r_sol_scaled_on_point+phi_i_sol_scaled_on_point*phi_i_sol_scaled_on_point)
#        else:
#            probe_z[inx,3]=np.sqrt(phi_r_sol(pnt)*phi_r_sol(pnt)+phi_i_sol(pnt)*phi_i_sol(pnt))
#    np.savetxt('Results_adaptive/Phi_Zprobe'+str(Field_calc_param.frequenc)+'.csv',  probe_z, delimiter=" ")

    #print("Tissue impedance: ", Z_tis)

#=============================================================================#
    if Field_calc_param.c_c==1 or Field_calc_param.CPE==1:

        Z_tissue = V_across/J_ground      # Tissue impedance
        logging.critical("Tissue impedance: {}".format(Z_tissue))

        if Field_calc_param.CPE==1:

            if len(Domains.Amp_vector)>2:
                logging.critical("Currently, CPE can be used only for simulations with two contacts. Please, assign the rest to 'None'")
                raise SystemExit

            from GUI_inp_dict import d as d_cpe
            CPE_param=[d_cpe["K_A"],d_cpe["beta"],d_cpe["K_A_ground"],d_cpe["beta_ground"]]

            from FEM_in_spectrum import get_CPE_corrected_Dirichlet_BC
            Dirichlet_bc_with_CPE,total_impedance=get_CPE_corrected_Dirichlet_BC(Field_calc_param.external_grounding,facets,boundaries_sol,CPE_param,Field_calc_param.EQS_mode,Field_calc_param.frequenc,Field_calc_param.frequenc,Domains.Active_contacts,Domains.Amp_vector,V_across,Z_tissue,V_space)

            logging.critical("Solving for an adjusted potential on contacts to account for CPE")
            start_math=tm.time()
            # to solve the Laplace equation for the adjusted Dirichlet
            phi_sol_CPE=define_variational_form_and_solve(V_space,Dirichlet_bc_with_CPE,kappa,Field_calc_param.EQS_mode,Cond_tensor,Solver_type)
            minutes=int((tm.time() - start_math)/60)
            secnds=int(tm.time() - start_math)-minutes*60
            logging.critical("--- assembled and solved in {} min {} sec ---".format(minutes, secnds))
            if Field_calc_param.EQS_mode=='EQS':
                (phi_r_CPE,phi_i_CPE)=phi_sol_CPE.split(deepcopy=True)
            else:
                phi_r_CPE=phi_sol_CPE
                phi_i_CPE=Function(V_space)
                phi_i_CPE.vector()[:] = 0.0

            # get current flowing through the grounded contact and the electric field in the whole domain
            J_ground_CPE,E_field_CPE,E_field_im_CPE = get_current(mesh_sol,facets,boundaries_sol,Field_calc_param.element_order,Field_calc_param.EQS_mode,Domains.Active_contacts,kappa,Cond_tensor,phi_r_CPE,phi_i_CPE,ground_index,get_E_field=True)
            # If EQS, J_ground is a complex number. If QS, E_field_CPE is a null function

            # to get current density function which is required for mesh refinement when checking current convergence
            j_dens_real_CPE,j_dens_im_CPE = get_current_density(mesh_sol,Field_calc_param.element_order,Field_calc_param.EQS_mode,kappa,Cond_tensor,E_field_CPE,E_field_im_CPE)
            # If QS, j_dens_im is null function

            # will be used for mesh refinement
            j_dens_real_unscaled=j_dens_real_CPE.copy(deepcopy=True)
            j_dens_im_unscaled=j_dens_im_CPE.copy(deepcopy=True)
            J_real_unscaled=copy.deepcopy(np.real(J_ground))
            J_im_unscaled=copy.deepcopy(np.imag(J_ground))

            E_norm=project(sqrt(inner(E_field_CPE,E_field_CPE)+inner(E_field_im_CPE,E_field_im_CPE)),V_normE,solver_type="cg", preconditioner_type="amg")
            max_E=E_norm.vector().max()

            file=File(os.environ['PATIENTDIR']+'/Results_adaptive/E_ampl_'+str(Field_calc_param.EQS_mode)+'.pvd')
            file<<E_norm,mesh_sol
            file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_Phi_r_field_'+str(Field_calc_param.EQS_mode)+'.pvd')
            file<<phi_r_CPE,mesh_sol
            if Field_calc_param.EQS_mode=='EQS':
                file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_Phi_im_field_'+str(Field_calc_param.EQS_mode)+'.pvd')
                file<<phi_i_CPE,mesh_sol

            return phi_r_CPE,phi_i_CPE,E_field_CPE,E_field_im_CPE,max_E,J_real_unscaled,J_im_unscaled,j_dens_real_unscaled,j_dens_im_unscaled

        if Field_calc_param.c_c==1:
            if Field_calc_param.EQS_mode=='EQS':     # For EQS, we need to scale the potential on boundaries (because the error is absolute) and recompute field, etc. Maybe we can scale them also directly?
                Dirichlet_bc_scaled=[]
                for bc_i in range(len(Domains.Active_contacts)):          #CPE estimation is valid only for one activa and one ground contact configuration
                    if Field_calc_param.EQS_mode=='EQS':
                        if Domains.Amp_vector[bc_i]!=0.0:
                            Active_with_CC=V_across*V_across/J_ground          #(impedance * current through the contact (V_across coincides with the assigned current magnitude))
                            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(0), np.real(Active_with_CC), boundaries_sol,Domains.Active_contacts[bc_i]))
                            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(1), np.imag(Active_with_CC), boundaries_sol,Domains.Active_contacts[bc_i]))
                        else:
                            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(0), Constant(0.0), boundaries_sol,Domains.Active_contacts[bc_i]))
                            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(1), Constant(0.0), boundaries_sol,Domains.Active_contacts[bc_i]))

                if Field_calc_param.external_grounding==True:
                    if Field_calc_param.EQS_mode == 'EQS':
                        Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(0),0.0,facets,1))
                        Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(1),0.0,facets,1))
                    else:
                        Dirichlet_bc_scaled.append(DirichletBC(V_space,0.0,facets,1))



                logging.critical("Solving for a scaled potential on contacts (to match the desired current)")
                start_math=tm.time()
                # to solve the Laplace equation for the adjusted Dirichlet
                phi_sol_scaled=define_variational_form_and_solve(V_space,Dirichlet_bc_scaled,kappa,Field_calc_param.EQS_mode,Cond_tensor,Solver_type)
                minutes=int((tm.time() - start_math)/60)
                secnds=int(tm.time() - start_math)-minutes*60
                logging.critical("--- assembled and solved in {} min {} sec ---".format(minutes, secnds))

                (phi_r_sol_scaled,phi_i_sol_scaled)=phi_sol_scaled.split(deepcopy=True)

                # get current flowing through the grounded contact and the electric field in the whole domain
                J_ground_scaled,E_field_scaled,E_field_im_scaled = get_current(mesh_sol,facets,boundaries_sol,Field_calc_param.element_order,Field_calc_param.EQS_mode,Domains.Active_contacts,kappa,Cond_tensor,phi_r_sol_scaled,phi_i_sol_scaled,ground_index,get_E_field=True)
                # If EQS, J_ground is a complex number. If QS, E_field_im is 0
            else:   # here we can simply scale the potential in the domain and recompute the E-field
                phi_i_sol_scaled=Function(V_space)
                phi_i_sol_scaled.vector()[:] = 0.0
                phi_r_sol_scaled=Function(V_space)
                phi_r_sol_scaled.vector()[:]=V_across*phi_r_sol.vector()[:]/J_ground

                J_ground_scaled,E_field_scaled,E_field_im_scaled = get_current(mesh_sol,facets,boundaries_sol,Field_calc_param.element_order,Field_calc_param.EQS_mode,Domains.Active_contacts,kappa,Cond_tensor,phi_r_sol_scaled,phi_i_sol_scaled,ground_index,get_E_field=True)
                #E_field_im_scale is a null function

            E_norm=project(sqrt(inner(E_field_scaled,E_field_scaled)+inner(E_field_im_scaled,E_field_im_scaled)),V_normE,solver_type="cg", preconditioner_type="amg")
            max_E=E_norm.vector().max()

            file=File(os.environ['PATIENTDIR']+'/Results_adaptive/E_ampl_'+str(Field_calc_param.EQS_mode)+'.pvd')
            file<<E_norm,mesh_sol
            file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_Phi_r_field_'+str(Field_calc_param.EQS_mode)+'.pvd')
            file<<phi_r_sol_scaled,mesh_sol
            if Field_calc_param.EQS_mode=='EQS':
                file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_Phi_im_field_'+str(Field_calc_param.EQS_mode)+'.pvd')
                file<<phi_i_sol_scaled,mesh_sol

            return phi_r_sol_scaled,phi_i_sol_scaled,E_field_scaled,E_field_im_scaled,max_E,J_real_unscaled,J_im_unscaled,j_dens_real_unscaled,j_dens_im_unscaled

    else:
        E_norm=project(sqrt(inner(E_field,E_field)+inner(E_field_im,E_field_im)),V_normE,solver_type="cg", preconditioner_type="amg")
        max_E=E_norm.vector().max()
        file=File(os.environ['PATIENTDIR']+'/Results_adaptive/E_ampl_'+str(Field_calc_param.EQS_mode)+'.pvd')
        file<<E_norm,mesh_sol
        file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_Phi_r_field_'+str(Field_calc_param.EQS_mode)+'.pvd')
        file<<phi_r_sol,mesh_sol
        if Field_calc_param.EQS_mode=='EQS':
            file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_Phi_im_field_'+str(Field_calc_param.EQS_mode)+'.pvd')
            file<<phi_i_sol,mesh_sol

        return phi_r_sol,phi_i_sol,E_field,E_field_im,max_E,J_real_unscaled,J_im_unscaled,j_dens_real_unscaled,j_dens_im_unscaled

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
