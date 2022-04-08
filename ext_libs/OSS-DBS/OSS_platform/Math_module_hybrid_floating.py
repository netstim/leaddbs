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

import warnings
warnings.filterwarnings("ignore")

from tissue_dielectrics import DielectricProperties

parameters['linear_algebra_backend']='PETSc'
set_log_active(False)   #turns off debugging info

def load_scaled_cond_tensor(xx,xy,xz,yy,yz,zz,mesh_tensor):

    # Code for C++ evaluation of conductivity
    conductivity_code = """

    #include <pybind11/pybind11.h>
    #include <pybind11/eigen.h>
    namespace py = pybind11;

    #include <dolfin/function/Expression.h>
    #include <dolfin/mesh/MeshFunction.h>

    class Conductivity : public dolfin::Expression
    {
    public:

      // Create expression with 6 components
      Conductivity() : dolfin::Expression(6) {}

      // Function for evaluating expression on each cell
      void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x, const ufc::cell& cell) const override
      {
        const uint cell_index = cell.index;
        values[0] = (*c00)[cell_index];
        values[1] = (*c01)[cell_index];
        values[2] = (*c02)[cell_index];
        values[3] = (*c11)[cell_index];
        values[4] = (*c12)[cell_index];
        values[5] = (*c22)[cell_index];
      }

      // The data stored in mesh functions
      std::shared_ptr<dolfin::MeshFunction<double>> c00;
      std::shared_ptr<dolfin::MeshFunction<double>> c01;
      std::shared_ptr<dolfin::MeshFunction<double>> c02;
      std::shared_ptr<dolfin::MeshFunction<double>> c11;
      std::shared_ptr<dolfin::MeshFunction<double>> c12;
      std::shared_ptr<dolfin::MeshFunction<double>> c22;

    };

    PYBIND11_MODULE(SIGNATURE, m)
    {
      py::class_<Conductivity, std::shared_ptr<Conductivity>, dolfin::Expression>
        (m, "Conductivity")
        .def(py::init<>())
        .def_readwrite("c00", &Conductivity::c00)
        .def_readwrite("c01", &Conductivity::c01)
        .def_readwrite("c02", &Conductivity::c02)
        .def_readwrite("c11", &Conductivity::c11)
        .def_readwrite("c12", &Conductivity::c12)
        .def_readwrite("c22", &Conductivity::c22);
    }

    """

    hdf = HDF5File(mesh_tensor.mpi_comm(), os.environ['PATIENTDIR']+"/Results_adaptive/Mesh_to_solve.h5", "r")
    hdf.read(xx, "/c00")
    hdf.read(xy, "/c01")
    hdf.read(xz, "/c02")
    hdf.read(yy, "/c11")
    hdf.read(yz, "/c12")
    hdf.read(zz, "/c22")
    hdf.close()

    c = CompiledExpression(compile_cpp_code(conductivity_code).Conductivity(),
                   c00=xx, c01=xy, c02=xz, c11=yy, c12=yz, c22=zz, degree=0)

    C_tensor = as_matrix(((c[0], c[1], c[2]), (c[1], c[3], c[4]),(c[2],c[4],c[5])))

    return C_tensor


#if calculating with MPI, the dielectic properties (kappa) and the scaled tensor were already prepared
def get_field_with_floats(external_grounding,mesh_sol,active_index,Domains,subdomains,boundaries_sol,default_material,element_order,anisotropy,frequenc,Laplace_mode,Solver_type,calc_with_MPI=False,kappa=False):


    set_log_active(False)   #turns off debugging info
    parameters['linear_algebra_backend']='PETSc'

    if Laplace_mode == 'EQS':
        Solver_type = 'MUMPS'  # always direct solver for EQS when multiple floats

    if calc_with_MPI==False:

        [cond_GM, perm_GM]=DielectricProperties(3).get_dielectrics(frequenc)        #3 for grey matter and so on (numeration as in voxel_data)
        [cond_WM, perm_WM]=DielectricProperties(2).get_dielectrics(frequenc)
        [cond_CSF, perm_CSF]=DielectricProperties(1).get_dielectrics(frequenc)

        [cond_default,perm_default]=DielectricProperties(default_material).get_dielectrics(frequenc)

        from GUI_inp_dict import d as d_encap
        [cond_encap, perm_encap]=DielectricProperties(d_encap['encap_tissue_type']).get_dielectrics(frequenc)
        cond_encap=cond_encap*d_encap['encap_scaling_cond']
        perm_encap=perm_encap*d_encap['encap_scaling_perm']

        conductivities=[cond_default,cond_GM,cond_WM,cond_CSF,cond_encap]
        rel_permittivities=[perm_default,perm_GM,perm_WM,perm_CSF,perm_encap]

        # to get conductivity (and permittivity if EQS formulation) mapped accrodingly to the subdomains. k_val_r is just a list of conductivities (S/mm!) in a specific order to scale the cond. tensor
        from FEM_in_spectrum import get_dielectric_properties_from_subdomains
        kappa,k_val_r=get_dielectric_properties_from_subdomains(mesh_sol,subdomains,Laplace_mode,Domains.Float_contacts,conductivities,rel_permittivities,frequenc)

    if anisotropy==1:
        # order xx,xy,xz,yy,yz,zz
        c00 = MeshFunction("double", mesh_sol, 3, 0.0)
        c01 = MeshFunction("double", mesh_sol, 3, 0.0)
        c02 = MeshFunction("double", mesh_sol, 3, 0.0)
        c11 = MeshFunction("double", mesh_sol, 3, 0.0)
        c12 = MeshFunction("double", mesh_sol, 3, 0.0)
        c22 = MeshFunction("double", mesh_sol, 3, 0.0)

        if calc_with_MPI==True:        #load predefined Tensor
            Cond_tensor = load_scaled_cond_tensor(c00,c01,c02,c11,c12,c22,mesh_sol)
        else:
            # load the unscaled diffusion tensor (should be normalized beforehand)
            hdf = HDF5File(mesh_sol.mpi_comm(), os.environ['PATIENTDIR']+"/Results_adaptive/Tensors_to_solve_num_el_"+str(mesh_sol.num_cells())+".h5", "r")
            hdf.read(c00, "/c00")
            hdf.read(c01, "/c01")
            hdf.read(c02, "/c02")
            hdf.read(c11, "/c11")
            hdf.read(c12, "/c12")
            hdf.read(c22, "/c22")
            hdf.close()

            unscaled_tensor=[c00,c01,c02,c11,c12,c22]

            # to get tensor scaled by the conductivity map (twice send frequenc to always get unscaled ellipsoid tensor for visualization)
            from FEM_in_spectrum import get_scaled_cond_tensor
            Cond_tensor=get_scaled_cond_tensor(mesh_sol,subdomains,frequenc,frequenc,unscaled_tensor,k_val_r)
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
            #print("one")
            if Laplace_mode == 'EQS':
                dirichlet_bc.append(DirichletBC(V_space.sub(0), Domains.Amp_vector[bc_i], boundaries_sol,Domains.Active_contacts[bc_i]))
                dirichlet_bc.append(DirichletBC(V_space.sub(1), Constant(0.0), boundaries_sol,Domains.Active_contacts[bc_i]))
            else:
                dirichlet_bc.append(DirichletBC(V_space, Domains.Amp_vector[bc_i], boundaries_sol,Domains.Active_contacts[bc_i]))

            if bc_i==active_index:
                facets_active.array()[boundaries_sol.array()==Domains.Active_contacts[bc_i]]=1
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
        return Float_potentials_real,Float_potentials_imag,J_real,J_im
    else:
        if Cond_tensor!=False:
            j_dens_real_contact = dot(Cond_tensor*E_field,-1*n)('-')*dsS_int(1)
        else:
            j_dens_real_contact = dot(kappa[0]*E_field,-1*n)('-')*dsS_int(1)

    J_real=assemble(j_dens_real_contact)

    #print("Shape float potentials: ",Float_potentials_real.shape[0])

    return Float_potentials_real,0.0,J_real,0.0

#if calculating with MPI, the dielectic properties (kappa) and the scaled tensor were already prepared
def get_field_with_scaled_BC(external_grounding,mesh_sol,Domains,Phi_scaled,subdomains,boundaries_sol,default_material,element_order,Laplace_mode,anisotropy,frequenc,Solver_type,calc_with_MPI=False,kappa=False):

    set_log_active(False)   #turns off debugging info
    parameters['linear_algebra_backend']='PETSc'

    if calc_with_MPI==False or MPI.comm_world.rank==1:
        logging.critical("Calculating field with scaled voltage on contacts")

    if calc_with_MPI==False:

        [cond_GM, perm_GM]=DielectricProperties(3).get_dielectrics(frequenc)        #3 for grey matter and so on (numeration as in voxel_data)
        [cond_WM, perm_WM]=DielectricProperties(2).get_dielectrics(frequenc)
        [cond_CSF, perm_CSF]=DielectricProperties(1).get_dielectrics(frequenc)

        [cond_default,perm_default]=DielectricProperties(default_material).get_dielectrics(frequenc)

        from GUI_inp_dict import d as d_encap
        [cond_encap, perm_encap]=DielectricProperties(d_encap['encap_tissue_type']).get_dielectrics(frequenc)
        cond_encap=cond_encap*d_encap['encap_scaling_cond']
        perm_encap=perm_encap*d_encap['encap_scaling_perm']

        conductivities=[cond_default,cond_GM,cond_WM,cond_CSF,cond_encap]
        rel_permittivities=[perm_default,perm_GM,perm_WM,perm_CSF,perm_encap]

        # to get conductivity (and permittivity if EQS formulation) mapped accrodingly to the subdomains. k_val_r is just a list of conductivities (S/mm!) in a specific order to scale the cond. tensor
        from FEM_in_spectrum import get_dielectric_properties_from_subdomains
        kappa,k_val_r=get_dielectric_properties_from_subdomains(mesh_sol,subdomains,Laplace_mode,Domains.Float_contacts,conductivities,rel_permittivities,frequenc)
        file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_subdomains_map.pvd')
        file<<subdomains
        file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_conductivity_map.pvd')
        file<<kappa[0]
        if Laplace_mode == 'EQS':
            file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_permittivity_map.pvd')
            file<<kappa[1]

    if anisotropy==1:
        # order xx,xy,xz,yy,yz,zz
        c00 = MeshFunction("double", mesh_sol, 3, 0.0)
        c01 = MeshFunction("double", mesh_sol, 3, 0.0)
        c02 = MeshFunction("double", mesh_sol, 3, 0.0)
        c11 = MeshFunction("double", mesh_sol, 3, 0.0)
        c12 = MeshFunction("double", mesh_sol, 3, 0.0)
        c22 = MeshFunction("double", mesh_sol, 3, 0.0)

        if calc_with_MPI==True:        #load predefined Tensor
            Cond_tensor = load_scaled_cond_tensor(c00,c01,c02,c11,c12,c22,mesh_sol)
        else:
            # load the unscaled diffusion tensor (should be normalized beforehand)
            hdf = HDF5File(mesh_sol.mpi_comm(), os.environ['PATIENTDIR']+"/Results_adaptive/Tensors_to_solve_num_el_"+str(mesh_sol.num_cells())+".h5", "r")
            hdf.read(c00, "/c00")
            hdf.read(c01, "/c01")
            hdf.read(c02, "/c02")
            hdf.read(c11, "/c11")
            hdf.read(c12, "/c12")
            hdf.read(c22, "/c22")
            hdf.close()

            unscaled_tensor=[c00,c01,c02,c11,c12,c22]

            # to get tensor scaled by the conductivity map (twice send frequenc to always get unscaled ellipsoid tensor for visualization)
            from FEM_in_spectrum import get_scaled_cond_tensor
            Cond_tensor=get_scaled_cond_tensor(mesh_sol,subdomains,frequenc,frequenc,unscaled_tensor,k_val_r,plot_tensors=True)
    else:
        Cond_tensor=False  #just to initialize

    from FEM_in_spectrum import get_solution_space_and_Dirichlet_BC
    V_space,facets=get_solution_space_and_Dirichlet_BC(external_grounding,1,mesh_sol,subdomains,boundaries_sol,element_order,Laplace_mode,Domains.Active_contacts,Phi_scaled,only_space=True)

    Dirichlet_bc_scaled=[]
    if calc_with_MPI==False or MPI.comm_world.rank==1:
        logging.critical("Scaled complex potential on contacts: ")
        for j in range(Phi_scaled.shape[0]):
            logging.critical("{}".format(Phi_scaled[j]))
    for bc_i in range(len(Domains.Active_contacts)):
        if Laplace_mode == 'EQS':
            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(0), np.real(Phi_scaled[bc_i]), boundaries_sol,Domains.Active_contacts[bc_i]))
            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(1), np.imag(Phi_scaled[bc_i]), boundaries_sol,Domains.Active_contacts[bc_i]))
        else:
            Dirichlet_bc_scaled.append(DirichletBC(V_space, Phi_scaled[bc_i], boundaries_sol,Domains.Active_contacts[bc_i]))

    if external_grounding==True:
        if Laplace_mode == 'EQS':
            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(0),0.0,facets,1))
            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(1),0.0,facets,1))
        else:
            Dirichlet_bc_scaled.append(DirichletBC(V_space,0.0,facets,1))


        #facets.array()[boundaries_sol.array()==Domains.Active_contacts[bc_i]]=bc_i+1

    # to solve the Laplace equation div(kappa*grad(phi))=0   (variational form: a(u,v)=L(v))
    from FEM_in_spectrum import define_variational_form_and_solve
    phi_sol=define_variational_form_and_solve(V_space,Dirichlet_bc_scaled,kappa,Laplace_mode,Cond_tensor,Solver_type)

    if Laplace_mode=='EQS':
        (phi_r_sol,phi_i_sol)=phi_sol.split(deepcopy=True)
    else:
        phi_r_sol=phi_sol
        phi_i_sol=Function(V_space)
        phi_i_sol.vector()[:] = 0.0

    # to get manual projections for the E-field
    from FEM_in_spectrum_multicontact import get_E_field
    E_field,E_field_im=get_E_field(mesh_sol,element_order,Laplace_mode,phi_r_sol,phi_i_sol)
    #if QS, E_field_im is a null function

    # to get current density function which is required for mesh refinement when checking current convergence
    from Math_module_hybrid import get_current_density
    j_dens_real,j_dens_im = get_current_density(mesh_sol,element_order,Laplace_mode,kappa,Cond_tensor,E_field,E_field_im)
    # If QS, j_dens_im is null function

    # to project the E-field magnitude
    if element_order>1:
        V_normE=FunctionSpace(mesh_sol,"CG",element_order-1)
    else:
        V_normE=FunctionSpace(mesh_sol,"CG",element_order)

    # to get current on the active contacts (inlcuding the ground)
    from FEM_in_spectrum_multicontact import get_current_on_multiple_contacts
    J_r_contacts,J_im_contacts = get_current_on_multiple_contacts(external_grounding,facets,mesh_sol,boundaries_sol,Laplace_mode,Domains.Active_contacts,Phi_scaled,E_field,E_field_im,kappa,Cond_tensor)
    # J_currents_imag is a zero array if 'QS' mode

    J_r_max=J_r_contacts.max()      # we need it to filter out contacts with very small currents (they might have very high quasi-impedance)
    ind_high_current=[]

    if calc_with_MPI==False:
        Vertices_get=read_csv(os.environ['PATIENTDIR']+'/Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', delimiter=' ', header=None)
        Vertices_array=Vertices_get.values

        Phi_ROI=np.zeros((Vertices_array.shape[0],4),float)

        for inx in range(Vertices_array.shape[0]):
            pnt=Point(Vertices_array[inx,0],Vertices_array[inx,1],Vertices_array[inx,2])

            Phi_ROI[inx,0]=Vertices_array[inx,0]
            Phi_ROI[inx,1]=Vertices_array[inx,1]
            Phi_ROI[inx,2]=Vertices_array[inx,2]
            Phi_ROI[inx,3]=np.sqrt(phi_r_sol(pnt)*phi_r_sol(pnt)+phi_i_sol(pnt)*phi_i_sol(pnt))

        np.savetxt(os.environ['PATIENTDIR']+'/Results_adaptive/Phi_'+str(frequenc)+'.csv',  Phi_ROI, delimiter=" ")      # this is amplitude, actually

    if external_grounding==True:
        Quasi_imp_real=np.zeros(len(Domains.Active_contacts)+1,float)       #not really, but gives an idea
        Quasi_imp_im=np.zeros(len(Domains.Active_contacts)+1,float)       #not really, but gives an idea
    else:
        Quasi_imp_real=np.zeros(len(Domains.Active_contacts),float)       #not really, but gives an idea
        Quasi_imp_im=np.zeros(len(Domains.Active_contacts),float)       #not really, but gives an idea

    for bc_i in range(len(Domains.Active_contacts)):
        if calc_with_MPI==False or MPI.comm_world.rank==1:
            logging.critical("J on contact {}: {} A".format(Domains.Active_on_lead[bc_i],J_r_contacts[bc_i]+1j*J_im_contacts[bc_i]))

        Quasi_imp_real[bc_i]=np.real(Phi_scaled[bc_i]/(J_r_contacts[bc_i]+1j*J_im_contacts[bc_i]))
        if Laplace_mode=='EQS':
            Quasi_imp_im[bc_i]=np.imag(Phi_scaled[bc_i]/(J_r_contacts[bc_i]+1j*J_im_contacts[bc_i]))

    if external_grounding==True:        # in the future I should choose a better metric
        if calc_with_MPI==False or MPI.comm_world.rank==1:
            logging.critical("J on external grounding :{} A".format(J_r_contacts[-1]+1j*J_im_contacts[-1]))
        Quasi_imp_real[-1]=np.real(1.0/(J_r_contacts[-1]+1j*J_im_contacts[-1]))
        if Laplace_mode=='EQS':
            Quasi_imp_im[-1]=np.imag(1.0/(J_r_contacts[-1]+1j*J_im_contacts[-1]))

        if J_r_contacts[-1]>=0.1*J_r_max:
            ind_high_current.append(len(Domains.Active_contacts))      #because it's the last

    for bc_i in range(len(Domains.Active_contacts)):
        if J_r_contacts[bc_i]>=0.1*J_r_max:
            ind_high_current.append(bc_i)



    Quasi_imp_real_total=np.sum(abs(Quasi_imp_real[ind_high_current]))
    Quasi_imp_im_total=np.sum(abs(Quasi_imp_im[ind_high_current]))

    if calc_with_MPI == True:
        J_Vector=Vector(MPI.comm_self,2)
        J_Vector.set_local(np.array([Quasi_imp_real_total,Quasi_imp_im_total],dtype=np.float64))
        Hdf=HDF5File(mesh_sol.mpi_comm(), os.environ['PATIENTDIR']+"/Results_adaptive/Solution_"+str(np.round(frequenc,6))+".h5", "w")
        Hdf.write(mesh_sol, "mesh_sol")
        Hdf.write(phi_sol, "solution_phi_full")
        Hdf.write(E_field, "solution_E_field")
        Hdf.write(E_field_im, "solution_E_field_im")
        Hdf.write(j_dens_real, "solution_j_real")
        Hdf.write(j_dens_im, "solution_j_im")
        Hdf.write(J_Vector, "/J_Vector")
        Hdf.close()

        return True
    else:
        logging.critical("Quasi_impedance (for current check): {}".format(sqrt(Quasi_imp_real_total**2+Quasi_imp_im_total**2)))
        E_norm=project(sqrt(inner(E_field,E_field)+inner(E_field_im,E_field_im)),V_normE,solver_type="cg", preconditioner_type="amg")
        max_E=E_norm.vector().max()
        if calc_with_MPI==False or MPI.comm_world.rank==1:
            file=File(os.environ['PATIENTDIR']+'/Results_adaptive/E_ampl_'+str(Laplace_mode)+'.pvd')
            file<<E_norm,mesh_sol
            file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_Phi_r_field_'+str(Laplace_mode)+'.pvd')
            file<<phi_r_sol,mesh_sol
            if Laplace_mode=='EQS':
                file=File(os.environ['PATIENTDIR']+'/Results_adaptive/Last_Phi_im_field_'+str(Laplace_mode)+'.pvd')
                file<<phi_i_sol,mesh_sol

        return (phi_r_sol,phi_i_sol,E_field,E_field_im,max_E,Quasi_imp_real_total,Quasi_imp_im_total,j_dens_real,j_dens_im)


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


def compute_field_with_superposition(mesh_sol,Domains,subdomains_assigned,subdomains,boundaries_sol,Field_calc_param):

    start_math=tm.time()

    from Math_module_hybrid import choose_solver_for_me
    if Field_calc_param.Solver_type=='Default':
        Solver_type=choose_solver_for_me(Field_calc_param.EQS_mode,Domains.Float_contacts)    #choses solver basing on the Laplace formulation and whether the floating conductors are used
    else:
        Solver_type=Field_calc_param.Solver_type      # just get the solver directly
    #IMPORTANT: for get_field_with_floats when solving EQS we always use direct solver MUMPS for stability issues (multiple floating conductors)

    if Field_calc_param.element_order==1:
        logging.critical("Selected element_order (1st) is too low for current-controlled stimulation, increasing to 2nd")
        Field_calc_param.element_order=2

    logging.critical("Computing field with superposition on mesh with {} elements".format(mesh_sol.num_cells()))
    logging.critical("{} computations are required for the iteration".format(len(Domains.Amp_vector)))

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
            if Domains.Active_on_lead[i] == Domains.Float_on_lead[j]:       #find the index of the floating conductor (in .med/.msh file) for the active contact (i)

                # to solve "one active contact (i) vs ground" system, get potentials on the rest of the active contacts (which are put to floating condcutors), get current on the active contact
                if Field_calc_param.EQS_mode == 'EQS':
                    phi_r_floating[glob_counter,:],phi_i_floating[glob_counter,:],J_real_current_contacts[glob_counter],J_im_current_contacts[glob_counter]=get_field_with_floats(Field_calc_param.external_grounding,mesh_sol,i,Domains,subdomains,boundaries_sol,Field_calc_param.default_material,Field_calc_param.element_order,Field_calc_param.anisotropy,Field_calc_param.frequenc,Field_calc_param.EQS_mode,Solver_type)
                else:
                    phi_r_floating[glob_counter,:],__,J_real_current_contacts[glob_counter],__=get_field_with_floats(Field_calc_param.external_grounding,mesh_sol,i,Domains,subdomains,boundaries_sol,Field_calc_param.default_material,Field_calc_param.element_order,Field_calc_param.anisotropy,Field_calc_param.frequenc,Field_calc_param.EQS_mode,Solver_type)

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

    # quasi_imp is a metric we use to assess the current convergence (an analog to unscaled current).
    phi_r_sol,phi_i_sol,Field_real,Field_imag,max_E,quasi_imp_real,quasi_imp_im,j_dens_real,j_dens_im=get_field_with_scaled_BC(Field_calc_param.external_grounding,mesh_sol,Domains,scaled_phi,subdomains,boundaries_sol,Field_calc_param.default_material,Field_calc_param.element_order,Field_calc_param.EQS_mode,Field_calc_param.anisotropy,Field_calc_param.frequenc,Solver_type)

    minutes=int((tm.time() - start_math)/60)
    secnds=int(tm.time() - start_math)-minutes*60

    logging.critical("----- Field with superposition was calculated in {} min {} sec -----\n".format(minutes, secnds))
    logging.critical("__________________________________")

    return phi_r_sol,phi_i_sol,Field_real,Field_imag,max_E,quasi_imp_real,quasi_imp_im,j_dens_real,j_dens_im
