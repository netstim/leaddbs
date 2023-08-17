# -*- coding: utf-8 -*-
"""
Created on Sun Aug 19 13:23:54 2018

@author: butenko
"""

from dolfin import *
from pandas import read_csv
from tissue_dielectrics import DielectricProperties
import numpy as np
import os
import subprocess
import pickle

import warnings
warnings.filterwarnings("ignore")

import importlib

# additions
#import CSF_refinement as CSF_imp
import time as tm

#def run_mpi(np):
#    print("Code runs in parallel, {} processors.".format(np))
#    subprocess.call(["mpirun", "-np", "{}".format(np), "python3", "Math_module_MPI_only_FEniCS.py"])

parameters["allow_extrapolation"]=True
parameters['linear_algebra_backend']='PETSc'
set_log_active(False)   #turns off debugging info


#def norm_E(E_r,E_im):
#    return sqrt(E_r[0] ** 2 + E_r[1] ** 2+ E_r[2] ** 2+E_im[0] ** 2 + E_im[1] ** 2+ E_im[2] ** 2)


#this function runs in series
#def get_solutions(EQS_form,Domains,boundaries_sol,frequency,el_order):         #know where to look for the solutions in h5py files
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

    if MPI.comm_world.rank==0:
        print("Unscaled current on the contact (to check current convergence): ",J_real,J_im)

        minutes=int((tm.time() - start_reassamble)/60)
        secnds=int(tm.time() - start_reassamble)-minutes*60
        print("--- solution reassambled in ",minutes," min ",secnds," s ---")

    return phi_r_sol,phi_i_sol,E_field,E_field_im,max_E,J_real,J_im,j_dens_real,j_dens_im

def get_field_on_points(phi_r,phi_i,c_c,J_r,J_i):

#================If we check on predefined vertices===========================#

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

def get_cond_tensor(mesh):

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

    # Define conductivity expression and matrix
    c00 = MeshFunction("double", mesh, 3, 0.0)
    c01 = MeshFunction("double", mesh, 3, 0.0)
    c02 = MeshFunction("double", mesh, 3, 0.0)
    c11 = MeshFunction("double", mesh, 3, 0.0)
    c12 = MeshFunction("double", mesh, 3, 0.0)
    c22 = MeshFunction("double", mesh, 3, 0.0)

    hdf = HDF5File(mesh.mpi_comm(), os.environ['PATIENTDIR']+"/Results_adaptive/Mesh_to_solve.h5", "r")
    hdf.read(c00, "/c00")
    hdf.read(c01, "/c01")
    hdf.read(c02, "/c02")
    hdf.read(c11, "/c11")
    hdf.read(c12, "/c12")
    hdf.read(c22, "/c22")
    hdf.close()

    c = CompiledExpression(compile_cpp_code(conductivity_code).Conductivity(),
                   c00=c00, c01=c01, c02=c02, c11=c11, c12=c12, c22=c22, degree=0)

    C = as_matrix(((c[0], c[1], c[2]), (c[1], c[3], c[4]),(c[2],c[4],c[5])))

    return C


def compute_field(mesh_sol,Domains,subdomains,boundaries_sol,kappa_r,Field_calc_param,kappa_i=False):

    set_log_active(False)   #turns off debugging info
    if MPI.comm_world.rank==1:
        print("_________________________")
    parameters['linear_algebra_backend']='PETSc'

    if Field_calc_param.EQS_mode == 'EQS':
        kappa = [kappa_r, kappa_i]
    else:
        kappa = [kappa_r]

    if Field_calc_param.anisotropy==1:
        Cond_tensor=get_cond_tensor(mesh_sol)       # unlike get_scaled_cond_tensor, this function does not scale tensor with conductivity (as it was already scaled)
    else:
        Cond_tensor=False       #just to initialize

    from Math_module_hybrid import choose_solver_for_me
    if Field_calc_param.Solver_type=='Default':
        Solver_type=choose_solver_for_me(Field_calc_param.EQS_mode,Domains.Float_contacts)    #choses solver basing on the Laplace formulation and whether the floating conductors are used
    else:
        Solver_type=Field_calc_param.Solver_type      # just get the solver directly

    #In case of current-controlled stimulation, Dirichlet_bc or the whole potential distribution will be scaled afterwards (due to the system's linearity)
    from FEM_in_spectrum import get_solution_space_and_Dirichlet_BC
    V_space,Dirichlet_bc,ground_index,facets=get_solution_space_and_Dirichlet_BC(Field_calc_param.external_grounding,Field_calc_param.c_c,mesh_sol,subdomains,boundaries_sol,Field_calc_param.element_order,Field_calc_param.EQS_mode,Domains.Contacts,Domains.Amp_vector)
    #ground index refers to the ground in .med/.msh file

    #facets = MeshFunction('size_t',mesh_sol,2)
    #facets.set_all(0)
    if Field_calc_param.external_grounding==False:
        facets.array()[boundaries_sol.array()==Domains.Contacts[ground_index]]=1
    dsS=Measure("ds",domain=mesh_sol,subdomain_data=facets)
    Ground_surface_size=assemble(1.0*dsS(1))
    dx = Measure("dx",domain=mesh_sol)

    # to solve the Laplace equation div(kappa*grad(phi))=0   (variational form: a(u,v)=L(v))
    start_math=tm.time()
    from FEM_in_spectrum import define_variational_form_and_solve
    phi_sol=define_variational_form_and_solve(V_space,Dirichlet_bc,kappa,Field_calc_param.EQS_mode,Cond_tensor,Solver_type)
    if MPI.comm_world.rank==1:
        minutes=int((tm.time() - start_math)/60)
        secnds=int(tm.time() - start_math)-minutes*60
        print("--- assembled and solved in ",minutes," min ",secnds," s ---")



    if Field_calc_param.EQS_mode == 'EQS':
        (phi_r_sol,phi_i_sol)=phi_sol.split(deepcopy=True)
    else:
        phi_r_sol=phi_sol
        phi_i_sol=Function(V_space)
        phi_i_sol.vector()[:] = 0.0

    if MPI.comm_world.rank==1:
        print("dofs on 2nd process: ",(max(V_space.dofmap().dofs())+1))

    # get current flowing through the grounded contact and the electric field in the whole domain
    from FEM_in_spectrum import get_current
    J_ground,E_field,E_field_im = get_current(mesh_sol,facets,boundaries_sol,Field_calc_param.element_order,Field_calc_param.EQS_mode,Domains.Contacts,kappa,Cond_tensor,phi_r_sol,phi_i_sol,ground_index,get_E_field=True)
    # If EQS, J_ground is a complex number. If QS, E_field_im is a null function

    # to get current density function which is required for mesh refinement when checking current convergence
    from Math_module_hybrid import get_current_density
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


### Do not probe anything when inside MPI process!
#    Vertices_get=read_csv('Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', delimiter=' ', header=None)
#    Vertices_array=Vertices_get.values

#    Phi_ROI=np.zeros((Vertices_array.shape[0],4),float)
#
#    for inx in range(Vertices_array.shape[0]):
#        pnt=Point(Vertices_array[inx,0],Vertices_array[inx,1],Vertices_array[inx,2])
#
#        Phi_ROI[inx,0]=Vertices_array[inx,0]
#        Phi_ROI[inx,1]=Vertices_array[inx,1]
#        Phi_ROI[inx,2]=Vertices_array[inx,2]
#        if Field_calc_param.c_c==1:
#            phi_r_sol_scaled_on_point=V_across*np.real((phi_r_sol(pnt)+1j*phi_i_sol(pnt))/(J_real+1j*J_im))
#            phi_i_sol_scaled_on_point=V_across*np.imag((phi_r_sol(pnt)+1j*phi_i_sol(pnt))/(J_real+1j*J_im))
#            Phi_ROI[inx,3]=np.sqrt(phi_r_sol_scaled_on_point*phi_r_sol_scaled_on_point+phi_i_sol_scaled_on_point*phi_i_sol_scaled_on_point)
#        else:
#            Phi_ROI[inx,3]=np.sqrt(phi_r_sol(pnt)*phi_r_sol(pnt)+phi_i_sol(pnt)*phi_i_sol(pnt))
#
#    np.savetxt('Results_adaptive/Phi_'+str(Field_calc_param.frequenc)+'.csv',  Phi_ROI, delimiter=" ")

    #save the results
    if Field_calc_param.c_c!=1 and Field_calc_param.CPE!=1:
        #Just to save total currunt through the ground in FEniCS hdf5
        J_Vector=Vector(MPI.comm_self,2)
        J_Vector.set_local(np.array([J_real_unscaled,J_im_unscaled],dtype=np.float64))
        Hdf=HDF5File(mesh_sol.mpi_comm(), os.environ['PATIENTDIR']+"/Results_adaptive/Solution_"+str(np.round(Field_calc_param.frequenc,6))+".h5", "w")
        Hdf.write(mesh_sol, "mesh_sol")
        Hdf.write(phi_sol, "solution_phi_full")
        Hdf.write(E_field, "solution_E_field")
        Hdf.write(E_field_im, "solution_E_field_im")
        Hdf.write(j_dens_real_unscaled, "solution_j_real")
        Hdf.write(j_dens_im_unscaled, "solution_j_im")
        Hdf.write(J_Vector, "/J_Vector")
        Hdf.close()
        return True

### Do not probe anything when inside MPI process!
#    #Probe_of_potential
#    probe_z=np.zeros((100,4),float)
#    for inx in range(100):
#        pnt=Point(75.5,78.5,27.865+inx/10.0)
#        probe_z[inx,0]=75.5
#        probe_z[inx,1]=78.5
#        probe_z[inx,2]=27.865+inx/10.0
#        if Field_calc_param.c_c==1:
#            phi_r_sol_scaled_on_point=V_across*np.real((phi_r_sol(pnt)+1j*phi_i_sol(pnt))/(J_real+1j*J_im))
#            phi_i_sol_scaled_on_point=V_across*np.imag((phi_r_sol(pnt)+1j*phi_i_sol(pnt))/(J_real+1j*J_im))
#            probe_z[inx,3]=np.sqrt(phi_r_sol_scaled_on_point*phi_r_sol_scaled_on_point+phi_i_sol_scaled_on_point*phi_i_sol_scaled_on_point)
#        else:
#            probe_z[inx,3]=np.sqrt(phi_r_sol(pnt)*phi_r_sol(pnt)+phi_i_sol(pnt)*phi_i_sol(pnt))
#    np.savetxt('Results_adaptive/Phi_Zprobe'+str(Field_calc_param.frequenc)+'.csv',  probe_z, delimiter=" ")

    #print("Tissue impedance: ", Z_tis)

#=============================================================================#
    if Field_calc_param.c_c==1 or Field_calc_param.CPE==1:

        Z_tissue = V_across/J_ground      # Tissue impedance
        if MPI.comm_world.rank==1:
            print("Tissue impedance: ", Z_tissue)

        if Field_calc_param.CPE==1:

            if len(Domains.Amp_vector)>2:
                print("Currently, CPE can be used only for simulations with two contacts. Please, assign the rest to 'None'")
                raise SystemExit

            from GUI_inp_dict import d as d_cpe
            CPE_param=[d_cpe["K_A"],d_cpe["beta"],d_cpe["K_A_ground"],d_cpe["beta_ground"]]

            from FEM_in_spectrum import get_CPE_corrected_Dirichlet_BC           #-1.0 to avoid printing
            Dirichlet_bc_with_CPE,total_impedance=get_CPE_corrected_Dirichlet_BC(Field_calc_param.external_grounding,facets,boundaries_sol,CPE_param,Field_calc_param.EQS_mode,Field_calc_param.frequenc,-1.0,Domains.Contacts,Domains.Amp_vector,V_across,Z_tissue,V_space)
            if MPI.comm_world.rank==1:
                print("Solving for an adjusted potential on contacts to account for CPE")
                start_math=tm.time()
            # to solve the Laplace equation for the adjusted Dirichlet
            phi_sol_CPE=define_variational_form_and_solve(V_space,Dirichlet_bc_with_CPE,kappa,Field_calc_param.EQS_mode,Cond_tensor,Solver_type)
            if MPI.comm_world.rank==1:
                minutes=int((tm.time() - start_math)/60)
                secnds=int(tm.time() - start_math)-minutes*60
                print("--- assembled and solved in ",minutes," min ",secnds," s ")
            if Field_calc_param.EQS_mode=='EQS':
                (phi_r_CPE,phi_i_CPE)=phi_sol_CPE.split(deepcopy=True)
            else:
                phi_r_CPE=phi_sol_CPE
                phi_i_CPE=Function(V_space)
                phi_i_CPE.vector()[:] = 0.0

            # get current flowing through the grounded contact and the electric field in the whole domain
            J_ground_CPE,E_field_CPE,E_field_im_CPE = get_current(mesh_sol,facets,boundaries_sol,Field_calc_param.element_order,Field_calc_param.EQS_mode,Domains.Contacts,kappa,Cond_tensor,phi_r_CPE,phi_i_CPE,ground_index,get_E_field=True)
            # If EQS, J_ground is a complex number. If QS, E_field_CPE is a null function

            # to get current density function which is required for mesh refinement when checking current convergence
            j_dens_real_CPE,j_dens_im_CPE = get_current_density(mesh_sol,Field_calc_param.element_order,Field_calc_param.EQS_mode,kappa,Cond_tensor,E_field_CPE,E_field_im_CPE)
            # If QS, j_dens_im is null function

            # will be used for mesh refinement
            j_dens_real_unscaled=j_dens_real_CPE.copy(deepcopy=True)
            j_dens_im_unscaled=j_dens_im_CPE.copy(deepcopy=True)
            J_real_unscaled=copy.deepcopy(np.real(J_ground))
            J_im_unscaled=copy.deepcopy(np.imag(J_ground))

            #Just to save total currunt through the ground in FEniCS hdf5
            J_Vector=Vector(MPI.comm_self,2)
            J_Vector.set_local(np.array([J_real_unscaled,J_im_unscaled],dtype=np.float64))
            Hdf=HDF5File(mesh_sol.mpi_comm(), os.environ['PATIENTDIR']+"/Results_adaptive/Solution_"+str(np.round(Field_calc_param.frequenc,6))+".h5", "w")
            Hdf.write(mesh_sol, "mesh_sol")
            Hdf.write(phi_sol_CPE, "solution_phi_full")
            Hdf.write(E_field_CPE, "solution_E_field")
            Hdf.write(E_field_im_CPE, "solution_E_field_im")
            Hdf.write(j_dens_real_unscaled, "solution_j_real")
            Hdf.write(j_dens_im_unscaled, "solution_j_im")
            Hdf.write(J_Vector, "/J_Vector")
            Hdf.close()

            return True

        if Field_calc_param.c_c==1:
            if Field_calc_param.EQS_mode=='EQS':     # For EQS, we need to scale the potential on boundaries (because the error is absolute) and recompute field, etc. Maybe we can scale them also directly?
                Dirichlet_bc_scaled=[]
                for bc_i in range(len(Domains.Contacts)):          #CPE estimation is valid only for one activa and one ground contact configuration
                    if Field_calc_param.EQS_mode=='EQS':
                        if Domains.Amp_vector[bc_i]!=0.0:
                            Active_with_CC=V_across*V_across/J_ground          #(impedance * current through the contact (V_across coincides with the assigned current magnitude))
                            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(0), np.real(Active_with_CC), boundaries_sol,Domains.Contacts[bc_i]))
                            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(1), np.imag(Active_with_CC), boundaries_sol,Domains.Contacts[bc_i]))
                        else:
                            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(0), Constant(0.0), boundaries_sol,Domains.Contacts[bc_i]))
                            Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(1), Constant(0.0), boundaries_sol,Domains.Contacts[bc_i]))

                if Field_calc_param.external_grounding==True:
                    if Sim_setup.Laplace_eq == 'EQS':
                        Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(0),0.0,facets,1))
                        Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(1),0.0,facets,1))
                    else:
                        Dirichlet_bc_scaled.append(DirichletBC(V_space,0.0,facets,1))

                    print("Solving for a scaled potential on contacts (to match the desired current)")
                    start_math=tm.time()
                    # to solve the Laplace equation for the adjusted Dirichlet
                    phi_sol_scaled=define_variational_form_and_solve(V_space,Dirichlet_bc_scaled,kappa,Field_calc_param.EQS_mode,Cond_tensor,Solver_type)
                    minutes=int((tm.time() - start_math)/60)
                    secnds=int(tm.time() - start_math)-minutes*60
                    print("--- assembled and solved in ",minutes," min ",secnds," s ---")

                    (phi_r_sol_scaled,phi_i_sol_scaled)=phi_sol_scaled.split(deepcopy=True)

                    # get current flowing through the grounded contact and the electric field in the whole domain
                    J_ground_scaled,E_field_scaled,E_field_im_scaled = get_current(mesh_sol,facets,boundaries_sol,Field_calc_param.element_order,Field_calc_param.EQS_mode,Domains.Contacts,kappa,Cond_tensor,phi_r_CPE,phi_i_CPE,ground_index,get_E_field=True)
                    # If EQS, J_ground is a complex number. If QS, E_field_im is 0

            else:   # here we can simply scale the potential in the domain and recompute the E-field
                phi_r_sol_scaled=Function(V_space)
                phi_i_sol_scaled=Function(V_space)
                phi_i_sol_scaled.vector()[:] = 0.0
                phi_r_sol_scaled.vector()[:]=V_across*phi_r_sol.vector()[:]/J_ground

                phi_sol_scaled=phi_r_sol_scaled

                J_ground_scaled,E_field_scaled,E_field_im_scaled = get_current(mesh_sol,facets,boundaries_sol,Field_calc_param.element_order,Field_calc_param.EQS_mode,Domains.Contacts,kappa,Cond_tensor,phi_r_sol_scaled,phi_i_sol_scaled,ground_index,get_E_field=True)
                #E_field_im_scale is a null function

            #Just to save total currunt through the ground in FEniCS hdf5 (we save unscaled currents!)
            J_Vector=Vector(MPI.comm_self,2)
            J_Vector.set_local(np.array([J_real_unscaled,J_im_unscaled],dtype=np.float64))
            Hdf=HDF5File(mesh_sol.mpi_comm(), os.environ['PATIENTDIR']+"/Results_adaptive/Solution_"+str(np.round(Field_calc_param.frequenc,6))+".h5", "w")
            Hdf.write(mesh_sol, "mesh_sol")
            Hdf.write(phi_sol_scaled, "solution_phi_full")
            Hdf.write(E_field_scaled, "solution_E_field")
            Hdf.write(E_field_im_scaled, "solution_E_field_im")
            Hdf.write(j_dens_real_unscaled, "solution_j_real")
            Hdf.write(j_dens_im_unscaled, "solution_j_im")
            Hdf.write(J_Vector, "/J_Vector")
            Hdf.close()
            return True

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

    if Field_calc_param.EQS_mode=='EQS':
        hdf.read(kappa_i, "/kappa_i")

    #anisotropy will be read in at the site

    hdf.close()

    if Field_calc_param.EQS_mode=='EQS':
        compute_field(mesh,Domains,subdomains,boundaries,kappa_r,Field_calc_param,kappa_i)
    elif Field_calc_param.EQS_mode=='QS':
        compute_field(mesh,Domains,subdomains,boundaries,kappa_r,Field_calc_param)
