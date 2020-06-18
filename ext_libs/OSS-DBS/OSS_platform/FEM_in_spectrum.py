# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 10:39:20 2017

@author: konstantin
"""
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=FutureWarning)
    import h5py
from dolfin import *
import numpy as np
import time
import pickle

from tissue_dielectrics import DielectricProperties

parameters['linear_algebra_backend']='PETSc'
set_log_active(False)   #turns off debugging info


def get_dielectric_properties_from_subdomains(mesh,subdomains,Laplace_formulation,float_conductors,conductivities,rel_permittivities,frequenc):
    
    cond_default,cond_GM,cond_WM,cond_CSF,cond_encap=conductivities[:]
    V0_r=FunctionSpace(mesh,'DG',0)   
    kappa_r=Function(V0_r)
    if float_conductors==-1:     #(-1 means no floating contacts)       #[Default,CSF,WM,GM,Encap,Floating]
        k_val_r=[cond_default*0.001,cond_CSF*0.001,cond_WM*0.001,cond_GM*0.001,cond_encap*0.001]
    else:
        k_val_r=[cond_default*0.001,cond_CSF*0.001,cond_WM*0.001,cond_GM*0.001,cond_encap*0.001,1000.0]
    
    help = np.asarray(subdomains.array(), dtype=np.int32)
    kappa_r.vector()[:] = np.choose(help, k_val_r)
    kappa=[kappa_r]
    
    if Laplace_formulation=='EQS':
        perm_default,perm_GM,perm_WM,perm_CSF,perm_encap=rel_permittivities
        V0_i=FunctionSpace(mesh,'DG',0)
        kappa_i=Function(V0_i)        
        omega_eps0=2*np.pi*frequenc*8.854e-12           #2*pi*f*eps0
        if float_conductors==-1:    #[Default,CSF,WM,GM,Encap,Floating]
            k_val_i=[omega_eps0*perm_default*0.001,omega_eps0*perm_CSF*0.001,omega_eps0*perm_WM*0.001,omega_eps0*perm_GM*0.001,1*omega_eps0*perm_encap*0.001]
        else:
            k_val_i=[omega_eps0*perm_default*0.001,omega_eps0*perm_CSF*0.001,omega_eps0*perm_WM*0.001,omega_eps0*perm_GM*0.001,1*omega_eps0*perm_encap*0.001,1000000000*omega_eps0]

        help = np.asarray(subdomains.array(), dtype=np.int32)
        kappa_i.vector()[:] = np.choose(help, k_val_i)        #because model is in mm
        kappa = [kappa_r, kappa_i]
        
    return kappa,k_val_r


def get_solution_space_and_Dirichlet_BC(mesh,boundaries,element_order,Laplace_eq,Contacts_indices,Phi_vector,only_space=False):
    # Function space where the electric potential will be computed
    if Laplace_eq=='EQS':     #complex numbers are not yet implemented in FEniCS, mixed space is used
        El_r = FiniteElement("Lagrange", mesh.ufl_cell(),element_order)
        El_i = FiniteElement("Lagrange", mesh.ufl_cell(),element_order)
        El_complex = El_r * El_i
        V = FunctionSpace(mesh, El_complex)
    else:
        V = FunctionSpace(mesh, "Lagrange",element_order)
        
    if only_space==True:        #for some setups we have custom assignment of boundaries
        return V

    # Dirichlet boundary condition (electric potenial on the contacts. In case of current-controlled stimulation, it will be scaled afterwards (due to the system linearity))
    bc=[]    
    for bc_i in range(len(Contacts_indices)):
        if Laplace_eq=='EQS':
            bc.append(DirichletBC(V.sub(0), Phi_vector[bc_i], boundaries,Contacts_indices[bc_i]))
            bc.append(DirichletBC(V.sub(1), Constant(0.0), boundaries,Contacts_indices[bc_i]))         # the imaginary part is set to 0 for the initial computation
        else:
            bc.append(DirichletBC(V, Phi_vector[bc_i], boundaries,Contacts_indices[bc_i]))
            
        if Phi_vector[bc_i]==0.0:       #we must have ground in every simulation
            ground_index=bc_i
    
    if not('ground_index' in locals()):
        print("No Dirichlet BC for grounding was found")
        raise SystemExit
        
    return V,bc,ground_index

def get_scaled_cond_tensor(mesh,subdomains,sine_freq,signal_freq,unscaled_tensor,cond_list,plot_tensors=False):

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
    
    #the Tensor order is xx,xy,xz,yy,yz,zz    
    # initializing here just to ensure stability
    c00 = MeshFunction("double", mesh, 3, 1.0)
    c01 = MeshFunction("double", mesh, 3, 0.0)
    c02 = MeshFunction("double", mesh, 3, 0.0)
    c11 = MeshFunction("double", mesh, 3, 1.0)
    c12 = MeshFunction("double", mesh, 3, 0.0)
    c22 = MeshFunction("double", mesh, 3, 1.0)

    if int(sine_freq)==int(signal_freq):  # will be used for visualization
        c_unscaled00 = unscaled_tensor[0]
        c_unscaled01 = unscaled_tensor[1]
        c_unscaled02 = unscaled_tensor[2]
        c_unscaled11 = unscaled_tensor[3]
        c_unscaled12 = unscaled_tensor[4]
        c_unscaled22 = unscaled_tensor[5]
    
    cell_Anis = MeshFunction('bool',mesh,3)
    cell_Anis.set_all(False)
  
    for cell in cells(mesh):
        scale_cond=cond_list[subdomains[cell]]      #check whick cond. value (S/mm) this cell was assigned
        
        if unscaled_tensor[0][cell]!=1.0 or unscaled_tensor[3][cell]!=1.0 or unscaled_tensor[5][cell]!=1.0:
            cell_Anis[cell]=True
    
        c00[cell]=unscaled_tensor[0][cell]*scale_cond
        c01[cell]=unscaled_tensor[1][cell]*scale_cond
        c02[cell]=unscaled_tensor[2][cell]*scale_cond
        c11[cell]=unscaled_tensor[3][cell]*scale_cond
        c12[cell]=unscaled_tensor[4][cell]*scale_cond
        c22[cell]=unscaled_tensor[5][cell]*scale_cond
        
        
    if plot_tensors==True:
        file=File('Tensors/c00_mapped.pvd')
        file<<c00,mesh
        file=File('Tensors/c11_mapped.pvd')
        file<<c11,mesh    
        file=File('Tensors/c22_mapped.pvd')
        file<<c22,mesh
        file=File('Tensors/c01_mapped.pvd')
        file<<c01,mesh
    
        file=File('Tensors/Anis_cells.pvd')
        file<<cell_Anis,mesh
    

    c = CompiledExpression(compile_cpp_code(conductivity_code).Conductivity(),
           c00=c00, c01=c01, c02=c02, c11=c11, c12=c12, c22=c22, degree=0)
    C_tensor = as_matrix(((c[0], c[1], c[2]), (c[1], c[3], c[4]),(c[2],c[4],c[5])))
    
    if int(sine_freq)==int(signal_freq):        
        c_unscaled = CompiledExpression(compile_cpp_code(conductivity_code).Conductivity(),
                       c00=c_unscaled00, c01=c_unscaled01, c02=c_unscaled02, c11=c_unscaled11, c12=c_unscaled12, c22=c_unscaled22, degree=0)
        tensor= as_tensor([[c_unscaled[0], c_unscaled[1], c_unscaled[2]], [c_unscaled[1], c_unscaled[3], c_unscaled[4]],[c_unscaled[2],c_unscaled[4],c_unscaled[5]]])
        f_vector_repr=project(tensor,TensorFunctionSpace(mesh, "Lagrange", 1),solver_type="cg", preconditioner_type="amg")
        file=File('Tensors/Ellipsoids_unscaled_at_'+str(signal_freq)+'_Hz.pvd')
        file<<f_vector_repr    
        
    return C_tensor

def define_variational_form_and_solve(V,dirichlet_bc,kappa,Laplace_eq,Cond_tensor,Solver_type):  # to solve the Laplace equation div(kappa*grad(phi))=0   (variational form: a(u,v)=L(v))    

    # to define the variational problem
    f = Constant(0.0)   #to keep the right side of Laplace equation 0    
    if Laplace_eq=='EQS':    
        u_r, u_i = TrialFunction(V)
        v_r, v_i = TestFunction(V)
        if Cond_tensor!=False:
            a = (inner(Cond_tensor*grad(u_r), grad(v_r))*dx
             -inner(kappa[1]*grad(u_i), grad(v_r))*dx
             -inner(kappa[1]*grad(u_r), grad(v_i))*dx
             -inner(Cond_tensor*grad(u_i), grad(v_i))*dx
             +inner(Cond_tensor*grad(u_r), grad(v_i))*dx
             -inner(kappa[1]*grad(u_i), grad(v_i))*dx
             +inner(kappa[1]*grad(u_r), grad(v_r))*dx
             +inner(Cond_tensor*grad(u_i), grad(v_r))*dx
             )
            
        else:
            a = (inner(kappa[0]*grad(u_r), grad(v_r))*dx
             -inner(kappa[1]*grad(u_i), grad(v_r))*dx
             -inner(kappa[1]*grad(u_r), grad(v_i))*dx
             -inner(kappa[0]*grad(u_i), grad(v_i))*dx
             +inner(kappa[0]*grad(u_r), grad(v_i))*dx
             -inner(kappa[1]*grad(u_i), grad(v_i))*dx
             +inner(kappa[1]*grad(u_r), grad(v_r))*dx
             +inner(kappa[0]*grad(u_i), grad(v_r))*dx
             )
            
        L = -(f*v_r+f*v_i)*dx 
    else:
        u = TrialFunction(V)
        v = TestFunction(V)
        if Cond_tensor!=False:
            a = inner(Cond_tensor*grad(u), grad(v))*dx        
        else:
            a = inner(kappa[0]*grad(u), grad(v))*dx
        
        L = f*v*dx

    u=Function(V)
    
    # to compute the solution
    if Solver_type=='MUMPS':        # slow but realiable, especially suitable with multiple floating conductors
        # solving the problem
        # solve iterative
        problem = LinearVariationalProblem(a, L, u, dirichlet_bc)
        solver = LinearVariationalSolver(problem)
        #solver.parameters.linear_solver = 'cg'
        solver.parameters["linear_solver"] = "mumps"
        solver.parameters['preconditioner'] = 'ilu'
        cg_prm = solver.parameters['krylov_solver']
        cg_prm['absolute_tolerance'] = 1E-7
        cg_prm['relative_tolerance'] = 1E-6
        solver.solve()
    elif Solver_type=='BiCGSTAB':       # efficient iterative solver suitable for problems with no more than one floating conductor
        A, b = assemble_system(a, L, dirichlet_bc)        
        solver = PETScKrylovSolver('bicgstab','petsc_amg')
        #solver.parameters['monitor_convergence'] = True
        #solver.parameters['report'] = True
        #print dir(self.solver.parameters)
        solver.set_operator(A)
        solver.solve(u.vector(), b)
    elif Solver_type=='GMRES':       # iterative solver suitable for QS, also suitable for QS problems with multiple floating conductors (check!)
        A, b = assemble_system(a, L, dirichlet_bc)
        solver = KrylovSolver(A,'gmres','hypre_amg')
        #solver.parameters['monitor_convergence'] = True
        #solver.parameters['report'] = True
        solver.solve(u.vector(), b)
    else:
        print("Solver was not found")
        raise SystemExit
    
    return u 

def get_current(mesh,boundaries,element_order,Laplace_eq,Contacts_indices,kappa,C_tensor,phi_real,phi_imag,ground_index,get_E_field=False):        
        
        if element_order>1:
            W = VectorFunctionSpace(mesh,'DG',element_order-1)
            W_i = VectorFunctionSpace(mesh,'DG',element_order-1)
        else:
            W = VectorFunctionSpace(mesh,'DG',element_order)
            W_i = VectorFunctionSpace(mesh,'DG',element_order)  
        
        facets=MeshFunction('size_t',mesh, 2)
        facets.set_all(0)
        facets.array()[boundaries.array()==Contacts_indices[ground_index]]=1
        ds=Measure("ds",domain=mesh,subdomain_data=facets)
        
        #Explicit E-field projection
        w = TestFunction(W)
        Pv = TrialFunction(W)
        E_field = Function(W)
        a_local = inner(w, Pv) * dx
        L_local = inner(w, -grad(phi_real)) * dx            
        A_local, b_local = assemble_system(a_local, L_local, bcs=[])            
        local_solver = PETScKrylovSolver('bicgstab')
        local_solver.solve(A_local,E_field.vector(),b_local)  
        
        n = FacetNormal(mesh)        
        if Laplace_eq == 'EQS':
            w_i = TestFunction(W_i)
            Pv_i = TrialFunction(W_i)
            E_field_im = Function(W_i)
            a_local = inner(w_i, Pv_i) * dx
            L_local = inner(w_i, -grad(phi_imag)) * dx                
            A_local, b_local = assemble_system(a_local, L_local, bcs=[])                
            local_solver = PETScKrylovSolver('bicgstab')
            local_solver.solve(A_local,E_field_im.vector(),b_local)   

            if C_tensor!=False:
                j_dens_real_ground = dot(C_tensor*E_field,-1*n)*ds(1)-dot(kappa[1]*E_field_im,-1*n)*ds(1)
                j_dens_im_ground= dot(C_tensor*E_field_im,-1*n)*ds(1)+dot(kappa[1]*E_field,-1*n)*ds(1)
            else:
                j_dens_real_ground = dot(kappa[0]*E_field,-1*n)*ds(1)-dot(kappa[1]*E_field_im,-1*n)*ds(1)
                j_dens_im_ground= dot(kappa[0]*E_field_im,-1*n)*ds(1)+dot(kappa[1]*E_field,-1*n)*ds(1)
                
            #we always assess current on the ground in 2 contact-case, so the sign should be flipped    
            J_real=-1*assemble(j_dens_real_ground)
            J_im=-1*assemble(j_dens_im_ground)
            J_complex_ground=J_real+1j*J_im 
            
            if get_E_field==True:
                return J_complex_ground,E_field,E_field_im
            else:                                      
                return J_complex_ground
          
        else:
            E_field_im = Function(W)
            E_field_im.vector()[:] = 0.0        #fake        
            if C_tensor!=False:
                j_dens_real_ground = dot(C_tensor*E_field,-1*n)*ds(1)
            else:
                j_dens_real_ground = dot(kappa[0]*E_field,-1*n)*ds(1)        
            #we always assess current on the ground in 2 contact-case, so the sign should be flipped
            J_real_ground=-1*assemble(j_dens_real_ground)

            if get_E_field==True:
                return J_real_ground,E_field,E_field_im
            else:                                      
                return J_real_ground

def get_CPE_corrected_Dirichlet_BC(boundaries,CPE_param,Laplace_eq,sine_freq,freq_signal,Contacts_indices,Phi_vector,Voltage_drop,Z_tissue,V_space):
            
    K_A,beta,K_A_ground,beta_ground = CPE_param

    if sine_freq==0.0:       #how to define this?
        Z_CPE=0.0
        Z_CPE_ground=0.0
    else:
        Z_CPE=K_A/((1j*2*np.pi*sine_freq)**(beta))
        Z_CPE_ground=K_A_ground/((1j*2*np.pi*sine_freq)**(beta_ground))

    if sine_freq==freq_signal:       
        print("Z_CPE_ground at "+str(sine_freq)+" Hz : ", Z_CPE_ground)
        print("Z_CPE at "+str(sine_freq)+" Hz : ", Z_CPE)

    bc_cpe=[]
    for bc_i in range(len(Contacts_indices)):          #CPE estimation is valid only for one activa and one ground contact configuration
        if Phi_vector[bc_i] == 0.0:
            if Laplace_eq=='EQS':
                Ground_with_CPE=0.0+(Voltage_drop/(Z_tissue+Z_CPE+Z_CPE_ground))*Z_CPE_ground
                bc_cpe.append(DirichletBC(V_space.sub(0), np.real(Ground_with_CPE), boundaries,Contacts_indices[bc_i]))
                bc_cpe.append(DirichletBC(V_space.sub(1), np.imag(Ground_with_CPE), boundaries,Contacts_indices[bc_i]))
            else:
                Ground_with_CPE=0.0+(Voltage_drop/(Z_tissue+Z_CPE+Z_CPE_ground))*Z_CPE_ground
                bc_cpe.append(DirichletBC(V_space, np.real(Ground_with_CPE), boundaries,Contacts_indices[bc_i]))
        else:
            if Laplace_eq=='EQS':
                Active_with_CPE=Phi_vector[bc_i]-(Voltage_drop/(Z_tissue+Z_CPE+Z_CPE_ground))*Z_CPE
                bc_cpe.append(DirichletBC(V_space.sub(0), np.real(Active_with_CPE), boundaries,Contacts_indices[bc_i]))
                bc_cpe.append(DirichletBC(V_space.sub(1), np.imag(Active_with_CPE), boundaries,Contacts_indices[bc_i]))
            else:
                Active_with_CPE=Phi_vector[bc_i]-(Voltage_drop/(Z_tissue+Z_CPE+Z_CPE_ground))*Z_CPE
                bc_cpe.append(DirichletBC(V_space, np.real(Active_with_CPE), boundaries,Contacts_indices[bc_i]))
                
                
    if sine_freq==freq_signal:
        print("Ground adjusted by CPE at "+str(sine_freq)+" Hz : ", Ground_with_CPE)
        print("Active contact adjusted by CPE at "+str(sine_freq)+" Hz : ", Active_with_CPE)

    if Laplace_eq=='EQS':
        comb_Z=np.vstack((np.real(Z_tissue+Z_CPE+Z_CPE_ground),np.imag(Z_tissue+Z_CPE+Z_CPE_ground),sine_freq)).T
    else:
        comb_Z=np.vstack((np.real(Z_tissue+Z_CPE+Z_CPE_ground),0.0,sine_freq)).T        #no imag. part for QS
        
    return bc_cpe,comb_Z


def solve_Laplace(Sim_setup,Solver_type,Vertices_array,Domains,core,Full_IFFT,output):
    set_log_active(False)   #turns off debugging info
    
    Sim_setup.mesh.coordinates()
    Sim_setup.mesh.init()
    
    # to get conductivity (and permittivity if EQS formulation) mapped accrodingly to the subdomains. k_val_r is just a list of conductivities (S/mm!) in a specific order to scale the cond. tensor
    kappa,k_val_r=get_dielectric_properties_from_subdomains(Sim_setup.mesh,Sim_setup.subdomains,Sim_setup.Laplace_eq,Domains.Float_contacts,Sim_setup.conductivities,Sim_setup.rel_permittivities,Sim_setup.sine_freq)
    if int(Sim_setup.sine_freq)==int(Sim_setup.signal_freq):
        file=File('Field_solutions/Conductivity_map_'+str(Sim_setup.signal_freq)+'Hz.pvd')
        file<<kappa[0]
        if Sim_setup.Laplace_eq == 'EQS':
            file=File('Field_solutions/Permittivity_map_'+str(Sim_setup.signal_freq)+'Hz.pvd')
            file<<kappa[1]

    # to get tensor scaled by the conductivity map
    if Sim_setup.anisotropy==1:
        Cond_tensor=get_scaled_cond_tensor(Sim_setup.mesh,Sim_setup.subdomains,Sim_setup.sine_freq,Sim_setup.signal_freq,Sim_setup.unscaled_tensor,k_val_r)
    else:
        Cond_tensor=False  #just to initialize

    #In case of current-controlled stimulation, Dirichlet_bc or the whole potential distribution will be scaled afterwards (due to the system's linearity)
    V_space,Dirichlet_bc,ground_index=get_solution_space_and_Dirichlet_BC(Sim_setup.mesh,Sim_setup.boundaries,Sim_setup.element_order,Sim_setup.Laplace_eq,Domains.Contacts,Domains.fi)
    #ground index refers to the ground in .med/.msh file
    
    # to solve the Laplace equation div(kappa*grad(phi))=0   (variational form: a(u,v)=L(v))    
    phi_sol=define_variational_form_and_solve(V_space,Dirichlet_bc,kappa,Sim_setup.Laplace_eq,Cond_tensor,Solver_type)        
    
    if Sim_setup.Laplace_eq=='EQS':
        (phi_r,phi_i)=phi_sol.split(deepcopy=True)
    else:
        phi_r=phi_sol
        phi_i=Function(V_space)
        phi_i.vector()[:] = 0.0

    #save unscaled real solution for plotting
    if int(Sim_setup.sine_freq)==int(Sim_setup.signal_freq):
        file=File('Field_solutions/Phi_real_unscaled_'+str(Sim_setup.signal_freq)+'Hz.pvd')
        file<<phi_r,Sim_setup.mesh
        print("DoFs on the mesh for "+Sim_setup.Laplace_eq+" : ", (max(V_space.dofmap().dofs())+1))
    

    if Sim_setup.c_c==1 or Sim_setup.CPE_status==1:     #we compute E-field, currents and impedances only for current-controlled or if CPE is used
        J_ground=get_current(Sim_setup.mesh,Sim_setup.boundaries,Sim_setup.element_order,Sim_setup.Laplace_eq,Domains.Contacts,kappa,Cond_tensor,phi_r,phi_i,ground_index)                
        #If EQS, J_ground is a complex number
                
        V_across=max(Domains.fi[:], key=abs)        # voltage drop in the system
        Z_tissue = V_across/J_ground                   # Tissue impedance
        if int(Sim_setup.sine_freq)==int(Sim_setup.signal_freq):
            print("Tissue impedance at the signal freq.: ",Z_tissue)
    
        if Sim_setup.CPE_status==1:   # in this case we need to estimate the voltage drop over the CPE and adjust the Derichlet BC accordingly 
            Dirichlet_bc_with_CPE,total_impedance=get_CPE_corrected_Dirichlet_BC(Sim_setup.boundaries,Sim_setup.CPE_param,Sim_setup.Laplace_eq,Sim_setup.sine_freq,Sim_setup.signal_freq,Domains.Contacts,Domains.fi,V_across,Z_tissue,V_space)

            f=open('Field_solutions/Impedance'+str(core)+'.csv','ab')
            np.savetxt(f, total_impedance, delimiter=" ")
            f.close()

            # to solve the Laplace equation for the adjusted Dirichlet   
            phi_sol_CPE=define_variational_form_and_solve(V_space,Dirichlet_bc_with_CPE,kappa,Sim_setup.Laplace_eq,Cond_tensor,Solver_type)   
            if Sim_setup.Laplace_eq=='EQS':
                (phi_r_CPE,phi_i_CPE)=phi_sol_CPE.split(deepcopy=True)
            else:
                phi_r_CPE=phi_sol_CPE
                phi_i_CPE=Function(V_space)
                phi_i_CPE.vector()[:] = 0.0

            J_ground_CPE=get_current(Sim_setup.mesh,Sim_setup.boundaries,Sim_setup.element_order,Sim_setup.Laplace_eq,Domains.Contacts,kappa,Cond_tensor,phi_r_CPE,phi_i_CPE,ground_index)                

            # just resaving
            phi_sol,phi_r,phi_i,J_ground=(phi_sol_CPE,phi_r_CPE,phi_i_CPE,J_ground_CPE)

    if Full_IFFT==1:
        Hdf=HDF5File(Sim_setup.mesh.mpi_comm(), "Field_solutions_functions/solution"+str(np.round(Sim_setup.sine_freq,6))+".h5", "w")
        Hdf.write(Sim_setup.mesh, "mesh")
        if Sim_setup.CPE_status!=1:
            Hdf.write(phi_sol, "solution_full")
            
        if Sim_setup.c_c==1:
            with open('Field_solutions_functions/current_scale'+str(np.round(Sim_setup.sine_freq,6))+'.file', 'wb') as f:
                pickle.dump(np.array([np.real(J_ground),np.imag(J_ground)]), f)
        Hdf.close()
        
    else:
        Phi_ROI=np.zeros((Vertices_array.shape[0],5),float)

        for inx in range(Vertices_array.shape[0]):
            pnt=Point(Vertices_array[inx,0],Vertices_array[inx,1],Vertices_array[inx,2])
            if Sim_setup.mesh.bounding_box_tree().compute_first_entity_collision(pnt)<Sim_setup.mesh.num_cells()*100:
    
                Phi_ROI[inx,0]=Vertices_array[inx,0]
                Phi_ROI[inx,1]=Vertices_array[inx,1]
                Phi_ROI[inx,2]=Vertices_array[inx,2]

                if Sim_setup.c_c==1:
                    Phi_ROI[inx,3]=np.real((phi_r(pnt)+1j*phi_i(pnt))/J_ground)    #*1A is left out here
                    Phi_ROI[inx,4]=np.imag((phi_r(pnt)+1j*phi_i(pnt))/J_ground)    #*1A is left out here
                else:
                    Phi_ROI[inx,3]=phi_r(pnt)
                    Phi_ROI[inx,4]=phi_i(pnt)
            else:
                print("Couldn't probe the potential at the point ",Vertices_array[inx,0],Vertices_array[inx,1],Vertices_array[inx,2])
                print("check the neuron array, exiting....")
                raise SystemExit
                
        fre_vector=[Sim_setup.sine_freq]*Phi_ROI.shape[0] 
        comb=np.vstack((Phi_ROI[:,0],Phi_ROI[:,1],Phi_ROI[:,2],Phi_ROI[:,3],Phi_ROI[:,4],fre_vector)).T    
               
        f = h5py.File('Field_solutions/sol_cor'+str(core)+'.h5','a')
        f.create_dataset(str(Sim_setup.sine_freq), data=comb)
        f.close()
                    
    if Sim_setup.c_c==1:
        comb_Z=np.vstack((np.real(Z_tissue),np.imag(Z_tissue),Sim_setup.sine_freq)).T
        f=open('Field_solutions/Impedance'+str(core)+'.csv','ab')
        np.savetxt(f, comb_Z, delimiter=" ")
        f.close()
        
        
    if Sim_setup.sine_freq==Sim_setup.signal_freq and Sim_setup.c_c==1:  # re-solve with the scaled potential just to check (to match 1 A)        
        Dirichlet_bc_scaled=[]       
        for bc_i in range(len(Domains.Contacts)):          #CPE estimation is valid only for one activa and one ground contact configuration
            if Sim_setup.Laplace_eq == 'EQS':
                if Domains.fi[bc_i]==0.0:
                    Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(0), 0.0, Sim_setup.boundaries,Domains.Contacts[bc_i]))
                    Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(1), 0.0, Sim_setup.boundaries,Domains.Contacts[bc_i]))
                else:
                    Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(0), np.real((Domains.fi[bc_i])/J_ground),Sim_setup.boundaries,Domains.Contacts[bc_i]))
                    Dirichlet_bc_scaled.append(DirichletBC(V_space.sub(1), np.imag((Domains.fi[bc_i])/J_ground),Sim_setup.boundaries,Domains.Contacts[bc_i]))
            else:
                if Domains.fi[bc_i]==0.0:
                    Dirichlet_bc_scaled.append(DirichletBC(V_space, 0.0, Sim_setup.boundaries,Domains.Contacts[bc_i]))               
                else:
                    Dirichlet_bc_scaled.append(DirichletBC(V_space, np.real((Domains.fi[bc_i])/J_ground),Sim_setup.boundaries,Domains.Contacts[bc_i]))

        phi_sol_check=define_variational_form_and_solve(V_space,Dirichlet_bc_scaled,kappa,Sim_setup.Laplace_eq,Cond_tensor,Solver_type)        

        if Sim_setup.Laplace_eq=='EQS':
            (phi_r_check,phi_i_check)=phi_sol_check.split(deepcopy=True)
        else:
            phi_r_check=phi_sol_check
            phi_i_check=Function(V_space)
            phi_i_check.vector()[:] = 0.0

        J_ground=get_current(Sim_setup.mesh,Sim_setup.boundaries,Sim_setup.element_order,Sim_setup.Laplace_eq,Domains.Contacts,kappa,Cond_tensor,phi_r_check,phi_i_check,ground_index)                        
        print("Current through the ground after normalizing to 1 A at the signal freq.: ",J_ground)
             
        file=File('Field_solutions/'+str(Sim_setup.Laplace_eq)+str(Sim_setup.signal_freq)+'_phi_r_1A.pvd')
        file<<phi_r_check
        
    output.put(1)

