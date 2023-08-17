"By Konstantin"

import multiprocessing as mp
import numpy as np
import time as tm
import subprocess
import matplotlib.pyplot as plt
import os
from pandas import read_csv

import logging
logging.getLogger('UFL').setLevel(logging.WARNING)
logging.getLogger('FFC').setLevel(logging.WARNING)

from dolfin import *
import h5py

from tissue_dielectrics import DielectricProperties
set_log_active(False)   #turns off debugging info

#Parallel_field_calc is the manager function (called in Launcher)

class Simulation_setup:
    def __init__(self,sine_freq,signal_freq,mesh,boundaries,subdomains,cond_vector,perm_vector,FEM_element_order,anisotropy,current_control_mode,unscaled_tensor,CPE_status,CPE_param,EQS_mode,external_grounding):
        self.mesh=mesh
        self.boundaries=boundaries
        self.subdomains=subdomains
        self.sine_freq=sine_freq
        self.signal_freq=signal_freq
        self.conductivities=cond_vector          #list
        self.rel_permittivities=perm_vector      #list
        self.anisotropy=anisotropy
        self.unscaled_tensor=unscaled_tensor        #list
        self.Laplace_eq=EQS_mode
        self.external_grounding=external_grounding
        self.element_order=FEM_element_order
        self.c_c=current_control_mode
        self.CPE_status=CPE_status
        self.CPE_param=CPE_param         #list

#solution should be re-sorted to the pointwise vectors, i.e. point1_0Hz,point1_130Hz,...,point2_0Hz...
def sort_full_solution(d,FR_vector,full_solution,number_of_points):
    FR_jump=FR_vector.shape[0]
    solution_sort=np.zeros((number_of_points*FR_vector.shape[0],6),float)

    logging.critical("--- Sorting the obtained solution")
    inx_sol_sort=0
    for point_number in np.arange(number_of_points):

        solution_point=[]
        counter_loc=0

        for i in np.arange(point_number,number_of_points*FR_jump,number_of_points):
            solution_point.append(full_solution[i,:])
            counter_loc=counter_loc+1

        solution_point=np.array(solution_point)
        solution_point=solution_point[solution_point[:,5].argsort(axis=0)]      # 5 is the index of the frequency column

        solution_sort[inx_sol_sort:inx_sol_sort+FR_jump,:]=solution_point
        inx_sol_sort=counter_loc+inx_sol_sort

    if d["spectrum_trunc_method"]=='No Truncation':

        hf = h5py.File(os.environ['PATIENTDIR']+'/Field_solutions/sorted_solution.h5', 'w')
        hf.create_dataset('dataset_1', data=solution_sort)
        hf.close()

    else:

        hf = h5py.File(os.environ['PATIENTDIR']+'/Field_solutions/sorted_solution_'+str(d["spectrum_trunc_method"])+'_'+str(d["trunc_param"])+'.h5', 'w')
        hf.create_dataset('dataset_1', data=solution_sort)
        hf.close()

    logging.critical("Saved sorted solution in Field_solutions/")

    del full_solution,solution_sort

    return True



def calculate_in_parallel(d,freq_list,Domains,MRI_param,DTI_param,anisotropy,number_of_points,cc_multicontact):

    start_paral=tm.time()

    Field_on_VTA=0          #temp solution
    if d["VTA_approx"]==1:
        Field_on_VTA=1
        d["VTA_approx"]=0

    #load adapted mesh
    mesh = Mesh(os.environ['PATIENTDIR']+'/Results_adaptive/mesh_adapt.xml.gz')
    boundaries = MeshFunction('size_t',mesh,os.environ['PATIENTDIR']+'/Results_adaptive/boundaries_adapt.xml')
    subdomains_assigned=MeshFunction('size_t',mesh,os.environ['PATIENTDIR']+'/Results_adaptive/subdomains_assigned_adapt.xml')

    #load the neuron array
    Vertices_get=read_csv(os.environ['PATIENTDIR']+'/Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', delimiter=' ', header=None)
    Vertices=Vertices_get.values

    #load CPE parameters if necessary
    CPE_param=[]  # just initialization
    if d["CPE_activ"]==1:
        CPE_param=[d["K_A"],d["beta"],d["K_A_ground"],d["beta_ground"]]



    #load unscaled tensors if DTI data are provided
    if anisotropy==1:

        from Tissue_marking_new import get_cellmap_tensors
        subdomains=get_cellmap_tensors(mesh,subdomains_assigned,Domains,MRI_param,DTI_param,d["default_material"])      #mapping of tissue onto the mesh

        #initiating with isotropic tensor
        c00 = MeshFunction("double", mesh, 3, 1.0)
        c01 = MeshFunction("double", mesh, 3, 0.0)
        c02 = MeshFunction("double", mesh, 3, 0.0)
        c11 = MeshFunction("double", mesh, 3, 1.0)
        c12 = MeshFunction("double", mesh, 3, 0.0)
        c22 = MeshFunction("double", mesh, 3, 1.0)

        hdf = HDF5File(mesh.mpi_comm(), os.environ['PATIENTDIR']+"/Results_adaptive/Tensors_to_solve_num_el_"+str(mesh.num_cells())+".h5", "r")
        hdf.read(c00, "/c00")
        hdf.read(c01, "/c01")
        hdf.read(c02, "/c02")
        hdf.read(c11, "/c11")
        hdf.read(c12, "/c12")
        hdf.read(c22, "/c22")
        hdf.close()

        DTI_tensor=[c00,c01,c02,c11,c12,c22]
    else:
        from Tissue_marking_new import get_cellmap
        subdomains=get_cellmap(mesh,subdomains_assigned,Domains,MRI_param,d["default_material"])    #mapping of tissue onto the mesh
        DTI_tensor=[0,0,0,0,0,0]       #initiating

    logging.critical("Subdomains file for parallel is saved in Field_solutions/parallel_Subdomains.pvd")
    file=File(os.environ['PATIENTDIR']+'/Field_solutions/parallel_Subdomains.pvd')
    file<<subdomains,mesh

    # choose solver
    if d['Solver_Type']=='Default':
        from Math_module_hybrid import choose_solver_for_me
        Solver_type=choose_solver_for_me(d["EQS_core"],Domains.Float_contacts)    #choses solver basing on the Laplace formulation and whether the floating conductors are used
    else:
        Solver_type=d['Solver_Type']      # just get the solver directly


    #with open(os.devnull, 'w') as FNULL: subprocess.call('python Paraview_adapted.py', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    i=0     #frequency index in the frequency list of the signal spectrum
    complete_solution=[]

    # this snippet will check at which frequency the FFEM computations were interrupted
    if d["Parallel_comp_interrupted"] == 1:
        pack_to_start_after=np.genfromtxt(os.environ['PATIENTDIR']+'/Field_solutions/last_completed_pack.csv', delimiter=' ')
        if pack_to_start_after.size == 1:
            rslt=np.where(freq_list == pack_to_start_after)
            i=rslt[0][0]+1
        elif pack_to_start_after[-1]==freq_list[-1]:
            i=freq_list.shape[0]
            logging.critical("All computations in frequency spectrum were already conducted")
        else:
            rslt=np.where(freq_list == pack_to_start_after[-1])
            i=rslt[0][0]+1

    #FFEM calculations are conducted in parallel
    while i<freq_list.shape[0]:
        proc=[]
        j=0  #counter for processes
        output=mp.Queue()
        freq_pack=[]
        while j<d["number_of_processors"] and i<freq_list.shape[0]:

            sine_freq=freq_list[i]
            freq_pack.append(sine_freq)
            i=i+1

            [cond_GM, perm_GM]=DielectricProperties(3).get_dielectrics(sine_freq)        #1 for grey matter and so on
            [cond_WM, perm_WM]=DielectricProperties(2).get_dielectrics(sine_freq)
            [cond_CSF, perm_CSF]=DielectricProperties(1).get_dielectrics(sine_freq)
            [cond_default,perm_default]=DielectricProperties(d["default_material"]).get_dielectrics(sine_freq)

            [cond_encap, perm_encap]=DielectricProperties(d["encap_tissue_type"]).get_dielectrics(sine_freq)
            cond_encap=cond_encap*d["encap_scaling_cond"]
            perm_encap=perm_encap*d["encap_scaling_perm"]

            cond_vector=[cond_default,cond_GM,cond_WM,cond_CSF,cond_encap]
            perm_vector=[perm_default,perm_GM,perm_WM,perm_CSF,perm_encap]

            Sim_setup=Simulation_setup(sine_freq,d["freq"],mesh,boundaries,subdomains,cond_vector,perm_vector,d["el_order"],anisotropy,d["current_control"],DTI_tensor,d["CPE_activ"],CPE_param,d["EQS_core"],d["external_grounding"])


            if cc_multicontact==True:
                import FEM_in_spectrum_multicontact
                processes=mp.Process(target=FEM_in_spectrum_multicontact.solve_Laplace_multicontact,args=(Sim_setup,Solver_type,Vertices,Domains,j,Field_on_VTA,output))
            else:
                import FEM_in_spectrum
                processes=mp.Process(target=FEM_in_spectrum.solve_Laplace,args=(Sim_setup,Solver_type,Vertices,Domains,j,Field_on_VTA,output))

            proc.append(processes)

            j=j+1

        for p in proc:
            p.start()
        for p in proc:
            p.join()

        # check if solutions on all cores were obtained (not a perfect check, works just for the first pack)
        for freq_i in range(j):
            if not os.path.isfile(os.environ['PATIENTDIR']+'/Field_solutions/sol_cor'+str(freq_i)+'.h5'):
                logging.critical('Not all cores returned results, check RAM consumption, exiting')
                raise SystemExit

        last_completed_pack=np.asarray(freq_pack)
        np.savetxt(os.environ['PATIENTDIR']+'/Field_solutions/last_completed_pack.csv', last_completed_pack, delimiter=" ")       #to recover the last frequency of FFEM was interrupted
        if d["freq"] in freq_pack:
            logging.critical("Processed frequencies: ")
        logging.critical("{}".format(' '.join(map(str, freq_pack))))
    if d["number_of_processors"]>freq_list.shape[0]:
        n_files=freq_list.shape[0]
    else:
        n_files=d["number_of_processors"]

    complete_impedance=np.zeros((freq_list.shape[0],3),float)

    if d["VTA_approx"]==1:     # for Full IFFT we have to re-sort only impedance results
        cnt_freq=0
        for core in range(n_files):
            if (d["CPE_activ"]==1 or d["current_control"]==1) and cc_multicontact==False:
                impedance_get=read_csv(os.environ['PATIENTDIR']+'/Field_solutions/Impedance'+str(core)+'.csv', delimiter=' ', header=None)
                impedance=impedance_get.values
                complete_impedance[cnt_freq:cnt_freq+impedance.shape[0],:]
                cnt_freq=cnt_freq+impedance.shape[0]

        if (d["CPE_activ"]==1 or d["current_control"]==1) and cc_multicontact==False:           # we calculate impedance with FFEM only for these cases
            sorted_impedance=complete_impedance[complete_impedance[:,2].argsort(axis=0)]      #sort by freq
            np.savetxt(os.environ['PATIENTDIR']+'/Field_solutions/sorted_impedance.csv', sorted_impedance, delimiter=" ")

        minutes = int((tm.time() - start_paral)/60)
        secnds = int(tm.time() - start_paral)-minutes*60
        logging.critical("----- Parallel calculations took {} min {} sec -----\n".format(minutes, secnds))

        return True
    else:
        inx_compl_sol=0
        complete_solution=np.zeros((freq_list.shape[0]*Vertices.shape[0],6),float)
        cnt_freq=0
        for core in range(n_files):
            hf = h5py.File(os.environ['PATIENTDIR']+'/Field_solutions/sol_cor'+str(core)+'.h5', 'r')
            lst=list(hf.keys())
            result_total=[]
            for i in lst:
                a=hf.get(i)
                a=np.array(a)
                result_total.append(a)
            result=np.concatenate(result_total)
            hf.close()

            complete_solution[inx_compl_sol:inx_compl_sol+result.shape[0],:]=result
            inx_compl_sol=inx_compl_sol+result.shape[0]
            if (d["CPE_activ"]==1 or d["current_control"]==1) and cc_multicontact==False:   # we calculate impedance with FFEM only for these cases
                impedance_get=read_csv(os.environ['PATIENTDIR']+'/Field_solutions/Impedance'+str(core)+'.csv', delimiter=' ', header=None)
                impedance=impedance_get.values
                complete_impedance[cnt_freq:cnt_freq+impedance.shape[0],:]=impedance
                cnt_freq=cnt_freq+impedance.shape[0]

        if (d["CPE_activ"]==1 or d["current_control"]==1) and cc_multicontact==False:   # we calculate impedance with FFEM only for these cases
            sorted_impedance=complete_impedance[complete_impedance[:,2].argsort(axis=0)]      #sort by freq
            np.savetxt(os.environ['PATIENTDIR']+'/Field_solutions/sorted_impedance.csv', sorted_impedance, delimiter=" ")
            plt.figure(101010)
            plt.plot(sorted_impedance[:,0],sorted_impedance[:,1],marker='o')
            #plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
            ###           ncol=2, mode="expand", borderaxespad=0.)
            plt.xlabel(r'Re(Z), $\Omega$')
            plt.ylabel(r'Im(Z), $\Omega$')
            plt.grid(True)
            plt.savefig(os.environ['PATIENTDIR']+'/Images/Imp_plot.eps', format='eps', dpi=1000)

            Ampl_imp=np.zeros(sorted_impedance.shape[0],dtype=float)
            for i_fr in range(sorted_impedance.shape[0]):
                Ampl_imp[i_fr]=np.sqrt(sorted_impedance[i_fr,0]*sorted_impedance[i_fr,0]+sorted_impedance[i_fr,1]*sorted_impedance[i_fr,1])
            plt.figure(2)
            plt.plot(sorted_impedance[:,2],Ampl_imp[:],marker='o')
            plt.xscale("log")
            plt.xlabel('f, Hz')
            plt.ylabel(r'Ampl(Z), $\Omega$')
            plt.grid(True)
            plt.savefig(os.environ['PATIENTDIR']+'/Images/Imp_Ampl_plot.eps', format='eps', dpi=1000)

        sort_full_solution(d,freq_list,complete_solution,number_of_points)
        del complete_solution

        minutes=int((tm.time() - start_paral)/60)
        secnds=int(tm.time() - start_paral)-minutes*60
        logging.critical("----- Parallel calculations took {} min {} sec -----\n".format(minutes, secnds))


        if Field_on_VTA==1:
            d["VTA_approx"]=1

        return True
