#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 13:31:19 2019

@author: butenko
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 15:41:38 2019

@author: butenko
"""

from dolfin import *
import numpy as np
import os
from pandas import read_csv
import pickle

import time as time_lib

def get_field_in_time(d,FR_vector_signal,Xs_signal_norm,t_vector):
    print("\n----- Conducting signal scaling and IFFT for the whole computational domain -----")

    Phi_vector_active_non_zero=[x for x in d["Phi_vector"] if (x is not None) and (x!=0.0)]
    cc_multicontact=False
    if d["current_control"]==1 and len(Phi_vector_active_non_zero)>1:       #multicontact current-controlled case
        cc_multicontact=True

    start_full_IFFT=time_lib.time()

    mesh = Mesh()
    f = HDF5File(mesh.mpi_comm(),os.environ['PATIENTDIR']+'/Field_solutions_functions/solution'+str(d["freq"]*1.0)+'.h5','r')
    f.read(mesh,"mesh", False)
    if d["EQS_core"]=='EQS':
        Er = FiniteElement("Lagrange", mesh.ufl_cell(),d["el_order"])
        Ei = FiniteElement("Lagrange", mesh.ufl_cell(),d["el_order"])
        Ec = Er * Ei
        V = FunctionSpace(mesh, Ec)
        phi_sol=Function(V)

        f.read(phi_sol,'solution_full')
        phi_r_sol1,phi_i_sol1=phi_sol.split(deepcopy=True)
        f.close()
    if d["EQS_core"]=='QS':
        V = FunctionSpace(mesh, "Lagrange",d["el_order"])
        phi_r_sol1=Function(V)
        f.read(phi_r_sol1,'solution_full')
        f.close()
        phi_i_sol1=Function(V)
        phi_i_sol1.vector()[:] = 0.0
        f.close()


    file=File(os.environ['PATIENTDIR']+'/Field_solutions_functions/phi_r_sol_unscaled_'+str(d["freq"]*1.0)+'.pvd')
    file<<phi_r_sol1

    #Xs_recovered = np.genfromtxt('Stim_Signal/Xs_storage_full.csv', delimiter=' ')
    #Xs_signal_norm=np.vectorize(complex)(Xs_recovered[:,0],Xs_recovered[:,1])
    #FR_vector_signal = np.genfromtxt('Stim_Signal/FR_vector_signal.csv', delimiter=' ')

    if d["spectrum_trunc_method"]=='Octave Band Method':
        FR_vector_signal_new = np.genfromtxt(os.environ['PATIENTDIR']+'/Stim_Signal/FR_vector_signal_octaves'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')

        Fr_corresp_arr = np.genfromtxt(os.environ['PATIENTDIR']+'/Stim_Signal/Fr_corresp_array'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')
        Fr_octave_vect = np.genfromtxt(os.environ['PATIENTDIR']+'/Stim_Signal/Fr_octave_vector_'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')

        Fr_corresp_arr=np.round(Fr_corresp_arr,6)
        Fr_octave_vect=np.round(Fr_octave_vect,6)

    #t_vector = np.genfromtxt('Stim_Signal/t_vector.csv', delimiter=' ')
    Phi_1=np.vectorize(complex)(phi_r_sol1.vector()[:],phi_i_sol1.vector()[:])
    Freq_number=FR_vector_signal.shape[0]

    Phi_Global=np.complex(1.0,1.0)*np.zeros((Freq_number,Phi_1.shape[0]),float)

    #phi_point=np.zeros((Freq_number,2),float)



    #==============================Load Full===================================#
    if d["spectrum_trunc_method"]=='No Truncation':
        i_fr=0
        for fr in FR_vector_signal:
            if os.path.isfile(os.environ['PATIENTDIR']+'/Field_solutions_functions/solution'+str(np.round(fr,6))+'.h5'):
                f = HDF5File(mesh.mpi_comm(),os.environ['PATIENTDIR']+'/Field_solutions_functions/solution'+str(np.round(fr,6))+'.h5','r')
                #print("Freq: ", fr)

                if d["EQS_core"]=='EQS':
                    Er = FiniteElement("Lagrange", mesh.ufl_cell(),d["el_order"])
                    Ei = FiniteElement("Lagrange", mesh.ufl_cell(),d["el_order"])
                    Ec = Er * Ei
                    V = FunctionSpace(mesh, Ec)
                    phi_sol=Function(V)

                    f.read(phi_sol,'solution_full')
                    phi_r_sol,phi_i_sol=phi_sol.split(deepcopy=True)
                    f.close()
                if d["EQS_core"]=='QS':
                    V = FunctionSpace(mesh, "Lagrange",d["el_order"])
                    phi_r_sol=Function(V)
                    f.read(phi_r_sol,'solution_full')
                    f.close()
                    phi_i_sol=Function(V)
                    phi_i_sol.vector()[:] = 0.0
                    f.close()
                if d["current_control"] == 1 and cc_multicontact==False:
                    with open(os.environ['PATIENTDIR']+'/Field_solutions_functions/current_scale'+str(np.round(fr,6))+'.file', "rb") as fpickle:
                        Currents = pickle.load(fpickle)

                    a=np.real((phi_r_sol.vector()[:]+1j*phi_i_sol.vector()[:])/(Currents[0]+1j*Currents[1]))
                    b=np.imag((phi_r_sol.vector()[:]+1j*phi_i_sol.vector()[:])/(Currents[0]+1j*Currents[1]))
                    Phi_complex=np.vectorize(complex)(a,b)
                else:
                    Phi_complex=np.vectorize(complex)(phi_r_sol.vector()[:],phi_i_sol.vector()[:])

                #Phi_Global_real[i,:]=phi_r_sol.vector()[:]
                #Phi_Global_imag[i,:]=phi_i_sol.vector()[:]
                Phi_Global[i_fr,:]=Phi_complex

            else:
                print("wrong frequency")
                raise SystemExit
                #Phi_complex=np.vectorize(complex)(phi_r_sol1.vector()[:]*0.0,phi_i_sol1.vector()[:]*0.0)
                #Phi_Global[i_fr,:]=Phi_complex

            i_fr=i_fr+1
        i_fr=0
#==============================Trunc1and2===================================#
    if d["spectrum_trunc_method"]=='High Amplitude Method' or d["spectrum_trunc_method"]=='Cutoff Method':
        i_fr=0
        for fr in FR_vector_signal:
            if os.path.isfile(os.environ['PATIENTDIR']+'/Field_solutions_functions/solution'+str(np.round(fr,6))+'.h5'):
                f = HDF5File(mesh.mpi_comm(),os.environ['PATIENTDIR']+'/Field_solutions_functions/solution'+str(np.round(fr,6))+'.h5','r')
                #print("Freq: ", fr)
                if d["EQS_core"]=='EQS':
                    Er = FiniteElement("Lagrange", mesh.ufl_cell(),d["el_order"])
                    Ei = FiniteElement("Lagrange", mesh.ufl_cell(),d["el_order"])
                    Ec = Er * Ei
                    V = FunctionSpace(mesh, Ec)
                    phi_sol=Function(V)

                    f.read(phi_sol,'solution_full')
                    phi_r_sol,phi_i_sol=phi_sol.split(deepcopy=True)
                    f.close()
                if d["EQS_core"]=='QS':
                    V = FunctionSpace(mesh, "Lagrange",d["el_order"])
                    phi_r_sol=Function(V)
                    f.read(phi_r_sol,'solution_full')
                    f.close()
                    phi_i_sol=Function(V)
                    phi_i_sol.vector()[:] = 0.0
                    f.close()

                if d["current_control"] == 1 and cc_multicontact==False:
                    with open(os.environ['PATIENTDIR']+'/Field_solutions_functions/current_scale'+str(np.round(fr,6))+'.file', "rb") as fpickle:
                        Currents = pickle.load(fpickle)

                    a=np.real((phi_r_sol.vector()[:]+1j*phi_i_sol.vector()[:])/(Currents[0]+1j*Currents[1]))
                    b=np.imag((phi_r_sol.vector()[:]+1j*phi_i_sol.vector()[:])/(Currents[0]+1j*Currents[1]))
                    Phi_complex=np.vectorize(complex)(a,b)
                else:
                    Phi_complex=np.vectorize(complex)(phi_r_sol.vector()[:],phi_i_sol.vector()[:])



                Phi_Global[i_fr,:]=Phi_complex

                #phi_point[i_fr,0]=phi_r_sol(Point(0.7071,0.982233,0.705))
                #phi_point[i_fr,1]=phi_i_sol(Point(0.7071,0.982233,0.705))
            else:
                #print("wrong frequency")
                #raise SystemExit
                Phi_complex=np.vectorize(complex)(phi_r_sol1.vector()[:]*0.0,phi_i_sol1.vector()[:]*0.0)
                Phi_Global[i_fr,:]=Phi_complex

                #phi_point[i_fr,0]=0.0
                #phi_point[i_fr,1]=0.0

            i_fr=i_fr+1
        i_fr=0
    #==============================Load Octaves===================================#
    if d["spectrum_trunc_method"]=='Octave Band Method':
        i_fr=0

        i_start_octv_rslt=np.where(np.round(FR_vector_signal_new,6)==Fr_octave_vect[0])
        i_start_octv=i_start_octv_rslt[0][0]

        stepper=0
        for fr in FR_vector_signal_new:
            if i_fr>=i_start_octv:

                rslt=np.where(Fr_corresp_arr[:,0]==np.round(fr,6))
                step_octv=rslt[0].shape[0]   #size of the freq. pack in the octave

                f = HDF5File(mesh.mpi_comm(),os.environ['PATIENTDIR']+'/Field_solutions_functions/solution'+str(np.round(fr,6))+'.h5','r')
                #f = HDF5File(mesh.mpi_comm(),'Field_solutions_functions/solution'+str(int(fr))+'.h5','r')
                #print("Freq: ", i_fr)
                if d["EQS_core"]=='EQS':
                    Er = FiniteElement("Lagrange", mesh.ufl_cell(),d["el_order"])
                    Ei = FiniteElement("Lagrange", mesh.ufl_cell(),d["el_order"])
                    Ec = Er * Ei
                    V = FunctionSpace(mesh, Ec)
                    phi_sol=Function(V)

                    f.read(phi_sol,'solution_full')
                    phi_r_sol,phi_i_sol=phi_sol.split(deepcopy=True)
                    f.close()
                if d["EQS_core"]=='QS':
                    V = FunctionSpace(mesh, "Lagrange",d["el_order"])
                    phi_r_sol=Function(V)
                    f.read(phi_r_sol,'solution_full')
                    f.close()
                    phi_i_sol=Function(V)
                    phi_i_sol.vector()[:] = 0.0
                    f.close()
                if d["current_control"] == 1 and cc_multicontact==False:
                    with open(os.environ['PATIENTDIR']+'/Field_solutions_functions/current_scale'+str(np.round(fr,6))+'.file', "rb") as fpickle:
                        Currents = pickle.load(fpickle)

                    a=np.real((phi_r_sol.vector()[:]+1j*phi_i_sol.vector()[:])/(Currents[0]+1j*Currents[1]))
                    b=np.imag((phi_r_sol.vector()[:]+1j*phi_i_sol.vector()[:])/(Currents[0]+1j*Currents[1]))
                    Phi_complex=np.vectorize(complex)(a,b)
                else:
                    Phi_complex=np.vectorize(complex)(phi_r_sol.vector()[:],phi_i_sol.vector()[:])

                Phi_Global[stepper:stepper+step_octv,:]=Phi_complex

                stepper=stepper+step_octv
                #if global_i_point==0:
                    #print "hrreee!"
                    #print step_octv
            else:
                f = HDF5File(mesh.mpi_comm(),os.environ['PATIENTDIR']+'/Field_solutions_functions/solution'+str(np.round(fr,6))+'.h5','r')
                #f = HDF5File(mesh.mpi_comm(),'Field_solutions_functions/solution'+str(int(fr))+'.h5','r')
                #print("Freq: ", i_fr)
                if d["EQS_core"]=='EQS':
                    Er = FiniteElement("Lagrange", mesh.ufl_cell(),d["el_order"])
                    Ei = FiniteElement("Lagrange", mesh.ufl_cell(),d["el_order"])
                    Ec = Er * Ei
                    V = FunctionSpace(mesh, Ec)
                    phi_sol=Function(V)

                    f.read(phi_sol,'solution_full')
                    phi_r_sol,phi_i_sol=phi_sol.split(deepcopy=True)
                    f.close()
                if d["EQS_core"]=='QS':
                    V = FunctionSpace(mesh, "Lagrange",d["el_order"])
                    phi_r_sol=Function(V)
                    f.read(phi_r_sol,'solution_full')
                    f.close()
                    phi_i_sol=Function(V)
                    phi_i_sol.vector()[:] = 0.0
                    f.close()

                if d["current_control"] == 1 and cc_multicontact==False:
                    with open(os.environ['PATIENTDIR']+'/Field_solutions_functions/current_scale'+str(np.round(fr,6))+'.file', "rb") as fpickle:
                        Currents = pickle.load(fpickle)

                    a=np.real((phi_r_sol.vector()[:]+1j*phi_i_sol.vector()[:])/(Currents[0]+1j*Currents[1]))
                    b=np.imag((phi_r_sol.vector()[:]+1j*phi_i_sol.vector()[:])/(Currents[0]+1j*Currents[1]))
                    Phi_complex=np.vectorize(complex)(a,b)
                else:
                    Phi_complex=np.vectorize(complex)(phi_r_sol.vector()[:],phi_i_sol.vector()[:])

                Phi_Global[stepper,:]=Phi_complex

                stepper=stepper+1

            i_fr=i_fr+1

        i_fr=0

    #np.save('Field_solutions/phi_point', phi_point)
    ##=============================================================================#
    #np.save('Stim_Signal_IFFT/Phi_Global_real', Phi_Global.real)
    #np.save('Stim_Signal_IFFT/Phi_Global_imag', Phi_Global.imag)
    #Phi_1_conv=np.multiply(Phi_1,Xs_signal_norm[1])
    #Phi_1_conv=Xs_signal_norm[1]*Phi_1
    #Phi_1_conv=np.ascontiguousarray(Phi_1_conv, dtype=np.float32)
    del f
    del Phi_complex
    #
    #print(Phi_Global.T[0,0:100])
    Phi_Global_conv_T=np.multiply(Phi_Global.T,Xs_signal_norm)
    #Phi_Global_conv_T=np.multiply(Xs_signal_norm.T*Phi_Global)
    #Phi_Global_conv_T=Phi_Global.T*Xs_signal_norm
    del Phi_Global

    Phi_Global_conv=Phi_Global_conv_T.T
    del Phi_Global_conv_T

    #Phi_Global_conv=np.ascontiguousarray(Phi_Global_conv, dtype=np.float32)

    t_vect=t_vector

    #Signal_t_conv=np.zeros((t_vect.shape[0],Phi_Global_conv.shape[1]),float)
    #if np.mod(t_vect.shape[0], 2):  # if the FT vector is odd
    #    fv_conj = np.conjugate(Phi_Global_conv[-1:0:-1,:])
    #else:  # if the FT vector is even
    #    fv_conj = np.conjugate(Phi_Global_conv[-2:0:-1,:])
    #
    #Y = np.concatenate((Phi_Global_conv, fv_conj), axis=0)
    #Signal_t_conv=np.fft.ifft(Y).real

    '''CAN BE parallelized'''
    Signal_t_conv=np.zeros((t_vect.shape[0],Phi_Global_conv.shape[1]),float)

    #print(Phi_Global_conv[0:100,0])

    dofs_quart=[int(Phi_Global_conv.shape[1]/4.0),int(2*Phi_Global_conv.shape[1]/4.0),int(3*Phi_Global_conv.shape[1]/4.0),int(4*Phi_Global_conv.shape[1]/4.0)]

    for j in range(Phi_Global_conv.shape[1]):
        if j in dofs_quart:
            print(int(j*100/Phi_Global_conv.shape[1])+1,"% of dofs were processed")
        #if j%10000==0:
            #print("dof number: ", j)
        Phi_Global_conv_vert=Phi_Global_conv[:,j]
    #    Phi_Global_conv_vert=np.vectorize(complex)(a[:,j].T,b[:,j].T)
        if np.mod(t_vect.shape[0], 2):  # if the FT vector is odd
            fv_conj = np.conjugate(Phi_Global_conv_vert[-1:0:-1])
        else:  # if the FT vector is even
            fv_conj = np.conjugate(Phi_Global_conv_vert[-2:0:-1])

        Y = np.concatenate((Phi_Global_conv_vert, fv_conj), axis=0)

        Signal_t_conv_vert=np.fft.ifft(Y).real
        Signal_t_conv[:,j]=Signal_t_conv_vert
        #if j>10000:
         #   break

    del Phi_Global_conv
    Signal_t_conv=np.ascontiguousarray(Signal_t_conv, dtype=np.float32)
    np.save(os.environ['PATIENTDIR']+'/Field_solutions/Signal_t_conv', Signal_t_conv.real)
    #print(Signal_t_conv[0:40,0])
    #print(Signal_t_conv[0:40,10000])

    #========================Parallel=============================================#

    #def conduct_parallel_IFFT(j_vert):
    #
    #    global Phi_Global_conv
    #    global t_vect
    #
    #    print("vertice number: ", j_vert)
    #    Phi_Global_conv_vert=Phi_Global_conv[:,j_vert]
    #    if np.mod(t_vect.shape[0], 2):  # if the FT vector is odd
    #        fv_conj = np.conjugate(Phi_Global_conv_vert[-1:0:-1])
    #    else:  # if the FT vector is even
    #        fv_conj = np.conjugate(Phi_Global_conv_vert[-2:0:-1])
    #
    #    Y = np.concatenate((Phi_Global_conv_vert, fv_conj), axis=0)
    #
    #    Signal_t_conv_vert=np.fft.ifft(Y).real
    #
    #    return Signal_t_conv_vert
    #
    #    #np.savetxt('Fields_in_time/Signal_t_conv'+str(j_vert)+'.csv', Signal_t_conv_vert.real, delimiter=" ")
    #
    #
    #    #Signal_t_conv[:,j_freq]=Signal_t_conv_vert
    #import multiprocessing as mp
    #i=0
    #number_of_processors=10
    #while i < Phi_Global_conv.shape[1]:     #number of points
    #    #output = mp.Queue()         #defines an output queue
    #    #proc=[]
    #
    #    if Phi_Global_conv.shape[1]-i>=number_of_processors:
    #        pool=mp.Pool(processes=number_of_processors)
    #        Signal_t_pack=[pool.apply(conduct_parallel_IFFT,args=(i_point,)) for i_point in range(i,i+number_of_processors)]
    #
    #        #Signal_t_pack_array=np.zeros((t_vect.shape[0],number_of_processors),float)
    #        #Signal_t_pack_array[:,:]=Signal_t_pack[:]
    #        #print(type(Signal_t_pack_array))
    #        Signal_t_conv[:,i:i+number_of_processors]=np.asarray(Signal_t_pack).T
    ##        loc_inx=0
    ##        for j in range(i,i+number_of_processors):
    ##            Signal_t_conv[:,j]=Signal_t_pack[loc_inx]
    ##            loc_inx=loc_inx+1
    #
    #        i=i+number_of_processors
    #    else:
    #        pool=mp.Pool(processes=Phi_Global_conv.shape[1]-i)
    #        Signal_t_pack=[pool.apply(conduct_parallel_IFFT,args=(i_point,)) for i_point in range(i,Phi_Global_conv.shape[1])]
    #        #Signal_t_conv[:,i:Phi_Global_conv.shape[1]]=Signal_t_pack[0:(Phi_Global_conv.shape[1]-i)]
    #        #Signal_t_pack_array=np.zeros((t_vect.shape[0],Phi_Global_conv.shape[1]-i),float)
    #        #Signal_t_pack_array[:,:]=Signal_t_pack[:]
    #        #print(type(Signal_t_pack_array))
    #        Signal_t_conv[:,i:Phi_Global_conv.shape[1]]=np.asarray(Signal_t_pack).T
    #
    ##        loc_inx=0
    ##        for j in range(i,Phi_Global_conv.shape[1]):
    ##            Signal_t_conv[:,j]=Signal_t_pack[0:loc_inx]
    ##            loc_inx=loc_inx+1
    #        i=Phi_Global_conv.shape[1]
    #
    #
    #
    #    if i<10:
    #        print(Signal_t_pack[0])
    #        print(Signal_t_pack)
    #
    #i=0

    #Signal_t_conv=np.ascontiguousarray(Signal_t_conv, dtype=np.float32)
    #np.savetxt('MPI_stuff/Phi_Global_conv.csv', Phi_Global_conv, delimiter=" ")


    minutes=int((time_lib.time() - start_full_IFFT)/60)
    secnds=int(time_lib.time() - start_full_IFFT)-minutes*60
    print("--- Signal scaling and IFFT for the whole domain took ",minutes," min ",secnds," s ")

    if d["VTA_from_NEURON"]==True:
        start_VTA_NEURON=time_lib.time()
        Vert_get=read_csv(os.environ['PATIENTDIR']+'/Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', delimiter=' ', header=None)
        Vert=Vert_get.values
        Signal_t_points=np.zeros((Vert.shape[0],t_vect.shape[0]),float)
    else:
        start_VTA_divE=time_lib.time()


    file=File(os.environ['PATIENTDIR']+'/Animation_Field_in_time/Field.pvd')
    for i in range(t_vect.shape[0]):
        t_step=int(i)

        if t_step>d["t_step_end"] and d["t_step_end"]!=0:
            break

        if t_step%100==0 and t_step!=0:
            print(t_step," time steps are completed")
        phi_sol_time=Function(V)
        if d["EQS_core"]=='EQS':
            phi_r_time,phi_i_time=phi_sol_time.split(deepcopy=True)

        if d["EQS_core"]=='QS':
            phi_r_time=phi_sol_time


        phi_r_time.vector()[:]=np.real(Signal_t_conv[t_step,:])*d["Ampl_scale"]
        phi_r_time.rename("phi_r_time", "phi_r_time")
        file<<phi_r_time,t_step

        if d["VTA_from_NEURON"]==True:
            for i_point in range(Vert.shape[0]):
                pnt=Point(Vert[i_point,0],Vert[i_point,1],Vert[i_point,2])
                Signal_t_points[i_point,i]=phi_r_time(pnt)

    if d["VTA_from_NEURON"]==True:
        for i_point in range(Vert.shape[0]):
            #np.savetxt('Points_in_time/Signal_t_conv'+str(i_point)+'.csv', Signal_t_points[i_point,:], delimiter=" ")
            np.save(os.environ['PATIENTDIR']+'/Points_in_time/Signal_t_conv'+str(i_point), Signal_t_points[i_point,:])

        minutes=int((time_lib.time() - start_VTA_NEURON)/60)
        secnds=int(time_lib.time() - start_VTA_NEURON)-minutes*60
        print("----- Extracting a pointwise solution from the whole domain IFFT took ",minutes," min ",secnds," s -----\n")

        return 0.0

    if d["VTA_from_divE"]==True or d["VTA_from_E"]==True:
        #bad way to define where the potential is max because we give the time step explicitly in the beginning. Instead, we can sum up potential on the nodes over time and choose a step with the highest
        ind_max=np.sort(np.argpartition(np.absolute(Signal_t_conv[10,:]),-1*1)[-1*1:])          #first check where the max is on time step 10
        ind_max_2=np.sort(np.argpartition(np.absolute(Signal_t_conv[:,ind_max[0]]),-1*1)[-1*1:])        #now check max on the particular comp
        print ("VTA is estimated at step ", ind_max_2)

        '''VTA: div E-field'''
        W=VectorFunctionSpace(mesh,'DG',d["el_order"]-1)
        #W=FunctionSpace(mesh,'RT',1)
        #E_field=project(grad(-1.0*phi_r_sol),W,solver_type="cg", preconditioner_type="jacobi")

        phi_sol_time=Function(V)
        if d["EQS_core"]=='EQS':
            phi_r_time,phi_i_time=phi_sol_time.split(deepcopy=True)

        if d["EQS_core"]=='QS':
            phi_r_time=phi_sol_time

        phi_r_time.vector()[:]=np.real(Signal_t_conv[ind_max_2[0],:])*d["Ampl_scale"]

        phi_i_sol=Function(V)
        phi_i_sol.vector()[:] = 0.0



        w = TestFunction(W)
        Pv = TrialFunction(W)
        E_field_real = Function(W)
        a_local = inner(w, Pv) * dx
        L_local = inner(w, -grad(phi_r_time)) * dx

        A_local, b_local = assemble_system(a_local, L_local, bcs=[])

        local_solver = PETScKrylovSolver('bicgstab')
        local_solver.solve(A_local,E_field_real.vector(),b_local)

        file=File(os.environ['PATIENTDIR']+'/Field_solutions_functions/E_field_at_stim_peak.pvd')
        file<<E_field_real

        #V_for_Enorm = FunctionSpace(mesh, "DG", 2)
        #E_norm = project(sqrt(inner(E_field_real, E_field_real)), V_for_Enorm)



        W_unit=FunctionSpace(mesh,'CG',1)
        Unit_function = Function(W_unit)
        Unit_function.vector()[:]=1.0
        VTA_size=0.0    #in mm3

        if d["VTA_from_E"]==True:
            W_amp=FunctionSpace(mesh,'DG',2)
            w_amp = TestFunction(W_amp)
            Pv_amp = TrialFunction(W_amp)
            E_norm = Function(W_amp)
            a_local = inner(w_amp, Pv_amp) * dx
            L_local = inner(w_amp, sqrt(dot(E_field_real,E_field_real))) * dx
            A_local, b_local = assemble_system(a_local, L_local, bcs=[])
            local_solver = PETScKrylovSolver('bicgstab')
            local_solver.solve(A_local,E_norm.vector(),b_local)

            file=File(os.environ['PATIENTDIR']+'/Field_solutions_functions/E_norm_at_stim_peak.pvd')
            file<<E_norm

            for cell in cells(mesh):
                cell_size=assemble_local(Unit_function*dx,cell)
                if abs(assemble_local(sqrt(dot(E_field_real,E_field_real))*dx,cell))/cell_size>d["Activation_threshold_VTA"]:           #0.4 is based on Astrom and visual estimation
                    VTA_size=VTA_size+abs(assemble_local(Unit_function*dx,cell))

        if d["VTA_from_divE"]==True:
###========================VTA with E-field=================================###

#        W_amp=FunctionSpace(mesh,'CG',1)
#        w_amp = TestFunction(W_amp)
#        Pv_amp = TrialFunction(W_amp)
#        E_amp_real = Function(W_amp)
#        a_local = inner(w_amp, Pv_amp) * dx
#        L_local = inner(w_amp, sqrt(inner(E_field_real,E_field_real))) * dx
#        A_local, b_local = assemble_system(a_local, L_local, bcs=[])
#
#        local_solver = PETScKrylovSolver('bicgstab')
#        local_solver.solve(A_local,E_amp_real.vector(),b_local)
#
#        file=File('Field_solutions_functions/E_field_at_stim_peak_ampl.pvd')
#        file<<E_amp_real
#
#        Unit_function = Function(W_amp)
#        Unit_function.vector()[:]=1.0
#        VTA_size=0.0    #in mm3
#        for cell in cells(mesh):
#            cell_size=assemble_local(Unit_function*dx,cell)
#            if abs(assemble_local(E_amp_real*dx,cell))/cell_size>0.255:           #should be 255 V/m
#                VTA_size=VTA_size+abs(assemble_local(Unit_function*dx,cell))
#        print("VTA size in mm3: ", VTA_size)

###========================VTA with div(E)==================================###
            W_amp=FunctionSpace(mesh,'CG',d["el_order"]-1)

            from ufl import nabla_div
            Second_deriv=project(nabla_div(E_field_real),W_amp,solver_type="cg", preconditioner_type="amg")

            Second_deriv_abs=Second_deriv
            a_abs=Second_deriv_abs.vector()[:]
            Second_deriv_abs.vector()[:]=abs(a_abs)

            file=File(os.environ['PATIENTDIR']+'/Field_solutions_functions/Second_derivative.pvd')
            file<<Second_deriv_abs

            '''This is only for the particular example! Check, where the marked cells are'''
            #mesh_SN=Mesh('SN_med_shifted.xml')
            #mesh_STN=Mesh('STN_med_shifted.xml')
            #mesh_EPN=Mesh('EPN_med_shifted.xml')


            VTA_SN=0.0
            VTA_STN=0.0
            VTA_EPN=0.0

            for cell in cells(mesh):
                cell_size=assemble_local(Unit_function*dx,cell)
                if abs(assemble_local(Second_deriv_abs*dx,cell))/cell_size>d["Activation_threshold_VTA"]:           #0.4 is based on Astrom and visual estimation
                    VTA_size=VTA_size+abs(assemble_local(Unit_function*dx,cell))

#                #For the particular study, we need to check the place of activation
#                x_coord=cell.midpoint().x()
#                y_coord=cell.midpoint().y()
#                z_coord=cell.midpoint().z()
#                pnt=Point(x_coord,y_coord,z_coord)
#
#                if mesh_SN.bounding_box_tree().compute_first_entity_collision(pnt)<mesh_SN.num_cells()*100:
#                    VTA_SN=VTA_SN+abs(assemble_local(Unit_function*dx,cell))
#                elif mesh_STN.bounding_box_tree().compute_first_entity_collision(pnt)<mesh_STN.num_cells()*100:
#                    VTA_STN=VTA_STN+abs(assemble_local(Unit_function*dx,cell))
#                elif mesh_EPN.bounding_box_tree().compute_first_entity_collision(pnt)<mesh_EPN.num_cells()*100:
#                    VTA_EPN=VTA_EPN+abs(assemble_local(Unit_function*dx,cell))
#
#
        print("VTA size in mm3: ", VTA_size)
#        print("VTA size in SN in mm3: ", VTA_SN)
#        print("VTA size in STN in mm3: ", VTA_STN)
#        print("VTA size in EPN in mm3: ", VTA_EPN)

        minutes=int((time_lib.time() - start_VTA_divE)/60)
        secnds=int(time_lib.time() - start_VTA_divE)-minutes*60
        print("----- Direct evaluation of VTA from Full IFFT took ",minutes," min ",secnds," s -----\n")

        return VTA_size#,VTA_SN,VTA_STN,VTA_EPN


