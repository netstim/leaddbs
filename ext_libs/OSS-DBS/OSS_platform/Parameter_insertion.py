# -*- coding: utf-8 -*-
"""
Created on Sun Oct 14 12:24:38 2018

@author: butenko
"""

import os
import fileinput

def paste_geom_dim(x_length,y_length,z_length,X_tip,Y_tip,Z_tip):
    #end=" "
    DX_MRI_line="DX_MRI"
    DX_MRI_input="DX_MRI={}\n".format(x_length)            #NEURON uses ms
    DY_MRI_line="DY_MRI"
    DY_MRI_input="DY_MRI={}\n".format(y_length)            #NEURON uses ms
    DZ_MRI_line="DZ_MRI"
    DZ_MRI_input="DZ_MRI={}\n".format(z_length)            #NEURON uses ms
    
    X_tip_line="X_tip"
    X_tip_input="X_tip={}\n".format(X_tip)            #NEURON uses ms
    Y_tip_line="Y_tip"
    Y_tip_input="Y_tip={}\n".format(Y_tip)            #NEURON uses ms
    Z_tip_line="Z_tip"
    Z_tip_input="Z_tip={}\n".format(Z_tip)            #NEURON uses ms

    fl = fileinput.input(files="Brain_substitute.py", inplace=1)
    for line in fl:
        if line.startswith(DX_MRI_line):
            line = DX_MRI_input
        if line.startswith(DY_MRI_line):
            line = DY_MRI_input
        if line.startswith(DZ_MRI_line):
            line = DZ_MRI_input
        if line.startswith(X_tip_line):
            line = X_tip_input
        if line.startswith(Y_tip_line):
            line = Y_tip_input
        if line.startswith(Z_tip_line):
            line = Z_tip_input
        print(line, end="")
    fl.close()
    
    return True
    
def paste_neuron_model_param(d,N_models,t_stps):
    
    #end=" "
    dt_line="dt"
    dt_input="dt={}\n".format(1000.0*d["t_step"])            #NEURON uses ms
    
    tstop_line="tstop"
    tstop_input="tstop={}\n".format(1000.0/d["freq"])            #NEURON uses ms
    
    nseg_line="n_Ranvier"
    nseg_input="n_Ranvier={}\n".format(d["n_Ranvier"])            #in one model
    
    nmod_line="N_models"
    nmod_input="N_models={}\n".format(N_models)            
    
    nv_init_line="v_init"
    nv_init_input="v_init={}\n".format(d["v_init"])            #normally, -80mv
    
    ndiam_line="'diameter':"
    ndiam_input="'diameter':{}\n".format(d["diam_fib"])            #fiber diameter
    
    ntsteps_line="t_steps"
    ntsteps_input="t_steps={}\n".format(t_stps)            #fiber diameter
    
    Ampl_scale_line="Ampl_scale"
    Ampl_scale_input="Ampl_scale={}\n".format(d["Ampl_scale"])            #fiber diameter

    n_proc_line="n_processors"
    n_proc_input="n_processors={}\n".format(d["number_of_processors"])            #for parallel
    
    #trial_line="Tis_max"
    #input_line="Tis_max={}\n".format(Tis_max)
    
    #x = fileinput.input(files="Axon_files/NEURON_sep_points_parallel.py", inplace=1)
    #for line in x:
    #    if line.startswith(dt_line):
    #        line = dt_input
    #    print(line,end="")
    #x.close()
    

 
    x = fileinput.input(files="Axon_files/NEURON_sep_points_parallel_python3.py", inplace=1)    
    for line in x:
        if line.startswith(dt_line):
            line = dt_input
        #print line,
        if line.startswith(tstop_line):
            line = tstop_input
        #print line,
        if line.startswith(nseg_line):
            line = nseg_input
        #print line,
        if line.startswith(nmod_line):
            line = nmod_input
        if line.startswith(nv_init_line):
            line = nv_init_input
        if line.startswith(ndiam_line):
            line = ndiam_input
        if line.startswith(ntsteps_line):
            line = ntsteps_input
        if line.startswith(Ampl_scale_line):
            line = Ampl_scale_input
        if line.startswith(n_proc_line):
            line = n_proc_input
        print(line,end="")
    x.close()
    
def paste_dif_neuron_model_param(d,N_models,t_stps,population_index,last_point):
    
    #end=" "
    dt_line="dt"
    dt_input="dt={}\n".format(1000.0*d["t_step"])            #NEURON uses ms
    
    tstop_line="tstop"
    tstop_input="tstop={}\n".format(1000.0/d["freq"])            #NEURON uses ms
    
    nseg_line="n_Ranvier"
    nseg_input="n_Ranvier={}\n".format(d["n_Ranvier"][population_index])            #in one model
    
    nmod_line="N_models"
    nmod_input="N_models={}\n".format(N_models)            
    
    nv_init_line="v_init"
    nv_init_input="v_init={}\n".format(d["v_init"])            #normally, -80mv
    
    ndiam_line="'diameter':"
    ndiam_input="'diameter':{}\n".format(d["diam_fib"])            #fiber diameter
    
    ntsteps_line="t_steps"
    ntsteps_input="t_steps={}\n".format(t_stps)            #fiber diameter
    
    Ampl_scale_line="Ampl_scale"
    Ampl_scale_input="Ampl_scale={}\n".format(d["Ampl_scale"])            #fiber diameter

    n_proc_line="n_processors"
    n_proc_input="n_processors={}\n".format(d["number_of_processors"])            #for parallel
    
    last_point_line="last_point"
    last_point_input="last_point={}\n".format(last_point)
    
    population_line="population_index"
    population_input="population_index={}\n".format(population_index)
    
    #trial_line="Tis_max"
    #input_line="Tis_max={}\n".format(Tis_max)
    
    #x = fileinput.input(files="Axon_files/NEURON_sep_points_parallel.py", inplace=1)
    #for line in x:
    #    if line.startswith(dt_line):
    #        line = dt_input
    #    print(line,end="")
    #x.close()
    
    
    x = fileinput.input(files="Axon_files/NEURON_sep_points_parallel_python3_dif_axons.py", inplace=1)
    for line in x:
        if line.startswith(dt_line):
            line = dt_input
        #print line,
        if line.startswith(tstop_line):
            line = tstop_input
        #print line,
        if line.startswith(nseg_line):
            line = nseg_input
        #print line,
        if line.startswith(nmod_line):
            line = nmod_input
        if line.startswith(nv_init_line):
            line = nv_init_input
        if line.startswith(ndiam_line):
            line = ndiam_input
        if line.startswith(ntsteps_line):
            line = ntsteps_input
        if line.startswith(Ampl_scale_line):
            line = Ampl_scale_input
        if line.startswith(n_proc_line):
            line = n_proc_input
        if line.startswith(last_point_line):
            line = last_point_input    
        if line.startswith(population_line):
            line = population_input             
        print(line,end="")
    x.close()    
    
    
def paste_to_hoc(axonnodes,paranodes1,paranodes2,axoninter,axontotal,v_init,fiberD,paralength1,paralength2,nodelength,nodeD,axonD,paraD1,paraD2,deltax,nl):
    #end=" "
    axonnodes_line="axonnodes="
    axonnodes_input="axonnodes={}\n".format(axonnodes)            
    
    paranodes1_line="paranodes1="
    paranodes1_input="paranodes1={}\n".format(paranodes1)            
    
    paranodes2_line="paranodes2="
    paranodes2_input="paranodes2={}\n".format(paranodes2)
    
    axoninter_line="axoninter="
    axoninter_input="axoninter={}\n".format(axoninter)

    axontotal_line="axontotal="
    axontotal_input="axontotal={}\n".format(axontotal)           
    
    nv_init_line="v_init="
    nv_init_input="v_init={}\n".format(v_init)            #normally, -80mv
    
    fiberD_line="fiberD="
    fiberD_input="fiberD={}\n".format(fiberD)            #fiber diameter
    
    paralength1_line="paralength1="
    paralength1_input="paralength1={}\n".format(paralength1)            
    
    paralength2_line="paralength2="
    paralength2_input="paralength2={}\n".format(paralength2)
    
    nodelength_line="nodelength="
    nodelength_input="nodelength={}\n".format(nodelength)
    
    nodeD_line="nodeD="
    nodeD_input="nodeD={}\n".format(nodeD)            
    
    axonD_line="axonD="
    axonD_input="axonD={}\n".format(axonD)            
    
    paraD1_line="paraD1="
    paraD1_input="paraD1={}\n".format(paraD1)            
    
    paraD2_line="paraD2="
    paraD2_input="paraD2={}\n".format(paraD2)            
    
    deltax_line="deltax="
    deltax_input="deltax={}\n".format(deltax)            
    
    nl_line="nl="
    nl_input="nl={}\n".format(nl)            
    
    #trial_line="Tis_max"
    #input_line="Tis_max={}\n".format(Tis_max)
    
    
    x = fileinput.input(files="Axon_files/axon4pyfull.hoc", inplace=1)
    for line in x:
        if line.startswith(axonnodes_line):
            line = axonnodes_input
        if line.startswith(paranodes1_line):
            line = paranodes1_input
        if line.startswith(paranodes2_line):
            line = paranodes2_input
        if line.startswith(axoninter_line):
            line = axoninter_input
        if line.startswith(axontotal_line):
            line = axontotal_input
        if line.startswith(nv_init_line):
            line = nv_init_input
        if line.startswith(fiberD_line):
            line = fiberD_input
        if line.startswith(paralength1_line):
            line = paralength1_input
        if line.startswith(paralength2_line):
            line = paralength2_input
        if line.startswith(nodelength_line):
            line = nodelength_input
            
        if line.startswith(nodeD_line):
            line = nodeD_input
        if line.startswith(axonD_line):
            line = axonD_input
        if line.startswith(paraD1_line):
            line = paraD1_input
        if line.startswith(paraD2_line):
            line = paraD2_input
        if line.startswith(deltax_line):
            line = deltax_input
        if line.startswith(nl_line):
            line = nl_input
        print(line,end="")
    x.close()
    
    return True

    
def paste_paraview_vis(Points_on_model,N_comp_in_between):
    #end=" "
    NPoints_line="Points_on_model"
    NPoints_input="Points_on_model={}\n".format(Points_on_model)            #NEURON uses ms
    N_comp_in_between_line="N_comp_in_between"
    N_comp_in_between_input="N_comp_in_between={}\n".format(N_comp_in_between)            #NEURON uses ms


    fl = fileinput.input(files="Visualization_files/Paraview_vis_axon.py", inplace=1)
    for line in fl:
        if line.startswith(NPoints_line):
            line = NPoints_input
        if line.startswith(N_comp_in_between_line):
            line = N_comp_in_between_input
        print(line,end="")
    fl.close()
    
    return True

def paste_paraview_connections_vis(list_of_connections):
    #end=" "
    NPoints_line="list_of_connections="
    NPoints_input="list_of_connections={}\n".format(list_of_connections)          

    fl = fileinput.input(files="Visualization_files/Paraview_connections_processed.py", inplace=1)
    for line in fl:
        if line.startswith(NPoints_line):
            line = NPoints_input
        print(line,end="")
    fl.close()
    
    fl = fileinput.input(files="Visualization_files/Paraview_connections_activation.py", inplace=1)
    for line in fl:
        if line.startswith(NPoints_line):
            line = NPoints_input
        print(line,end="")
    fl.close()
    
    return True


def paste_paraview_clipping(x_t,y_t,z_t):
    #end=" "
    #clip1.ClipType.Origin = [1.0532820108523, -0.793227985190315, 2.0314999999999994]
    clip_line="clip1.ClipType.Origin"
    clip_input="clip1.ClipType.Origin=[{},{},{}]\n".format(x_t,y_t,z_t)
    
    cell_clip_line="extractCellsByRegion1.IntersectWith.Origin"
    cell_clip_input="extractCellsByRegion1.IntersectWith.Origin=[{},{},{}]\n".format(x_t,y_t,z_t)

#    fl = fileinput.input(files="Paraview_adapted_field.py", inplace=1)
#    for line in fl:
#        if line.startswith(clip_line):
#            line = clip_input
#        print line,
#    fl.close()
    
    fl = fileinput.input(files="Visualization_files/Paraview_adapted.py", inplace=1)
    for line in fl:
        if line.startswith(cell_clip_line):
            line = cell_clip_input
        print(line,end="")
    fl.close()
    
    fl = fileinput.input(files="Visualization_files/Paraview_CSFref.py", inplace=1)
    for line in fl:
        if line.startswith(cell_clip_line):
            line = cell_clip_input
        print(line,end="")
    fl.close()
    
    fl = fileinput.input(files="Visualization_files/Paraview_CSFunref.py", inplace=1)
    for line in fl:
        if line.startswith(cell_clip_line):
            line = cell_clip_input
        print(line,end="")
    fl.close()
    
    fl = fileinput.input(files="Visualization_files/Paraview_field_CSF_ref.py", inplace=1)
    for line in fl:
        if line.startswith(clip_line):
            line = clip_input
        print(line,end="")
    fl.close()
    
    return True
    
def insert_f_name(f_nam,file_name):
    #end=" "
    f_line="f_name="
    f_input='f_name={}\n'.format(f_nam)
    
    fl = fileinput.input(files=file_name, inplace=1)
    for line in fl:
        if line.startswith(f_line):
            line = f_input
        print(line,end="")
    fl.close()
    
    return True
    
#    fl = fileinput.input(files="Paraview_adapted_field.py", inplace=1)
#    for line in fl:
#        if line.startswith(clip_line):
#            line = clip_input
#        print line,
#    fl.close()
    

def get_Para_Array_name (vtu_file_name):            #this will not insert anything, but will find the name of the vector in .vtu file
    from Profile_Process_V6 import words_detect
    f=open(vtu_file_name,'r');  
    for index,line in enumerate(f):      
        var_list = words_detect('Scalars=',line);
        if (var_list[0]):
            word=line[var_list[1]+8::];###  <PointData  Scalars="f_127"> 
            return (str(word[:len(word)-3]));
            break;
    f.close();
