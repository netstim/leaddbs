#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 11:20:57 2019

@author: butenko
"""
import os
import fileinput

def paste_to_hoc_python3(axonnodes,paranodes1,paranodes2,axoninter,axontotal,v_init,fiberD,paralength1,paralength2,nodelength,nodeD,axonD,paraD1,paraD2,deltax,nl,steps_per_ms):
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

    steps_per_ms_line="steps_per_ms="
    steps_per_ms_input="steps_per_ms={}\n".format(steps_per_ms)       

    #trial_line="Tis_max"
    #input_line="Tis_max={}\n".format(Tis_max)
    
    
    #x = fileinput.input(files="Axon_files/axon4pyfull.hoc", inplace=1)
    x = fileinput.input(files="axon4pyfull.hoc", inplace=1)
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
        if line.startswith(steps_per_ms_line):
            line = steps_per_ms_input            
        print(line,end="")
    x.close()
    
    return True

def paste_paraview_vis_python3(Points_on_model,N_comp_in_between):
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
