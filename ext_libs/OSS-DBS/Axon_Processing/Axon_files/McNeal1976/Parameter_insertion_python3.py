#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 11:20:57 2019

@author: butenko
"""
import os
import fileinput

def paste_to_hoc_python3(axonnodes,axoninter,axontotal,v_init,steps_per_ms):
    #end=" "
    #print(os.path.dirname(os.path.realpath(__file__)))

    NNODES_line="NNODES ="
    NNODES_input="NNODES = {}\n".format(axonnodes)

    axonnodes_line="axonnodes="
    axonnodes_input="axonnodes={}\n".format(axonnodes)

    nv_init_line="v_init="
    nv_init_input="v_init={}\n".format(v_init)            #normally, -80mv

    steps_per_ms_line="steps_per_ms="
    steps_per_ms_input="steps_per_ms={}\n".format(steps_per_ms)

    #trial_line="Tis_max"
    #input_line="Tis_max={}\n".format(Tis_max)



    #x = fileinput.input(files="axon4pyfull.hoc", inplace=1)
    x = fileinput.input(files="init_B5_extracellular.hoc", inplace=1)
    for line in x:
        if line.startswith(axonnodes_line):
            line = axonnodes_input
        if line.startswith(nv_init_line):
            line = nv_init_input
        if line.startswith(steps_per_ms_line):
            line = steps_per_ms_input
        print(line,end="")
    x.close()

    x = fileinput.input(files="axon5.hoc", inplace=1)
    for line in x:
        if line.startswith(NNODES_line):
            line = NNODES_input
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
