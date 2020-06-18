# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 23:15:08 2018

@author: trieu
"""
from Electrode_files.Profile_Process_V6 import words_detect

def get_Para_Array_name (vtu_file_name):
    f=open(vtu_file_name,'r');  
    for index,line in enumerate(f):      
        var_list = words_detect('Scalars=',line);
        if (var_list[0]):
            word=line[var_list[1]+8::];###  <PointData  Scalars="f_127"> 
            return (str(word[:len(word)-3]));
            break;
    f.close();
############################### Test ##############
#vtu_file = "Results_adaptive/Phi_r_field_EQS000000.vtu";
#name = Para_Array_name(vtu_file);
#print name
