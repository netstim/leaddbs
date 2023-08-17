# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 10:55:09 2018

@author: trieu
"""

import os

from Electrode_files.Profile_Process_V6 import words_detect

def create_geometry_script (d,Brain_map,electrode_profile,Impl_coords,Second_point_coords,ROI_radial,shift_to_PO,Vertice_enable):

   electrode_profile=electrode_profile
   # check parameters inputs
   name_idx=len(electrode_profile)
   check_profile_name = words_detect ('_profile.py', "Electrode_files/"+electrode_profile);
   if ( check_profile_name[0] == False):
       print ("ERROR: DBS lead profile name should be a string ends with _profile.py")
   else:
       f3=open("Electrode_files/"+electrode_profile,'r')
       f2=open(os.environ['PATIENTDIR']+"/"+electrode_profile[:name_idx-11] + '_position.py','w+') # new file with new position
       #print(electrode_profile[:name_idx-11] + '_position.py')
       for index,line in enumerate(f3):
               line_replace = False;
               var_list = words_detect("##### VARIABLE LIST #####",line)
               if (var_list[0]):
                   line_replace =True; #replace the code
                   f2.write("##### VARIABLE LIST #####\n"
                       +'Lead2nd_Enable = False\n'
                       +'Xm = {}\n'.format(shift_to_PO[0])
                       +'Ym = {}\n'.format(shift_to_PO[1])
                       +'Zm = {}\n'.format(shift_to_PO[2])
                       +'Xt = {}\n'.format(Impl_coords[0])
                       +'Yt = {}\n'.format(Impl_coords[1])
                       +'Zt = {}\n'.format(Impl_coords[2])
                       +'X_2nd = {}\n'.format(Second_point_coords[0])
                       +'Y_2nd = {}\n'.format(Second_point_coords[1])
                       +'Z_2nd = {}\n'.format(Second_point_coords[2])
                       +'OZ_angle = {}\n'.format(d["Rotation_Z"])
                       +'encap_thickness = {}\n'.format(d["encap_thickness"])
                       +'ROI_radial = {}\n'.format(ROI_radial)
                       +'Vertice_enable = {}\n'.format(Vertice_enable)
                       +"Brain_map = '{}'\n".format(Brain_map)
		       +"Phi_vector = {}\n".format(d["Pulse_amp"])
		       +"stretch = {}\n".format(d["stretch"])
                         );

               if(line_replace == False): f2.write(line);
       f2.close();
       f3.close();

###############################################################################
