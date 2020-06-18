# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 10:55:09 2018

@author: trieu
"""


from Electrode_files.Profile_Process_V6 import words_detect

def create_geometry_script (Phi_vector,Brain_map,electrode_profile,Xt,Yt,Zt,X_2nd,Y_2nd,Z_2nd,OZ_angle,Xt2,Yt2,Zt2,OX_angle2,OY_angle2,OZ_angle2,encap_thickness,ROI_radial,Xm,Ym,Zm,Vertice_enable,Lead2nd_Enable):
    
   #electrode_position=electrode_profile
   electrode_profile=electrode_profile	 
   # check parameters inputs
   name_idx=len(electrode_profile)
   check_profile_name = words_detect ('_profile.py', "Electrode_files/"+electrode_profile);
   if ( check_profile_name[0] == False):
       print ("ERROR: DBS lead profile name should be a string ends with _profile.py")
   else:
       f3=open("Electrode_files/"+electrode_profile,'r')                        
       f2=open(electrode_profile[:name_idx-11] + '_position.py','w+') # new file with new position
       #print(electrode_profile[:name_idx-11] + '_position.py')
       for index,line in enumerate(f3):  
               line_replace = False;
               var_list = words_detect("##### VARIABLE LIST #####",line)
               if (var_list[0]):
                   line_replace =True; #replace the code
                   f2.write("##### VARIABLE LIST #####\n"
                       +'Lead2nd_Enable = {}\n'.format(Lead2nd_Enable)
                       +'Xm = {}\n'.format(Xm)
                       +'Ym = {}\n'.format(Ym)
                       +'Zm = {}\n'.format(Zm)
                       +'Xt = {}\n'.format(Xt)
                       +'Yt = {}\n'.format(Yt)
                       +'Zt = {}\n'.format(Zt)
                       +'X_2nd = {}\n'.format(X_2nd)
                       +'Y_2nd = {}\n'.format(Y_2nd)
                       +'Z_2nd = {}\n'.format(Z_2nd)
                       +'OZ_angle = {}\n'.format(OZ_angle)
                       +'encap_thickness = {}\n'.format(encap_thickness)
                       +'ROI_radial = {}\n'.format(ROI_radial)
                       +'Vertice_enable = {}\n'.format(Vertice_enable)
                       +"Brain_map = '{}'\n".format(Brain_map)
		       +"Phi_vector = {}\n".format(Phi_vector)
                       +'if(Lead2nd_Enable):\n'
                       +'   Xt2 = {}\n'.format(Xt2)
                       +'   Yt2 = {}\n'.format(Yt2)
                       +'   Zt2 = {}\n'.format(Zt2)
                       +'   OX_angle2 = {}\n'.format(OX_angle2)
                       +'   OY_angle2 = {}\n'.format(OY_angle2)
                       +'   OZ_angle2 = {}\n'.format(OZ_angle2)
                         );
            
               if(line_replace == False): f2.write(line);      
       f2.close();
       f3.close();

###############################################################################
