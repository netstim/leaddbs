#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 26 20:54:56 2019

@author: butenko
"""

import os
import subprocess
import pickle
import numpy as np
import time as time_lib

class Mesh_ind:
    def __init__(self,Tis_indx,ROI_indx,Contact_indx,Rest_indx,Flt_cnt_indx,Encup_indx,Contacts_indx,fi_vector,Active_contacts_on_lead=-1,Float_contacts_on_lead=-1):
        self.Tis_index=Tis_indx             #list: stores indices of the rest of the tissue and the outer encap. layey from Meshes/Mesh_unref.msh
        self.ROI_index=ROI_indx             # stores index of the Region of Interest from Meshes/Mesh_unref.msh
        self.Contact_index=Contact_indx         # stores index of the inner encap layer (in the vicinity of the electrode contacts)
        self.Rest_index=Rest_indx               # stores index of the Rest of the tissue from Meshes/Mesh_unref.msh
        self.Float_contacts=Flt_cnt_indx        # stores index of unactive contacts (uninsulated metals in general) from Meshes/Mesh_unref.msh. List if multicontact current-controlled.
        self.Encup_index=Encup_indx      #list:  # indices for inner and out encap. layer
        self.Contacts=Contacts_indx        #list: # indices of the active contacts 
        self.fi=fi_vector                 #list:  # amplitudes assign to the active contacts (corresponds to self.Contacts)
        self.Active_on_lead=Active_contacts_on_lead  #list    #needed for multicontact current-controlled. Stores simply number of the active contacts in a vector. E.g. if 3 active, then [1,2,3]
        self.Float_on_lead=Float_contacts_on_lead #list    #needed for multicontact current-controlled. Stores simply number of the floating contacts + contacts with assigned currents in a vector. E.g. if 3 floating, 1 ground and 2 active, then [1,2,3,4,5].

def kill_SALOME_port():         #to ensure that Salome processes were terminated
    import os
    import subprocess
    direct = os.getcwd()
    port_file = open(direct+'/salomePort.txt','r')
    killPort = int(port_file.readline())
    port_file.close()
    #Kill the session with the specified port:
    subprocess.call('/opt/SALOME-8.3.0-UB16.04/BINARIES-UB16.04/KERNEL/bin/salome/killSalomeWithPort.py %s' % killPort,shell=True)   # if you run without docker, this command should be changing according to the setup on your machine
    return True    
       
def build_brain_approx(d,MRI_param):
    
    if d['Approximating_Dimensions'][0]==0:    #build box or ellipsoid using MRI dimensions (starting in 0,0,0)
        x_length=abs(MRI_param.x_max-MRI_param.x_min)+MRI_param.x_vox_size
        y_length=abs(MRI_param.y_max-MRI_param.y_min)+MRI_param.y_vox_size
        z_length=abs(MRI_param.z_max-MRI_param.z_min)+MRI_param.z_vox_size
    else:                   #build box or ellipsoid using given dimensions
        x_length,y_length,z_length=d['Approximating_Dimensions'][:]
            
    if d["Aprox_geometry_center"]==0: #"Centering approximation on the MRI data"           
        Geom_center_x=(MRI_param.x_max+MRI_param.x_min)/2
        Geom_center_y=(MRI_param.y_max+MRI_param.y_min)/2
        Geom_center_z=(MRI_param.z_max+MRI_param.z_min)/2
    else:   #centering on the given coordinates
        Geom_center_x,Geom_center_y,Geom_center_z=(d["Aprox_geometry_center"][0],d["Aprox_geometry_center"][1],d["Aprox_geometry_center"][2])          #this will shift only the approximating geometry, not the MRI data set!
    
    from Parameter_insertion import paste_geom_dim
    paste_geom_dim(x_length,y_length,z_length,Geom_center_x,Geom_center_y,Geom_center_z)        #directly inserts parameters to Brain_substitute.py
    direct = os.getcwd()
    print("----- Creating brain approximation in SALOME -----")            
    with open(os.devnull, 'w') as FNULL: subprocess.call('salome -t python '+ 'Brain_substitute.py' +' --ns-port-log='+direct+'/salomePort.txt', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    kill_SALOME_port()

    print("Brain_substitute.brep was created\n")
    with open(os.devnull, 'w') as FNULL: subprocess.call('gmsh Meshes/Mesh_brain_substitute_max_ROI.med -3 -v 0 -o Meshes/Mesh_brain_substitute_max_ROI.msh2 && mv Meshes/Mesh_brain_substitute_max_ROI.msh2 Meshes/Mesh_brain_substitute_max_ROI.msh',shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    with open(os.devnull, 'w') as FNULL: subprocess.call('dolfin-convert Meshes/Mesh_brain_substitute_max_ROI.msh Meshes/Mesh_brain_substitute_max_ROI.xml',shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    
    return x_length,y_length,z_length

def build_final_geometry(d,MRI_param,Brain_shape_name,ROI_radius,cc_multicontact):
    
    start_final_geom=time_lib.time()
    
    from Electrode_files.DBS_lead_position_V10 import create_geometry_script
    Brain_link = str(Brain_shape_name)
    
    if cc_multicontact==True:           #here we will also create floating volumes for active contacts with assigned currents
        Electrode_profile=d["Electrode_type"]+'_floating_profile.py'
        position_script_name=d["Electrode_type"]+"_floating_position.py"
    else:
        Electrode_profile=d["Electrode_type"]+'_profile.py'
        position_script_name=d["Electrode_type"]+"_position.py"
    
    create_geometry_script(d["Phi_vector"],Brain_link,Electrode_profile,d["Implantation_coordinate_X"],d["Implantation_coordinate_Y"],d["Implantation_coordinate_Z"],d["Second_coordinate_X"],d["Second_coordinate_Y"],d["Second_coordinate_Z"],d["Rotation_Z"],0.0,0.0,0.0,0.0,0.0,0.0,d["encap_thickness"],ROI_radius,MRI_param.x_shift,MRI_param.y_shift,MRI_param.z_shift,False,False)
    
    direct = os.getcwd()
    # Make a subprocess call to the salome executable and store the used port in a text file:
    
    print("----- Creating final geometry in SALOME -----")

    with open(os.devnull, 'w') as FNULL: subprocess.call('salome -t python '+ position_script_name +' --ns-port-log='+direct+'/salomePort.txt', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    kill_SALOME_port()

    with open(os.devnull, 'w') as FNULL: subprocess.call('gmsh Meshes/Mesh_unref.med -3 -v 0 -o Meshes/Mesh_unref.msh2 && mv Meshes/Mesh_unref.msh2 Meshes/Mesh_unref.msh',shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    
    Phi_vector=[x for x in d["Phi_vector"] if x is not None] #now we don't need None values, Contacts point to the active ones
    
    if cc_multicontact==True:               #if multicontact with current-controlled, then we need another dictionary
        from MeshTransfer import read_mesh_indicies_extended
        dict_ind = {
        'Tis_ind'                    :  [0,0] ,
        'ROI_ind'                    :  0 ,
        'Contact_ind'                :  0,       #subdomain of encap. layer where contacts are
        'Rest_ind'                   :  0 ,
        'Flt_contacts'               :  -1 ,     #here it will become a list! contains indices of floatin in gmsh
        'Encup_ind'                  :  [0,0] ,
        'Contacts'                   :  [0,0,0,0],   #contains indices of contacts in gmsh
        'Active_contacts_on_lead'    : [0,0],    #contains indixes of the contacts on the lead that are active   (starting from 1)
        'Float_contacts_on_lead'     : [0,0]         #contains indices of the floats on the lead  (starting from 1)
        }
        
        read_mesh_indicies_extended(dict_ind)       # to get indices from Meshes/Mesh_unref.med
        Tis_ind,ROI_ind,Contact_ind,Rest_ind,Flt_cnt,Encup_ind,Contacts,Active_on_lead,Float_on_lead=(dict_ind["Tis_ind"],dict_ind["ROI_ind"],dict_ind["Contact_ind"],dict_ind["Rest_ind"],dict_ind["Flt_contacts"],dict_ind["Encup_ind"],dict_ind["Contacts"],dict_ind["Active_contacts_on_lead"],dict_ind["Float_contacts_on_lead"])
        Domains=Mesh_ind(Tis_ind,ROI_ind,Contact_ind,Rest_ind,Flt_cnt,Encup_ind,Contacts,Phi_vector,Active_on_lead,Float_on_lead)
    else:
        from MeshTransfer import read_mesh_indicies
        dict_ind = {
        'Tis_ind'                    :  [0,0] ,
        'ROI_ind'                    :  0 ,
        'Contact_ind'                :  0,       #subdomain of encap. layer where contacts are
        'Rest_ind'                   :  0 ,
        'Flt_cnt'                    :  -1 ,     #only one value, because they are collected in one mesh group. -1 to avoid confusion
        'Encup_ind'                  :  [0,0] ,
        'Contacts'                   :  [0,0,0,0],
        }
    
        read_mesh_indicies(dict_ind)            # to get indices from Meshes/Mesh_unref.med
        Tis_ind,ROI_ind,Contact_ind,Rest_ind,Flt_cnt,Encup_ind,Contacts=(dict_ind["Tis_ind"],dict_ind["ROI_ind"],dict_ind["Contact_ind"],dict_ind["Rest_ind"],dict_ind["Flt_cnt"],dict_ind["Encup_ind"],dict_ind["Contacts"])
        Domains=Mesh_ind(Tis_ind,ROI_ind,Contact_ind,Rest_ind,Flt_cnt,Encup_ind,Contacts,Phi_vector)
    
    if Domains.Tis_index==-1:
        print("ROI is the whole computational domain! Employing a bigger geometrical domain is necessary")
    
    with open('Meshes/Mesh_ind.file', "wb") as f:
        pickle.dump(Domains, f, pickle.HIGHEST_PROTOCOL)
        
    with open(os.devnull, 'w') as FNULL: subprocess.call('dolfin-convert Meshes/Mesh_unref.msh Meshes/Mesh_unref.xml',shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    minutes=int((time_lib.time() - start_final_geom)/60)
    secnds=int(time_lib.time() - start_final_geom)-minutes*60

    print("----- Final geometry was created and meshed in ",minutes," min ",secnds," s , the files are stored in Meshes/ -----")
    
    x_imp=MRI_param.x_shift+d["Implantation_coordinate_X"]
    y_imp=MRI_param.y_shift+d["Implantation_coordinate_Y"]
    z_imp=MRI_param.z_shift+d["Implantation_coordinate_Z"]    
    #print("Loading clipping")
    from Parameter_insertion import paste_paraview_clipping
    paste_paraview_clipping(x_imp,y_imp,z_imp)      #clipping on this point in x-direction    
    print("Coordinates of the electrode tip in the positive octant coordinates: ",x_imp,y_imp,z_imp,"\n")

    return Domains

