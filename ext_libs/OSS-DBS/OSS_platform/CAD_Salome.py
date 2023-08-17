#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 26 20:54:56 2019
@author: Konstantin Butenko
Routines for geometry generation in SALOME 8.3.0 (for later versions, .med export has to be changed)
We do not define a geometry in SALOME, but only adjust the existing (see scripts in Eletrodes_files/)
Therefore, the commands here are not flexible.

class Mesh_ind contains indices of mesh entities in .xml files (taken from .msh)
"""

import os
import subprocess
import pickle
import time as time_lib
import logging

class Mesh_ind:
    def __init__(self,Tis_indx,ROI_indx,Contact_indx,Rest_indx,Flt_cnt_indx,Encup_indx,Contacts_indx,pulse_amp_vector,Active_contacts_on_lead=-1,Float_contacts_on_lead=-1):
        self.Tis_index = Tis_indx             # list: stores indices of the rest of the tissue and the outer encap. layey from Meshes/Mesh_unref.msh
        self.ROI_index = ROI_indx             # stores index of the Region of Interest from Meshes/Mesh_unref.msh
        self.Contact_index = Contact_indx         # stores index of the inner encap layer (in the vicinity of the electrode contacts)
        self.Rest_index = Rest_indx               # stores index of the Rest of the tissue from Meshes/Mesh_unref.msh
        self.Float_contacts = Flt_cnt_indx        # stores index of unactive contacts (uninsulated metals in general) from Meshes/Mesh_unref.msh. List if multicontact current-controlled.
        self.Encup_index = Encup_indx      # list:  # indices for inner and out encap. layer
        self.Active_contacts = Contacts_indx        # list: # indices of the active contacts
        self.Amp_vector = pulse_amp_vector                 # list:  # amplitudes (negative for cathods) assign to the active contacts (corresponds to self.Contacts)
        self.Active_on_lead = Active_contacts_on_lead  # list    #needed for multicontact current-controlled. Stores simply number of the active contacts in a vector. E.g. if 3 active, then [1,2,3]
        self.Float_on_lead = Float_contacts_on_lead # list    #needed for multicontact current-controlled. Stores simply number of the floating contacts + contacts with assigned currents in a vector. E.g. if 3 floating, 1 ground and 2 active, then [1,2,3,4,5].

def kill_SALOME_port():         #to ensure that Salome processes were terminated
    """ ensures that SALOME is terminated """

    import os
    import subprocess
    direct = os.environ['PATIENTDIR']
    port_file = open(direct+'/salomePort.txt','r')
    killPort = int(port_file.readline())
    port_file.close()
    #Kill the session with the specified port:
    subprocess.call('/opt/SALOME-9.7.0-UB16.04-SRC/BINARIES-UB16.04/KERNEL/bin/salome/killSalomeWithPort.py %s' % killPort,shell=True)
    # if you run without docker, this command should be changing according to the setup on your machine
    return True

def check_approx_dimensions(approx_dimensions, approx_geom_center, ROI_radius):

    """ Checks if brain approx. dimensions encompass all neurons (defined by ROI_radius).
        If not increases the dimensions """

    needs_a_rebuilt = 0
    for i in range(3):
        if ROI_radius > (approx_dimensions[i] / 2):
            approx_dimensions[i] = ROI_radius * 2 + 0.1
            logging.critical("increasing brain approx. along axis {} to encompass the neuron array\n".format(i))
            needs_a_rebuilt = 1

    if needs_a_rebuilt == 1:
        logging.critical("Increasing the dimensions of brain approximation\n")
        approx_dimensions = build_brain_approx(approx_dimensions, approx_geom_center)  # returns the same, but creates enlarged 'Brain_substitute.brep'

    # senseless check
    if ROI_radius > min((approx_dimensions[0] / 2), (approx_dimensions[1] / 2),
                                (approx_dimensions[2] / 2)):
        logging.critical("ROI_radius: ", ROI_radius)
        logging.critical("ROI is still bigger than the computational domain, check settings")
        raise SystemExit

    return True


def build_brain_approx(approx_dimensions, approx_geom_center, MRI_param = 0):

    """ Builds brain approximation
        if approx_geom_center == 0 or approx_dimensions[0]=0 (i.e. default), they will be estimated
        from the provided MRI dimensions """

    #  define dimensions of elliptic brain approximation
    if approx_dimensions[0] == 0:    #  build box or ellipsoid using MRI dimensions (starting in 0,0,0)
        x_length = abs(MRI_param.last_vox_coords[0] - MRI_param.first_vox_coords[0]) + MRI_param.voxel_dims[0]
        y_length = abs(MRI_param.last_vox_coords[1] - MRI_param.first_vox_coords[1]) + MRI_param.voxel_dims[1]
        z_length = abs(MRI_param.last_vox_coords[2] - MRI_param.first_vox_coords[2]) + MRI_param.voxel_dims[2]
    else:                   #build box or ellipsoid using given dimensions
        x_length, y_length, z_length = approx_dimensions[:]

    #  define center coordinates
    if approx_geom_center == 0: #  Centering approximation on the MRI data
        geom_center_x = (MRI_param.last_vox_coords[0] + MRI_param.first_vox_coords[0]) / 2
        geom_center_y = (MRI_param.last_vox_coords[1] + MRI_param.first_vox_coords[1]) / 2
        geom_center_z = (MRI_param.last_vox_coords[2] + MRI_param.first_vox_coords[2]) / 2
    else:   #centering on the given coordinates
        geom_center_x, geom_center_y, geom_center_z = (approx_geom_center[:])          #this will shift only the approximating geometry, not the MRI data set!

    # copy a template Salome file to the stim folder
    from shutil import copy2
    oss_dbs_folder = '/opt/OSS-DBS/OSS_platform'
    copy2(oss_dbs_folder +'/Brain_substitute.py', os.environ['PATIENTDIR'])
    brain_substitute_path = os.environ['PATIENTDIR'] + "/Brain_substitute.py"

    #  directly inserts parameters to Brain_substitute.py
    from Parameter_insertion import paste_geom_dim
    paste_geom_dim(x_length, y_length, z_length, geom_center_x, geom_center_y, geom_center_z)
    direct = os.environ['PATIENTDIR']
    logging.critical("----- Creating brain approximation in SALOME -----")

    #  run Brain_substitute.py using Salome
    with open(os.devnull, 'w') as FNULL: subprocess.call('salome -t python3 '+ brain_substitute_path +' --ns-port-log='+direct+'/salomePort.txt', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    #kill_SALOME_port()
    logging.critical("Brain_substitute.brep was created\n")

    #  run mesh converters (med -> msh2 -> msh -> xml)
    #  the last conversion should create 3 files (Mesh_unref.xml, ..._facet_region.xml, _physical_region.xml)
    with open(os.devnull, 'w') as FNULL: subprocess.call('gmsh ' + os.environ['PATIENTDIR']+'/Meshes/Mesh_brain_substitute_max_ROI.med -3 -v 0 -o ' + os.environ['PATIENTDIR']+'/Meshes/Mesh_brain_substitute_max_ROI.msh2 && mv ' + os.environ['PATIENTDIR']+'/Meshes/Mesh_brain_substitute_max_ROI.msh2 ' + os.environ['PATIENTDIR']+'/Meshes/Mesh_brain_substitute_max_ROI.msh',shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    with open(os.devnull, 'w') as FNULL: subprocess.call('dolfin-convert ' + os.environ['PATIENTDIR']+'/Meshes/Mesh_brain_substitute_max_ROI.msh ' + os.environ['PATIENTDIR']+'/Meshes/Mesh_brain_substitute_max_ROI.xml',shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    #  make sure Salome processes are terminated
    pro = subprocess.Popen(['salome','killall'])  
        
    # the dimensions might need to be adjusted later
    return [x_length, y_length, z_length]

def build_final_geometry(d, MRI_param, ROI_radius, cc_multicontact):

    """ Builds the final geometry (brain, electrode, encapsulation)
        if cc_multicontact = True, builds active contacts as conductive volumes """

    start_final_geom = time_lib.time()

    Brain_link = os.environ['PATIENTDIR']+ '/' +str(d["Brain_shape_name"])

    if cc_multicontact == True:      #  here also creates floating volumes for active contacts (with assigned currents)
        Electrode_profile = d["Electrode_type"] + '_floating_profile.py'
        position_script_name = os.environ['PATIENTDIR'] + "/" + d["Electrode_type"] + "_floating_position.py"
    else:
        Electrode_profile = d["Electrode_type"]+'_profile.py'
        position_script_name=os.environ['PATIENTDIR']+"/"+d["Electrode_type"]+"_position.py"

    # correction for human electrodes. Otherwise, SALOME might fail to find encap_outer
    if ROI_radius < 6.5 and d['rodent_electrode'] == False:
        ROI_radius = 6.5

    impl_coords = [d["Implantation_coordinate_X"], d["Implantation_coordinate_Y"], d["Implantation_coordinate_Z"]]
    second_coords = d["Second_coordinate_X"], d["Second_coordinate_Y"], d["Second_coordinate_Z"]  # second point on the lead to define a line
    Shift_PO = MRI_param.MRI_shift  # shift all coordinates so that MRI starts at 0,0,0 (Positive Octant)

    # copies the electrode profile file (from Electrode_files/) and inserts there the specified geom. parameters
    from Electrode_files.DBS_lead_position_V10 import create_geometry_script
    create_geometry_script(d, Brain_link, Electrode_profile, impl_coords, second_coords, ROI_radius, Shift_PO, False)

    direct = os.environ['PATIENTDIR'] # there was a glitch with paths, so I assign it to a variable

    logging.critical("\n----- Creating final geometry in SALOME -----")
    # Make a subprocess call to the salome executable and store the used port in a text file:
    with open(os.devnull, 'w') as FNULL: subprocess.call('salome -t python3 '+ position_script_name +' --ns-port-log='+direct+'/salomePort.txt', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    #kill_SALOME_port()  # terminate SALOME

    #  run mesh converters (med -> msh2 -> msh)
    with open(os.devnull, 'w') as FNULL: subprocess.call('gmsh ' + os.environ['PATIENTDIR']+'/Meshes/Mesh_unref.med -3 -v 0 -o ' + os.environ['PATIENTDIR']+'/Meshes/Mesh_unref.msh2 && mv ' + os.environ['PATIENTDIR']+'/Meshes/Mesh_unref.msh2 ' + os.environ['PATIENTDIR']+'/Meshes/Mesh_unref.msh',shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    # now we don't need None values, Contacts point to the active only
    Pulse_amp_active = [x for x in d["Pulse_amp"] if x is not None]

    if cc_multicontact == True:               #if multicontact with current-controlled, then we need addit.entries

        dict_ind = {             # I have separate dictionaries for clarity
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

        from MeshTransfer import read_mesh_indicies_extended
        read_mesh_indicies_extended(dict_ind)       # to get entities' indices from Meshes/Mesh_unref.msh

        # store in class
        Domains = Mesh_ind(dict_ind["Tis_ind"],dict_ind["ROI_ind"],dict_ind["Contact_ind"],dict_ind["Rest_ind"],dict_ind["Flt_contacts"],dict_ind["Encup_ind"],dict_ind["Contacts"],Pulse_amp_active,dict_ind["Active_contacts_on_lead"],dict_ind["Float_contacts_on_lead"])
        Domains.Active_on_lead.sort()
    else:

        dict_ind = {
        'Tis_ind'                    :  [0,0] ,
        'ROI_ind'                    :  0 ,
        'Contact_ind'                :  0,       #subdomain of encap. layer where contacts are
        'Rest_ind'                   :  0 ,
        'Flt_cnt'                    :  -1 ,     #only one value, because they are collected in one mesh group. -1 to avoid confusion
        'Encup_ind'                  :  [0,0] ,
        'Contacts'                   :  [0,0,0,0],
        }

        from MeshTransfer import read_mesh_indicies
        read_mesh_indicies(dict_ind)            # to get entities' indices from Meshes/Mesh_unref.msh

        # store in class
        Domains = Mesh_ind(dict_ind["Tis_ind"],dict_ind["ROI_ind"],dict_ind["Contact_ind"],dict_ind["Rest_ind"],dict_ind["Flt_cnt"],dict_ind["Encup_ind"],dict_ind["Contacts"],Pulse_amp_active)

    if Domains.Tis_index == -1:
        logging.critical("ROI is the whole computational domain! Employing a bigger geometrical domain is necessary")
        raise SystemExit

    # save the class instance
    with open(os.environ['PATIENTDIR']+'/Meshes/Mesh_ind.file', "wb") as f:
        pickle.dump(Domains, f, pickle.HIGHEST_PROTOCOL)

    # run mesh converters (msh -> xml)
    with open(os.devnull, 'w') as FNULL: subprocess.call('dolfin-convert ' + os.environ['PATIENTDIR']+'/Meshes/Mesh_unref.msh ' + os.environ['PATIENTDIR']+'/Meshes/Mesh_unref.xml',shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    x_imp = MRI_param.MRI_shift[0]+d["Implantation_coordinate_X"]
    y_imp = MRI_param.MRI_shift[1]+d["Implantation_coordinate_Y"]
    z_imp = MRI_param.MRI_shift[2]+d["Implantation_coordinate_Z"]
    #print("Loading clipping")
    #from Parameter_insertion import paste_paraview_clipping
    #paste_paraview_clipping(x_imp,y_imp,z_imp)      #clipping on this point in x-direction
    logging.critical("Coordinates of the implantation in the positive octant: {} {} {}".format(x_imp,y_imp,z_imp))

    minutes = int((time_lib.time() - start_final_geom) / 60)
    secnds = int(time_lib.time() - start_final_geom) - minutes * 60
    logging.critical("----- Final geometry was created and meshed in {} min {} sec, the files are stored in Meshes/ -----\n".format(minutes, secnds))

    return Domains

