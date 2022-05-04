# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.3.0 with dump python functionality
###### Run with DPS_lead_position_V9.py

import sys
import salome

salome.salome_init()

###
### GEOM component
###
########################################### extra code 1 V10 15/12/18#############################################
###### This file runs with DBS_lead_position_V10.py
import os
sys.path.insert( 0, r'{}'.format(os.getcwd()))
sys.path.append('/usr/local/lib/python2.7/dist-packages')
#from pandas import read_csv

##### DEFAULT LIST #####

#Lead2nd_Enable = True
#Xt = 0
#Yt = 5
#Zt = 0
#X_2nd = 0
#Y_2nd = 5
#Z_2nd = 0
#OZ_angle = 0
#Xm = 0
#Ym = 0
#Zm = 0
#encap_thickness = 0.1
#ROI_radial = 13
#Vertice_enable = False
#Brain_map = '/home/trieu/electrode_dir/brain_elipse.brep'
#if(Lead2nd_Enable):
#   Xt2 = 0
#   Yt2 = -5
#   Zt2 = 0
#   OX_angle2 = 0
#   OY_angle2 = 0
#   OZ_angle2 = 0

##### VARIABLE LIST #####


########## End of variable list#############


if Z_2nd == Zt:
    Z_2nd_artif = Zt+1.0 # just to ensure the rotation is possible
else:
    Z_2nd_artif=Z_2nd

#for Lead-DBS, the tip point should be shifted down (they use the middle of the lowest contact as the reference point)
Zt_tip=Zt-2.25		#for Medtronic3389

Vert_array =[0];
number_vertex = len(Vert_array)
Vert = []
VolumeObject1 = []
ContactObject1 = []
VolumeObject2 = []
ContactObject2 = []
print( " DBS_lead's Geometry buid\n")
######################################### end of extra code 1 ########################################
######################################################################################################
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Circle_1 = geompy.MakeCircle(O, OZ, 0.635)
Contact_1 = geompy.MakePrismVecH(Circle_1, OZ, 0.75*stretch+0.75)
geompy.TranslateDXDYDZ(Contact_1, 0, 0, 0.865)
Contact_1_full = geompy.MakePrismVecH(Circle_1, OZ, 1.5*stretch)
geompy.TranslateDXDYDZ(Contact_1_full, 0, 0, 0.865)

Contact_2 = geompy.MakeTranslation(Contact_1_full, 0, 0, 1.25*stretch+0.75)
Contact_3 = geompy.MakeTranslation(Contact_1_full, 0, 0, 3.25*stretch+0.75)
Contact_4 = geompy.MakeTranslation(Contact_1, 0, 0, 5.25*stretch+0.75)
Cylinder_1 = geompy.MakeCylinderRH(0.635, 149.365)
Sphere_1 = geompy.MakeSphereR(0.635)
Fuse_1 = geompy.MakeFuseList([Cylinder_1, Sphere_1], True, True)
Cylinder_2 = geompy.MakeCylinderRH(encap_thickness+0.635, 149.365)
Sphere_2 = geompy.MakeSphereR(encap_thickness+0.635)
Fuse_2 = geompy.MakeFuseList([Cylinder_2, Sphere_2], True, True)
encap_layer = geompy.MakeCutList(Fuse_2, [Fuse_1], True)
geompy.TranslateDXDYDZ(Circle_1, 0, 0, 0.635)
geompy.TranslateDXDYDZ(Contact_1, 0, 0, 0.635)
geompy.TranslateDXDYDZ(Contact_2, 0, 0, 0.635)
geompy.TranslateDXDYDZ(Contact_3, 0, 0, 0.635)
geompy.TranslateDXDYDZ(Contact_4, 0, 0, 0.635)
geompy.TranslateDXDYDZ(Cylinder_1, 0, 0, 0.635)
geompy.TranslateDXDYDZ(Sphere_1, 0, 0, 0.635)
geompy.TranslateDXDYDZ(Fuse_1, 0, 0, 0.635)
geompy.TranslateDXDYDZ(Cylinder_2, 0, 0, 0.635)
geompy.TranslateDXDYDZ(Sphere_2, 0, 0, 0.635)
geompy.TranslateDXDYDZ(Fuse_2, 0, 0, 0.635)
geompy.TranslateDXDYDZ(encap_layer, 0, 0, 0.635)
Sphere_ROI = geompy.MakeSphereR(ROI_radial)
encap_outer_ROI = geompy.MakeCutList(encap_layer, [Sphere_ROI], True)
encap_inner_ROI = geompy.MakeCutList(encap_layer, [encap_outer_ROI], True)
Fuse_all_lead_encap_ROI = geompy.MakeFuseList([Sphere_ROI, Fuse_2], True, True)
ROI = geompy.MakeCutList(Sphere_ROI, [Fuse_2], True)
CV1 = geompy.MakeCylinderRH(0.635, 0.75*stretch+0.75)
geompy.TranslateDXDYDZ(CV1, 0, 0, 1.5)

CV1_full = geompy.MakeCylinderRH(0.635, 1.5*stretch)
geompy.TranslateDXDYDZ(CV1_full, 0, 0, 1.5)

CV2 = geompy.MakeTranslation(CV1_full, 0, 0, 0.75+1.25*stretch)
CV3 = geompy.MakeTranslation(CV1_full, 0, 0, 3.25*stretch+0.75)
CV4 = geompy.MakeTranslation(CV1, 0, 0, 5.25*stretch+0.75)
##################################################################################################################
########################################### extra code 2 V10 15/12/18#############################################
print( " Load brain image \n")
if (Brain_map[-4:] == 'brep'):
    brain_solid = geompy.ImportBREP( Brain_map )
elif (Brain_map[-4:] == 'step'):
    brain_solid = geompy.ImportSTEP( Brain_map )
elif (Brain_map[-4:] == 'iges'):
    brain_solid = geompy.ImportIGES( Brain_map )
elif (Brain_map[-4:] == '.stl'):
    brain_solid = geompy.ImportSTL( Brain_map )
else:
    print( " unknow imported file format")
Fuse_all_lead_encap_ROI_no_internal_face = geompy.RemoveInternalFaces(Fuse_all_lead_encap_ROI)

#################################################### Geometry and extra code interface ##############################################################
VolumeObject1 = [ encap_outer_ROI,ROI,encap_inner_ROI,CV1,CV2,CV3,CV4]         # Declare objects included to partition, encap_outer_ROI always @1st position
Volume_name1  = ['encap_outer_ROI1','ROI1','encap_inner_ROI1','CV1_1','CV1_2','CV1_3','CV1_4'] # Declare name of the group in the partition for volume
ContactObject1 = [Contact_1,Contact_2,Contact_3,Contact_4]
Contact_name1 = ['Contact1_1','Contact1_2','Contact1_3','Contact1_4']

if(Lead2nd_Enable): ##################  2nd LEAD ###############################################
  VolumeObject2 = [ROI]*len(VolumeObject1)
  ContactObject2 = [Contact_1]*len(ContactObject1)
  Volume_name2  = [ 'encap_outer_ROI2','ROI2','encap_inner_ROI2','CV2_1','CV2_2','CV2_3','CV2_4']
  Contact_name2 = ['Contact2_1','Contact2_2','Contact2_3','Contact2_4']
##############################################################################################################################################
  print( "Position 2nd Fuse all object at [{},{},{}], [{}',{}',{}']\n".format(Xt2,Yt2,Zt2,OX_angle2,OY_angle2,OZ_angle2))
  Fuse_all_lead_encap_ROI_no_internal_face2 = geompy.MakeTranslation(Fuse_all_lead_encap_ROI_no_internal_face,Xt2,Yt2,Zt2)
  OX2 = geompy.MakeTranslation(OX,Xt2,Yt2,Zt2)
  OY2 = geompy.MakeTranslation(OY,Xt2,Yt2,Zt2)
  OZ2 = geompy.MakeTranslation(OZ,Xt2,Yt2,Zt2)
  geompy.Rotate(Fuse_all_lead_encap_ROI_no_internal_face2, OX2,OX_angle2*math.pi/180.0)
  geompy.Rotate(Fuse_all_lead_encap_ROI_no_internal_face2, OY2,OY_angle2*math.pi/180.0)
  geompy.Rotate(Fuse_all_lead_encap_ROI_no_internal_face2, OZ2,OZ_angle2*math.pi/180.0)
  print( "Position 2nd Lead at [{},{},{}], [{}',{}',{}']\n".format(Xt2,Yt2,Zt2,OX_angle2,OY_angle2,OZ_angle2))
  for i in range(0,len(VolumeObject1)):
    VolumeObject2[i] = geompy.MakeTranslation(VolumeObject1[i],Xt2,Yt2,Zt2)
    geompy.Rotate(VolumeObject2[i], OX2,OX_angle2*math.pi/180.0)
    geompy.Rotate(VolumeObject2[i], OY2,OY_angle2*math.pi/180.0)
    geompy.Rotate(VolumeObject2[i], OZ2,OZ_angle2*math.pi/180.0)
  for i in range(0,len(ContactObject1)):
    ContactObject2[i] = geompy.MakeTranslation(ContactObject1[i],Xt2,Yt2,Zt2)
    geompy.Rotate(ContactObject2[i], OX2,OX_angle2*math.pi/180.0)
    geompy.Rotate(ContactObject2[i], OY2,OY_angle2*math.pi/180.0)
    geompy.Rotate(ContactObject2[i], OZ2,OZ_angle2*math.pi/180.0)
  print( "Cut outer ROI2 with brain\n")
  cut_outer_ROI = geompy.MakeCutList(VolumeObject2[0], [brain_solid], True)
  VolumeObject2[0] = geompy.MakeCutList(VolumeObject2[0], [cut_outer_ROI], True)
  print( "Cut ROI2 with brain\n")
  VolumeObject2[1] = geompy.MakeCommonList([VolumeObject2[1], brain_solid], True)
  print( "Group 2nd:volume and area extraction for group ID identification process\n")
  Volume2_Pro = [geompy.BasicProperties( VolumeObject2[0])]*len(VolumeObject2)
  Contact2_Pro = [geompy.BasicProperties( ContactObject2[0])]*len(ContactObject2)
  for i in range(0,len(VolumeObject2)):
    Volume2_Pro[i] = geompy.BasicProperties( VolumeObject2[i])
  for i in range(0,len(ContactObject2)):
    Contact2_Pro[i] = geompy.BasicProperties( ContactObject2[i])

################## LEAD 1st #############################################################
#print( "Position 1st Fuse all object at [{},{},{}], [{}',{}',{}']\n".format(Xt,Yt,Zt,OX_angle,OY_angle,OZ_angle))
geompy.TranslateDXDYDZ(Fuse_all_lead_encap_ROI_no_internal_face,Xt,Yt,Zt_tip)

OX1 = geompy.MakeTranslation(OX,Xt,Yt,Zt_tip)
OY1 = geompy.MakeTranslation(OY,Xt,Yt,Zt_tip)
OZ1 = geompy.MakeTranslation(OZ,Xt,Yt,Zt_tip)

geompy.Rotate(Fuse_all_lead_encap_ROI_no_internal_face, OZ1,OZ_angle*math.pi/180.0)

Vertex_1 = geompy.MakeVertex(X_2nd,Y_2nd,Z_2nd)
Vertex_O = geompy.MakeVertex(Xt,Yt,Zt)
Vertex_3 = geompy.MakeVertex(Xt,Yt,Z_2nd_artif)

if X_2nd!=Xt or Y_2nd!=Yt:
        Fuse_all_lead_encap_ROI_no_internal_face=geompy.MakeRotationThreePoints(Fuse_all_lead_encap_ROI_no_internal_face, Vertex_O, Vertex_3, Vertex_1)


#print( "Position 1st Lead at [{},{},{}], [{}',{}',{}']\n".format(Xt,Yt,Zt,OX_angle,OY_angle,OZ_angle))
for i in range(0,len(VolumeObject1)):
    geompy.TranslateDXDYDZ(VolumeObject1[i],Xt,Yt,Zt_tip)
    geompy.Rotate(VolumeObject1[i], OZ1,OZ_angle*math.pi/180.0)
    if X_2nd!=Xt or Y_2nd!=Yt:
        VolumeObject1[i]=geompy.MakeRotationThreePoints(VolumeObject1[i], Vertex_O, Vertex_3, Vertex_1)

for i in range(0,len(ContactObject1)):
    geompy.TranslateDXDYDZ(ContactObject1[i],Xt,Yt,Zt_tip)
    geompy.Rotate(ContactObject1[i], OZ1,OZ_angle*math.pi/180.0)
    if X_2nd!=Xt or Y_2nd!=Yt:
        ContactObject1[i]=geompy.MakeRotationThreePoints(ContactObject1[i], Vertex_O, Vertex_3, Vertex_1)






print( "Cut outer ROI1 with brain\n")
cut_outer_ROI = geompy.MakeCutList(VolumeObject1[0], [brain_solid], True)
VolumeObject1[0] = geompy.MakeCutList(VolumeObject1[0], [cut_outer_ROI], True)
print( "Cut ROI1 with brain\n")
VolumeObject1[1] = geompy.MakeCommonList([VolumeObject1[1], brain_solid], True)
print( "Group 1st:volume and area extraction for group ID identification process\n")

Volume1_Pro = [geompy.BasicProperties( VolumeObject1[0])]*len(VolumeObject1)
Contact1_Pro = [geompy.BasicProperties( ContactObject1[0])]*len(ContactObject1)
for i in range(0,len(VolumeObject1)):
    Volume1_Pro[i] = geompy.BasicProperties( VolumeObject1[i])
for i in range(0,len(ContactObject1)):
    Contact1_Pro[i] = geompy.BasicProperties( ContactObject1[i])


print( "Create reference groups for ID identification process\n")

if(Lead2nd_Enable):

  Rest = geompy.MakeCutList(brain_solid, [Fuse_all_lead_encap_ROI_no_internal_face,Fuse_all_lead_encap_ROI_no_internal_face2], True)
  Partition_profile = geompy.MakePartition(VolumeObject1+VolumeObject2+ContactObject1+ContactObject2, [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)

###reference_volume
  reference_volume = VolumeObject1 + VolumeObject2
  reference_volume_Pro = Volume1_Pro + Volume2_Pro
  Volume_name = Volume_name1+Volume_name2
### reference_area
  reference_surface = ContactObject1 + ContactObject2
  reference_surface_Pro = Contact1_Pro + Contact2_Pro
  Contact_name = Contact_name1+Contact_name2
  Group_volume  = [geompy.CreateGroup(Partition_profile, geompy.ShapeType["SOLID"])] * (len(VolumeObject1)+len(VolumeObject2)+1) # +1 is Rest Group
  Group_surface = [geompy.CreateGroup(Partition_profile, geompy.ShapeType["FACE"])] * (len(ContactObject1)+len(ContactObject2))
else:

  Rest = geompy.MakeCutList(brain_solid, [Fuse_all_lead_encap_ROI_no_internal_face], True)
  Partition_profile = geompy.MakePartition(VolumeObject1+ContactObject1, [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
###reference_volume
  reference_volume = VolumeObject1
  reference_volume_Pro = Volume1_Pro
  Volume_name = Volume_name1
### reference_area
  reference_surface = ContactObject1
  reference_surface_Pro = Contact1_Pro
  Contact_name = Contact_name1
  Group_volume  = [geompy.CreateGroup(Partition_profile, geompy.ShapeType["SOLID"])] * (len(VolumeObject1)+1) # +1 is Rest Group
  Group_surface = [geompy.CreateGroup(Partition_profile, geompy.ShapeType["FACE"])] * len(ContactObject1)
### find out subshape and subshape ID

Group_surface_ListIDs =[]
Group_volume_ListIDs =[]
Group_partition_volume = []
Group_partition_surface = []

### find group volume ID ######################################################################
Partition_volume_IDsList = geompy.SubShapeAllIDs(Partition_profile, geompy.ShapeType["SOLID"]) # list all sub shape volume in Partition
print( "Partition_volume_IDsList",Partition_volume_IDsList, '\n')

for ref_ind in range (0, len(reference_volume)):
    temp_volume = []
    for sub_ind in range (0, len (Partition_volume_IDsList)):
        subshape = geompy.GetSubShape(Partition_profile, [Partition_volume_IDsList[sub_ind]]) # get subshape
        subshape_Pro = geompy.BasicProperties(subshape)       # extract volume of subshape
        Common_volume = geompy.MakeCommonList([subshape, reference_volume[ref_ind]], True) # check common intersection
        Common_volume_Pro = geompy.BasicProperties(Common_volume)
        print( "volume difference",abs(Common_volume_Pro[2]-subshape_Pro[2]),"/",abs(Common_volume_Pro[2]-reference_volume_Pro[ref_ind][2]))
        # if ( common volume = subshape) and (common volume = ref volume) => ref volume = sub shape
        if (abs(Common_volume_Pro[2]-subshape_Pro[2])< 0.0003) and (abs(Common_volume_Pro[2]-reference_volume_Pro[ref_ind][2])<0.0003):

            Group_partition_volume.append([Volume_name[ref_ind],Partition_volume_IDsList[sub_ind]])
        # if ( common volume = subshape) and (common volume < ref volume) => sub shape belong to ref volume
        elif (abs(Common_volume_Pro[2]-subshape_Pro[2])< 0.0003) and ((Common_volume_Pro[2] - reference_volume_Pro[ref_ind][2])<-0.0003):
            temp_volume.append( Partition_volume_IDsList[sub_ind] )
    if len(temp_volume) >1 : # the volume is devided
        Group_partition_volume.append([Volume_name[ref_ind],temp_volume ])
        print( Volume_name[ref_ind]," is devided and has sub IDs:{}\n".format(temp_volume))
if len(reference_volume) != len(Group_partition_volume):
    print( "Geometry-volume error please check ROI diameter and DBS lead Position ",len(reference_volume),len(Group_partition_volume))
print( 'Group_partition_volume',Group_partition_volume,'\n')

### find group surface ID ######################################################################
Partition_surface_IDsList = geompy.SubShapeAllIDs(Partition_profile, geompy.ShapeType["FACE"]) # list all sub shape face in Partition
print( 'Partition_surface_IDsList',Partition_surface_IDsList,'\n')
sub_face = [] ## store devided faces
for reff_ind in range (0, len (reference_surface)):
    temp_surface = []
    for subf_ind in range (0, len(Partition_surface_IDsList)):
        subshapef = geompy.GetSubShape(Partition_profile, [Partition_surface_IDsList[subf_ind]]) # get subshape
        Common_face = geompy.MakeCommonList([subshapef, reference_surface[reff_ind]], True) # check common intersection
        Common_face_Pro = geompy.BasicProperties(Common_face)
        subshapef_Pro = geompy.BasicProperties(subshapef) # extract volume of subshape
        print( "area difference",abs(Common_face_Pro[1]-subshapef_Pro[1]),"/",abs(Common_face_Pro[1]-reference_surface_Pro[reff_ind][1]))
        # if ( common face = subface) and (common face = ref face) => ref face = sub face
        if (abs(Common_face_Pro[1]-subshapef_Pro[1])<0.000001 )and (abs(Common_face_Pro[1]-reference_surface_Pro[reff_ind][1])<0.000001):
            Group_partition_surface.append([ Contact_name[reff_ind],Partition_surface_IDsList[subf_ind] ])
        # if ( common face = subface) and (common face < ref face) => sub face belong to ref face
        elif (abs(Common_face_Pro[1]-subshapef_Pro[1])<0.000001 ) and ((Common_face_Pro[1] - reference_surface_Pro[reff_ind][1])<-0.000001):
            temp_surface.append(Partition_surface_IDsList[subf_ind])
    if len(temp_surface) >1 : # the face is devided
        Group_partition_surface.append( [Contact_name[reff_ind],temp_surface ])
        print( Contact_name[reff_ind]," is devided and has sub IDs:{}\n".format(temp_surface))
if len(reference_surface) != len(Group_partition_surface): #+len(Group_partition_Multi_surface):
    print( "Geometry-Surface error please check ROI diameter and DBS lead Position ",len(reference_surface),len(Group_partition_surface),'\n')

print( 'Group_partition_surface',Group_partition_surface,'\n')

if(Lead2nd_Enable):
   Partition_profile = geompy.MakePartition(VolumeObject1+VolumeObject2+ContactObject1+ContactObject2+[Rest], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
else:
   Partition_profile = geompy.MakePartition(VolumeObject1+ContactObject1+[Rest], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)

new_volume_ID= geompy.SubShapeAllIDs(Partition_profile, geompy.ShapeType["SOLID"])
ID= list(set(Partition_volume_IDsList) ^ set (new_volume_ID))
Group_partition_volume.append(['Rest_1',ID[0]])
print( "REST ID:",ID)
print( 'Group_partition_volume',Group_partition_volume,'\n')
print("Create volume and surface group under partition_profile\n")

for i_solid in range (0,len (Group_partition_volume)):
    Group_volume[i_solid] = geompy.CreateGroup(Partition_profile, geompy.ShapeType["SOLID"])
    if (isinstance (Group_partition_volume[i_solid][1],list) == False):
        geompy.UnionIDs(Group_volume[i_solid], [Group_partition_volume[i_solid][1]])
    if (isinstance (Group_partition_volume[i_solid][1],list) == True):
        geompy.UnionIDs(Group_volume[i_solid], Group_partition_volume[i_solid][1])

#############################################

for i_surface in range (0,len (Group_partition_surface)):
    Group_surface[i_surface] = geompy.CreateGroup(Partition_profile, geompy.ShapeType["FACE"])
    if (isinstance (Group_partition_surface[i_surface][1],list) == False): # not a list
        geompy.UnionIDs(Group_surface[i_surface], [Group_partition_surface[i_surface][1]])
    if (isinstance (Group_partition_surface[i_surface][1],list) == True): #  it is a list
        geompy.UnionIDs(Group_surface[i_surface], Group_partition_surface[i_surface][1])
print( "Translate whole partition to Xm,Ym,Zm\n")
geompy.TranslateDXDYDZ(Partition_profile, Xm, Ym, Zm)
### add Vertices to geometry
if(Vertice_enable):
   for ver_ind in range (0,number_vertex):
       print("Add vertices to model\n")
       Vert.append(geompy.MakeVertex(Vert_array[ver_ind][0],Vert_array[ver_ind][1],Vert_array[ver_ind][2]))
       geompy.TranslateDXDYDZ(Vert[ver_ind], Xm, Ym, Zm)       ###Translate vertices to Xm,Ym,Zm
       geompy.addToStudy( Vert[ver_ind], 'Vert_{}'.format(ver_ind))
print("add to study\n")
############################################ end of extra code 2 ############################################
#############################################################################################################
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
#geompy.addToStudy( Circle_1, 'Circle_1' )
geompy.addToStudy( Contact_1, 'Contact_1' )
geompy.addToStudy( Contact_2, 'Contact_2' )
geompy.addToStudy( Contact_3, 'Contact_3' )
geompy.addToStudy( Contact_4, 'Contact_4' )
geompy.addToStudy( CV1, 'CV1' )
geompy.addToStudy( CV2, 'CV2' )
geompy.addToStudy( CV3, 'CV3' )
geompy.addToStudy( CV4, 'CV4' )
#geompy.addToStudy( Cylinder_1, 'Cylinder_1' )
#geompy.addToStudy( Sphere_1, 'Sphere_1' )
#geompy.addToStudy( Fuse_1, 'Fuse_1' )
#geompy.addToStudy( Cylinder_2, 'Cylinder_2' )
#geompy.addToStudy( Sphere_2, 'Sphere_2' )
#geompy.addToStudy( Fuse_2, 'Fuse_2' )
#geompy.addToStudy( encap_layer, 'encap_layer' )
#geompy.addToStudy( Sphere_ROI, 'Sphere_ROI' )
geompy.addToStudy( ROI, 'ROI' )
geompy.addToStudy( encap_outer_ROI, 'encap_outer_ROI' )
geompy.addToStudy( encap_inner_ROI, 'encap_inner_ROI' )
geompy.addToStudy( Fuse_all_lead_encap_ROI, 'Fuse_all_lead_encap_ROI' )
################################################################################################################
####################################### extra code 3 V10  15/12/18##############################################/
#for i in range(0,len(VolumeObject2)):/
#	geompy.addToStudy( VolumeObject2[i], 'VolumeObject2_{}'.format(i) )
#for i in range(0,len(ContactObject2)):
#	geompy.addToStudy( ContactObject2[i], 'ContactObject2_{}'.format(i) )
#for i in range(0,len(VolumeObject1)):
#	geompy.addToStudy( VolumeObject1[i], 'VolumeObject1_{}'.format(i) )
#for i in range(0,len(ContactObject1)):
#	geompy.addToStudy( ContactObject1[i], 'ContactObject1_{}'.format(i) )

geompy.addToStudy( Partition_profile, 'Partition_profile' )
for i_solid1 in range (0,len (Group_partition_volume)):
    geompy.addToStudyInFather( Partition_profile, Group_volume [i_solid1], Group_partition_volume[i_solid1][0])

for i_surface1 in range (0,len (Group_partition_surface)):
    geompy.addToStudyInFather( Partition_profile, Group_surface [i_surface1], Group_partition_surface[i_surface1][0])

##################################### end of extra code 3##########################################
###################################################################################################
Contact1_1=Group_surface[0]
Contact1_2=Group_surface[1]
Contact1_3=Group_surface[2]
Contact1_4=Group_surface[3]

encap_inner_ROI1=Group_volume[2]
encap_outer_ROI1=Group_volume[0]
ROI1=Group_volume[1]
Rest_1=Group_volume[7]

Floating_contacts=[]
float_indices=[]
for i in range(len(Phi_vector)):
    Floating_contacts.append(Group_volume[i+3])     #because the first contact is Group_volume[3]
    float_indices.append(i+3)

Auto_group_for_floating = geompy.CreateGroup(Partition_profile, geompy.ShapeType["SOLID"])
geompy.UnionList(Auto_group_for_floating, Floating_contacts[:])
geompy.addToStudyInFather( Partition_profile, Auto_group_for_floating, 'Auto_group_for_floating' )
###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
Mesh_1 = smesh.Mesh(Partition_profile)
NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
NETGEN_3D_Parameters_1.SetMaxSize( 25.4615 )
NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 0 )
NETGEN_3D_Parameters_1.SetMinSize( 0.000374134 )
NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )

NETGEN_1D_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=Contact1_1)
Sub_mesh_1 = NETGEN_1D_2D.GetSubMesh()
NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize( 0.05 )
NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 4 )
NETGEN_2D_Parameters_1.SetMinSize( 0.0001 )
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )

NETGEN_1D_2D_1 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=Contact1_2)
Sub_mesh_2 = NETGEN_1D_2D_1.GetSubMesh()
status = Mesh_1.AddHypothesis(NETGEN_2D_Parameters_1,Contact1_2)
NETGEN_1D_2D_2 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=Contact1_3)
Sub_mesh_3 = NETGEN_1D_2D_2.GetSubMesh()
status = Mesh_1.AddHypothesis(NETGEN_2D_Parameters_1,Contact1_3)
NETGEN_1D_2D_3 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=Contact1_4)
Sub_mesh_4 = NETGEN_1D_2D_3.GetSubMesh()
status = Mesh_1.AddHypothesis(NETGEN_2D_Parameters_1,Contact1_4)


NETGEN_1D_2D_3D_1 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=encap_inner_ROI1)
Sub_mesh_5 = NETGEN_1D_2D_3D_1.GetSubMesh()
NETGEN_3D_Parameters_2 = NETGEN_1D_2D_3D_1.Parameters()
NETGEN_3D_Parameters_2.SetMaxSize( encap_thickness )
NETGEN_3D_Parameters_2.SetSecondOrder( 0 )
NETGEN_3D_Parameters_2.SetOptimize( 1 )
NETGEN_3D_Parameters_2.SetFineness( 2 )
NETGEN_3D_Parameters_2.SetMinSize( 0.00283583 )
NETGEN_3D_Parameters_2.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_2.SetFuseEdges( 1 )
NETGEN_3D_Parameters_2.SetQuadAllowed( 0 )
isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_4, Sub_mesh_3, Sub_mesh_2, Sub_mesh_1, Sub_mesh_5 ] ])
NETGEN_1D_2D_3D_2 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=encap_outer_ROI1)
Sub_mesh_6 = NETGEN_1D_2D_3D_2.GetSubMesh()
NETGEN_3D_Parameters_3 = NETGEN_1D_2D_3D_2.Parameters()
NETGEN_3D_Parameters_3.SetMaxSize( encap_thickness )
NETGEN_3D_Parameters_3.SetSecondOrder( 0 )
NETGEN_3D_Parameters_3.SetOptimize( 1 )
NETGEN_3D_Parameters_3.SetFineness( 2 )
NETGEN_3D_Parameters_3.SetMinSize( 0.0333798 )
NETGEN_3D_Parameters_3.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_3.SetFuseEdges( 1 )
NETGEN_3D_Parameters_3.SetQuadAllowed( 0 )
isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_4, Sub_mesh_3, Sub_mesh_2, Sub_mesh_1, Sub_mesh_5, Sub_mesh_6 ] ])
NETGEN_1D_2D_3D_3 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=ROI1)
Sub_mesh_7 = NETGEN_1D_2D_3D_3.GetSubMesh()
NETGEN_3D_Parameters_4 = NETGEN_1D_2D_3D_3.Parameters()
NETGEN_3D_Parameters_4.SetMaxSize( 25.4615 )
NETGEN_3D_Parameters_4.SetSecondOrder( 0 )
NETGEN_3D_Parameters_4.SetOptimize( 1 )
NETGEN_3D_Parameters_4.SetFineness( 2 )
NETGEN_3D_Parameters_4.SetMinSize( 0.00328242 )
NETGEN_3D_Parameters_4.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_4.SetFuseEdges( 1 )
NETGEN_3D_Parameters_4.SetQuadAllowed( 0 )
isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_4, Sub_mesh_3, Sub_mesh_2, Sub_mesh_1, Sub_mesh_5, Sub_mesh_6, Sub_mesh_7 ] ])
NETGEN_1D_2D_3D_4 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=Rest_1)
Sub_mesh_8 = NETGEN_1D_2D_3D_4.GetSubMesh()
NETGEN_3D_Parameters_5 = NETGEN_1D_2D_3D_4.Parameters()
NETGEN_3D_Parameters_5.SetMaxSize( 2.5 )
NETGEN_3D_Parameters_5.SetSecondOrder( 0 )
NETGEN_3D_Parameters_5.SetOptimize( 1 )
NETGEN_3D_Parameters_5.SetFineness( 2 )
NETGEN_3D_Parameters_5.SetMinSize( 0.000374134 )
NETGEN_3D_Parameters_5.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_5.SetFuseEdges( 1 )
NETGEN_3D_Parameters_5.SetQuadAllowed( 0 )


NETGEN_1D_2D_3D_5 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=Group_volume[3])
Sub_mesh_9 = NETGEN_1D_2D_3D_5.GetSubMesh()
NETGEN_3D_Parameters_6 = NETGEN_1D_2D_3D_5.Parameters()
NETGEN_3D_Parameters_6.SetMaxSize( 25.4615 )
NETGEN_3D_Parameters_6.SetSecondOrder( 0 )
NETGEN_3D_Parameters_6.SetOptimize( 1 )
NETGEN_3D_Parameters_6.SetFineness( 2 )
NETGEN_3D_Parameters_6.SetMinSize( 0.000374134 )
NETGEN_3D_Parameters_6.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_6.SetFuseEdges( 1 )

NETGEN_1D_2D_3D_6 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=Group_volume[4])
Sub_mesh_10 = NETGEN_1D_2D_3D_6.GetSubMesh()
status = Mesh_1.AddHypothesis(NETGEN_3D_Parameters_6,Group_volume[4])

NETGEN_1D_2D_3D_7 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=Group_volume[5])
Sub_mesh_11 = NETGEN_1D_2D_3D_7.GetSubMesh()
status = Mesh_1.AddHypothesis(NETGEN_3D_Parameters_6,Group_volume[5])

NETGEN_1D_2D_3D_8 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=Group_volume[6])
Sub_mesh_12 = NETGEN_1D_2D_3D_8.GetSubMesh()
status = Mesh_1.AddHypothesis(NETGEN_3D_Parameters_6,Group_volume[6])

NETGEN_3D_Parameters_6.SetQuadAllowed( 0 )
#isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_4, Sub_mesh_3, Sub_mesh_2, Sub_mesh_1, Sub_mesh_5,Sub_mesh_9,Sub_mesh_6, Sub_mesh_7, Sub_mesh_8 ] ])
isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_4, Sub_mesh_3, Sub_mesh_2, Sub_mesh_1, Sub_mesh_5,Sub_mesh_9,Sub_mesh_10,Sub_mesh_11,Sub_mesh_12,Sub_mesh_6, Sub_mesh_7, Sub_mesh_8 ] ])

#if Phi_vector[0]==None:
#    Mesh_1.GetMesh().RemoveSubMesh( Sub_mesh_1 )
#if Phi_vector[1]==None:
#    Mesh_1.GetMesh().RemoveSubMesh( Sub_mesh_2 )
#if Phi_vector[2]==None:
#    Mesh_1.GetMesh().RemoveSubMesh( Sub_mesh_3 )
#if Phi_vector[3]==None:
#    Mesh_1.GetMesh().RemoveSubMesh( Sub_mesh_4 )

if Phi_vector[0]==0.0:
    Mesh_1.GetMesh().RemoveSubMesh( Sub_mesh_9 )
if Phi_vector[1]==0.0:
    Mesh_1.GetMesh().RemoveSubMesh( Sub_mesh_10 )
if Phi_vector[2]==0.0:
    Mesh_1.GetMesh().RemoveSubMesh( Sub_mesh_11 )
if Phi_vector[3]==0.0:
    Mesh_1.GetMesh().RemoveSubMesh( Sub_mesh_12 )

isDone = Mesh_1.Compute()

if Phi_vector[0]!=0.0:
    Mesh_1.GroupOnGeom(Group_volume[3],'Flt_cnt1',SMESH.VOLUME)
if Phi_vector[1]!=0.0:
    Mesh_1.GroupOnGeom(Group_volume[4],'Flt_cnt2',SMESH.VOLUME)
if Phi_vector[2]!=0.0:
    Mesh_1.GroupOnGeom(Group_volume[5],'Flt_cnt3',SMESH.VOLUME)
if Phi_vector[3]!=0.0:
    Mesh_1.GroupOnGeom(Group_volume[6],'Flt_cnt4',SMESH.VOLUME)

if Phi_vector[0]!=None:
    Mesh_1.GroupOnGeom(Contact1_1,'C1_1',SMESH.FACE)
if Phi_vector[1]!=None:
    Mesh_1.GroupOnGeom(Contact1_2,'C1_2',SMESH.FACE)
if Phi_vector[2]!=None:
    Mesh_1.GroupOnGeom(Contact1_3,'C1_3',SMESH.FACE)
if Phi_vector[3]!=None:
    Mesh_1.GroupOnGeom(Contact1_4,'C1_4',SMESH.FACE)

Encap_contact = Mesh_1.GroupOnGeom(encap_inner_ROI1,'Encap_contact',SMESH.VOLUME)
Encap_rest = Mesh_1.GroupOnGeom(encap_outer_ROI1,'Encap_rest',SMESH.VOLUME)
RegOfInt = Mesh_1.GroupOnGeom(ROI1,'RegOfInt',SMESH.VOLUME)
Rst = Mesh_1.GroupOnGeom(Rest_1,'Rst',SMESH.VOLUME)


## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')
smesh.SetName(NETGEN_3D_Parameters_2, 'NETGEN 3D Parameters_2')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
smesh.SetName(NETGEN_3D_Parameters_5, 'NETGEN 3D Parameters_5')
smesh.SetName(NETGEN_3D_Parameters_6, 'NETGEN 3D Parameters_6')

smesh.SetName(NETGEN_3D_Parameters_3, 'NETGEN 3D Parameters_3')

smesh.SetName(NETGEN_3D_Parameters_4, 'NETGEN 3D Parameters_4')
smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')

smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Rst, 'Rst')
smesh.SetName(Sub_mesh_12, 'Sub_mesh_12')
smesh.SetName(Sub_mesh_11, 'Sub_mesh_11')
smesh.SetName(Sub_mesh_10, 'Sub_mesh_10')
smesh.SetName(RegOfInt, 'RegOfInt')
smesh.SetName(Encap_rest, 'Encap_rest')
smesh.SetName(Encap_contact, 'Encap_contact')
smesh.SetName(Sub_mesh_7, 'Sub-mesh_7')
smesh.SetName(Sub_mesh_6, 'Sub-mesh_6')
smesh.SetName(Sub_mesh_5, 'Sub-mesh_5')
smesh.SetName(Sub_mesh_8, 'Sub-mesh_8')
smesh.SetName(Sub_mesh_9, 'Sub-mesh_9')

#if Phi_vector[0]!=None:
#
#    smesh.SetName(C1_1, 'C1_1')
#if Phi_vector[1]!=None:
#
#    smesh.SetName(C1_2, 'C1_2')
#if Phi_vector[2]!=None:
#
#    smesh.SetName(C1_3, 'C1_3')
#if Phi_vector[3]!=None:
#
#    smesh.SetName(C1_4, 'C1_4')

Mesh_1.ExportMED(os.environ['PATIENTDIR']+'/Meshes/Mesh_unref.med', 0, 33)

#if salome.sg.hasDesktop():
#  salome.sg.updateObjBrowser(True)

import killSalome
killSalome.killAllPorts()
