# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.3.0 with dump python functionality
###### Run with DPS_lead_position_V9.py 

import sys
import salome
import os

salome.salome_init()
sys.path.insert( 0, "r'"+os.getcwd())

###
### GEOM component
###
######################################################################################################
########################################### extra code 1 V10 15/12/18#############################################
###### This file runs with DBS_lead_position_V10.py

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

Zt_tip=Zt

#for Lead-DBS, the tip point should be shifted down (they use the middle of the lowest contact as the reference point)

Vert_array =[0];
number_vertex = len(Vert_array)
Vert = []
VolumeObject1 = []
ContactObject1 = []
VolumeObject2 = []
ContactObject2 = []
print(" DBS_lead's Geometry\n")
######################################### end of extra code 1 ########################################
######################################################################################################
import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS

#Monopolar Electrode:
print('Monopolar Electrode C04')
a_semi=8.3805
b_semi=14.23
c_semi=5.62
r_lead=0.1125
####################################################
##       Begin of NoteBook variables section      ##
####################################################
import salome_notebook
notebook = salome_notebook.notebook
notebook.set("r_lead", r_lead)
notebook.set("encap_thickness", encap_thickness)
notebook.set("ROI_radial", ROI_radial)
notebook.set("Xm", Xm)
notebook.set("Ym", Ym)
notebook.set("Zm", Zm)
notebook.set("Xt", Xt)
notebook.set("Yt", Yt)
notebook.set("Zt", Zt)
notebook.set("OZ_angle", OZ_angle)
notebook.set("rPlusEncaps", "r_lead+encap_thickness")
notebook.set("a_semi", a_semi)
notebook.set("b_semi", b_semi)
notebook.set("c_semi", c_semi)
####################################################
##        End of NoteBook variables section       ##
####################################################
geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
#lead:
Sphere_1 = geompy.MakeSphereR("r_lead")
Cylinder_1 = geompy.MakeCylinderRH("r_lead", 15)
DBS_lead_preTrans1 = geompy.MakeFuseList([Sphere_1, Cylinder_1], True, True)
#encaps
Sphere_2 = geompy.MakeSphereR("rPlusEncaps")
Cylinder_2 = geompy.MakeCylinderRH("rPlusEncaps", 15)
encaps_fused = geompy.MakeFuseList([Sphere_2, Cylinder_2], True, True)
encaps = geompy.MakeCutList(encaps_fused, [DBS_lead_preTrans1], True)
Sphere_3 = geompy.MakeSphereR("ROI_radial")
Common_1 = geompy.MakeCommonList([encaps_fused, Sphere_3], True)
encap_outer_ROI_preTrans1 = geompy.MakeCutList(encaps, [Common_1], True)
encap_inner_ROI_preTrans1 = geompy.MakeCommonList([encaps, Sphere_3], True)
#ROI
ROI_preTrans1 = geompy.MakeCutList(Sphere_3, [encaps_fused], True)
#Translate so that tip at coordinate of DBS Stim.
Translation_1_lead = geompy.MakeTranslation(DBS_lead_preTrans1, 0, 0, "r_lead")
Translation_1_encapO = geompy.MakeTranslation(encap_outer_ROI_preTrans1, 0, 0, "r_lead")
Translation_1_encapI = geompy.MakeTranslation(encap_inner_ROI_preTrans1, 0, 0, "r_lead")
Translation_1_ROI = geompy.MakeTranslation(ROI_preTrans1, 0, 0, "r_lead")
DBS_lead_preCut = geompy.MakeTranslation(Translation_1_lead, "Xt", "Yt", "Zt")
encap_outer_ROI_preCut = geompy.MakeTranslation(Translation_1_encapO, "Xt", "Yt", "Zt")
encap_inner_ROI_preCut = geompy.MakeTranslation(Translation_1_encapI, "Xt", "Yt", "Zt")
ROI_preCut = geompy.MakeTranslation(Translation_1_ROI, "Xt", "Yt", "Zt")
#Import Brain Geometry and define modelling domain
Fuse_1 = geompy.MakeFuseList([DBS_lead_preCut, encap_outer_ROI_preCut, encap_inner_ROI_preCut, ROI_preCut], True, True)
#    print (" Load brain image \n")
#    if (Brain_map[-4:] == 'brep'):
#    	brain_solid = geompy.ImportBREP( Brain_map )
#    elif (Brain_map[-4:] == 'step'):
#    	brain_solid = geompy.ImportSTEP( Brain_map )
#    elif (Brain_map[-4:] == 'iges'):
#    	brain_solid = geompy.ImportIGES( Brain_map )
#    elif (Brain_map[-4:] == '.stl'):
#    	brain_solid = geompy.ImportSTL( Brain_map )
#    else:
#    	print (" unknow imported file format")
Sphere_forEllipsoid=geompy.MakeSphereR(1)
Scale_forEllipsoid=geompy.MakeScaleAlongAxes(Sphere_forEllipsoid,None,a_semi,b_semi,c_semi)
brain_solid = geompy.MakeTranslation(Scale_forEllipsoid, 0, -2, 2.2)
ImportBrain = brain_solid# geompy.ImportBREP("/data/andrea/OSS-DBS/OSS_platform/Brain_substitute.brep" )
Cut_1 = geompy.MakeCutList(Fuse_1, [ImportBrain], True)
encap_outer_ROI_preTrans2 = geompy.MakeCutList(encap_outer_ROI_preCut, [Cut_1], True)
DBS_lead_preTrans2 = geompy.MakeCutList(DBS_lead_preCut, [Cut_1], True)
ROI_preTrans2 = geompy.MakeCutList(ROI_preCut, [Cut_1], True)
encap_inner_ROI_preTrans2 = geompy.MakeCutList(encap_inner_ROI_preCut, [Cut_1], True)
Rest_preTrans2 = geompy.MakeCutList(ImportBrain, [Fuse_1], True)
#Translate to fit MRI
DBS_lead = geompy.MakeTranslation(DBS_lead_preTrans2, "Xm", "Ym", "Zm")
encap_outer_ROI = geompy.MakeTranslation(encap_outer_ROI_preTrans2, "Xm", "Ym", "Zm")
encap_inner_ROI = geompy.MakeTranslation(encap_inner_ROI_preTrans2, "Xm", "Ym", "Zm")
ROI = geompy.MakeTranslation(ROI_preTrans2, "Xm", "Ym", "Zm")
Rest = geompy.MakeTranslation(Rest_preTrans2, "Xm", "Ym", "Zm")#rest is brain solid and surface
#define Contact_1 and Contact_2:
[Contact_1,LeadRest_1,LeadRest_2] = geompy.ExtractShapes(DBS_lead, geompy.ShapeType["FACE"], True)
[Contact_2,Rest_1,Rest_2] = geompy.ExtractShapes(Rest, geompy.ShapeType["FACE"], True)
#Make Partition:
Partition_profile = geompy.MakePartition([ROI, Rest, DBS_lead, encap_outer_ROI, encap_inner_ROI], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
#[Contact_2,encapsTip,Contact_1,roi,encapsInner,insulationInner,connector,encapsOuter,insulationOuter,insulationTop,encapsulationTop] = geompy.ExtractShapes(Partition_profile, geompy.ShapeType["FACE"], True)
#[f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11] = geompy.ExtractShapes(Partition_profile, geompy.ShapeType["FACE"], True)

###################################################### Geometry and extra code interface ###############################################################
VolumeObject1 = [ encap_outer_ROI,ROI,encap_inner_ROI,Rest]         # Declare objects included to partition, encap_outer_ROI always @1st position
Volume_name1  = ['encap_outer_ROI1','ROI1','encap_inner_ROI1','Rest'] # Declare name of the group in the partition for volume
ContactObject1 = [Contact_1,Contact_2]   
Contact_name1 = ['Contact1_1','Contact1_2']


Volume1_Pro = [geompy.BasicProperties( VolumeObject1[0])]*len(VolumeObject1)
Contact1_Pro = [geompy.BasicProperties( ContactObject1[0])]*len(ContactObject1)
for i in range(0,len(VolumeObject1)):
	Volume1_Pro[i] = geompy.BasicProperties( VolumeObject1[i])
for i in range(0,len(ContactObject1)):
	Contact1_Pro[i] = geompy.BasicProperties( ContactObject1[i])


print("Create reference groups for ID identification process\n")
  
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

#### find group volume ID ######################################################################
Partition_volume_IDsList = geompy.SubShapeAllIDs(Partition_profile, geompy.ShapeType["SOLID"]) # list all sub shape volume in Partition
print("Partition_volume_IDsList",Partition_volume_IDsList, '\n')

for ref_ind in range (0, len(reference_volume)):
	temp_volume = []
	for sub_ind in range (0, len (Partition_volume_IDsList)):
		subshape = geompy.GetSubShape(Partition_profile, [Partition_volume_IDsList[sub_ind]]) # get subshape
		subshape_Pro = geompy.BasicProperties(subshape)       # extract volume of subshape
		Common_volume = geompy.MakeCommonList([subshape, reference_volume[ref_ind]], True) # check common intersection
		Common_volume_Pro = geompy.BasicProperties(Common_volume) 
		print("volume difference",abs(Common_volume_Pro[2]-subshape_Pro[2]),"/",abs(Common_volume_Pro[2]-reference_volume_Pro[ref_ind][2]))
		# if ( common volume = subshape) and (common volume = ref volume) => ref volume = sub shape
		if (abs(Common_volume_Pro[2]-subshape_Pro[2])< 0.00001) and (abs(Common_volume_Pro[2]-reference_volume_Pro[ref_ind][2])<0.00001):
		
			Group_partition_volume.append([Volume_name[ref_ind],Partition_volume_IDsList[sub_ind]])
		# if ( common volume = subshape) and (common volume < ref volume) => sub shape belong to ref volume
		elif (abs(Common_volume_Pro[2]-subshape_Pro[2])< 0.00001) and ((Common_volume_Pro[2] - reference_volume_Pro[ref_ind][2])<-0.00001):
			temp_volume.append( Partition_volume_IDsList[sub_ind] )
	if len(temp_volume) >1 : # the volume is devided
		Group_partition_volume.append([Volume_name[ref_ind],temp_volume ])
		print(Volume_name[ref_ind]," is devided and has sub IDs:{}\n".format(temp_volume))
if len(reference_volume) != len(Group_partition_volume):
	print("Geometry-volume error please check ROI diameter and DBS lead Position ",len(reference_volume),len(Group_partition_volume))
print('Group_partition_volume',Group_partition_volume,'\n')

#### find group surface ID ######################################################################
Partition_surface_IDsList = geompy.SubShapeAllIDs(Partition_profile, geompy.ShapeType["FACE"]) # list all sub shape face in Partition
print('Partition_surface_IDsList',Partition_surface_IDsList,'\n')
sub_face = [] ## store devided faces
for reff_ind in range (0, len (reference_surface)):
	temp_surface = []
	for subf_ind in range (0, len(Partition_surface_IDsList)):
		subshapef = geompy.GetSubShape(Partition_profile, [Partition_surface_IDsList[subf_ind]]) # get subshape
		Common_face = geompy.MakeCommonList([subshapef, reference_surface[reff_ind]], True) # check common intersection
		Common_face_Pro = geompy.BasicProperties(Common_face) 
		subshapef_Pro = geompy.BasicProperties(subshapef) # extract volume of subshape
		print("area difference",abs(Common_face_Pro[1]-subshapef_Pro[1]),"/",abs(Common_face_Pro[1]-reference_surface_Pro[reff_ind][1]))
		# if ( common face = subface) and (common face = ref face) => ref face = sub face
		if (abs(Common_face_Pro[1]-subshapef_Pro[1])<0.000001 )and (abs(Common_face_Pro[1]-reference_surface_Pro[reff_ind][1])<0.000001): 
			Group_partition_surface.append([ Contact_name[reff_ind],Partition_surface_IDsList[subf_ind] ])
		# if ( common face = subface) and (common face < ref face) => sub face belong to ref face 
		elif (abs(Common_face_Pro[1]-subshapef_Pro[1])<0.000001 ) and ((Common_face_Pro[1] - reference_surface_Pro[reff_ind][1])<-0.000001):
			temp_surface.append(Partition_surface_IDsList[subf_ind])
	if len(temp_surface) >1 : # the face is devided
		Group_partition_surface.append( [Contact_name[reff_ind],temp_surface ])		
		print(Contact_name[reff_ind]," is devided and has sub IDs:{}\n".format(temp_surface))
if len(reference_surface) != len(Group_partition_surface): #+len(Group_partition_Multi_surface):
	print("Geometry-Surface error please check ROI diameter and DBS lead Position ",len(reference_surface),len(Group_partition_surface),'\n')

print('Group_partition_surface',Group_partition_surface,'\n')


for i_solid in range (0,len(Group_partition_volume)):
    Group_volume[i_solid] = geompy.CreateGroup(Partition_profile, geompy.ShapeType["SOLID"])
    #print(Group_partition_volume[i_solid][0],Group_partition_volume[i_solid][1])
    if (isinstance (Group_partition_volume[i_solid][1],list) == False):
        geompy.UnionIDs(Group_volume[i_solid], [Group_partition_volume[i_solid][1]])
    if (isinstance (Group_partition_volume[i_solid][1],list) == True):
        geompy.UnionIDs(Group_volume[i_solid], Group_partition_volume[i_solid][1])

#############################################

for i_surface in range (0,len(Group_partition_surface)):
	Group_surface[i_surface] = geompy.CreateGroup(Partition_profile, geompy.ShapeType["FACE"])
	if (isinstance (Group_partition_surface[i_surface][1],list) == False): # not a list
		geompy.UnionIDs(Group_surface[i_surface], [Group_partition_surface[i_surface][1]])
	if (isinstance (Group_partition_surface[i_surface][1],list) == True): #  it is a list
		geompy.UnionIDs(Group_surface[i_surface], Group_partition_surface[i_surface][1])

print("add to study\n")
############################################## end of extra code 2 ##########################
#geompy.addToStudyInFather( Partition_profile, insulationInner, 'insulationInner' )
#geompy.addToStudyInFather( Partition_profile, connector, 'connector' )
#geompy.addToStudyInFather( Partition_profile, encapsOuter, 'encapsOuter' )
#geompy.addToStudyInFather( Partition_profile, insulationOuter, 'insulationOuter' )
#geompy.addToStudyInFather( Partition_profile, insulationTop, 'insulationTop' )##################
###############################################################################################################
#############################################################Diesen Teil weiter nach unten schieben XD und ab Zeile 235 weiter machen
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( ROI, 'ROI' )
geompy.addToStudy( Rest, 'Rest' )
geompy.addToStudy( DBS_lead, 'DBS_lead' )
geompy.addToStudy( encap_outer_ROI, 'encap_outer_ROI' )
geompy.addToStudy( encap_inner_ROI, 'encap_inner_ROI' )
#geompy.addToStudy( Partition_profile, 'Partition_profile' )
#geompy.addToStudyInFather( Partition_profile, Contact_2, 'Contact_2' )
#geompy.addToStudyInFather( Partition_profile, Contact_1, 'Contact_1' )
#geompy.addToStudyInFather( Partition_profile, roi, 'roi' )
#geompy.addToStudyInFather( Partition_profile, encapsInner, 'encapsInner' )
#geompy.addToStudyInFather( Partition_profile, encapsulationTop, 'encapsulationTop' )

##################################################################################################################
######################################### extra code 3 V10  15/12/18##############################################/
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

####################################### end of extra code 3##########################################
#####################################################################################################
anode_surf=Contact1_Pro[0][1]
cath_surf=Contact1_Pro[1][1]


Contact1_1=Group_surface[0]
Contact1_2=Group_surface[1]

encap_inner_ROI1=Group_volume[2]
encap_outer_ROI1=Group_volume[0]
ROI1=Group_volume[1]
Rest_1=Group_volume[3]

########################################Meshing################################
import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

anode_mesh_max=0.005
cathode_mesh_max=1.0

smesh = smeshBuilder.New()
Mesh_1 = smesh.Mesh(Partition_profile)
NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
NETGEN_3D_Parameters_1.SetMaxSize( 2.25 )
NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 2 )
NETGEN_3D_Parameters_1.SetMinSize( 0.1 )
NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_1D_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=Contact1_1)
Sub_mesh_1 = NETGEN_1D_2D.GetSubMesh()
NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize( anode_mesh_max )
NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 3 )
NETGEN_2D_Parameters_1.SetMinSize( 1e-05 )
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_1D_2D_3D_1 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=encap_inner_ROI1)
Sub_mesh_2 = NETGEN_1D_2D_3D_1.GetSubMesh()
NETGEN_3D_Parameters_2 = NETGEN_1D_2D_3D_1.Parameters()
NETGEN_3D_Parameters_2.SetMaxSize( encap_thickness )
NETGEN_3D_Parameters_2.SetSecondOrder( 0 )
NETGEN_3D_Parameters_2.SetOptimize( 2 )
NETGEN_3D_Parameters_2.SetMinSize( 0.0001 )
NETGEN_3D_Parameters_2.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_2.SetFuseEdges( 1 )
NETGEN_3D_Parameters_2.SetQuadAllowed( 0 )
isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_1, Sub_mesh_2 ] ])
NETGEN_1D_2D_3D_2 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=encap_outer_ROI1)
Sub_mesh_3 = NETGEN_1D_2D_3D_2.GetSubMesh()
NETGEN_3D_Parameters_3 = NETGEN_1D_2D_3D_2.Parameters()
NETGEN_3D_Parameters_3.SetMaxSize( encap_thickness )
NETGEN_3D_Parameters_3.SetSecondOrder( 0 )
NETGEN_3D_Parameters_3.SetOptimize( 1 )
NETGEN_3D_Parameters_3.SetFineness( 3 )
NETGEN_3D_Parameters_3.SetMinSize( 0.01 )
NETGEN_3D_Parameters_3.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_3.SetFuseEdges( 1 )
NETGEN_3D_Parameters_3.SetQuadAllowed( 0 )
isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_1, Sub_mesh_2, Sub_mesh_3 ] ])
NETGEN_3D_Parameters_2.SetFineness( 3 )
NETGEN_1D_2D_3D_3 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=ROI1)
Sub_mesh_4 = NETGEN_1D_2D_3D_3.GetSubMesh()
NETGEN_3D_Parameters_4 = NETGEN_1D_2D_3D_3.Parameters()
NETGEN_3D_Parameters_4.SetMaxSize( 0.35 )
NETGEN_3D_Parameters_4.SetSecondOrder( 0 )
NETGEN_3D_Parameters_4.SetOptimize( 1 )
NETGEN_3D_Parameters_4.SetFineness( 2 )
NETGEN_3D_Parameters_4.SetMinSize( 0.001 )
NETGEN_3D_Parameters_4.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_4.SetFuseEdges( 1 )
NETGEN_3D_Parameters_4.SetQuadAllowed( 0 )
isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_1, Sub_mesh_2, Sub_mesh_3, Sub_mesh_4 ] ])
NETGEN_1D_2D_3D_4 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=Rest_1)
Sub_mesh_5 = NETGEN_1D_2D_3D_4.GetSubMesh()
NETGEN_3D_Parameters_5 = NETGEN_1D_2D_3D_4.Parameters()
NETGEN_3D_Parameters_5.SetMaxSize( cathode_mesh_max )
NETGEN_3D_Parameters_5.SetSecondOrder( 0 )
NETGEN_3D_Parameters_5.SetOptimize( 1 )
NETGEN_3D_Parameters_5.SetFineness( 2 )
NETGEN_3D_Parameters_5.SetMinSize( 0.1 )
NETGEN_3D_Parameters_5.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_5.SetFuseEdges( 1 )
NETGEN_3D_Parameters_5.SetQuadAllowed( 0 )
isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_1, Sub_mesh_2, Sub_mesh_3, Sub_mesh_4, Sub_mesh_5 ] ])
isDone = Mesh_1.Compute()

C1_1 = Mesh_1.GroupOnGeom(Contact1_1,'C1_1',SMESH.FACE)
C1_2 = Mesh_1.GroupOnGeom(Contact1_2,'C1_2',SMESH.FACE)
Encap_rest = Mesh_1.GroupOnGeom(encap_outer_ROI1,'Encap_rest',SMESH.VOLUME)
Encap_contact = Mesh_1.GroupOnGeom(encap_inner_ROI1,'Encap_contact',SMESH.VOLUME)
RegOfInt = Mesh_1.GroupOnGeom(ROI1,'RegOfInt',SMESH.VOLUME)
Rst = Mesh_1.GroupOnGeom(Rest_1,'Rst',SMESH.VOLUME)
measure = smesh.CreateMeasurements()
An_meshed_surf=measure.Area(C1_1)
Cat_meshed_surf=measure.Area(C1_2)


print("Core_div: ")
print(abs(An_meshed_surf-anode_surf)/anode_surf)

print("Outer_div: ")
print(abs(Cat_meshed_surf-cath_surf)/cath_surf)

smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
smesh.SetName(NETGEN_3D_Parameters_2, 'NETGEN 3D Parameters_2')
smesh.SetName(NETGEN_3D_Parameters_3, 'NETGEN 3D Parameters_3')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')
#smesh.SetName(NETGEN_2D_Parameters_2, 'NETGEN 2D Parameters_2')
smesh.SetName(C1_1, 'C1_1')
smesh.SetName(NETGEN_3D_Parameters_4, 'NETGEN 3D Parameters_4')
smesh.SetName(C1_2, 'C1_2')
smesh.SetName(NETGEN_3D_Parameters_5, 'NETGEN 3D Parameters_5')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
#smesh.SetName(Sub_mesh_6, 'Sub-mesh_6')
smesh.SetName(Sub_mesh_5, 'Sub-mesh_5')
smesh.SetName(Rst, 'Rst')
smesh.SetName(RegOfInt, 'RegOfInt')
smesh.SetName(Encap_rest, 'Encap_rest')
smesh.SetName(Encap_contact, 'Encap_contact')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')

#Mesh_1.ExportMED('opt/patient/Meshes/Mesh_unref.med')
#Mesh_1.ExportMED('Mesh_unref.med')
Mesh_1.ExportMED(os.environ['PATIENTDIR']+'/Meshes/Mesh_unref.med', 0, 33)

print("Mesh was saved\n")
print(os.getcwd())
#if salome.sg.hasDesktop():#
#  salome.sg.updateObjBrowser(True)

import killSalome
killSalome.killAllPorts()
