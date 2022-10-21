# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.3.0 with dump python functionality
###

import sys
import salome
import os
salome.salome_init()
sys.path.insert( 0, "r'"+os.getcwd())

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS

#input variables
DX_MRI=80.0
DY_MRI=80.0
DZ_MRI=80.0

X_tip=-11.1660453318873
Y_tip=-14.7114212603611
Z_tip=-9.62969801365702

##################
geompy = geomBuilder.New()
print("brain model file is saved at"+ os.getcwd()+"/Brain_substitute.brep")
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Box1_1 = geompy.MakeBoxDXDYDZ(DX_MRI, DY_MRI, DZ_MRI)
Sphere1_1 = geompy.MakeSphereR(DX_MRI/2)
Brain_model = geompy.MakeScaleAlongAxes(Sphere1_1, None, 1, DY_MRI/DX_MRI, DZ_MRI/DX_MRI)
geompy.TranslateDXDYDZ(Brain_model, X_tip, Y_tip, Z_tip)

MRI_DX_max_ROI=DX_MRI*0.9
MRI_DY_max_ROI=DY_MRI*0.9
MRI_DZ_max_ROI=DZ_MRI*0.9
Sphere1_2 = geompy.MakeSphereR(MRI_DX_max_ROI/2)
Brain_model_ROI = geompy.MakeScaleAlongAxes(Sphere1_2, None, 1, MRI_DY_max_ROI/MRI_DX_max_ROI, MRI_DZ_max_ROI/MRI_DX_max_ROI)
geompy.TranslateDXDYDZ(Brain_model_ROI, X_tip, Y_tip, Z_tip)


geompy.ExportBREP(Brain_model, os.environ['PATIENTDIR']+"/Brain_substitute.brep" )
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
#geompy.addToStudy( Box_1, 'Box_1' )
#geompy.addToStudy( Sphere_1, 'Sphere_1' )
#geompy.addToStudy( Point_1, 'Point_1' )
geompy.addToStudy( Brain_model, 'Brain_model' )
geompy.addToStudy( Brain_model_ROI, 'Brain_model_ROI' )

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
Mesh_1 = smesh.Mesh(Brain_model_ROI)
NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
NETGEN_3D_Parameters_1.SetMaxSize( 100.732 )
NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 4 )
NETGEN_3D_Parameters_1.SetMinSize( 0.00533808 )
NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )
isDone = Mesh_1.Compute()
smesh.SetName(Mesh_1, 'Mesh_1')

Mesh_1.ExportMED(os.environ['PATIENTDIR']+'/Meshes/Mesh_brain_substitute_max_ROI.med', 0, 33)

smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')

#if salome.sg.hasDesktop():
#  salome.sg.updateObjBrowser(True)

import killSalome
killSalome.killAllPorts()
