# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.3.0 with dump python functionality
###

import sys
import salome

salome.salome_init()

import os
sys.path.insert( 0, r'{}'.format(os.getcwd()))
sys.path.append('/usr/local/lib/python2.7/dist-packages')

#sys.path.insert( 0, r'/data/butenko/Salome_GMSH_gen/OS_DBS_rat/OSSA_core3')

###
### PARAVIS component
###

def make_a_screenshot(path_to_insert):
	import pvsimple
	pvsimple.ShowParaviewView()
	#### import the simple module from the paraview
	from pvsimple import *
	#### disable automatic camera reset on 'Show'
	pvsimple._DisableFirstRenderCameraReset()

	# create a new 'MED Reader'
	home_dir=os.path.expanduser("~")
	mesh_unrefmed = MEDReader(FileName=home_dir+path_to_insert+'/Meshes/Mesh_unref.med')

	# get active view
	renderView1 = GetActiveViewOrCreate('RenderView')
	# uncomment following to set a specific view size
	# renderView1.ViewSize = [1470, 639]

	# show data in view
	mesh_unrefmedDisplay = Show(mesh_unrefmed, renderView1)

	# reset view to fit data
	renderView1.ResetCamera()

	# create a new 'CSV Reader'
	vert_of_Neural_model_NEURONcsv = CSVReader(FileName=[home_dir+path_to_insert+'/Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv'])

	# Properties modified on vert_of_Neural_model_NEURONcsv
	vert_of_Neural_model_NEURONcsv.HaveHeaders = 0
	vert_of_Neural_model_NEURONcsv.FieldDelimiterCharacters = ' '

	# Create a new 'SpreadSheet View'
	spreadSheetView1 = CreateView('SpreadSheetView')
	spreadSheetView1.ColumnToSort = ''
	spreadSheetView1.BlockSize = 1024
	# uncomment following to set a specific view size
	# spreadSheetView1.ViewSize = [400, 400]

	# get layout
	layout1 = GetLayout()

	# place view in the layout
	layout1.AssignView(2, spreadSheetView1)

	# show data in view
	vert_of_Neural_model_NEURONcsvDisplay = Show(vert_of_Neural_model_NEURONcsv, spreadSheetView1)

	# create a new 'Table To Points'
	tableToPoints1 = TableToPoints(Input=vert_of_Neural_model_NEURONcsv)

	# Properties modified on tableToPoints1
	tableToPoints1.YColumn = 'Field 1'
	tableToPoints1.ZColumn = 'Field 2'

	# show data in view
	tableToPoints1Display = Show(tableToPoints1, spreadSheetView1)

	# hide data in view
	Hide(vert_of_Neural_model_NEURONcsv, spreadSheetView1)

	# set active view
	SetActiveView(renderView1)

	# set active source
	SetActiveSource(tableToPoints1)

	# show data in view
	tableToPoints1Display_1 = Show(tableToPoints1, renderView1)

	# hide data in view
	Hide(mesh_unrefmed, renderView1)

	# set active source
	SetActiveSource(mesh_unrefmed)

	# show data in view
	mesh_unrefmedDisplay = Show(mesh_unrefmed, renderView1)

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.98

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.8

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.78

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.75

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.64

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.63

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.48

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.36

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.35

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.25

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.23

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.22

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.21

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.22

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.24

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.25

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.26

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.27

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.26

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.25

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.24

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.23

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.22

	# Properties modified on mesh_unrefmedDisplay
	mesh_unrefmedDisplay.Opacity = 0.21

	# set active source
	SetActiveSource(tableToPoints1)

	# change solid color
	tableToPoints1Display_1.DiffuseColor = [0.00392156862745098, 1.0, 0.00392156862745098]

	#### saving camera placements for all active views

	# current camera placement for renderView1
	renderView1.CameraPosition = [-11.97298107774897, 3.0091410210384057, 9.63659807262625]
	renderView1.CameraFocalPoint = [1.024347367682281, 1.1129612638758126, 1.0326026334725262]
	renderView1.CameraViewUp = [0.5421585701688229, -0.07947879959151795, 0.8365089391082374]
	renderView1.CameraParallelScale = 7.199611409666049

	renderView1.ResetCamera()
	# save screenshot
	renderView1.ViewSize = [1600, 1000]
	renderView1.ResetCamera()
	SaveScreenshot(home_dir+path_to_insert+'/Images/InitMesh_and_Neurons.png', magnification=1, quality=100, view=renderView1)


	#if salome.sg.hasDesktop():
	#  salome.sg.updateObjBrowser(True)

	import killSalome
	killSalome.killAllPorts()

if __name__ == '__main__':
	make_a_screenshot(*sys.argv[1:])
