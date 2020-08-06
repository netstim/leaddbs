#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 14:19:55 2020

@author: butenko
"""

#IMPORTANT: if you want to load this script direclty to Paraview, copy it to OSS_Platform/ and load from there


#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

list_of_connections=['Direct_upsampled_GPi_Str', 'Direct_upsampled_SN_Str', 'HDP', 'HDP_STN_GPi_upsampled_30', 'HDP_STN_SN_10', 'Indirect_upsampled_STN_GPe', 'Indirect_upsampled_STN_SN_10']

import os
def make_a_screenshot(path_to_insert):
	# create a new 'PVD Reader'
	home_dir=os.path.expanduser("~")
	for connection in list_of_connections:
	    
	    connection_csv = CSVReader(FileName=[home_dir+path_to_insert+'/Field_solutions/Activation/Neuron_model_results_'+connection+'.csv'])

	    #gPi_mask_GPe_mask_5_Indirect_Branch_MC_Str_GPe_GPicsv = CSVReader(FileName=['/home/konstantin/Documents/brukerMRI-master/Classified_by_Sosoho/Segment/GPi_mask_GPe_mask_5_Indirect_Branch_MC_Str_GPe_GPi.csv'])
	    
	    # Properties modified on gPi_mask_GPe_mask_5_Indirect_Branch_MC_Str_GPe_GPicsv
	    connection_csv.HaveHeaders = 0
	    connection_csv.FieldDelimiterCharacters = ' '
	    
	    # get active view
	    renderView1 = GetActiveViewOrCreate('RenderView')
	    viewLayout1 = GetLayout()
	    
	    # Create a new 'SpreadSheet View'
	    spreadSheetView1 = CreateView('SpreadSheetView')
	    spreadSheetView1.BlockSize = 1024L
	    spreadSheetView1.ColumnToSort = ''
	    viewLayout1.AssignView(2, spreadSheetView1)
	    
	    # show data in view
	    connection_csvDisplay = Show(connection_csv, spreadSheetView1)
	    # trace defaults for the display properties.
	    connection_csvDisplay.FieldAssociation = 'Row Data'
	    
	    # create a new 'Table To Points'
	    #tableToPoints1 = TableToPoints(Input=neuron_model_results_STN_mask_GPe_mask_20_Excitatory_stn2gpe_asscsv)   
	    ## rename source object
	    #RenameSource('Nmodels', tableToPoints1)   

	    globals()[connection+'_table']= TableToPoints(Input=connection_csv)        #very bad way, but it is just a standalone script    
	    globals()[connection+'_table'].XColumn = 'Field 0'
	    globals()[connection+'_table'].YColumn = 'Field 1'
	    globals()[connection+'_table'].ZColumn = 'Field 2'
		
	#    tableToPoints1 = TableToPoints(Input=connection_csv)
	#    tableToPoints1.XColumn = 'Field 0'
	#    tableToPoints1.YColumn = 'Field 0'
	#    tableToPoints1.ZColumn = 'Field 0'
	#    
	#    # Properties modified on tableToPoints1
	#    tableToPoints1.YColumn = 'Field 1'
	#    tableToPoints1.ZColumn = 'Field 2'
	    
	    # show data in view
	    #tableToPoints1Display = Show(tableToPoints1, spreadSheetView1)  
	    globals()[connection+'_tableDisplay']= Show(globals()[connection+'_table'], spreadSheetView1)
	    RenameSource(connection,globals()[connection+'_table'])  
	    # hide data in view
	    Hide(connection_csv, spreadSheetView1)
	    
	    # set active view
	    SetActiveView(renderView1)
	    
	    # set active source
	    #SetActiveSource(tableToPoints1)
	    SetActiveSource(globals()[connection+'_table'])
	    
	    ## show data in view
	    #tableToPoints1Display_1 = Show(tableToPoints1, renderView1)
	    ## trace defaults for the display properties.
	    #tableToPoints1Display_1.ColorArrayName = [None, '']
	    #tableToPoints1Display_1.GlyphType = 'Arrow'

	    globals()[connection+'_tableDisplay_1'] = Show(globals()[connection+'_table'], renderView1)
	    # trace defaults for the display properties.
	    globals()[connection+'_tableDisplay_1'].ColorArrayName = [None, '']
	    globals()[connection+'_tableDisplay_1'].GlyphType = 'Arrow'
	    
	    # reset view to fit data
	    renderView1.ResetCamera()
	    
	    # set active view
	    SetActiveView(spreadSheetView1)
	    
	    # destroy spreadSheetView1
	    Delete(spreadSheetView1)
	    del spreadSheetView1
	    
	    # close an empty frame
	    viewLayout1.Collapse(2)
	    
	    # set active view
	    SetActiveView(renderView1)

	    # set active source
	    SetActiveSource(globals()[connection+'_table'])
	    
	    # show data in view
	    # trace defaults for the display properties.
	    globals()[connection+'_tableDisplay_1'].Representation = 'Surface'
	    globals()[connection+'_tableDisplay_1'].OSPRayScaleArray = 'Field 3'
	    globals()[connection+'_tableDisplay_1'].OSPRayScaleFunction = 'PiecewiseFunction'
	    globals()[connection+'_tableDisplay_1'].SelectOrientationVectors = 'Field 3'
	    globals()[connection+'_tableDisplay_1'].ScaleFactor = 1.506762191
	    globals()[connection+'_tableDisplay_1'].SelectScaleArray = 'Field 3'
	    #tableToPoints1Display.GlyphType = 'Arrow'
	    globals()[connection+'_tableDisplay_1'].GlyphTableIndexArray = 'Field 3'
	    globals()[connection+'_tableDisplay_1'].DataAxesGrid = 'GridAxesRepresentation'
	    globals()[connection+'_tableDisplay_1'].PolarAxes = 'PolarAxesRepresentation'
	    
	    # reset view to fit data
	    renderView1.ResetCamera()
	    
	    # set scalar coloring
	    ColorBy(globals()[connection+'_tableDisplay_1'], ('POINTS', 'Field 3'))
	    
	    # rescale color and/or opacity maps used to include current data range
	    globals()[connection+'_tableDisplay_1'].RescaleTransferFunctionToDataRange(True, False)
	    
	    # show color bar/color legend
	    globals()[connection+'_tableDisplay_1'].SetScalarBarVisibility(renderView1, True)
	    
	    # get color transfer function/color map for 'Field3'
	    field3LUT = GetColorTransferFunction('Field3')


	# reset view to fit data
	renderView1.ResetCamera()

	#### saving camera placements for all active views

	# current camera placement for renderView1
	# current camera placement for renderView1
	renderView1.CameraPosition = [-5.943765067238889, 30.063310553609174, 19.481178764885886]
	renderView1.CameraFocalPoint = [20.037838396112274, -8.840344561377364, -4.671394748401762]
	renderView1.CameraViewUp = [-0.4484871730366652, 0.23572557572149833, -0.8621442504432473]
	renderView1.CameraParallelScale = 24.899825806468712
	renderView1.ViewSize = [1600, 1000]
	renderView1.ResetCamera()
	SaveScreenshot(home_dir+path_to_insert+'/Images/Axon_activation.png', magnification=1, quality=100, view=renderView1)

	#### uncomment the following to render all views
	# RenderAllViews()
	# alternatively, if you want to write images, you can use SaveScreenshot(...).

if __name__ == '__main__':
	make_a_screenshot(*sys.argv[1:])
