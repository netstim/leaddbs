
#IMPORTANT: if you want to load this script direclty to Paraview, copy it to OSS_Platform/ and load from there


#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
#get this list using arr_srt = [str(r) for r in list(hf.keys())] 

list_of_connections=['GPi_mask_GPe_mask_5_Indirect_gpe2stn_sm', 'GPi_mask_GPe_mask_6_Indirect_gpe2stn_ass', 'GPi_mask_Th_mask_40_Inhibitory_ansa_lenticularis', 'GPi_mask_Th_mask_40_Inhibitory_lenticular_fasciculus', 'STN_mask_GPe_mask_20_Excitatory_stn2gpe_ass', 'STN_mask_GPe_mask_20_Indirect_gpe2stn_ass', 'STN_mask_GPe_mask_8_Excitatory_stn2gpe_sm', 'STN_mask_GPe_mask_8_Indirect_gpe2stn_sm', 'STN_mask_GPi_mask_11_Excitatory_stn2gpi_ass', 'STN_mask_GPi_mask_9_Excitatory_stn2gpi_sm', 'STN_mask_MC_mask_35_HDP_Branch_face', 'STN_mask_MC_mask_35_HDP_Branch_lowerex', 'STN_mask_MC_mask_35_HDP_Branch_upperex', 'STN_mask_PMC_mask_35_HDP_Branch_Premotor']

#prepare a color map
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

viridis = cm.get_cmap('viridis', len(list_of_connections))
index_population=0

for connection in list_of_connections:
    
    connection_csv = CSVReader(FileName=['Neuron_model_arrays/'+connection+'.csv'])


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
    
    ## change solid color
    ##tableToPoints1Display_1.DiffuseColor = [0.6196078431372549, 1.0, 0.12549019607843137]
    #if "Indirect" in connection:
    #    globals()[connection+'_tableDisplay_1'].DiffuseColor = [0.6196078431372549, 1.0, 0.12549019607843137]
    #elif "Direct" in connection:
    #    globals()[connection+'_tableDisplay_1'].DiffuseColor = [1.0, 0.0, 0.0]
    #elif "HDP" in connection:
    #    globals()[connection+'_tableDisplay_1'].DiffuseColor = [0.12156862745098039, 0.5843137254901961, 0.8941176470588236]
    #else:
    #    globals()[connection+'_tableDisplay_1'].DiffuseColor = [1.0,1.0,1.0]
    
    globals()[connection+'_tableDisplay_1'].DiffuseColor = list(viridis(index_population))[:3]
    index_population=index_population+1

#
## create a new 'CSV Reader'
#gPi_mask_Str_mask_5_Direct_Branch_MC_Str_GPicsv = CSVReader(FileName=['/home/konstantin/Documents/brukerMRI-master/Classified_by_Sosoho/Segment/GPi_mask_Str_mask_5_Direct_Branch_MC_Str_GPi.csv'])
#
## Properties modified on gPi_mask_Str_mask_5_Direct_Branch_MC_Str_GPicsv
#gPi_mask_Str_mask_5_Direct_Branch_MC_Str_GPicsv.HaveHeaders = 0
#gPi_mask_Str_mask_5_Direct_Branch_MC_Str_GPicsv.FieldDelimiterCharacters = ' '
#
## Create a new 'SpreadSheet View'
#spreadSheetView1 = CreateView('SpreadSheetView')
#spreadSheetView1.BlockSize = 1024L
## uncomment following to set a specific view size
## spreadSheetView1.ViewSize = [400, 400]
#
## place view in the layout
#viewLayout1.AssignView(2, spreadSheetView1)
#
## show data in view
#gPi_mask_Str_mask_5_Direct_Branch_MC_Str_GPicsvDisplay = Show(gPi_mask_Str_mask_5_Direct_Branch_MC_Str_GPicsv, spreadSheetView1)
## trace defaults for the display properties.
#gPi_mask_Str_mask_5_Direct_Branch_MC_Str_GPicsvDisplay.FieldAssociation = 'Row Data'
#
## create a new 'Table To Points'
#tableToPoints2 = TableToPoints(Input=gPi_mask_Str_mask_5_Direct_Branch_MC_Str_GPicsv)
#tableToPoints2.XColumn = 'Field 0'
#tableToPoints2.YColumn = 'Field 0'
#tableToPoints2.ZColumn = 'Field 0'
#
## Properties modified on tableToPoints2
#tableToPoints2.YColumn = 'Field 1'
#tableToPoints2.ZColumn = 'Field 2'
#
## show data in view
#tableToPoints2Display = Show(tableToPoints2, spreadSheetView1)
#
## hide data in view
#Hide(gPi_mask_Str_mask_5_Direct_Branch_MC_Str_GPicsv, spreadSheetView1)
#
## set active view
#SetActiveView(renderView1)
#
## set active source
#SetActiveSource(tableToPoints2)
#
## show data in view
#tableToPoints2Display_1 = Show(tableToPoints2, renderView1)
## trace defaults for the display properties.
#tableToPoints2Display_1.ColorArrayName = [None, '']
#tableToPoints2Display_1.GlyphType = 'Arrow'
#
## set active view
#SetActiveView(spreadSheetView1)
#
## destroy spreadSheetView1
#Delete(spreadSheetView1)
#del spreadSheetView1
#
## close an empty frame
#viewLayout1.Collapse(2)
#
## set active view
#SetActiveView(renderView1)
#
## change solid color
#tableToPoints2Display_1.DiffuseColor = [1.0, 0.0, 0.0]
#
## create a new 'CSV Reader'
#sTN_mask_MC_mask_35_HDP_Branch_STNcsv = CSVReader(FileName=['/home/konstantin/Documents/brukerMRI-master/Classified_by_Sosoho/Segment/STN_mask_MC_mask_35_HDP_Branch_STN.csv'])
#
## Properties modified on sTN_mask_MC_mask_35_HDP_Branch_STNcsv
#sTN_mask_MC_mask_35_HDP_Branch_STNcsv.HaveHeaders = 0
#sTN_mask_MC_mask_35_HDP_Branch_STNcsv.FieldDelimiterCharacters = ' '
#
## Create a new 'SpreadSheet View'
#spreadSheetView1 = CreateView('SpreadSheetView')
#spreadSheetView1.BlockSize = 1024L
## uncomment following to set a specific view size
## spreadSheetView1.ViewSize = [400, 400]
#
## place view in the layout
#viewLayout1.AssignView(2, spreadSheetView1)
#
## show data in view
#sTN_mask_MC_mask_35_HDP_Branch_STNcsvDisplay = Show(sTN_mask_MC_mask_35_HDP_Branch_STNcsv, spreadSheetView1)
## trace defaults for the display properties.
#sTN_mask_MC_mask_35_HDP_Branch_STNcsvDisplay.FieldAssociation = 'Row Data'
#
## create a new 'Table To Points'
#tableToPoints3 = TableToPoints(Input=sTN_mask_MC_mask_35_HDP_Branch_STNcsv)
#tableToPoints3.XColumn = 'Field 0'
#tableToPoints3.YColumn = 'Field 0'
#tableToPoints3.ZColumn = 'Field 0'
#
## Properties modified on tableToPoints3
#tableToPoints3.YColumn = 'Field 1'
#tableToPoints3.ZColumn = 'Field 2'
#
## show data in view
#tableToPoints3Display = Show(tableToPoints3, spreadSheetView1)
#
## hide data in view
#Hide(sTN_mask_MC_mask_35_HDP_Branch_STNcsv, spreadSheetView1)
#
## destroy spreadSheetView1
#Delete(spreadSheetView1)
#del spreadSheetView1
#
## close an empty frame
#viewLayout1.Collapse(2)
#
## set active view
#SetActiveView(renderView1)
#
## set active source
#SetActiveSource(tableToPoints3)
#
## show data in view
#tableToPoints3Display = Show(tableToPoints3, renderView1)
## trace defaults for the display properties.
#tableToPoints3Display.ColorArrayName = [None, '']
#tableToPoints3Display.GlyphType = 'Arrow'
#
## change solid color
#tableToPoints3Display.DiffuseColor = [0.12156862745098039, 0.5843137254901961, 0.8941176470588236]




#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [-5.943765067238889, 30.063310553609174, 19.481178764885886]
renderView1.CameraFocalPoint = [20.037838396112274, -8.840344561377364, -4.671394748401762]
renderView1.CameraViewUp = [-0.4484871730366652, 0.23572557572149833, -0.8621442504432473]
renderView1.CameraParallelScale = 24.899825806468712
renderView1.ViewSize = [1600, 1000]
renderView1.ResetCamera()
SaveScreenshot('Images/Axon_connections.png', magnification=1, quality=100, view=renderView1)
    
    #### uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).