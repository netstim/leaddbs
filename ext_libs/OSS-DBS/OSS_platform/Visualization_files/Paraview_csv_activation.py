#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
neuron_model_resultscsv = CSVReader(FileName=['Field_solutions/Activation/Neuron_model_results.csv'])

# Properties modified on neuron_model_resultscsv
neuron_model_resultscsv.HaveHeaders = 0
neuron_model_resultscsv.FieldDelimiterCharacters = ' '

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# get layout
layout1 = GetLayout()

# place view in the layout
layout1.AssignView(2, spreadSheetView1)

# show data in view
neuron_model_resultscsvDisplay = Show(neuron_model_resultscsv, spreadSheetView1)
# trace defaults for the display properties.
neuron_model_resultscsvDisplay.FieldAssociation = 'Row Data'

# find view
renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [518, 724]

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
spreadSheetView1.Update()

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=neuron_model_resultscsv)
tableToPoints1.XColumn = 'Field 0'
tableToPoints1.YColumn = 'Field 0'
tableToPoints1.ZColumn = 'Field 0'

# Properties modified on tableToPoints1
tableToPoints1.YColumn = 'Field 1'
tableToPoints1.ZColumn = 'Field 2'

# show data in view
tableToPoints1Display = Show(tableToPoints1, spreadSheetView1)

# hide data in view
Hide(neuron_model_resultscsv, spreadSheetView1)

# update the view to ensure updated data information
spreadSheetView1.Update()

# set active view
SetActiveView(renderView1)

# set active view
SetActiveView(spreadSheetView1)

# destroy spreadSheetView1
Delete(spreadSheetView1)
del spreadSheetView1

# close an empty frame
layout1.Collapse(2)

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(tableToPoints1)

# show data in view
tableToPoints1Display = Show(tableToPoints1, renderView1)
# trace defaults for the display properties.
tableToPoints1Display.Representation = 'Surface'
tableToPoints1Display.ColorArrayName = [None, '']
tableToPoints1Display.OSPRayScaleArray = 'Field 3'
tableToPoints1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tableToPoints1Display.SelectOrientationVectors = 'Field 3'
tableToPoints1Display.ScaleFactor = 0.3399610000000003
tableToPoints1Display.SelectScaleArray = 'Field 3'
tableToPoints1Display.GlyphType = 'Arrow'
tableToPoints1Display.GlyphTableIndexArray = 'Field 3'
tableToPoints1Display.DataAxesGrid = 'GridAxesRepresentation'
tableToPoints1Display.PolarAxes = 'PolarAxesRepresentation'

# reset view to fit data
renderView1.ResetCamera()

# set scalar coloring
ColorBy(tableToPoints1Display, ('POINTS', 'Field 3'))

# rescale color and/or opacity maps used to include current data range
tableToPoints1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
tableToPoints1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Field3'
field3LUT = GetColorTransferFunction('Field3')

# reset view to fit data
renderView1.ResetCamera()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [11.8661065, 11.502742601993386, 8.924483]
renderView1.CameraFocalPoint = [11.8661065, 21.668008, 8.924483]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 2.6309642835257674

renderView1.ViewSize = [1600, 1000]
renderView1.ResetCamera()
SaveScreenshot('Images/Activated_neurons.png', magnification=1, quality=100, view=renderView1)

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).