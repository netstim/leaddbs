#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
cSF_Subdomains_unrefpvd = PVDReader(FileName='CSF_ref/CSF_Subdomains_full_ref.pvd')

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1164, 808]

# get color transfer function/color map for 'f'
fLUT = GetColorTransferFunction('f')

# show data in view
cSF_Subdomains_unrefpvdDisplay = Show(cSF_Subdomains_unrefpvd, renderView1)
# trace defaults for the display properties.
cSF_Subdomains_unrefpvdDisplay.ColorArrayName = ['CELLS', 'f']
cSF_Subdomains_unrefpvdDisplay.LookupTable = fLUT
cSF_Subdomains_unrefpvdDisplay.GlyphType = 'Arrow'
cSF_Subdomains_unrefpvdDisplay.ScalarOpacityUnitDistance = 0.42132744935143623

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
cSF_Subdomains_unrefpvdDisplay.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'f'
fPWF = GetOpacityTransferFunction('f')

# create a new 'Extract Cells By Region'
extractCellsByRegion1 = ExtractCellsByRegion(Input=cSF_Subdomains_unrefpvd)
extractCellsByRegion1.IntersectWith = 'Plane'

# init the 'Plane' selected for 'IntersectWith'
extractCellsByRegion1.IntersectWith.Origin=[106.929,119.883,66.303]

# show data in view
extractCellsByRegion1Display = Show(extractCellsByRegion1, renderView1)
# trace defaults for the display properties.
extractCellsByRegion1Display.ColorArrayName = ['CELLS', 'f']
extractCellsByRegion1Display.LookupTable = fLUT
extractCellsByRegion1Display.GlyphType = 'Arrow'
extractCellsByRegion1Display.ScalarOpacityUnitDistance = 0.5033337160022229

# hide data in view
Hide(cSF_Subdomains_unrefpvd, renderView1)

# show color bar/color legend
extractCellsByRegion1Display.SetScalarBarVisibility(renderView1, True)

# Properties modified on extractCellsByRegion1.IntersectWith
extractCellsByRegion1.IntersectWith.Normal = [-1.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on fLUT
fLUT.RGBPoints = [1.0, 0.06274509803921569, 0.6039215686274509, 0.9607843137254902, 2.0, 1.0, 1.0, 1.0, 3.0, 0.4980392156862745, 0.5019607843137255, 0.5725490196078431, 4.0, 0.9098039215686274, 0.0, 0.0, 5.0, 0.10588235294117647, 0.7058823529411765, 0.0]


# reset view to fit data
renderView1.ResetCamera()

# current camera placement for renderView1
renderView1.CameraPosition = [-35.61998003697748, 1.1129615306854248, 1.0326027870178223]
renderView1.CameraFocalPoint = [1.0243473052978516, 1.1129615306854248, 1.0326027870178223]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 9.484249811151892

# save screenshot
renderView1.ViewSize = [1600, 1000]
renderView1.ResetCamera()
SaveScreenshot('Images/CSF_full_refinement.png', magnification=1, quality=100, view=renderView1)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [-35.61998003697748, 1.1129615306854248, 1.0326027870178223]
renderView1.CameraFocalPoint = [1.0243473052978516, 1.1129615306854248, 1.0326027870178223]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 9.484249811151892

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
