#### import the simple module from the paraview
from paraview.simple import *
#from paraview_find_arrayname import get_Para_Array_name 
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
#import os
#terminal_path=os.getcwd() + "\n"
#cSF_Subdomains_unrefpvd = PVDReader(FileName=terminal_path+'/Results_adaptive/parallel_Subdomains.pvd')

import os
def make_a_screenshot(path_to_insert):
	# create a new 'PVD Reader'
	home_dir=os.path.expanduser("~")

	cSF_Subdomains_adapted = PVDReader(FileName=[home_dir+path_to_insert+'/Field_solutions/parallel_Subdomains.pvd'])


	# get active view
	renderView1 = GetActiveViewOrCreate('RenderView')
	# uncomment following to set a specific view size
	# renderView1.ViewSize = [1164, 808]

	# get color transfer function/color map for 'f'
	fLUT = GetColorTransferFunction('f')

	# show data in view
	cSF_Subdomains_unrefpvdDisplay = Show(cSF_Subdomains_adapted, renderView1)
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
	extractCellsByRegion1 = ExtractCellsByRegion(Input=cSF_Subdomains_adapted)
	extractCellsByRegion1.IntersectWith = 'Plane'

	# init the 'Plane' selected for 'IntersectWith'
	extractCellsByRegion1.IntersectWith.Origin=[100.226087,135.58128000000002,22.247337]

	# show data in view
	extractCellsByRegion1Display = Show(extractCellsByRegion1, renderView1)
	# trace defaults for the display properties.
	extractCellsByRegion1Display.ColorArrayName = ['CELLS', 'f']
	extractCellsByRegion1Display.LookupTable = fLUT
	extractCellsByRegion1Display.GlyphType = 'Arrow'
	extractCellsByRegion1Display.ScalarOpacityUnitDistance = 0.5033337160022229

	# hide data in view
	Hide(cSF_Subdomains_adapted, renderView1)

	# show color bar/color legend
	extractCellsByRegion1Display.SetScalarBarVisibility(renderView1, True)

	# Properties modified on extractCellsByRegion1.IntersectWith
	extractCellsByRegion1.IntersectWith.Normal = [-1.0, 0.0, 0.0]


	# show color bar/color legend
	extractCellsByRegion1Display.SetScalarBarVisibility(renderView1, True)

	# update the view to ensure updated data information
	renderView1.Update()

	# Properties modified on fLUT
	fLUT.RGBPoints = [1.0, 0.231373, 0.298039, 0.752941, 1.4484847784042358, 0.36470588235294116, 0.48627450980392156, 0.9019607843137255, 3.0, 0.865003, 0.865003, 0.865003, 5.0, 0.705882, 0.0156863, 0.14902]

	# Properties modified on fLUT
	fLUT.RGBPoints = [1.0, 0.2823529411764706, 0.6039215686274509, 0.8509803921568627, 1.4484847784042358, 0.36470588235294116, 0.48627450980392156, 0.9019607843137255, 3.0, 0.865003, 0.865003, 0.865003, 5.0, 0.705882, 0.0156863, 0.14902]

	# Properties modified on fLUT
	fLUT.RGBPoints = [1.0, 0.06274509803921569, 0.6039215686274509, 0.9607843137254902, 1.4484847784042358, 0.36470588235294116, 0.48627450980392156, 0.9019607843137255, 3.0, 0.865003, 0.865003, 0.865003, 5.0, 0.705882, 0.0156863, 0.14902]

	# Properties modified on fLUT
	fLUT.RGBPoints = [1.0, 0.06274509803921569, 0.6039215686274509, 0.9607843137254902, 3.0, 0.865003, 0.865003, 0.865003, 5.0, 0.705882, 0.0156863, 0.14902]

	# Properties modified on fLUT
	fLUT.RGBPoints = [1.0, 0.06274509803921569, 0.6039215686274509, 0.9607843137254902, 3.0, 0.4980392156862745, 0.5019607843137255, 0.5725490196078431, 5.0, 0.705882, 0.0156863, 0.14902]

	# Properties modified on fLUT
	fLUT.RGBPoints = [1.0, 0.06274509803921569, 0.6039215686274509, 0.9607843137254902, 1.860606074333191, 0.4745098039215686, 0.5882352941176471, 0.8156862745098039, 3.0, 0.4980392156862745, 0.5019607843137255, 0.5725490196078431, 5.0, 0.705882, 0.0156863, 0.14902]

	# Properties modified on fLUT
	fLUT.RGBPoints = [1.0, 0.06274509803921569, 0.6039215686274509, 0.9607843137254902, 2.0, 0.4745098039215686, 0.5882352941176471, 0.8156862745098039, 3.0, 0.4980392156862745, 0.5019607843137255, 0.5725490196078431, 5.0, 0.705882, 0.0156863, 0.14902]

	# Properties modified on fLUT
	fLUT.RGBPoints = [1.0, 0.06274509803921569, 0.6039215686274509, 0.9607843137254902, 2.0, 1.0, 1.0, 1.0, 3.0, 0.4980392156862745, 0.5019607843137255, 0.5725490196078431, 5.0, 0.705882, 0.0156863, 0.14902]

	# Properties modified on fLUT
	fLUT.RGBPoints = [1.0, 0.06274509803921569, 0.6039215686274509, 0.9607843137254902, 2.0, 1.0, 1.0, 1.0, 3.0, 0.4980392156862745, 0.5019607843137255, 0.5725490196078431, 4.090909004211426, 0.9098039215686274, 0.8392156862745098, 0.8, 5.0, 0.705882, 0.0156863, 0.14902]

	# Properties modified on fLUT
	fLUT.RGBPoints = [1.0, 0.06274509803921569, 0.6039215686274509, 0.9607843137254902, 2.0, 1.0, 1.0, 1.0, 3.0, 0.4980392156862745, 0.5019607843137255, 0.5725490196078431, 4.0, 0.9098039215686274, 0.8392156862745098, 0.8, 5.0, 0.705882, 0.0156863, 0.14902]

	# Properties modified on fLUT
	fLUT.RGBPoints = [1.0, 0.06274509803921569, 0.6039215686274509, 0.9607843137254902, 2.0, 1.0, 1.0, 1.0, 3.0, 0.4980392156862745, 0.5019607843137255, 0.5725490196078431, 4.0, 0.9098039215686274, 0.0, 0.0, 5.0, 0.705882, 0.0156863, 0.14902]

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
	SaveScreenshot(home_dir+path_to_insert+'/Images/Adapted_mesh.png', magnification=1, quality=100, view=renderView1)

	#### saving camera placements for all active views

	# current camera placement for renderView1
	renderView1.CameraPosition = [-35.61998003697748, 1.1129615306854248, 1.0326027870178223]
	renderView1.CameraFocalPoint = [1.0243473052978516, 1.1129615306854248, 1.0326027870178223]
	renderView1.CameraViewUp = [0.0, 0.0, 1.0]
	renderView1.CameraParallelScale = 9.484249811151892

	#### uncomment the following to render all views
	# RenderAllViews()
	# alternatively, if you want to write images, you can use SaveScreenshot(...).

if __name__ == '__main__':
	make_a_screenshot(*sys.argv[1:])
