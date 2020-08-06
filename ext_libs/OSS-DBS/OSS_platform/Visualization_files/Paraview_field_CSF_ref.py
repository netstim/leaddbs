#### import the simple module from the paraview

f_name="f_234"

from paraview.simple import *


# create a new 'PVD Reader'
import os
def make_a_screenshot(path_to_insert):
	# create a new 'PVD Reader'
	home_dir=os.path.expanduser("~")

	field_r_10pvd = PVDReader(FileName=[home_dir+path_to_insert+'/CSF_ref/Field_r_1.0.pvd'])
	#### disable automatic camera reset on 'Show'
	paraview.simple._DisableFirstRenderCameraReset()

	renderView1 = GetActiveViewOrCreate('RenderView')
	# uncomment following to set a specific view size
	# renderView1.ViewSize = [1164, 727]

	# get color transfer function/color map for f_name
	f344LUT = GetColorTransferFunction(f_name)

	# show data in view
	field_r_10pvdDisplay = Show(field_r_10pvd, renderView1)
	# trace defaults for the display properties.
	field_r_10pvdDisplay.ColorArrayName = ['POINTS', f_name]
	field_r_10pvdDisplay.LookupTable = f344LUT
	field_r_10pvdDisplay.GlyphType = 'Arrow'
	field_r_10pvdDisplay.ScalarOpacityUnitDistance = 0.40071001113096655

	# reset view to fit data
	renderView1.ResetCamera()

	# show color bar/color legend
	field_r_10pvdDisplay.SetScalarBarVisibility(renderView1, True)

	# get opacity transfer function/opacity map for 'f344'
	f344PWF = GetOpacityTransferFunction(f_name)

	# create a new 'Clip'
	clip1 = Clip(Input=field_r_10pvd)
	clip1.ClipType = 'Plane'
	clip1.Scalars = ['POINTS', f_name]
	clip1.Value = 0.5000000000000001

	# init the 'Plane' selected for 'ClipType'
	clip1.ClipType.Origin=[100.226087,135.58128000000002,22.247337]

	# Properties modified on clip1.ClipType
	clip1.ClipType.Origin=[100.226087,135.58128000000002,22.247337]

	# show data in view
	clip1Display = Show(clip1, renderView1)
	# trace defaults for the display properties.
	clip1Display.ColorArrayName = ['POINTS', f_name]
	clip1Display.LookupTable = f344LUT
	clip1Display.GlyphType = 'Arrow'
	clip1Display.ScalarOpacityUnitDistance = 0.4049629627814168

	# hide data in view
	Hide(field_r_10pvd, renderView1)

	# show color bar/color legend
	clip1Display.SetScalarBarVisibility(renderView1, True)

	# reset view to fit data
	renderView1.ResetCamera()

	# reset view to fit data
	renderView1.ResetCamera()

	# reset view to fit data
	renderView1.ResetCamera()

	# toggle 3D widget visibility (only when running from the GUI)
	Hide3DWidgets(proxy=clip1)


	# current camera placement for renderView1
	renderView1.CameraPosition = [-53.534594294245835, -0.7932281494140625, 2.0314996242523193]
	renderView1.CameraFocalPoint = [0.5645360946655273, -0.7932281494140625, 2.0314996242523193]
	renderView1.CameraViewUp = [0.0, 0.0, 1.0]
	renderView1.CameraParallelScale = 24.805193867502165

	# save screenshot
	renderView1.ViewSize = [1600, 1000]
	renderView1.ResetCamera()
	SaveScreenshot(home_dir+path_to_insert+'/Images/Field_Adapted_CSF.png', magnification=1, quality=100, view=renderView1)

	#### saving camera placements for all active views

	# current camera placement for renderView1
	renderView1.CameraPosition = [-53.534594294245835, -0.7932281494140625, 2.0314996242523193]
	renderView1.CameraFocalPoint = [0.5645360946655273, -0.7932281494140625, 2.0314996242523193]
	renderView1.CameraViewUp = [0.0, 0.0, 1.0]
	renderView1.CameraParallelScale = 24.805193867502165

	#### uncomment the following to render all views
	# RenderAllViews()
	# alternatively, if you want to write images, you can use SaveScreenshot(...).

if __name__ == '__main__':
	make_a_screenshot(*sys.argv[1:])
