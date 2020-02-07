import vtk
import os
import slicer
import SimpleITK as sitk
import sitkUtils
import numpy as np
import glob

try:
  import h5py
except:
  from pip._internal import main as pipmain
  failed = pipmain(['install', 'h5py'])
  if not failed:
    import h5py


def getTransformRASToIJK(transformNode):
  origin = transformNode.GetTransformFromParent().GetDisplacementGrid().GetOrigin()
  spacing = transformNode.GetTransformFromParent().GetDisplacementGrid().GetSpacing()
  IJKToRAS = [ 
              [spacing[0],          0 ,         0 ,  origin[0] ],
              [        0 ,  spacing[1],         0 ,  origin[1] ],
              [        0 ,          0 , spacing[2],  origin[2] ],
              [        0 ,          0 ,         0 ,         1  ],
              ] 
  m = vtk.vtkMatrix4x4()
  for i in range(4):
    for j in range(4):
      m.SetElement(i,j,IJKToRAS[i][j])
  m.Invert()
  return m


def createSegmentationFromAtlas(atlasPath):
  """
  Used to create Slicer segmentation from a Lead-DBS atlas folder
  """
  try:
    h5py
  except:
    msg=qt.QMessageBox().warning(qt.QWidget(),'','h5py module is not installed')
    return None
  
  # get variables stored in atlas_index.mat
  with h5py.File(os.path.join(atlasPath,'atlas_index.mat'),'r') as f:
    # get threshold
    threshold = f['atlases']['threshold']['value'][()][0]
    # get color map
    try:
      colormap = f['atlases']['colormap'][()].transpose()
    except: # colormap not present
      colormap = np.array([[ 0, 0,0.5156],[ 0, 0,0.5312],[ 0, 0,0.5469],[ 0, 0,0.5625],[ 0, 0,0.5781],[ 0, 0,0.5938],[ 0, 0,0.6094],[ 0, 0,0.6250],[ 0, 0,0.6406],[ 0, 0,0.6562],[ 0, 0,0.6719],[ 0, 0,0.6875],[ 0, 0,0.7031],[ 0, 0,0.7188],[ 0, 0,0.7344],[ 0, 0,0.7500],[ 0, 0,0.7656],[ 0, 0,0.7812],[ 0, 0,0.7969],[ 0, 0,0.8125],[ 0, 0,0.8281],[ 0, 0,0.8438],[ 0, 0,0.8594],[ 0, 0,0.8750],[ 0, 0,0.8906],[ 0, 0,0.9062],[ 0, 0,0.9219],[ 0, 0,0.9375],[ 0, 0,0.9531],[ 0, 0,0.9688],[ 0, 0,0.9844],[ 0, 0,1.0000],[ 0,0.0156,1.0000],[ 0,0.0312,1.0000],[ 0,0.0469,1.0000],[ 0,0.0625,1.0000],[ 0,0.0781,1.0000],[ 0,0.0938,1.0000],[ 0,0.1094,1.0000],[ 0,0.1250,1.0000],[ 0,0.1406,1.0000],[ 0,0.1562,1.0000],[ 0,0.1719,1.0000],[ 0,0.1875,1.0000],[ 0,0.2031,1.0000],[ 0,0.2188,1.0000],[ 0,0.2344,1.0000],[ 0,0.2500,1.0000],[ 0,0.2656,1.0000],[ 0,0.2812,1.0000],[ 0,0.2969,1.0000],[ 0,0.3125,1.0000],[ 0,0.3281,1.0000],[ 0,0.3438,1.0000],[ 0,0.3594,1.0000],[ 0,0.3750,1.0000],[ 0,0.3906,1.0000],[ 0,0.4062,1.0000],[ 0,0.4219,1.0000],[ 0,0.4375,1.0000],[ 0,0.4531,1.0000],[ 0,0.4688,1.0000],[ 0,0.4844,1.0000],[ 0,0.5000,1.0000],[ 0,0.5156,1.0000],[ 0,0.5312,1.0000],[ 0,0.5469,1.0000],[ 0,0.5625,1.0000],[ 0,0.5781,1.0000],[ 0,0.5938,1.0000],[ 0,0.6094,1.0000],[ 0,0.6250,1.0000],[ 0,0.6406,1.0000],[ 0,0.6562,1.0000],[ 0,0.6719,1.0000],[ 0,0.6875,1.0000],[ 0,0.7031,1.0000],[ 0,0.7188,1.0000],[ 0,0.7344,1.0000],[ 0,0.7500,1.0000],[ 0,0.7656,1.0000],[ 0,0.7812,1.0000],[ 0,0.7969,1.0000],[ 0,0.8125,1.0000],[ 0,0.8281,1.0000],[ 0,0.8438,1.0000],[ 0,0.8594,1.0000],[ 0,0.8750,1.0000],[ 0,0.8906,1.0000],[ 0,0.9062,1.0000],[ 0,0.9219,1.0000],[ 0,0.9375,1.0000],[ 0,0.9531,1.0000],[ 0,0.9688,1.0000],[ 0,0.9844,1.0000],[ 0,1.0000,1.0000],[0.0156,1.0000,0.9844],[0.0312,1.0000,0.9688],[0.0469,1.0000,0.9531],[0.0625,1.0000,0.9375],[0.0781,1.0000,0.9219],[0.0938,1.0000,0.9062],[0.1094,1.0000,0.8906],[0.1250,1.0000,0.8750],[0.1406,1.0000,0.8594],[0.1562,1.0000,0.8438],[0.1719,1.0000,0.8281],[0.1875,1.0000,0.8125],[0.2031,1.0000,0.7969],[0.2188,1.0000,0.7812],[0.2344,1.0000,0.7656],[0.2500,1.0000,0.7500],[0.2656,1.0000,0.7344],[0.2812,1.0000,0.7188],[0.2969,1.0000,0.7031],[0.3125,1.0000,0.6875],[0.3281,1.0000,0.6719],[0.3438,1.0000,0.6562],[0.3594,1.0000,0.6406],[0.3750,1.0000,0.6250],[0.3906,1.0000,0.6094],[0.4062,1.0000,0.5938],[0.4219,1.0000,0.5781],[0.4375,1.0000,0.5625],[0.4531,1.0000,0.5469],[0.4688,1.0000,0.5312],[0.4844,1.0000,0.5156],[0.5000,1.0000,0.5000],[0.5156,1.0000,0.4844],[0.5312,1.0000,0.4688],[0.5469,1.0000,0.4531],[0.5625,1.0000,0.4375],[0.5781,1.0000,0.4219],[0.5938,1.0000,0.4062],[0.6094,1.0000,0.3906],[0.6250,1.0000,0.3750],[0.6406,1.0000,0.3594],[0.6562,1.0000,0.3438],[0.6719,1.0000,0.3281],[0.6875,1.0000,0.3125],[0.7031,1.0000,0.2969],[0.7188,1.0000,0.2812],[0.7344,1.0000,0.2656],[0.7500,1.0000,0.2500],[0.7656,1.0000,0.2344],[0.7812,1.0000,0.2188],[0.7969,1.0000,0.2031],[0.8125,1.0000,0.1875],[0.8281,1.0000,0.1719],[0.8438,1.0000,0.1562],[0.8594,1.0000,0.1406],[0.8750,1.0000,0.1250],[0.8906,1.0000,0.1094],[0.9062,1.0000,0.0938],[0.9219,1.0000,0.0781],[0.9375,1.0000,0.0625],[0.9531,1.0000,0.0469],[0.9688,1.0000,0.0312],[0.9844,1.0000,0.0156],[1.0000,1.0000, 0],[1.0000,0.9844, 0],[1.0000,0.9688, 0],[1.0000,0.9531, 0],[1.0000,0.9375, 0],[1.0000,0.9219, 0],[1.0000,0.9062, 0],[1.0000,0.8906, 0],[1.0000,0.8750, 0],[1.0000,0.8594, 0],[1.0000,0.8438, 0],[1.0000,0.8281, 0],[1.0000,0.8125, 0],[1.0000,0.7969, 0],[1.0000,0.7812, 0],[1.0000,0.7656, 0],[1.0000,0.7500, 0],[1.0000,0.7344, 0],[1.0000,0.7188, 0],[1.0000,0.7031, 0],[1.0000,0.6875, 0],[1.0000,0.6719, 0],[1.0000,0.6562, 0],[1.0000,0.6406, 0],[1.0000,0.6250, 0],[1.0000,0.6094, 0],[1.0000,0.5938, 0],[1.0000,0.5781, 0],[1.0000,0.5625, 0],[1.0000,0.5469, 0],[1.0000,0.5312, 0],[1.0000,0.5156, 0],[1.0000,0.5000, 0],[1.0000,0.4844, 0],[1.0000,0.4688, 0],[1.0000,0.4531, 0],[1.0000,0.4375, 0],[1.0000,0.4219, 0],[1.0000,0.4062, 0],[1.0000,0.3906, 0],[1.0000,0.3750, 0],[1.0000,0.3594, 0],[1.0000,0.3438, 0],[1.0000,0.3281, 0],[1.0000,0.3125, 0],[1.0000,0.2969, 0],[1.0000,0.2812, 0],[1.0000,0.2656, 0],[1.0000,0.2500, 0],[1.0000,0.2344, 0],[1.0000,0.2188, 0],[1.0000,0.2031, 0],[1.0000,0.1875, 0],[1.0000,0.1719, 0],[1.0000,0.1562, 0],[1.0000,0.1406, 0],[1.0000,0.1250, 0],[1.0000,0.1094, 0],[1.0000,0.0938, 0],[1.0000,0.0781, 0],[1.0000,0.0625, 0],[1.0000,0.0469, 0],[1.0000,0.0312, 0],[1.0000,0.0156, 0],[1.0000, 0, 0],[0.9844, 0, 0],[0.9688, 0, 0],[0.9531, 0, 0],[0.9375, 0, 0],[0.9219, 0, 0],[0.9062, 0, 0],[0.8906, 0, 0],[0.8750, 0, 0],[0.8594, 0, 0],[0.8438, 0, 0],[0.8281, 0, 0],[0.8125, 0, 0],[0.7969, 0, 0],[0.7812, 0, 0],[0.7656, 0, 0],[0.7500, 0, 0],[0.7344, 0, 0],[0.7188, 0, 0],[0.7031, 0, 0],[0.6875, 0, 0],[0.6719, 0, 0],[0.6562, 0, 0],[0.6406, 0, 0],[0.6250, 0, 0],[0.6094, 0, 0],[0.5938, 0, 0],[0.5781, 0, 0],[0.5625, 0, 0],[0.5469, 0, 0],[0.5312, 0, 0],[0.5156, 0, 0],[0.5000, 0, 0]])
    # get colors
    colors = f['atlases']['colors'][()]
    # get names
    names = []
    atlases = f['atlases']
    for column in atlases['names']:
      row_data = []
      for row_number in range(len(column)):      
        row_data.append(''.join(map(unichr, atlases[column[row_number]][:])))   
      names.append(row_data)
  
  
  # init aux nodes
  labelNode = slicer.mrmlScene.AddNode(slicer.vtkMRMLLabelMapVolumeNode())
  segmentationNode = slicer.vtkMRMLSegmentationNode()
  slicer.mrmlScene.AddNode(segmentationNode)
  segmentationLogic = slicer.modules.segmentations.logic()
  
  # get nifti structures
  listing = sorted(glob.glob(os.path.join(atlasPath,'*','*.nii*')))
  if len(listing) == 0:
    return None
  
  for filename in listing:
    path, structure = os.path.split(filename)
    base = os.path.basename(path) # side: lh, rh, mixed
    structureName = os.path.splitext(os.path.splitext(structure)[0])[0]
    structureFullName = structureName + '_' + base if base != "mixed" else structureName
    # load volume
    loadSucces, volumeNode = slicer.util.loadVolume(filename, properties = {'name':'tmpVol', 'show':False}, returnNode = True)
    labelNode.SetName(structureFullName) # so the segment id is the same as name
    # cast volume to float
    if volumeNode.GetImageData().GetScalarTypeAsString() not in ['float','float']:
      sitkFilter = sitk.CastImageFilter()
      inputImage = sitk.ReadImage(sitkUtils.GetSlicerITKReadWriteAddress(volumeNode.GetName()))
      outputImage = sitkFilter.Execute(inputImage)
      nodeWriteAddress = sitkUtils.GetSlicerITKReadWriteAddress(volumeNode.GetName())
      sitk.WriteImage(outputImage, nodeWriteAddress)
    # run simpleITK threshold 
    sitkFilter = sitk.BinaryThresholdImageFilter()
    sitkFilter.SetLowerThreshold(float(threshold))
    inputImage = sitk.ReadImage(sitkUtils.GetSlicerITKReadWriteAddress(volumeNode.GetName()))
    outputImage = sitkFilter.Execute(inputImage)
    nodeWriteAddress = sitkUtils.GetSlicerITKReadWriteAddress(labelNode.GetName())
    sitk.WriteImage(outputImage, nodeWriteAddress)
    # add to segmentation
    index = int(colors[names.index([structure])] - 1) # index of colormap. -1 because its matlab index
    segmentationLogic.ImportLabelmapToSegmentationNode(labelNode, segmentationNode)
    segmentationNode.GetSegmentation().GetNthSegment(segmentationNode.GetSegmentation().GetNumberOfSegments()-1).SetColor(*colormap[index])
    segmentationNode.GetSegmentation().GetNthSegment(segmentationNode.GetSegmentation().GetNumberOfSegments()-1).SetName(structureFullName)
    slicer.mrmlScene.RemoveNode(volumeNode)
  
  # reorder
  count=0
  for n in names:
    struture = os.path.splitext(os.path.splitext(n[0])[0])[0]
    for base in ['lh','rh']:
      n = struture + '_' + base
      segIndex = segmentationNode.GetSegmentation().GetSegmentIndex(n)
      if segIndex > -1:
        segmentationNode.GetSegmentation().SetSegmentIndex(n, count)
        count += 1
  

  slicer.mrmlScene.RemoveNode(labelNode)
  
  return segmentationNode


def emptyTransform(transformNode):
  # init
  transformSize = [193,229,193]
  voxelType=vtk.VTK_FLOAT
  transformOrigin = [-96.0, -132.0, -78.0]
  transformSpacing = [1.0, 1.0, 1.0]
  transformPathections = [[1,0,0], [0,1,0], [0,0,1]]
  fillVoxelValue = 0
  # Create an empty image volume, filled with fillVoxelValue
  imageData = vtk.vtkImageData()
  imageData.SetDimensions(transformSize)
  imageData.AllocateScalars(voxelType, 3)
  imageData.GetPointData().GetScalars().Fill(fillVoxelValue)
  # Create transform
  transform = slicer.vtkOrientedGridTransform()
  transform.SetInterpolationModeToCubic()
  transform.SetDisplacementGridData(imageData)
  # Create transform node
  transformNode.SetAndObserveTransformFromParent(transform)
  transformNode.GetTransformFromParent().GetDisplacementGrid().SetOrigin(transformOrigin)
  transformNode.GetTransformFromParent().GetDisplacementGrid().SetSpacing(transformSpacing)
  #transformNode.CreateDefaultDisplayNodes()
  transformNode.CreateDefaultStorageNode()  