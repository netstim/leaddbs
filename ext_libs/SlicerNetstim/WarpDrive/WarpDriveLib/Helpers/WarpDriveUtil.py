import qt, vtk, slicer
import SimpleITK as sitk
import sitkUtils

import numpy as np
from scipy import ndimage


from . import GridNodeHelper



def delay(msecs = 500):
  # delay application
  dieTime = qt.QTime().currentTime().addMSecs(msecs)
  while qt.QTime().currentTime() < dieTime:
    qt.QCoreApplication.processEvents(qt.QEventLoop.AllEvents, 100)

def previewWarp(sourceNode, targetNode, outNode):
  # points
  sourcePoints = vtk.vtkPoints()
  sourceNode.GetControlPointPositionsWorld(sourcePoints)
  targetPoints = vtk.vtkPoints()
  targetNode.GetControlPointPositionsWorld(targetPoints)
  # thin plate
  transform=vtk.vtkThinPlateSplineTransform()
  transform.SetSourceLandmarks(sourcePoints)
  transform.SetTargetLandmarks(targetPoints)
  transform.SetBasisToR()
  transform.Inverse()
  transformNode=slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
  transformNode.SetAndObserveTransformFromParent(transform)
  # display
  transformNode.CreateDefaultDisplayNodes()
  transformNode.GetDisplayNode().SetVisibility(1)
  transformNode.GetDisplayNode().SetVisibility3D(0)
  transformNode.GetDisplayNode().SetAndObserveGlyphPointsNode(sourceNode)
  transformNode.GetDisplayNode().SetVisibility2D(1)
  delay() # update display
  return transformNode


def createFolderDisplayNode(folderID, color=[0.66,0.66,0.66]):
  # from qSlicerSubjectHierarchyFolderPlugin.cxx
  shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
  displayNode = slicer.vtkMRMLFolderDisplayNode()
  displayNode.SetName(shNode.GetItemName(folderID))
  displayNode.SetHideFromEditors(0)
  displayNode.SetAttribute('SubjectHierarchy.Folder', "1")
  displayNode.SetColor(*color)
  shNode.GetScene().AddNode(displayNode)
  shNode.SetItemDataNode(folderID, displayNode)
  shNode.ItemModified(folderID)


def getPointsFromAttribute(attributeName):
  # iterate over markups fiducial and get points from nodes with specified attribute
  points = vtk.vtkPoints() # init
  shNode = slicer.mrmlScene.GetSubjectHierarchyNode() # subject hierarchy
  nMarkups = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLMarkupsFiducialNode')
  for i in range(nMarkups): # iterate
    markupNode = slicer.mrmlScene.GetNthNodeByClass(i, 'vtkMRMLMarkupsFiducialNode')
    if attributeName in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(markupNode)) and shNode.GetItemDataNode(shNode.GetItemParent(shNode.GetItemByDataNode(markupNode))).GetDescription() == 'Enabled':
      # get points of current node
      p = vtk.vtkPoints()
      markupNode.GetControlPointPositionsWorld(p)
      # add to output points
      points.InsertPoints(points.GetNumberOfPoints(), p.GetNumberOfPoints(), 0, p)
  # return
  return points

def getMaxSpread():
  # iterate over labelmaps created in the module and get max spread
  maxSpread = 0
  shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
  nLabelMaps = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLLabelMapVolumeNode') 
  for i in range(nLabelMaps):
    labelMapNode = slicer.mrmlScene.GetNthNodeByClass(i, 'vtkMRMLLabelMapVolumeNode')
    if 'correction' in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(labelMapNode)) and shNode.GetItemDataNode(shNode.GetItemParent(shNode.GetItemByDataNode(labelMapNode))).GetDescription() == 'Enabled':
      if int(labelMapNode.GetDescription()) > maxSpread:
        maxSpread = int(labelMapNode.GetDescription())
  
  # set default if none
  if maxSpread == 0:
    maxSpread = 30 

  return maxSpread 


def getDistanceMap(labelMapNode):
  # labelmap to scalar
  outNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLScalarVolumeNode')
  slicer.modules.volumes.logic().CreateScalarVolumeFromVolume(slicer.mrmlScene, outNode, labelMapNode)
  # cast type to int
  imageCast = vtk.vtkImageCast()
  imageCast.SetInputData(outNode.GetImageData())
  imageCast.SetOutputScalarTypeToInt()
  imageCast.Update()
  outNode.SetAndObserveImageData(imageCast.GetOutput())
  # replace values above 1 to 1
  a = slicer.util.array(outNode.GetID())
  a[a>0] = 1
  outNode.Modified()
  # filter
  myFilter = sitk.ApproximateSignedDistanceMapImageFilter()
  myFilter.SetDebug(False)
  myFilter.SetInsideValue(1.0)
  myFilter.SetOutsideValue(0.0)
  inputImage = sitkUtils.PullVolumeFromSlicer(outNode)
  outputImage = myFilter.Execute(inputImage)
  sitkUtils.PushVolumeToSlicer(outputImage, outNode)
  return outNode

def getMaskVolume(referenceVolume):
  shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
  # create segmentation node
  segmentationNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode")
  segmentationNode.CreateDefaultDisplayNodes()
  segmentationNode.GetDisplayNode().SetVisibility(False)
  # add masks
  nLabelMaps = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLLabelMapVolumeNode') 
  for i in range(nLabelMaps):
    labelMapNode = slicer.mrmlScene.GetNthNodeByClass(i, 'vtkMRMLLabelMapVolumeNode')
    if 'correction' in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(labelMapNode)):
      slicer.modules.segmentations.logic().ImportLabelmapToSegmentationNode(labelMapNode, segmentationNode)
  # if no segment return zeros
  if not segmentationNode.GetSegmentation().GetNumberOfSegments():
    slicer.mrmlScene.RemoveNode(segmentationNode)
    scalarVolumeNode = GridNodeHelper.emptyVolume(referenceVolume.GetImageData().GetDimensions(), referenceVolume.GetOrigin(), referenceVolume.GetSpacing())
    return scalarVolumeNode
  # segmentation to label
  labelMapNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLabelMapVolumeNode")
  slicer.modules.segmentations.logic().ExportSegmentsToLabelmapNode(segmentationNode, vtk.vtkStringArray(), labelMapNode, referenceVolume)
  slicer.mrmlScene.RemoveNode(segmentationNode)
  # distance filter
  distanceFilterNode = getDistanceMap(labelMapNode)
  slicer.mrmlScene.RemoveNode(labelMapNode)
  distanceFilterArray = slicer.util.array(distanceFilterNode.GetID())
  distanceFilterArray[:] = ndimage.gaussian_filter( 1 - (np.tanh(distanceFilterArray)+1)/2  , 1)
  distanceFilterNode.Modified()
  return distanceFilterNode



def generateSphereModel(center, radius):
  # vtk sphere
  sphere = vtk.vtkSphereSource()
  sphere.SetCenter(center)
  sphere.SetRadius(radius)
  sphere.SetPhiResolution(30)
  sphere.SetThetaResolution(30)
  sphere.Update()
  # to model
  modelsLogic = slicer.modules.models.logic()
  model = modelsLogic.AddModel(sphere.GetOutput())
  return model

def modelToLabelmap(modelNode, referenceNode):
  # reference
  size,origin,spacing = GridNodeHelper.getGridDefinition(referenceNode)
  auxVolume = GridNodeHelper.emptyVolume(size, origin, spacing)
  # output
  outNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLabelMapVolumeNode")
  # cli
  cliParams = {
    "labelValue" : 1,
    "InputVolume" : auxVolume.GetID(),
    "surface" : modelNode.GetID(),
    "OutputVolume" : outNode.GetID(),
    } 
  cliNode = slicer.cli.run(slicer.modules.modeltolabelmap, None, cliParams, wait_for_completion=True, update_display=False)
  # delete
  slicer.mrmlScene.RemoveNode(auxVolume)
  # return
  return outNode

def fiducialToModel(fiducialNode):
  modelNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode')
  radius = 1
  # see function defenition for inputs
  slicer.modules.markupstomodel.logic().UpdateOutputCurveModel(fiducialNode, modelNode, \
          slicer.vtkMRMLMarkupsToModelNode.CardinalSpline,
          False, radius, 8, 5, True, 3, \
          slicer.vtkMRMLMarkupsToModelNode.RawIndices, \
          None, \
          slicer.vtkMRMLMarkupsToModelNode.GlobalLeastSquares, \
          0.5, \
          slicer.vtkMRMLMarkupsToModelNode.Rectangular)
  return modelNode

def dilate(labelmapNode, spread):
  outNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLabelMapVolumeNode")
  # kernel radius
  kernel_radius = int(spread / labelmapNode.GetSpacing()[0]) # TODO: anisotropic
  # sitk filter
  myFilter = sitk.BinaryDilateImageFilter()
  myFilter.SetBackgroundValue(0.0)
  myFilter.SetBoundaryToForeground(False)
  myFilter.SetDebug(False)
  myFilter.SetForegroundValue(1.0)
  myFilter.SetKernelRadius((kernel_radius, kernel_radius, kernel_radius))
  myFilter.SetKernelType(1)
  inputImage = sitkUtils.PullVolumeFromSlicer(labelmapNode)
  outputImage = myFilter.Execute(inputImage)
  sitkUtils.PushVolumeToSlicer(outputImage, outNode)

  return outNode


def generateMaskFromPointsSpread(fiducialNode, spread, referenceNode):

  if fiducialNode.GetNumberOfControlPoints() == 1: # one point -> create sphere
    center = [0] * 3
    fiducialNode.GetNthControlPointPositionWorld(0,center)
    modelNode = generateSphereModel(center, spread)
    labelmapNode = modelToLabelmap(modelNode, referenceNode)
    # remove
    slicer.mrmlScene.RemoveNode(modelNode)
  else:
    # create skeleton model from fiducials
    skeletonModel = fiducialToModel(fiducialNode)
    # model to label
    skeletonLabelmap = modelToLabelmap(skeletonModel, referenceNode)
    # dilate
    labelmapNode = dilate(skeletonLabelmap, spread)
    # remove
    slicer.mrmlScene.RemoveNode(skeletonModel)
    slicer.mrmlScene.RemoveNode(skeletonLabelmap)
  
  return labelmapNode


def addCorrection(sourceNode, targetNode, spread, referenceNode):

  # create folder display
  shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
  folderID = shNode.CreateFolderItem(shNode.GetSceneItemID(), targetNode.GetName())
  createFolderDisplayNode(folderID)
  shNode.SetItemAttribute(folderID, 'correction', '1')
  shNode.GetItemDataNode(folderID).SetDescription('Enabled')

  # add fiducials
  for fiducial,attribute in zip([sourceNode, targetNode],['source','target']):
    fiducial.SetName(attribute)
    fiducial.SetLocked(1)
    fiducial.GetDisplayNode().SetVisibility(attribute == 'target') # visualize target fiducials
    fiducial.GetDisplayNode().SetTextScale(0)
    shNode.SetItemAttribute(shNode.GetItemByDataNode(fiducial), 'correction', '1')
    shNode.SetItemAttribute(shNode.GetItemByDataNode(fiducial), attribute, '1')
    shNode.SetItemParent(shNode.GetItemByDataNode(fiducial), folderID)

  # add mask spread
  labelmapNode = generateMaskFromPointsSpread(sourceNode, spread, referenceNode)
  labelmapNode.SetDescription(str(spread))
  labelmapNode.SetName('spread')
  slicer.util.setSliceViewerLayers(labelOpacity=0.4)
  shNode.SetItemAttribute(shNode.GetItemByDataNode(labelmapNode), 'correction', '1')
  shNode.SetItemParent(shNode.GetItemByDataNode(labelmapNode), folderID)


def addFixedPointsFolder():
  # create folder dislpay node to store fixed points
  shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
  folderID = shNode.CreateFolderItem(shNode.GetSceneItemID(), '')
  createFolderDisplayNode(folderID)
  shNode.SetItemAttribute(folderID, 'correction', '1')
  shNode.GetItemDataNode(folderID).SetDescription('Enabled')
  folderNode = shNode.GetItemDataNode(folderID)
  folderNode.SetName('Fixed Points')
  return folderNode

def addFixedPoint(node):

  # get fixed points folder
  try:
    folderNode = slicer.util.getNode('Fixed Points')
  except:
    folderNode = addFixedPointsFolder()

  # add node to the folder and set attributes
  shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
  shNode.SetItemAttribute(shNode.GetItemByDataNode(node), 'correction', '1')
  shNode.SetItemAttribute(shNode.GetItemByDataNode(node), 'fixed', '1')
  shNode.SetItemParent(shNode.GetItemByDataNode(node), shNode.GetItemByDataNode(folderNode))
