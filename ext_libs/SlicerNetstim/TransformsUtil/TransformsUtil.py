import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging

import numpy as np

#
# TransformsUtil
#

class TransformsUtil(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "TransformsUtil" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Netstim"]
    self.parent.dependencies = []
    self.parent.contributors = ["John Doe (AnyWare Corp.)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
This is an example of scripted loadable module bundled in an extension.
It performs a simple thresholding on the input volume and optionally captures a screenshot.
"""
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# TransformsUtilWidget
#

class TransformsUtilWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...


    #
    # input warp selector
    #

    inputFrame = qt.QFrame()
    inputFrame.setLayout(qt.QHBoxLayout())

    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ["vtkMRMLGridTransformNode"]
    self.inputSelector.selectNodeUponCreation = False
    self.inputSelector.addEnabled = True
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = True
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the input to the algorithm." )
    

    inputFrame.layout().addWidget(qt.QLabel('Input Transform: '))
    inputFrame.layout().addWidget(self.inputSelector)

    self.layout.addWidget(inputFrame)

    #
    # Grid reference
    #
    emptyTransformCollapsibleButton = ctk.ctkCollapsibleButton()
    emptyTransformCollapsibleButton.text = "Grid Reference"
    self.layout.addWidget(emptyTransformCollapsibleButton)

    # Layout within the dummy collapsible button
    emptyTransformFormLayout = qt.QFormLayout(emptyTransformCollapsibleButton)

    self.transformSizeSelect = ctk.ctkCoordinatesWidget()
    self.transformSizeSelect.coordinates = '193,229,193'
    self.transformSizeSelect.setDecimals(0)
    emptyTransformFormLayout.addRow("Size", self.transformSizeSelect)

    self.transformOriginSelect = ctk.ctkCoordinatesWidget()
    self.transformOriginSelect.coordinates = '-96.0, -132.0, -78.0'
    emptyTransformFormLayout.addRow("Origin", self.transformOriginSelect)

    self.transformSpacingSelect = ctk.ctkCoordinatesWidget()
    self.transformSpacingSelect.coordinates = '1.0, 1.0, 1.0'
    emptyTransformFormLayout.addRow("Spacing", self.transformSpacingSelect)

    #
    # Operations
    #
    operationsCollapsibleButton = ctk.ctkCollapsibleButton()
    operationsCollapsibleButton.text = "Operations"
    self.layout.addWidget(operationsCollapsibleButton)
    # Layout within the dummy collapsible button
    operationsFormLayout = qt.QFormLayout(operationsCollapsibleButton)

    # Empty Button
    self.emptyButton = qt.QPushButton("Empty Transform")
    self.emptyButton.enabled = False
    operationsFormLayout.addRow(self.emptyButton)

    # Flatten Button & opts

    flattenFrame = qt.QFrame()
    flattenFrame.setLayout(qt.QHBoxLayout())

    self.includeFirstLayerCB = qt.QCheckBox('Include First Layer')
    self.includeFirstLayerCB.toolTip = "If selected all the layers of the composite transform will be flattened taking the first grid as reference. If not, the layers starting from the second one will be flattened and append to the first one. In this case the reference grid is used to flatten the transform"
    flattenFrame.layout().addWidget(self.includeFirstLayerCB)

    self.flattenButton = qt.QPushButton("Flatten Transform")
    self.flattenButton.enabled = False
    flattenFrame.layout().addWidget(self.flattenButton,2)

    operationsFormLayout.addRow(flattenFrame)

    # Remove last layer
    self.removeLastLayerButton = qt.QPushButton("Remove Last Layer")
    self.removeLastLayerButton.enabled = False
    operationsFormLayout.addRow(self.removeLastLayerButton)

    # connections
    self.emptyButton.connect('clicked(bool)', self.onEmptyButton)
    self.flattenButton.connect('clicked(bool)', self.onFlattenButton)
    self.removeLastLayerButton.connect('clicked(bool)', self.onRemoveLastLayerButton)
    self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.emptyButton.enabled = self.inputSelector.currentNode()
    self.flattenButton.enabled = self.inputSelector.currentNode()
    self.removeLastLayerButton.enabled = self.inputSelector.currentNode()


  def onEmptyButton(self):
    logic = TransformsUtilLogic()
    size =    np.fromstring(self.transformSizeSelect.coordinates, dtype=int, sep=',')
    origin =  np.fromstring(self.transformOriginSelect.coordinates, sep=',')
    spacing = np.fromstring(self.transformSpacingSelect.coordinates, sep=',')
    logic.emptyGridTransfrom(size, origin, spacing, self.inputSelector.currentNode())

  def onFlattenButton(self):
    logic = TransformsUtilLogic()
    size =    np.fromstring(self.transformSizeSelect.coordinates, dtype=int, sep=',')
    origin =  np.fromstring(self.transformOriginSelect.coordinates, sep=',')
    spacing = np.fromstring(self.transformSpacingSelect.coordinates, sep=',')
    logic.flattenTransform(self.inputSelector.currentNode(), self.includeFirstLayerCB.checked)

  def onRemoveLastLayerButton(self):
    logic = TransformsUtilLogic()
    logic.removeLastLayer(self.inputSelector.currentNode())

#
# TransformsUtilLogic
#

class TransformsUtilLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """


  def getMNIGrid(self, resolution):
    size = [394, 466, 378]
    origin = [-98, -134, -72]
    spacing = [0.5, 0.5, 0.5]
    if resolution != 0.5:
      size = [int(round(s * spacing[0] / resolution)) for s in size]
      spacing = [resolution] * 3

    return size,origin,spacing

  def emptyGridTransfrom(self, transformSize = [193,229,193], transformOrigin = [-96.0, -132.0, -78.0], transformSpacing = [1.0, 1.0, 1.0], transformNode = None):
    """
    Run the actual algorithm
    """
    if not transformNode:
      transformNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLGridTransformNode')

    voxelType = vtk.VTK_FLOAT
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

    return transformNode


  def emptySplineTransfrom(self, transformSize = [19,23,19], transformOrigin = [-96.0, -132.0, -78.0], transformSpacing = [10, 10, 10], transformNode = None):
    """
    Run the actual algorithm
    """
    if not transformNode:
      transformNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLBSplineTransformNode')

    voxelType = vtk.VTK_FLOAT
    fillVoxelValue = 0
    # Create an empty image volume, filled with fillVoxelValue
    imageData = vtk.vtkImageData()
    imageData.SetDimensions(transformSize)
    imageData.AllocateScalars(voxelType, 3)
    imageData.GetPointData().GetScalars().Fill(fillVoxelValue)
    # Create transform
    transform = slicer.vtkOrientedBSplineTransform()
    transform.SetCoefficientData(imageData)
    # Create transform node
    transformNode.SetAndObserveTransformFromParent(transform)
    transformNode.GetTransformFromParent().GetCoefficientData().SetOrigin(transformOrigin)
    transformNode.GetTransformFromParent().GetCoefficientData().SetSpacing(transformSpacing)
    #transformNode.CreateDefaultDisplayNodes()
    transformNode.CreateDefaultStorageNode()  

    return transformNode

  def createEmpyVolume(self, imageSize, imageOrigin, imageSpacing):
    voxelType = vtk.VTK_UNSIGNED_CHAR
    imageDirections = [[1,0,0], [0,1,0], [0,0,1]]
    fillVoxelValue = 0
    # Create an empty image volume, filled with fillVoxelValue
    imageData = vtk.vtkImageData()
    imageData.SetDimensions(imageSize)
    imageData.AllocateScalars(voxelType, 1)
    imageData.GetPointData().GetScalars().Fill(fillVoxelValue)
    # Create volume node
    volumeNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
    volumeNode.SetOrigin(imageOrigin)
    volumeNode.SetSpacing(imageSpacing)
    volumeNode.SetIJKToRASDirections(imageDirections)
    volumeNode.SetAndObserveImageData(imageData)
    volumeNode.CreateDefaultDisplayNodes()
    volumeNode.CreateDefaultStorageNode()

    return volumeNode

  def getTransformNodesInScene(self):
    transformNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLTransformNode')
    return set([transformNodes.GetItemAsObject(i).GetName() for i in range(transformNodes.GetNumberOfItems())])

  def getGridDefinition(self, transformNode):
    if not transformNode:
      return 0,0,[1,1,1]
    fp = transformNode.GetTransformFromParent()
    tp = transformNode.GetTransformToParent()
    if isinstance(fp, slicer.vtkOrientedGridTransform) and fp.GetDisplacementGrid():
      grid = fp.GetDisplacementGrid()
    elif isinstance(tp, slicer.vtkOrientedGridTransform) and tp.GetDisplacementGrid():
      grid = tp.GetDisplacementGrid()
    elif isinstance(fp, vtk.vtkGeneralTransform) and fp.GetConcatenatedTransform(fp.GetNumberOfConcatenatedTransforms()-1).GetDisplacementGrid():
      grid = fp.GetConcatenatedTransform(fp.GetNumberOfConcatenatedTransforms()-1).GetDisplacementGrid()
    elif isinstance(tp, vtk.vtkGeneralTransform) and tp.GetConcatenatedTransform(tp.GetNumberOfConcatenatedTransforms()-1).GetDisplacementGrid():
      grid = tp.GetConcatenatedTransform(tp.GetNumberOfConcatenatedTransforms()-1).GetDisplacementGrid()
    elif isinstance(fp, slicer.vtkOrientedBSplineTransform) and fp.GetCoefficientData():
      grid = fp.GetCoefficientData()
    else:
      return False
    
    size = grid.GetDimensions()
    origin = grid.GetOrigin()
    spacing = grid.GetSpacing()

    return size,origin,spacing


  def getNumberOfLayers(self, transformNode):
    if not transformNode:
      return 0
    if isinstance(transformNode.GetTransformFromParent(), vtk.vtkGeneralTransform):
      return max(transformNode.GetTransformFromParent().GetNumberOfConcatenatedTransforms(), transformNode.GetTransformToParent().GetNumberOfConcatenatedTransforms())
    else:
      return 1

  def hasMinimumNumberOfLayers(self, transformNode, N):
    return self.getNumberOfLayers(transformNode) >= N


  def splitAndGetNodeNames(self, transformNode):
    # split transform and get the created node names
    preSplitNames = self.getTransformNodesInScene()
    transformNode.Split()
    postSplitNames = self.getTransformNodesInScene()
    newNodeNames = list(postSplitNames - preSplitNames)
    # sort transforms accordig to the last character in the name
    newNodeNames.sort(key=lambda e: int(e[-1]) if e[-1].isdigit() else 0)    

    return newNodeNames

  def flattenTransform(self, transformNode, includeFirstLayer):

    # check that there are at least a number of layers to flatten the transform
    minimumNumberOfLayers = 2 if includeFirstLayer else 3
    if not self.hasMinimumNumberOfLayers(transformNode, minimumNumberOfLayers):
      print('already flat')
      return

    # get current transform ID applied and save for after flatten
    lastLayerID = transformNode.GetTransformNodeID()
    transformNode.SetAndObserveTransformNodeID('')

    size, origin, spacing = self.getGridDefinition(transformNode)
    newNodeNames = self.splitAndGetNodeNames(transformNode)

    if includeFirstLayer:
      newNodeNames.append(transformNode.GetName()) # add first layer in the end

    for nodeName in newNodeNames[1:]:
      node = slicer.util.getNode(nodeName)
      node.HardenTransform()

    outNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
    referenceVolume = self.createEmpyVolume(size, origin, spacing)    
    transformsLogic = slicer.modules.transforms.logic()
    transformsLogic.ConvertToGridTransform(node, referenceVolume, outNode)

    if includeFirstLayer:
      transformNode.SetAndObserveTransformFromParent(outNode.GetTransformFromParent())
      newNodeNames.pop() # remove so that is not deleted later
    else:
      transformNode.SetAndObserveTransformNodeID(outNode.GetID())
      transformNode.HardenTransform()
    
    # cleanup
    slicer.mrmlScene.RemoveNode(referenceVolume)
    slicer.mrmlScene.RemoveNode(outNode)
    for nodeName in newNodeNames:
      slicer.mrmlScene.RemoveNode(slicer.util.getNode(nodeName))

    # re-apply ID
    transformNode.SetAndObserveTransformNodeID(lastLayerID)

    return True

  def removeLastLayer(self, transformNode):

    if not self.hasMinimumNumberOfLayers(transformNode, 2):
      return None

    newNodeNames = self.splitAndGetNodeNames(transformNode)
    newNodeNames.append(transformNode.GetName())

    # remove last transform
    node = slicer.util.getNode(newNodeNames[1])
    node.SetAndObserveTransformNodeID('')

    lastLayer = slicer.util.getNode(newNodeNames.pop(0)) # get last layer and remove from list

    for nodeName in newNodeNames[1:]:
      slicer.util.getNode(nodeName).HardenTransform()

    newNodeNames.pop() # so as to not delete it

    for nodeName in newNodeNames:
      slicer.mrmlScene.RemoveNode(slicer.util.getNode(nodeName))
    
    return lastLayer.GetID()


  def arrayFromGeneralTransform(self, transformNode, componentNumber):
    # https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/util.py
    transformGrid = transformNode.GetTransformFromParent().GetConcatenatedTransform(componentNumber)
    if not isinstance(transformGrid, slicer.vtkOrientedGridTransform):
      return False
    displacementGrid = transformGrid.GetDisplacementGrid()
    nshape = tuple(reversed(displacementGrid.GetDimensions()))
    import vtk.util.numpy_support
    nshape = nshape + (3,)
    narray = vtk.util.numpy_support.vtk_to_numpy(displacementGrid.GetPointData().GetScalars()).reshape(nshape)
    return narray

  def getTransformRASToIJK(self, transformNode):
    size,origin,spacing = self.getGridDefinition(transformNode)
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
    


class TransformsUtilTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_TransformsUtil1()

  def test_TransformsUtil1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import SampleData
    SampleData.downloadFromURL(
      nodeNames='FA',
      fileNames='FA.nrrd',
      uris='http://slicer.kitware.com/midas3/download?items=5767',
      checksums='SHA256:12d17fba4f2e1f1a843f0757366f28c3f3e1a8bb38836f0de2a32bb1cd476560')
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = TransformsUtilLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
