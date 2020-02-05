import os
import unittest
import vtk, qt, ctk, slicer
import logging
import numpy as np
from subprocess import call
import shutil
from math import sqrt
from Helpers import CircleEffect
from Helpers.TransformUtil import getTransformRASToIJK
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin

#
# SmudgeModule
#

class SmudgeModule(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "SmudgeModule" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Examples"]
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
# SmudgeModuleWidget
#

class SmudgeModuleWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  def __init__(self, parent=None):
    ScriptedLoadableModuleWidget.__init__(self, parent)
    VTKObservationMixin.__init__(self)

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # hide reload and test in slicelet mode
    self.reloadCollapsibleButton.setVisible(not __name__ == "__main__")

    # Init parameter node
    self.parameterNode = SmudgeModuleLogic().getParameterNode()
    self.addObserver(self.parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGuiFromMRML)

    # Instantiate and connect widgets ...

    #
    # Inputs Area
    #
    inputsCollapsibleButton = ctk.ctkCollapsibleButton()
    inputsCollapsibleButton.text = "Inputs"
    self.layout.addWidget(inputsCollapsibleButton)

    # hide inputs and test in slicelet mode
    inputsCollapsibleButton.setVisible(not __name__ == "__main__")

    # Layout within the dummy collapsible button
    inputsFormLayout = qt.QFormLayout(inputsCollapsibleButton)

    #
    # input volume selector
    #
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.inputSelector.selectNodeUponCreation = False
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the input to the algorithm." )
    inputsFormLayout.addRow("Image Volume: ", self.inputSelector)


    #
    # affine selector
    #
    self.affineSelector = slicer.qMRMLNodeComboBox()
    self.affineSelector.nodeTypes = ["vtkMRMLLinearTransformNode"]
    self.affineSelector.selectNodeUponCreation = False
    self.affineSelector.addEnabled = True
    self.affineSelector.removeEnabled = False
    self.affineSelector.noneEnabled = False
    self.affineSelector.showHidden = False
    self.affineSelector.showChildNodeTypes = False
    self.affineSelector.setMRMLScene( slicer.mrmlScene )
    self.affineSelector.setToolTip( "Pick the affine transform (applied before the warp)." )
    inputsFormLayout.addRow("Affine: ", self.affineSelector)


    #
    # transform selector
    #
    self.transformSelector = slicer.qMRMLNodeComboBox()
    self.transformSelector.nodeTypes = ["vtkMRMLGridTransformNode"]
    self.transformSelector.selectNodeUponCreation = False
    self.transformSelector.addEnabled = False
    self.transformSelector.removeEnabled = False
    self.transformSelector.noneEnabled = False
    self.transformSelector.showHidden = False
    self.transformSelector.showChildNodeTypes = False
    self.transformSelector.setMRMLScene( slicer.mrmlScene )
    self.transformSelector.setToolTip( "Pick the output to the algorithm." )
    inputsFormLayout.addRow("Warp: ", self.transformSelector)


    #
    # segmentation selector
    #
    self.segmentationSelector = slicer.qMRMLNodeComboBox()
    self.segmentationSelector.nodeTypes = ["vtkMRMLSegmentationNode"]
    self.segmentationSelector.selectNodeUponCreation = False
    self.segmentationSelector.addEnabled = False
    self.segmentationSelector.removeEnabled = False
    self.segmentationSelector.noneEnabled = False
    self.segmentationSelector.showHidden = False
    self.segmentationSelector.showChildNodeTypes = False
    self.segmentationSelector.setMRMLScene( slicer.mrmlScene )
    self.segmentationSelector.setToolTip( "Pick the segmentation to display." )
    inputsFormLayout.addRow("Segmentation: ", self.segmentationSelector)

    self.initializeButton = qt.QPushButton("Initialize")
    self.initializeButton.setEnabled(False)
    inputsFormLayout.addRow(self.initializeButton)

    #
    # Tools Area
    #
    toolsCollapsibleButton = ctk.ctkCollapsibleButton()
    toolsCollapsibleButton.text = "Tools"
    self.layout.addWidget(toolsCollapsibleButton)

    # Layout within the dummy collapsible button
    toolsFormLayout = qt.QFormLayout(toolsCollapsibleButton)


    #
    # Edit options pushbuttons
    #
    optionsFrame = qt.QFrame()
    optionsFrame.setLayout(qt.QHBoxLayout())
    optionsFrame.setMinimumHeight(50)
    
    smudgePixmap = qt.QPixmap(self.resourcePath(os.path.join('Icons','smudgeIcon.png')))
    smudgeIcon = qt.QIcon(smudgePixmap)
    self.smudgeButton = qt.QPushButton()
    self.smudgeButton.setIcon(smudgeIcon)
    self.smudgeButton.setIconSize(smudgePixmap.rect().size())
    self.smudgeButton.setCheckable(True)
    self.smudgeButton.setEnabled(False)

    optionsFrame.layout().addWidget(self.smudgeButton)
    optionsFrame.layout().addWidget(qt.QPushButton("Other"))
    optionsFrame.layout().addWidget(qt.QPushButton("Other"))
    toolsFormLayout.addRow(optionsFrame)


    #
    # radius value
    #
    self.radiusSlider = ctk.ctkSliderWidget()
    self.radiusSlider.singleStep = 0.1
    self.radiusSlider.minimum = 1
    self.radiusSlider.maximum = 50
    self.radiusSlider.decimals = 1
    self.radiusSlider.value = float(self.parameterNode.GetParameter("radius"))
    toolsFormLayout.addRow("Radius (mm):", self.radiusSlider)

    #
    # blurr
    #
    self.blurrSlider = ctk.ctkSliderWidget()
    self.blurrSlider.singleStep = 1
    self.blurrSlider.minimum = 0
    self.blurrSlider.maximum = 100
    self.blurrSlider.decimals = 0
    self.blurrSlider.value = float(self.parameterNode.GetParameter("blurr"))
    toolsFormLayout.addRow("Blurr (%):", self.blurrSlider)

    #
    # blurr
    #
    self.hardnessSlider = ctk.ctkSliderWidget()
    self.hardnessSlider.singleStep = 1
    self.hardnessSlider.minimum = 0
    self.hardnessSlider.maximum = 100
    self.hardnessSlider.decimals = 0
    self.hardnessSlider.value = float(self.parameterNode.GetParameter("hardness"))
    toolsFormLayout.addRow("Hardness (%):", self.hardnessSlider)

    #
    # Flatten pb
    #

    self.flattenButton = qt.QPushButton("Flatten Transform")
    self.flattenButton.setEnabled(False)
    toolsFormLayout.addRow(self.flattenButton)
    
    #
    # Undo Redo
    #
   
    undoredoFrame = qt.QFrame()
    undoredoFrame.setLayout(qt.QHBoxLayout())

    self.undoButton = qt.QPushButton("<-")
    self.undoButton.setEnabled(False)
    self.redoButton = qt.QPushButton("->")

    undoredoFrame.layout().addWidget(self.undoButton)
    undoredoFrame.layout().addWidget(self.redoButton)
    
    toolsFormLayout.addRow(undoredoFrame)


    #
    # View Area
    #
    viewCollapsibleButton = ctk.ctkCollapsibleButton()
    viewCollapsibleButton.text = "View"
    self.layout.addWidget(viewCollapsibleButton)

    # Layout within the dummy collapsible button
    viewFormLayout = qt.QFormLayout(viewCollapsibleButton)

    #
    # Display
    #
    self.layoutComboBox = qt.QComboBox()
    self.layoutComboBox.addItems(['Four Up','Tabbed Slices'])
    viewFormLayout.addRow("Layout:", self.layoutComboBox)



    #
    # Modality
    #
    self.modalityComboBox = qt.QComboBox()
    self.modalityComboBox.addItems(['T1','T2'])
    viewFormLayout.addRow("Modality:", self.modalityComboBox)

    #
    # Template
    #
    self.templateSlider = ctk.ctkSliderWidget()
    self.templateSlider.singleStep = 0.1
    self.templateSlider.minimum = 0
    self.templateSlider.maximum = 1
    self.templateSlider.decimals = 1
    self.templateSlider.value = 0
    viewFormLayout.addRow("Template Image:", self.templateSlider)

    #
    # Warp
    #

    self.warpViewCheckBox = qt.QCheckBox('Slice View')
    viewFormLayout.addRow("Warp Field:", self.warpViewCheckBox)

    #
    # Segmentation
    #

    self.segmentationOutlineSlider = ctk.ctkSliderWidget()
    self.segmentationOutlineSlider.singleStep = 0.1
    self.segmentationOutlineSlider.minimum = 0
    self.segmentationOutlineSlider.maximum = 1
    self.segmentationOutlineSlider.decimals = 1
    self.segmentationOutlineSlider.value = 1
    viewFormLayout.addRow("Segmentation Outline:", self.segmentationOutlineSlider)

    self.segmentationFillSlider = ctk.ctkSliderWidget()
    self.segmentationFillSlider.singleStep = 0.1
    self.segmentationFillSlider.minimum = 0
    self.segmentationFillSlider.maximum = 1
    self.segmentationFillSlider.decimals = 1
    self.segmentationFillSlider.value = 0.5
    viewFormLayout.addRow("Segmentation Fill:", self.segmentationFillSlider)

    self.segmentationList = qt.QListWidget()
    self.segmentationList.setSelectionMode(qt.QAbstractItemView().ExtendedSelection)
    viewFormLayout.addRow("",self.segmentationList)

    # Add vertical spacer
    self.layout.addStretch(1)

    #
    # Save pushbutton
    #

    self.saveButton = qt.QPushButton("Save")
    #self.flattenButton.setEnabled(False)
    #parametersFormLayout.addRow(self.saveButton)
    self.layout.addWidget(self.saveButton)

    # connections
    self.smudgeButton.connect('clicked(bool)', self.onSmudgeButton)
    self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.exit)
    self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateMRMLFromGUI)
    self.transformSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.exit)
    self.transformSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateMRMLFromGUI)
    self.affineSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateMRMLFromGUI)
    #self.segmentationSelector.connect("currentNodeChanged(vtkMRMLNode*)", )
    self.segmentationSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.segmentationNodeChanged)
    self.initializeButton.connect("clicked(bool)", self.onInitializeButton)
    self.radiusSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.blurrSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.hardnessSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.flattenButton.connect("clicked(bool)", self.onFlattenButton)
    self.saveButton.connect("clicked(bool)", self.onSaveButton)
    self.undoButton.connect("clicked(bool)", self.onUndoButton)
    self.redoButton.connect("clicked(bool)", self.onRedoButton)

    self.layoutComboBox.connect('currentIndexChanged(int)', self.onLayoutChanged)
    self.modalityComboBox.connect('currentIndexChanged(int)', self.onModalityChanged)
    self.templateSlider.connect('valueChanged(double)', self.updateTemplateView)
    self.warpViewCheckBox.connect('stateChanged(int)', self.onWarpCheckBoxChange)
    self.segmentationOutlineSlider.connect('valueChanged(double)', self.updateSegmentationOutline)
    self.segmentationFillSlider.connect('valueChanged(double)', self.updateSegmentationFill)
    self.segmentationList.itemSelectionChanged.connect(self.segmentationItemChanged)
    self.segmentationList.itemDoubleClicked.connect(self.segmentationDoubleClicked)

    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)    

    # Refresh Apply button state
    self.updateGuiFromMRML()
    self.segmentationNodeChanged()
  
  

  def cleanup(self):
    self.exit()

  def enter(self):
    pass

  def updateTemplateView(self,value):
    layoutManager = slicer.app.layoutManager()
    for sliceViewName in layoutManager.sliceViewNames():
      compositeNode = layoutManager.sliceWidget(sliceViewName).sliceLogic().GetSliceCompositeNode()
      compositeNode.SetForegroundOpacity(value)

  def onLayoutChanged(self,index):
    layoutManager = slicer.app.layoutManager()
    if index == 0:
      layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpView)
    elif index == 1:
      layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutTabbedSliceView)

  def onModalityChanged(self,index):
    SmudgeModuleLogic().undoInitialize()
    modality = self.modalityComboBox.itemText(index).lower() # get modality
    imageWinLev = {'t1':{'Win':760,'Lev':420}, 't2':{'Win':1787,'Lev':893}}
    # subject
    imageID = self.parameterNode.GetParameter("imageID")
    if imageID != "":
      subjectPath = self.parameterNode.GetParameter("subjectPath")
      if os.path.isfile(os.path.join(subjectPath,"anat_"+modality+".nii")):
        slicer.mrmlScene.RemoveNode(slicer.util.getNode(imageID))
        ls, imageNode = slicer.util.loadVolume(os.path.join(subjectPath,"anat_"+modality+".nii"), properties = {'name':"anat_"+modality, 'show':False}, returnNode = True)
        self.parameterNode.SetParameter("imageID", imageNode.GetID())
        slicer.util.setSliceViewerLayers(background=imageNode)
        print(imageWinLev[modality]['Win'])
        imageNode.GetDisplayNode().AutoWindowLevelOff()
        imageNode.GetDisplayNode().SetWindow(imageWinLev[modality]['Win'])
        imageNode.GetDisplayNode().SetLevel(imageWinLev[modality]['Lev'])
    # template
    MNIPath = self.parameterNode.GetParameter("MNIPath")
    if MNIPath != "":
      slicer.mrmlScene.RemoveNode(slicer.util.getNode(self.parameterNode.GetParameter("templateID")))
      ls, templateNode = slicer.util.loadVolume(os.path.join(MNIPath,modality+".nii"), properties = {'name':"template_"+modality, 'show':False}, returnNode = True)
      self.parameterNode.SetParameter("templateID", templateNode.GetID())
      slicer.util.setSliceViewerLayers(foreground=templateNode)
      templateNode.GetDisplayNode().AutoWindowLevelOff()
      templateNode.GetDisplayNode().SetWindow(100)
      templateNode.GetDisplayNode().SetLevel(70)
    SmudgeModuleLogic().initialize()

  def onWarpCheckBoxChange(self,state):
    affineTransformID = self.parameterNode.GetParameter( "affineTransformID")
    if affineTransformID != "":
      affineTransformNode = slicer.util.getNode(affineTransformID)
      if not affineTransformNode.GetDisplayNode():
        affineTransformNode.CreateDefaultDisplayNodes()
      affineTransformNode.GetDisplayNode().SetSliceIntersectionVisibility(state)

  def onInitializeButton(self):
    if bool(int(self.parameterNode.GetParameter("initialized"))):
      SmudgeModuleLogic().undoInitialize()    
    SmudgeModuleLogic().initialize()

  def segmentationDoubleClicked(self,item):
    segmentationNode = slicer.util.getNode(self.parameterNode.GetParameter("segmentationID"))
    item.setSelected(1)
    seg = segmentationNode.GetSegmentation().GetSegment(item.text())
    pd = seg.GetRepresentation('Closed surface')
    center = vtk.vtkCenterOfMass()
    center.SetInputData(pd)
    center.Update()
    segCenter = center.GetCenter()
    # create markups node, add center as fiducial and jump and center slices
    markupsNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    markupsNode.AddFiducialFromArray(np.array(segCenter),'')
    markupsLogic = slicer.modules.markups.logic()
    markupsLogic.JumpSlicesToNthPointInMarkup(markupsNode.GetID(),0,True)
    slicer.mrmlScene.RemoveNode(markupsNode)

  def segmentationItemChanged(self):
    segmentationNode = slicer.util.getNode(self.parameterNode.GetParameter("segmentationID"))
    for i in range(segmentationNode.GetSegmentation().GetNumberOfSegments()):
      segmentationNode.GetDisplayNode().SetSegmentVisibility(self.segmentationList.item(i).text(), self.segmentationList.item(i).isSelected())

  def updateSegmentationOutline(self,value):
    segmentationNode = slicer.util.getNode(self.parameterNode.GetParameter("segmentationID"))
    segmentationNode.GetDisplayNode().SetAllSegmentsOpacity2DOutline(value)

  def updateSegmentationFill(self,value):
    segmentationNode = slicer.util.getNode(self.parameterNode.GetParameter("segmentationID"))
    segmentationNode.GetDisplayNode().SetAllSegmentsOpacity2DFill(value)

  def segmentationNodeChanged(self):
    self.segmentationList.clear()
    self.updateMRMLFromGUI()
    
  def onUndoButton(self):
    SmudgeModuleLogic().undoOperation()

  def onRedoButton(self):
    SmudgeModuleLogic().redoOperation()
   
  def onFlattenButton(self):
    SmudgeModuleLogic().flattenTransform()
  
  def onSaveButton(self):
    SmudgeModuleLogic().flattenTransform(overwriteBool=1)

  def updateGuiFromMRML(self, caller="", event=""):
    radius = float(self.parameterNode.GetParameter("radius"))
    self.radiusSlider.setValue( radius )
    if radius < self.radiusSlider.minimum or radius > self.radiusSlider.maximum:
      self.updateMRMLFromGUI()
    self.blurrSlider.setValue(float(self.parameterNode.GetParameter("blurr")))
    self.hardnessSlider.setValue(float(self.parameterNode.GetParameter("hardness")))
    redoTransformID = self.parameterNode.GetParameter("redoTransformID")
    self.redoButton.setEnabled(not redoTransformID == "")
    if not self.segmentationList.item(0): # empty list
      segmentationID = self.parameterNode.GetParameter("segmentationID")
      if segmentationID != "":
        segmentationNode = slicer.util.getNode(segmentationID)
        segmentationNode.CreateDefaultDisplayNodes()
        for i in range(segmentationNode.GetSegmentation().GetNumberOfSegments()):
          self.segmentationList.addItem(segmentationNode.GetSegmentation().GetNthSegment(i).GetName())
        segmentationNode.CreateClosedSurfaceRepresentation()
        segmentationNode.GetDisplayNode().SetAllSegmentsVisibility(0)
        segmentationNode.GetDisplayNode().SetAllSegmentsVisibility3D(0)
    self.smudgeButton.setEnabled(bool(int(self.parameterNode.GetParameter( "initialized"))))
    self.flattenButton.setEnabled(bool(int(self.parameterNode.GetParameter( "initialized"))))
    self.undoButton.setEnabled(bool(int(self.parameterNode.GetParameter( "initialized"))))

  def updateMRMLFromGUI(self):
    self.parameterNode.SetParameter( "radius", str(self.radiusSlider.value) )
    self.parameterNode.SetParameter( "blurr", str(self.blurrSlider.value) )
    self.parameterNode.SetParameter( "hardness", str(self.hardnessSlider.value) )
    if self.transformSelector.currentNode():
      self.parameterNode.SetParameter( "transformID", self.transformSelector.currentNode().GetID())
    if self.inputSelector.currentNode():
      self.parameterNode.SetParameter( "imageID", self.inputSelector.currentNode().GetID())
    if self.segmentationSelector.currentNode():
      self.parameterNode.SetParameter( "segmentationID", self.segmentationSelector.currentNode().GetID())
    if self.affineSelector.currentNode():
      self.parameterNode.SetParameter( "affineTransformID", self.affineSelector.currentNode().GetID())
    self.initializeButton.setEnabled(self.inputSelector.currentNode() and self.transformSelector.currentNode() and self.affineSelector.currentNode())
    
    

  def onSmudgeButton(self, buttonDown):
    if buttonDown:
      SmudgeModuleLogic().smudgeOn()
    else:
      SmudgeModuleLogic().smudgeOff()
      
  def exit(self):
    SmudgeModuleLogic().smudgeOff()
    SmudgeModuleLogic().undoInitialize()    
    imageID = self.parameterNode.GetParameter("imageID")
    if imageID != "":
      imageNode = slicer.util.getNode(imageID)
      imageNode.SetAndObserveTransformNodeID("")
    self.smudgeButton.setChecked(False)

  def onSceneStartClose(self, caller, event):
    self.parameterNode.RemoveAllObservers()
    self.exit()


#
# SmudgeEffectTool
#

class SmudgeEffectTool(CircleEffect.CircleEffectTool):

  # aux variable used to store all instances so they are cleaned later
  _instances = set()

  def __init__(self,sliceWidget):
    self.parameterNode = SmudgeModuleLogic().getParameterNode()
    super(SmudgeEffectTool,self).__init__(sliceWidget)

    self.transformNode = slicer.util.getNode(self.parameterNode.GetParameter("transformID"))

    # transform data
    self.auxTransformNode = slicer.util.getNode(self.parameterNode.GetParameter("auxTransformID"))
    self.auxTransformSpacing = self.auxTransformNode.GetTransformFromParent().GetDisplacementGrid().GetSpacing()[0] # Asume isotropic!
    self.auxTransfromRASToIJK = getTransformRASToIJK(self.auxTransformNode)  
    self.auxTransformArray = slicer.util.array(self.auxTransformNode.GetID())

    self.previousPoint = [0,0,0]   
    self.smudging = False
    # deactivate window level adjust
    interactorStyle = self.sliceView.sliceViewInteractorStyle()
    interactorStyle.SetActionEnabled(interactorStyle.AdjustWindowLevelBackground, False)

    self._instances.add(self)


  def processEvent(self, caller=None, event=None):
    super(SmudgeEffectTool,self).processEvent(caller, event)

    if event == 'LeftButtonPressEvent':
      # get aux transform array in case it has been modified externalyy
      self.auxTransformArray = slicer.util.array(self.auxTransformNode.GetID())
      # clean redo transform
      redoTransformID = self.parameterNode.GetParameter("redoTransformID")
      if redoTransformID != "":
        slicer.mrmlScene.RemoveNode(slicer.util.getNode(self.parameterNode.GetParameter("redoTransformID")))
        self.parameterNode.SetParameter("redoTransformID","")
      self.smudging = True
      xy = self.interactor.GetEventPosition()
      xyToRAS = self.sliceLogic.GetSliceNode().GetXYToRAS()
      self.previousPoint = xyToRAS.MultiplyDoublePoint( (xy[0], xy[1], 0, 1) )[0:3]
    elif event == 'LeftButtonReleaseEvent':
      self.smudging = False
      self.transformNode.HardenTransform()
      SmudgeModuleLogic().createEmptyTransform(transformNode=self.auxTransformNode)
      self.transformNode.SetAndObserveTransformNodeID(self.parameterNode.GetParameter("auxTransformID"))
      
    elif event == 'MouseMoveEvent':
      if self.smudging:
        # get current IJK coord
        xy = self.interactor.GetEventPosition()
        xyToRAS = self.sliceLogic.GetSliceNode().GetXYToRAS()
        currentPoint = xyToRAS.MultiplyDoublePoint( (xy[0], xy[1], 0, 1) )[0:3]
        pos_i,pos_j,pos_k,aux = self.auxTransfromRASToIJK.MultiplyDoublePoint(currentPoint + (1,))
        k,j,i = int(round(pos_k)),int(round(pos_j)),int(round(pos_i))
        # create a sphere with redius
        r = int(round(float(self.parameterNode.GetParameter("radius")) / self.auxTransformSpacing))
        xx, yy, zz = np.mgrid[:2*r+1, :2*r+1, :2*r+1]
        sphereResult = (xx-r) ** 2 + (yy-r) ** 2 + (zz-r) ** 2
        sphereResult[r][r][r] = 1 # replace 0 by 1
        sphereLarge = sphereResult <= (r**2+1) # sphere that the mouse shows
        sphereSmall = sphereResult <= ((r * (1-float(self.parameterNode.GetParameter("blurr")) / 100.0)) **2 + 1 ) # Blurr amount
        sphereResult = 1.0 / sphereResult # invert
        # get value in the edge of the small sphere
        i1,i2,i3 = np.nonzero(sphereSmall)
        newMaxValue = sphereResult[i1[0]][i2[0]][i3[0]]
        # set same value inside the small sphere
        sphereResult[sphereSmall] = sphereSmall[sphereSmall] * newMaxValue
        # delete outside values 
        sphereResult = sphereResult * sphereLarge
        # set range to [0-1]
        newMinValue = sphereResult.min()
        sphereResult = (sphereResult - newMinValue) / (newMaxValue - newMinValue)
        # set hardness
        sphereResult = sphereResult * float(self.parameterNode.GetParameter("hardness")) / 100.0
        # apply to transform array
        index = slice(k-r,k+r+1), slice(j-r,j+r+1), slice(i-r,i+r+1)
        self.auxTransformArray[index] += np.stack([(sphereResult) * i for i in (np.array(self.previousPoint) - np.array(currentPoint))],3)
        # update view
        self.auxTransformNode.Modified()
        # update previous point
        self.previousPoint = currentPoint

  def cleanup(self):
    interactorStyle = self.sliceView.sliceViewInteractorStyle()
    interactorStyle.SetActionEnabled(interactorStyle.AdjustWindowLevelBackground, True)
    super(SmudgeEffectTool,self).cleanup()

  @classmethod
  def empty(cls):
    # clean instances and reset
    for inst in cls._instances:
      inst.cleanup()
    cls._instances = set()
  

#
# SmudgeModuleLogic
#

class SmudgeModuleLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self):
    ScriptedLoadableModuleLogic.__init__(self)
    #super(SmudgeModuleLogic, self).__init__()
    self.isSingletonParameterNode = False

  def initialize(self):
    parameterNode = self.getParameterNode()
    if not bool(int(parameterNode.GetParameter("initialized"))):
      if parameterNode.GetParameter("auxTransformID") == "":
        self.createEmptyTransform(parameterName="auxTransformID")
      # apply affine transform to image
      affineNode = slicer.util.getNode(parameterNode.GetParameter("affineTransformID"))
      imageNode = slicer.util.getNode(parameterNode.GetParameter("imageID"))
      imageNode.ApplyTransform(affineNode.GetTransformToParent())
      # invert affine and set to image 
      affineNode.Inverse()
      imageNode.SetAndObserveTransformNodeID(parameterNode.GetParameter("affineTransformID"))
      # set transform to affine
      affineNode.SetAndObserveTransformNodeID(parameterNode.GetParameter("transformID"))
      # set aux transform to transform
      transformNode = slicer.util.getNode(parameterNode.GetParameter("transformID"))
      transformNode.SetAndObserveTransformNodeID(parameterNode.GetParameter("auxTransformID"))

      parameterNode.SetParameter("initialized", "1")

  def undoInitialize(self):
    parameterNode = self.getParameterNode()

    if bool(int(parameterNode.GetParameter("initialized"))):
    
      imageNode = slicer.util.getNode(parameterNode.GetParameter("imageID"))
      affineNode = slicer.util.getNode(parameterNode.GetParameter("affineTransformID"))
      transformNode = slicer.util.getNode(parameterNode.GetParameter("transformID"))

      transformNode.SetAndObserveTransformNodeID("")
      affineNode.SetAndObserveTransformNodeID("")
      imageNode.SetAndObserveTransformNodeID("")
      affineNode.Inverse()
      imageNode.ApplyTransform(affineNode.GetTransformFromParent())

      parameterNode.SetParameter("initialized", "0")

  def createParameterNode(self):
    node = ScriptedLoadableModuleLogic.createParameterNode(self)
    node.SetParameter("transformID", "")
    node.SetParameter("auxTransformID", "")
    node.SetParameter("redoTransformID", "")
    node.SetParameter("affineTransformID", "")
    node.SetParameter("imageID", "")
    node.SetParameter("templateID", "")
    node.SetParameter("segmentationID","")
    node.SetParameter("radius", "10")
    node.SetParameter("blurr", "40")
    node.SetParameter("hardness", "100")
    node.SetParameter("initialized", "0")
    node.SetParameter("subjectPath", "")
    node.SetParameter("MNIPath", "")
    node.SetParameter("atnsApplyTransformsPath", "")
    return node

  def smudgeOn(self):

    parameterNode = self.getParameterNode()
    
    if parameterNode.GetParameter("auxTransformID") == "":
      self.createEmptyTransform(parameterName="auxTransformID")

    for color in ['Red','Green','Yellow']:
      sliceWidget = slicer.app.layoutManager().sliceWidget(color)
      SmudgeEffectTool(sliceWidget)

  
  def smudgeOff(self):
    parameterNode = self.getParameterNode()
    SmudgeEffectTool.empty()
    auxTransformID = parameterNode.GetParameter("auxTransformID")
    if auxTransformID != "":
      slicer.mrmlScene.RemoveNode(slicer.util.getNode(auxTransformID))
      parameterNode.SetParameter("auxTransformID", "")
    redoTransformID = parameterNode.GetParameter("redoTransformID")
    if redoTransformID != "":
      slicer.mrmlScene.RemoveNode(slicer.util.getNode(redoTransformID))
      parameterNode.SetParameter("redoTransformID", "")
    
  def createEmptyTransform(self, transformNode = None, parameterName = ""):
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
    if not transformNode:
      transformNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLGridTransformNode")
      parameterNode = self.getParameterNode()
      parameterNode.SetParameter(parameterName, transformNode.GetID())
    transformNode.SetAndObserveTransformFromParent(transform)
    transformNode.GetTransformFromParent().GetDisplacementGrid().SetOrigin(transformOrigin)
    transformNode.GetTransformFromParent().GetDisplacementGrid().SetSpacing(transformSpacing)
    #transformNode.CreateDefaultDisplayNodes()
    transformNode.CreateDefaultStorageNode()    

  def flattenTransform(self, overwriteBool = 0):
    """
    Flat all manual transformations
    """
    # Get nodes
    parameterNode = self.getParameterNode()
    transformNode = slicer.util.getNode(parameterNode.GetParameter("transformID"))
    auxTransformNode = slicer.util.getNode(parameterNode.GetParameter("auxTransformID"))  
    
    # get all transforms in a collection
    col = vtk.vtkCollection()
    transformNode.FlattenGeneralTransform(col,transformNode.GetTransformToParent())

    if (col.GetNumberOfItems() < 3 and not overwriteBool): # already flat!
      return

    path,transformFilename = os.path.split(transformNode.GetStorageNode().GetFileName()) # get working dir
    os.mkdir(os.path.join(path,"tmp")) # create temp folder
    
    if overwriteBool:
      transformFullFilename = os.path.join(path,transformFilename)
      backupFile = os.path.join(path,"glanatComposite_backup.nii.gz")
      if not os.path.isfile(backupFile):
        shutil.copyfile(transformFullFilename,backupFile)
      command = "-t " + transformFullFilename + " "
    else:
      command = ""
    
    outname = os.path.join(path,"tmp","out.nii.gz") # flat transform name
    command = command + "-o [" + outname + ",1]"

    for i in range(1,col.GetNumberOfItems()): # save each transform and append name to command
    	con = vtk.vtkWarpTransform.SafeDownCast(col.GetItemAsObject(i))
    	auxTransformNode.SetAndObserveTransformToParent(con)
    	fullname = os.path.join(path, "tmp", "test" + str(i) + ".nii.gz")
    	slicer.util.saveNode(auxTransformNode,fullname)
    	command = "-t " + fullname + " " + command

    self.createEmptyTransform(transformNode=auxTransformNode) # reset aux transform

    if overwriteBool:
      referenceFile = os.path.join(path,"glanat.nii")
    else:
      referenceFile = fullname

    command = parameterNode.GetParameter("atnsApplyTransformsPath") + " -r " + referenceFile + " " + command
    commandOut = call(command, env=slicer.util.startupEnvironment(), shell=True) # run antsApplyTransforms

    if overwriteBool:
      os.remove(transformFullFilename)
      shutil.copyfile(outname,transformFullFilename)
      print('Saved!')
    else:
      ls,flatTransformNode = slicer.util.loadTransform(outname, returnNode = True) # load result
      # reset original transform
      con = vtk.vtkWarpTransform.SafeDownCast(col.GetItemAsObject(0))
      transformNode.SetAndObserveTransformToParent(con)
      # apply the flat transform
      transformNode.SetAndObserveTransformNodeID(flatTransformNode.GetID())
      transformNode.HardenTransform()
      #remove node
      slicer.mrmlScene.RemoveNode(flatTransformNode)
      # reset actions
      transformNode.SetAndObserveTransformNodeID(parameterNode.GetParameter("auxTransformID"))
      print('Flattened!')
    
    # cleanup
    shutil.rmtree(os.path.join(path,"tmp"))

  def undoOperation(self):

    # Get nodes
    parameterNode = self.getParameterNode()
    transformNode = slicer.util.getNode(parameterNode.GetParameter("transformID"))
    auxTransformNode = slicer.util.getNode(parameterNode.GetParameter("auxTransformID")) 

    col = vtk.vtkCollection()
    transformNode.FlattenGeneralTransform(col,transformNode.GetTransformToParent())

    if col.GetNumberOfItems() < 2:
      return
    
    con = vtk.vtkWarpTransform.SafeDownCast(col.GetItemAsObject(0))
    transformNode.SetAndObserveTransformToParent(con)
    
    for i in range(1,col.GetNumberOfItems()-1): # save each transform and append name to command
      con = vtk.vtkWarpTransform.SafeDownCast(col.GetItemAsObject(i))
      auxTransformNode.SetAndObserveTransformToParent(con)
      transformNode.SetAndObserveTransformNodeID(auxTransformNode.GetID())
      transformNode.HardenTransform()

    self.createEmptyTransform(transformNode=auxTransformNode) # reset aux transform

    # save last op for redo
    redoTransformID = parameterNode.GetParameter("redoTransformID")
    if redoTransformID == "":
      self.createEmptyTransform(parameterName="redoTransformID")
      redoTransformID = parameterNode.GetParameter("redoTransformID")
    redoTransformNode =  slicer.util.getNode(redoTransformID)
    con = vtk.vtkWarpTransform.SafeDownCast(col.GetItemAsObject(col.GetNumberOfItems()-1))
    redoTransformNode.SetAndObserveTransformToParent(con)

    transformNode.SetAndObserveTransformNodeID(parameterNode.GetParameter("auxTransformID")) 
  
  def redoOperation(self):
    # Get nodes
    parameterNode = self.getParameterNode()
    transformNode = slicer.util.getNode(parameterNode.GetParameter("transformID"))
    redoTransformNode = slicer.util.getNode(parameterNode.GetParameter("redoTransformID")) 
    transformNode.SetAndObserveTransformNodeID(redoTransformNode.GetID())
    transformNode.HardenTransform()
    slicer.mrmlScene.RemoveNode(redoTransformNode)
    parameterNode.SetParameter("redoTransformID","")
    transformNode.SetAndObserveTransformNodeID(parameterNode.GetParameter("auxTransformID"))


class SmudgeModuleTest(ScriptedLoadableModuleTest):
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
    self.test_SmudgeModule1()

  def test_SmudgeModule1(self):
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

    imagePath = '/Users/netstim/Desktop/smudge_test/anat_t1_affine.nii'
    transformPath = '/Users/netstim/Desktop/smudge_test/glanatComposite_downsampleModified.nii.gz'

    ls,imageNode = slicer.util.loadVolume(imagePath, properties = {'show':True}, returnNode = True)
    ls,transformNode = slicer.util.loadTransform(transformPath, returnNode = True)

    #node = SmudgeModuleLogic().getParameterNode()
    #node.SetParameter( "transformID", transformNode.GetID())
    #node.SetParameter( "imageID", imageNode.GetID())

    #SmudgeModuleLogic().smudgeOn()




## Slicelet


class SmudgeModuleSlicelet(object):
  """A slicer slicelet is a module widget that comes up in stand alone mode
  implemented as a python class.
  This class provides common wrapper functionality used by all slicer modlets.
  """

  def __init__(self, widgetClass=None):
    self.parent = qt.QFrame()
    self.parent.setLayout( qt.QHBoxLayout() )

    # Left Frame
    self.leftFrame = qt.QFrame(self.parent)
    self.leftFrame.setLayout(qt.QVBoxLayout())
    self.parent.layout().addWidget(self.leftFrame, 1)

    # Rigth Frame
    layoutWidget = slicer.qMRMLLayoutWidget()
    layoutManager = slicer.qSlicerLayoutManager()
    layoutManager.setMRMLScene(slicer.mrmlScene)
    #layoutManager.setScriptedDisplayableManagerDirectory(slicer.app.slicerHome + "/bin/Python/mrmlDisplayableManager")
    layoutWidget.setLayoutManager(layoutManager)
    slicer.app.setLayoutManager(layoutManager)
    layoutWidget.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpView)
    self.parent.layout().addWidget(layoutWidget, 3)

    # add module
    self.widget = SmudgeModuleWidget(self.leftFrame)
    self.widget.setup()
    self.parent.showMaximized()



if __name__ == "__main__":

  import sys
  print( sys.argv )

  slicelet = SmudgeModuleSlicelet()

  # load subject data
  if len(sys.argv) > 1:
    subjectPath = sys.argv[1]
    leadPath = sys.argv[2]

    MNIPath = os.path.join(leadPath,'templates','space','MNI_ICBM_2009b_NLIN_ASYM')

    from sys import platform
    if platform == "darwin":
      ext = "maci64"

    # TODO: other platforms

    atnsApplyTransformsPath = os.path.join(leadPath,'ext_libs','ANTs','antsApplyTransforms.' + ext)

    imagePath = os.path.join(subjectPath,'anat_t1.nii')
    templatePath = os.path.join(MNIPath,'t1.nii')
    transformPath = os.path.join(subjectPath,'glanatComposite.nii.gz')
    segmentationPath = os.path.join(leadPath,'ext_libs','SmudgeModule','DISTAL Minimal (Ewert 2017).seg.nrrd')
    affinePath = os.path.join(subjectPath,'glanat0GenericAffine_backup.mat')

    ls, segmentationNode = slicer.util.loadSegmentation(segmentationPath, returnNode=True)
    ls, imageNode = slicer.util.loadVolume(imagePath, properties = {'name':'anat_t1', 'show':False}, returnNode = True)
    ls, templateNode = slicer.util.loadVolume(templatePath, properties = {'name':'template_t1', 'show':False}, returnNode = True)
    ls, transformNode = slicer.util.loadTransform(transformPath, returnNode = True)
    
    if not os.path.isfile(affinePath):
      msg=qt.QMessageBox().warning(qt.QWidget(),'','Affine transform not available. The Displayed warp is affine+deformable')
      affineNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLinearTransformNode')
    else:
      ls, affineNode = slicer.util.loadTransform(affinePath, returnNode = True)
    
    parameterNode = SmudgeModuleLogic().getParameterNode()
    
    parameterNode.SetParameter("transformID", transformNode.GetID())
    parameterNode.SetParameter("affineTransformID", affineNode.GetID())
    parameterNode.SetParameter("imageID", imageNode.GetID())
    parameterNode.SetParameter("templateID", templateNode.GetID())
    parameterNode.SetParameter("segmentationID", segmentationNode.GetID())
    parameterNode.SetParameter("subjectPath", subjectPath)
    parameterNode.SetParameter("MNIPath", MNIPath)
    parameterNode.SetParameter("atnsApplyTransformsPath", atnsApplyTransformsPath)
    
    SmudgeModuleLogic().initialize()
    slicer.util.setSliceViewerLayers(background=imageNode, foreground=templateNode)

    templateNode.GetDisplayNode().AutoWindowLevelOff()
    templateNode.GetDisplayNode().SetWindow(100)
    templateNode.GetDisplayNode().SetLevel(70)
    imageNode.GetDisplayNode().AutoWindowLevelOff()
    imageNode.GetDisplayNode().SetWindow(760)
    imageNode.GetDisplayNode().SetLevel(420)

    # Link views
    #for sliceCompositeNode in slicer.util.getNodesByClass('vtkMRMLSliceCompositeNode'):
    #  sliceCompositeNode.SetLinkedControl(True)






