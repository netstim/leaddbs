import os, sys
import unittest
import vtk, qt, ctk, slicer
import logging
import numpy as np
from subprocess import call
import shutil
from math import sqrt
from Helpers import PointerEffect, WarpEffect
#sys.path.append(os.path.join(os.path.dirname(__file__), '..'))  # somehow needs this to find module
from Helpers import FunctionsUtil
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin
import glob
import SimpleITK as sitk
import h5py # should already be installed

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

    # hide reload and test in slicelet (when in developer mode)
    if slicer.util.settingsValue('Developer/DeveloperMode', False, converter=slicer.util.toBool):
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
    self.segmentationSelector.nodeTypes = ["vtkMRMLModelHierarchyNode"]
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
    # Subject
    #

    subjectFrame = qt.QFrame()
    subjectFrame.setLayout(qt.QHBoxLayout())
    self.subjectNameLabel = qt.QLabel('')

    subjectFrame.layout().addWidget(qt.QLabel('Subject: '))
    subjectFrame.layout().addWidget(self.subjectNameLabel)
    subjectFrame.layout().addStretch()

    self.layout.addWidget(subjectFrame)

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
    #self.smudgeButton.setAutoExclusive(True)
    optionsFrame.layout().addWidget(self.smudgeButton)

    pencilPixmap = qt.QPixmap(self.resourcePath(os.path.join('Icons','pencilIcon.png')))
    pencilIcon = qt.QIcon(pencilPixmap)
    self.snapButton = qt.QPushButton()
    self.snapButton.setIcon(pencilIcon)
    self.snapButton.setIconSize(pencilPixmap.rect().size())
    self.snapButton.setCheckable(True)
    self.snapButton.setEnabled(False)
    #self.snapButton.setAutoExclusive(True)
    optionsFrame.layout().addWidget(self.snapButton)

    #optionsFrame.layout().addWidget(qt.QPushButton("Other"))
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
    # hardness
    #
    self.hardnessSlider = ctk.ctkSliderWidget()
    self.hardnessSlider.singleStep = 1
    self.hardnessSlider.minimum = 0
    self.hardnessSlider.maximum = 100
    self.hardnessSlider.decimals = 0
    self.hardnessSlider.value = float(self.parameterNode.GetParameter("hardness"))
    toolsFormLayout.addRow("Hardness (%):", self.hardnessSlider)


    #
    # History Area
    #
    historyCollapsibleButton = ctk.ctkCollapsibleButton()
    historyCollapsibleButton.text = "History"
    self.layout.addWidget(historyCollapsibleButton)

    historyGridLayout = qt.QGridLayout(historyCollapsibleButton)

    #
    # Undo Redo Flatten
    #   
    undoredoFrame = qt.QFrame()
    undoredoFrame.setLayout(qt.QVBoxLayout())

    flattenPixmap = qt.QPixmap(self.resourcePath(os.path.join('Icons','flattenIcon.png')))
    flattenIcon = qt.QIcon(flattenPixmap)
    self.flattenButton = qt.QPushButton()
    self.flattenButton.setIcon(flattenIcon)
    self.flattenButton.setIconSize(flattenPixmap.rect().size())
    self.flattenButton.setEnabled(False)    
    self.flattenButton.toolTip = 'Flatten'

    undoPixmap = qt.QPixmap(self.resourcePath(os.path.join('Icons','undoIcon.png')))
    undoIcon = qt.QIcon(undoPixmap)
    self.undoButton = qt.QPushButton()
    self.undoButton.setIcon(undoIcon)
    self.undoButton.setIconSize(undoPixmap.rect().size())
    self.undoButton.setEnabled(False)
    self.undoButton.toolTip = 'Undo'

    redoPixmap = qt.QPixmap(self.resourcePath(os.path.join('Icons','redoIcon.png')))
    redoIcon = qt.QIcon(redoPixmap)
    self.redoButton = qt.QPushButton()
    self.redoButton.setIcon(redoIcon)
    self.redoButton.setIconSize(redoPixmap.rect().size())
    self.redoButton.setEnabled(False)
    self.redoButton.toolTip = 'Redo'

    undoredoFrame.layout().addWidget(self.undoButton)
    undoredoFrame.layout().addWidget(self.redoButton)
    undoredoFrame.layout().addWidget(self.flattenButton)
    
    historyGridLayout.addWidget(undoredoFrame,0,0)

    # History Stack

    self.historyList = qt.QListWidget()
    self.historyList.addItem('glanatComposite')
    self.historyList.setCurrentRow(0)

    historyGridLayout.addWidget(self.historyList,0,1)


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

    displayFrame = qt.QFrame()
    displayFrame.setLayout(qt.QHBoxLayout())

    self.layoutComboBox = qt.QComboBox()
    self.layoutComboBox.addItems(['Four Up','Tabbed Slices'])

    self.sliceIntersectionVisibilityCheckBox = qt.QCheckBox('Slice Intersection')

    displayFrame.layout().addWidget(self.layoutComboBox)
    displayFrame.layout().addWidget(self.sliceIntersectionVisibilityCheckBox)

    viewFormLayout.addRow("Layout:", displayFrame)

    #
    # Modality
    #
    self.modalityComboBox = qt.QComboBox()
    self.modalityComboBox.addItems(['T1','T2'])
    viewFormLayout.addRow("Modality:", self.modalityComboBox)
    self.modalityComboBox.setEnabled(__name__=="__main__")

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

    self.segmentationComboBox = qt.QComboBox()
    self.segmentationComboBox.addItem('Choose Atlas')
    viewFormLayout.addRow("Atlas:", self.segmentationComboBox)

    self.segmentationOutlineSlider = ctk.ctkSliderWidget()
    self.segmentationOutlineSlider.singleStep = 0.1
    self.segmentationOutlineSlider.minimum = 0
    self.segmentationOutlineSlider.maximum = 1
    self.segmentationOutlineSlider.decimals = 1
    self.segmentationOutlineSlider.value = 1
    viewFormLayout.addRow("Segmentation Outline:", self.segmentationOutlineSlider)


    self.segmentationList = qt.QListWidget()
    self.segmentationList.setSelectionMode(qt.QAbstractItemView().ExtendedSelection)
    viewFormLayout.addRow("",self.segmentationList)

    # Add vertical spacer
    self.layout.addStretch(1)

    #
    # Save Area
    #

    self.saveButton = qt.QPushButton("Save and Next")
    self.saveButton.setMinimumHeight(30)
    self.saveButton.setStyleSheet("background-color: green")
    self.layout.addWidget(self.saveButton)

    # connections
    self.smudgeButton.connect('clicked(bool)', self.onSmudgeButton)
    self.smudgeButton.toggled.connect(self.toogleTools)
    self.snapButton.connect('clicked(bool)', self.onSnapButton)
    self.snapButton.toggled.connect(self.toogleTools)
    self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.exit)
    self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateMRMLFromGUI)
    self.transformSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.exit)
    self.transformSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateMRMLFromGUI)
    self.affineSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateMRMLFromGUI)
    self.segmentationSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.segmentationNodeChanged)
    self.initializeButton.connect("clicked(bool)", self.onInitializeButton)
    self.radiusSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.blurrSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.hardnessSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.flattenButton.connect("clicked(bool)", self.onFlattenButton)
    self.saveButton.connect("clicked(bool)", self.onSaveButton)
    self.undoButton.connect("clicked(bool)", self.onUndoButton)
    self.redoButton.connect("clicked(bool)", self.onRedoButton)
    self.historyList.itemSelectionChanged.connect(self.historyItemChanged)

    self.layoutComboBox.connect('currentIndexChanged(int)', self.onLayoutChanged)
    self.sliceIntersectionVisibilityCheckBox.connect('stateChanged(int)', self.onSliceIntersectionVisibilityChanged)
    self.modalityComboBox.connect('currentIndexChanged(int)', self.onModalityChanged)
    self.templateSlider.connect('valueChanged(double)', self.updateTemplateView)
    self.warpViewCheckBox.connect('stateChanged(int)', self.onWarpCheckBoxChange)
    self.segmentationComboBox.connect('currentIndexChanged(int)', self.onSegmentationChanged)
    self.segmentationOutlineSlider.connect('valueChanged(double)', self.updateSegmentationOutline)
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

  def historyItemChanged(self):
    self.historyList.setCurrentRow(int(self.parameterNode.GetParameter("currentLayer"))) # keep same value

  def onSliceIntersectionVisibilityChanged(self, state):
    viewNodes = slicer.util.getNodesByClass('vtkMRMLSliceCompositeNode')
    for viewNode in viewNodes:
      viewNode.SetSliceIntersectionVisibility(state)


  def onSegmentationChanged(self,index):
    if index != 0:
      MNIPath = self.parameterNode.GetParameter("MNIPath")
      atlasPath = os.path.join(MNIPath,'atlases',self.segmentationComboBox.itemText(index))
      modelID = self.parameterNode.GetParameter("modelID")
      if modelID != "":
        modelHierarchyNode = slicer.util.getNode(modelID)
        FunctionsUtil.loadAtlas(atlasPath, modelHierarchyNode)
        self.segmentationNodeChanged()
  
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
    modelNode = slicer.util.getNode(self.parameterNode.GetParameter("modelID"))
    item.setSelected(1)
    seg = slicer.util.getNode(item.text())
    pd = seg.GetPolyData()
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
    modelNode = slicer.util.getNode(self.parameterNode.GetParameter("modelID"))
    col = vtk.vtkCollection()
    modelNode.GetChildrenModelNodes(col)
    for i in range(col.GetNumberOfItems()):
      col.GetItemAsObject(i).GetDisplayNode().SetSliceIntersectionVisibility(self.segmentationList.item(i).isSelected())


  def updateSegmentationOutline(self,value):
    modelID = self.parameterNode.GetParameter("modelID")
    if modelID != "":
      modelNode = slicer.util.getNode(modelID)
      col = vtk.vtkCollection()
      modelNode.GetChildrenModelNodes(col)
      for i in range(col.GetNumberOfItems()):
        col.GetItemAsObject(i).GetDisplayNode().SetSliceIntersectionOpacity(value)


  def segmentationNodeChanged(self):
    self.segmentationList.clear()
    self.updateGuiFromMRML()
    modelID = self.parameterNode.GetParameter("modelID")
    if modelID != "":
      self.updateSegmentationOutline(self.segmentationOutlineSlider.value)
    
  def onUndoButton(self):
    SmudgeModuleLogic().undoOperation()

  def onRedoButton(self):
    SmudgeModuleLogic().redoOperation()
   
  def onFlattenButton(self):
    SmudgeModuleLogic().flattenTransform()
  
  def onSaveButton(self):
    SmudgeModuleLogic().applyChanges()
    nextSubjectN = int(self.parameterNode.GetParameter("subjectN"))+1
    subjectPaths = self.parameterNode.GetParameter("subjectPaths").split(' ')
    
    # initialize multiple old lead folders
    #for nextSubjectN in range(1,len(subjectPaths)):
    #  self.parameterNode.SetParameter("subjectN", str(nextSubjectN))
    #  self.parameterNode.SetParameter("subjectPath", subjectPaths[nextSubjectN])
    #  SmudgeModuleLogic().resetSubjectData()
    #slicer.util.exit() 
    
    if nextSubjectN < len(subjectPaths):
      self.parameterNode.SetParameter("subjectN", str(nextSubjectN))
      self.parameterNode.SetParameter("subjectPath", subjectPaths[nextSubjectN])
      self.exit()
      SmudgeModuleLogic().resetSubjectData()
      # reset history
      self.parameterNode.SetParameter("currentLayer","-1")
      self.parameterNode.SetParameter("currentLayer","0")
      SmudgeModuleLogic().initialize()
      self.updateGuiFromMRML()
    else:
      slicer.util.exit()

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
      # initialize segmentaion
      modelID = self.parameterNode.GetParameter("modelID")
      if modelID != "":
        modelNode = slicer.util.getNode(modelID)
        childrenCollection = vtk.vtkCollection()
        modelNode.GetChildrenModelNodes(childrenCollection)
        #modelNode.CreateDefaultDisplayNodes()
        for i in range(childrenCollection.GetNumberOfItems()):
          self.segmentationList.addItem(childrenCollection.GetItemAsObject(i).GetName())
          childrenCollection.GetItemAsObject(i).GetDisplayNode().SetOpacity(0) # hide from 3d viewer
    if not self.segmentationComboBox.itemText(1):
      # initialize combo box
      MNIPath = self.parameterNode.GetParameter("MNIPath")
      if MNIPath != "":
        atlasesPath = os.path.join(MNIPath,'atlases')
        atlases = sorted(os.listdir(atlasesPath))
        for atlas in atlases:
          if atlas[0] != '.':
            self.segmentationComboBox.addItem(atlas)
    self.smudgeButton.setEnabled(bool(int(self.parameterNode.GetParameter("initialized"))))
    self.snapButton.setEnabled(bool(int(self.parameterNode.GetParameter("initialized"))))
    self.flattenButton.setEnabled(bool(int(self.parameterNode.GetParameter("initialized"))))
    self.undoButton.setEnabled(bool(int(self.parameterNode.GetParameter("initialized"))))
    self.initializeButton.setEnabled(self.inputSelector.currentNode() and self.transformSelector.currentNode() and self.affineSelector.currentNode())
    # add history layer
    if self.historyList.model().rowCount()-1 < int(self.parameterNode.GetParameter("currentLayer")):
      self.historyList.addItem("Transform Layer " + self.parameterNode.GetParameter("currentLayer"))
    # remove history layers when current layer is -1 (flatten transform)
    if int(self.parameterNode.GetParameter("currentLayer")) == -1:
      while self.historyList.model().rowCount() > 1:
        self.historyList.item(self.historyList.model().rowCount()-1).delete()
    # set current row as current layer
    self.historyList.setCurrentRow(int(self.parameterNode.GetParameter("currentLayer")))
    # remove history layer when undoing more than once
    if self.historyList.currentRow == self.historyList.model().rowCount()-3:
      self.historyList.item(self.historyList.model().rowCount()-1).delete()
    
    # subject text
    subjectN = int(self.parameterNode.GetParameter("subjectN"))
    subjectPaths = self.parameterNode.GetParameter("subjectPaths").split(' ')
    self.subjectNameLabel.text = os.path.split(subjectPaths[subjectN])[-1]
    self.saveButton.text = 'Save and Exit' if subjectN == len(subjectPaths)-1 else 'Save and Next'

  def updateMRMLFromGUI(self):
    self.parameterNode.SetParameter("radius", str(self.radiusSlider.value) )
    self.parameterNode.SetParameter("blurr", str(self.blurrSlider.value) )
    self.parameterNode.SetParameter("hardness", str(self.hardnessSlider.value) )
    self.parameterNode.SetParameter("modality", self.modalityComboBox.currentText.lower() )
    if not __name__ == "__main__":
      if self.transformSelector.currentNode():
        self.parameterNode.SetParameter("transformID", self.transformSelector.currentNode().GetID())
      if self.inputSelector.currentNode():
        self.parameterNode.SetParameter("imageID", self.inputSelector.currentNode().GetID())
      if self.segmentationSelector.currentNode():
        self.parameterNode.SetParameter("segmentationID", self.segmentationSelector.currentNode().GetID())
      if self.affineSelector.currentNode():
        self.parameterNode.SetParameter("affineTransformID", self.affineSelector.currentNode().GetID())
      
    
  def toogleTools(self):
    WarpEffect.WarpEffectTool.empty()

  def onSmudgeButton(self, buttonDown):
    if buttonDown:
      self.snapButton.setChecked(False)
      SmudgeModuleLogic().smudgeOn()
    else:
      SmudgeModuleLogic().smudgeOff()

  def onSnapButton(self, buttonDown):
    if buttonDown:
      self.smudgeButton.setChecked(False)
      SmudgeModuleLogic().snapOn()
    else:
      SmudgeModuleLogic().snapOff()
    
      
  def exit(self):
    SmudgeModuleLogic().smudgeOff()
    SmudgeModuleLogic().snapOff()
    SmudgeModuleLogic().undoInitialize()
    self.smudgeButton.setChecked(False)
    self.snapButton.setChecked(False)

  def onSceneStartClose(self, caller, event):
    self.parameterNode.RemoveAllObservers()
    self.exit()

  


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
      for auxID in ["auxTransformID"]:
        if parameterNode.GetParameter(auxID) == "":
          transformNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLGridTransformNode')
          FunctionsUtil.emptyTransform(transformNode)
          parameterNode.SetParameter(auxID, transformNode.GetID())
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
    node.SetParameter("modelID","")
    node.SetParameter("radius", "10")
    node.SetParameter("blurr", "60")
    node.SetParameter("hardness", "100")
    node.SetParameter("modality", "t1")
    node.SetParameter("initialized", "0")
    node.SetParameter("subjectPath", "")
    node.SetParameter("subjectPath", "")
    node.SetParameter("subjectN", "0")
    node.SetParameter("MNIPath", "")
    node.SetParameter("antsApplyTransformsPath", "")
    node.SetParameter("currentLayer", "0")
    return node

  def smudgeOn(self):

    parameterNode = self.getParameterNode()
    
    if parameterNode.GetParameter("auxTransformID") == "":
      transformNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLGridTransformNode')
      FunctionsUtil.emptyTransform(transformNode)
      parameterNode.SetParameter("auxTransformID", transformNode.GetID())

    for color in ['Red','Green','Yellow']:
      sliceWidget = slicer.app.layoutManager().sliceWidget(color)
      WarpEffect.SmudgeEffectTool(parameterNode,sliceWidget)
  
  def smudgeOff(self):
    parameterNode = self.getParameterNode()
    WarpEffect.WarpEffectTool.empty()
    auxTransformID = parameterNode.GetParameter("auxTransformID")
    if auxTransformID != "":
      slicer.mrmlScene.RemoveNode(slicer.util.getNode(auxTransformID))
      parameterNode.SetParameter("auxTransformID", "")
    redoTransformID = parameterNode.GetParameter("redoTransformID")
    if redoTransformID != "":
      slicer.mrmlScene.RemoveNode(slicer.util.getNode(redoTransformID))
      parameterNode.SetParameter("redoTransformID", "")

  def snapOn(self):
    parameterNode = self.getParameterNode()
    for color in ['Red','Green','Yellow']:
      sliceWidget = slicer.app.layoutManager().sliceWidget(color)
      WarpEffect.SnapEffectTool(parameterNode, sliceWidget)
  
  def snapOff(self):
    WarpEffect.WarpEffectTool.empty()
    
  
  def applyChanges(self):
    parameterNode = self.getParameterNode()
    if int(parameterNode.GetParameter("currentLayer")) == 0: # no changes
      msgBox = qt.QMessageBox()
      msgBox.setText('No modifications in warp')
      msgBox.setInformativeText('Save subject as approved?')
      msgBox.setStandardButtons(qt.QMessageBox().Save | qt.QMessageBox().Discard)
      ret = msgBox.exec_()
      if ret == qt.QMessageBox().Save:
        FunctionsUtil.saveApprovedData(parameterNode.GetParameter("subjectPath"))
      return
    tmpFolder, displacementFileName = self.flattenTransform(cleanup=False)
    inverseDisplacementFileName = self.invertTransform(displacementFileName)
    self.overwriteTransform(displacementFileName, inverse = False)
    self.overwriteTransform(inverseDisplacementFileName, inverse = True)

    shutil.rmtree(tmpFolder)

  def overwriteTransform(self, transformFileName, inverse):
    parameterNode = self.getParameterNode()
    subjectPath = parameterNode.GetParameter("subjectPath")
    if not inverse:
      originalTransformFileName = os.path.join(subjectPath,'glanatComposite.nii.gz')
      referenceFileName = os.path.join(subjectPath, 'glanat.nii')
    else:
      originalTransformFileName = os.path.join(subjectPath, 'glanatInverseComposite.nii.gz')
      referenceFileName = os.path.join(subjectPath, 'anat_t1.nii')
    
    command = parameterNode.GetParameter("antsApplyTransformsPath") + " -r " + referenceFileName + " "

    if not inverse:
      command = command + "-t " + transformFileName + " -t " + originalTransformFileName
    else:
      command = command + "-t " + originalTransformFileName + " -t " + transformFileName

    command = command + " -o [" + originalTransformFileName + ",1] -v 1"

    # buckup
    backupFile = os.path.join(subjectPath, "bu_" + str(int(inverse)) + "_composite.nii.gz")
    if not os.path.isfile(backupFile) and False: # buckup transform before overwrite
      shutil.copyfile(originalTransformFileName,backupFile)

    commandOut = call(command, env=slicer.util.startupEnvironment(), shell=True) # run antsApplyTransforms


  def invertTransform(self, transformFileName):
    reader = sitk.ImageFileReader()
    reader.SetOutputPixelType(sitk.sitkVectorFloat64)
    reader.SetFileName(transformFileName)
    displacementImage = reader.Execute()
        
    inverseDisplacement = sitk.InvertDisplacementField(displacementImage,\
                                                        maximumNumberOfIterations = 20,\
                                                        maxErrorToleranceThreshold = 0.01,\
                                                        meanErrorToleranceThreshold = 0.0001,\
                                                        enforceBoundaryCondition = True)                
        
    inverseDisplacementFileName = os.path.join(os.path.dirname(transformFileName), 'inverse.nii.gz')
    
    writer = sitk.ImageFileWriter()
    writer.SetFileName(inverseDisplacementFileName)
    writer.Execute(inverseDisplacement)

    return inverseDisplacementFileName

  def flattenTransform(self, cleanup=True):
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

    if int(parameterNode.GetParameter("currentLayer")) == 0: # already flat!
      return

    path,transformFilename = os.path.split(transformNode.GetStorageNode().GetFileName()) # get working dir
    tmpFolder = os.path.join(path,"tmp")
    os.mkdir(tmpFolder) # create temp folder

    command = ""

    for i in range(1,int(parameterNode.GetParameter("currentLayer"))+1): # save each transform and append name to command
      con = vtk.vtkWarpTransform.SafeDownCast(col.GetItemAsObject(i))
      auxTransformNode.SetAndObserveTransformToParent(con)
      fullname = os.path.join(tmpFolder, "layer" + str(i) + ".nii.gz")
      slicer.util.saveNode(auxTransformNode,fullname)
      command = "-t " + fullname + " " + command

    FunctionsUtil.emptyTransform(auxTransformNode) # reset aux transform

    command = parameterNode.GetParameter("antsApplyTransformsPath") + " -r " + fullname + " " + command

    outname = os.path.join(tmpFolder, "foreward.nii.gz") # flat transform name
    command = command + "-o [" + outname + ",1] --verbose 1"

    commandOut = call(command, env=slicer.util.startupEnvironment(), shell=True) # run antsApplyTransforms

    if cleanup:
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
      shutil.rmtree(tmpFolder)
      # reset history
      parameterNode.SetParameter("currentLayer","-1")
      parameterNode.SetParameter("currentLayer","1")
      print('Flattened!')
    else:
      return tmpFolder,outname
    

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

    FunctionsUtil.emptyTransform(auxTransformNode) # reset aux transform

    # save last op for redo
    redoTransformID = parameterNode.GetParameter("redoTransformID")
    if redoTransformID == "":
      redoTransformNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLGridTransformNode')
      FunctionsUtil.emptyTransform(redoTransformNode)
      parameterNode.SetParameter("redoTransformID", redoTransformNode.GetID())
    redoTransformNode =  slicer.util.getNode(parameterNode.GetParameter("redoTransformID"))
    con = vtk.vtkWarpTransform.SafeDownCast(col.GetItemAsObject(col.GetNumberOfItems()-1))
    redoTransformNode.SetAndObserveTransformToParent(con)

    transformNode.SetAndObserveTransformNodeID(parameterNode.GetParameter("auxTransformID")) 

    parameterNode.SetParameter("currentLayer",str(int(parameterNode.GetParameter("currentLayer"))-1))
  
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
    parameterNode.SetParameter("currentLayer",str(int(parameterNode.GetParameter("currentLayer"))+1))


  def resetSubjectData(self):
    parameterNode = self.getParameterNode()
    #parameterNode.RemoveAllObservers()
    # delete nodes
    imageID = parameterNode.GetParameter("imageID")
    transformID = parameterNode.GetParameter("transformID")
    affineTransformID = parameterNode.GetParameter("affineTransformID")

    removeIDs = [imageID,transformID,affineTransformID] 
   
    subjectPath = parameterNode.GetParameter("subjectPath")
    # check for old version of transforms
    if os.path.isfile(os.path.join(subjectPath,'glanatComposite.h5')):
      self.updateTransforms(subjectPath)

    # load nodes
    imagePath = os.path.join(subjectPath,'anat_' + parameterNode.GetParameter("modality") + '.nii')
    if not os.path.isfile(imagePath):
      imagePath = os.path.join(subjectPath,'anat_t1.nii')
    transformPath = os.path.join(subjectPath,'glanatComposite.nii.gz')
    affinePath = os.path.join(subjectPath,'glanat0GenericAffine_backup.mat')

    ls, transformNode = slicer.util.loadTransform(transformPath, returnNode = True)
    ls, imageNode = slicer.util.loadVolume(imagePath, properties = {'name':'anat_t1', 'show':False}, returnNode = True)
    if not os.path.isfile(affinePath):
      msg=qt.QMessageBox().warning(qt.QWidget(),'','Affine transform not available. The Displayed warp is affine+deformable')
      affineNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLinearTransformNode')
    else:
      ls, affineNode = slicer.util.loadTransform(affinePath, returnNode = True)

    # set parameters    
    parameterNode.SetParameter("transformID", transformNode.GetID())
    parameterNode.SetParameter("affineTransformID", affineNode.GetID())
    parameterNode.SetParameter("imageID", imageNode.GetID())

    # init
    slicer.util.setSliceViewerLayers(background=imageNode)
    imageNode.GetDisplayNode().AutoWindowLevelOff()
    imageNode.GetDisplayNode().SetWindow(760)
    imageNode.GetDisplayNode().SetLevel(420)

    if int(parameterNode.GetParameter("subjectN")) != 0:
      for ID in removeIDs:
        if ID != "":
          slicer.mrmlScene.RemoveNode(slicer.util.getNode(ID))


  def updateTransforms(self,subjectPath):
    # get affine parameters from h5
    with h5py.File(os.path.join(subjectPath,'glanatComposite.h5'),'r') as f:
      try:
        parameters  = f['TransformGroup/1']['TransformParameters'][()]
        fixedParameters = f['TransformGroup/1']['TransformFixedParameters'][()]
      except:
        parameters  = f['TransformGroup/1']['TranformParameters'][()]
        fixedParameters = f['TransformGroup/1']['TranformFixedParameters'][()]
    # save parameters to .txt transform
    tmpTransform = os.path.join(subjectPath,'tmpTransform.txt')
    file1 = open(tmpTransform,"w") 
    L = ["#Insight Transform File V1.0\n",\
    	 "#Transform 0\n",\
    	 "Transform: AffineTransform_double_3_3\n",\
    	 "Parameters: {} {} {} {} {} {} {} {} {} {} {} {}\n".format(*parameters),\
    	 "FixedParameters: {} {} {}\n".format(*fixedParameters)]
    file1.writelines(L) 
    file1.close()
    # load transform and save as a .mat
    ls, tmpTransformNode = slicer.util.loadTransform(tmpTransform, returnNode = True)
    slicer.util.saveNode(tmpTransformNode, os.path.join(subjectPath,'glanat0GenericAffine_backup.mat'))
    slicer.mrmlScene.RemoveNode(tmpTransformNode)
    os.remove(tmpTransform)
    # use ants apply transforms to change .h5 ext to .nii.gz
    for transform,reference in zip(['glanatComposite','glanatInverseComposite'],['glanat','anat_t1']):
      command = antsApplyTransformsPath + " -r " + os.path.join(subjectPath,reference + '.nii') + " -t " + os.path.join(subjectPath,transform + '.h5') + " -o [" + os.path.join(subjectPath,transform + '.nii.gz') + ",1] -v 1"
      commandOut = call(command, env=slicer.util.startupEnvironment(), shell=True) # run antsApplyTransforms
      os.remove(os.path.join(subjectPath,transform + '.h5'))


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
    pass




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

    self.parent.setGeometry(51, 45, 1869, 849)
    #self.parent.showMaximized()
    self.parent.show()



if __name__ == "__main__":

  import sys
  print( sys.argv )

  slicelet = SmudgeModuleSlicelet()

  # load subject data
  if len(sys.argv) > 1:
    
    leadPath = sys.argv[1]
    subjectPaths = ' '.join(sys.argv[2:])

    MNIPath = os.path.join(leadPath,'templates','space','MNI_ICBM_2009b_NLIN_ASYM')

    if sys.platform == "darwin":
      ext = "maci64"

    # TODO: other platforms

    antsApplyTransformsPath = os.path.join(leadPath,'ext_libs','ANTs','antsApplyTransforms.' + ext)

    # template data
    templatePath = os.path.join(MNIPath,'t1.nii')
    ls, templateNode = slicer.util.loadVolume(templatePath, properties = {'name':'template_t1', 'show':False}, returnNode = True)
    slicer.util.setSliceViewerLayers(foreground=templateNode)
    templateNode.GetDisplayNode().AutoWindowLevelOff()
    templateNode.GetDisplayNode().SetWindow(100)
    templateNode.GetDisplayNode().SetLevel(70)

    # load model
    modelHierarchyNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelHierarchyNode')
    FunctionsUtil.loadAtlas(os.path.join(MNIPath,'atlases','DISTAL Minimal (Ewert 2017)'), modelHierarchyNode)

    parameterNode = SmudgeModuleLogic().getParameterNode()

    parameterNode.SetParameter("subjectPaths", subjectPaths)
    parameterNode.SetParameter("subjectN", "0")
    parameterNode.SetParameter("subjectPath", subjectPaths.split(' ')[0])
    parameterNode.SetParameter("MNIPath", MNIPath)
    parameterNode.SetParameter("antsApplyTransformsPath", antsApplyTransformsPath)
    parameterNode.SetParameter("modelID", modelHierarchyNode.GetID())
    parameterNode.SetParameter("templateID", templateNode.GetID())

    SmudgeModuleLogic().resetSubjectData()
    SmudgeModuleLogic().initialize()

    


    





