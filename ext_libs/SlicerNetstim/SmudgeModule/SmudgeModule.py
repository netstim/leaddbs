import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging

import sys
import numpy as np
import shutil
from math import sqrt
from slicer.util import VTKObservationMixin
import glob
import SimpleITK as sitk

# import other netstim modules
from Helpers import WarpEffect, FunctionsUtil, Toolbar

# other netstim modules
import TransformsUtil
import ImportSubject
import ImportAtlas

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

    # Init parameter node
    self.parameterNode = SmudgeModuleLogic().getParameterNode()
    self.addObserver(self.parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGuiFromMRML)


    # Instantiate and connect widgets ...


    #
    # Inputs Area
    #
    self.inputsCollapsibleButton = ctk.ctkCollapsibleButton()
    self.inputsCollapsibleButton.text = "Inputs"
    self.layout.addWidget(self.inputsCollapsibleButton)

    # Layout within the dummy collapsible button
    inputsFormLayout = qt.QFormLayout(self.inputsCollapsibleButton)

    #
    # warp selector
    #
    self.warpSelector = slicer.qMRMLNodeComboBox()
    self.warpSelector.nodeTypes = ["vtkMRMLGridTransformNode"]
    self.warpSelector.selectNodeUponCreation = False
    self.warpSelector.addEnabled = False
    self.warpSelector.removeEnabled = False
    self.warpSelector.noneEnabled = True
    self.warpSelector.showHidden = False
    self.warpSelector.showChildNodeTypes = False
    self.warpSelector.setMRMLScene( slicer.mrmlScene )
    self.warpSelector.setToolTip( "Pick the warp to refine." )
    inputsFormLayout.addRow("Warp: ", self.warpSelector)


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

    nonePixmap = qt.QPixmap(self.resourcePath(os.path.join('Icons','cursor.png')))
    noneIcon = qt.QIcon(nonePixmap)
    self.noneButton = qt.QPushButton()
    self.noneButton.setIcon(noneIcon)
    self.noneButton.setIconSize(nonePixmap.rect().size())
    self.noneButton.setAutoExclusive(True)
    self.noneButton.setCheckable(True)
    self.noneButton.setChecked(True)
    optionsFrame.layout().addWidget(self.noneButton)
    
    smudgePixmap = qt.QPixmap(self.resourcePath(os.path.join('Icons','smudgeIcon.png')))
    smudgeIcon = qt.QIcon(smudgePixmap)
    self.smudgeButton = qt.QPushButton()
    self.smudgeButton.setIcon(smudgeIcon)
    self.smudgeButton.setIconSize(smudgePixmap.rect().size())
    self.smudgeButton.setCheckable(True)
    self.smudgeButton.setEnabled(False)
    self.smudgeButton.setAutoExclusive(True)
    optionsFrame.layout().addWidget(self.smudgeButton)

    pencilPixmap = qt.QPixmap(self.resourcePath(os.path.join('Icons','pencilIcon.png')))
    pencilIcon = qt.QIcon(pencilPixmap)
    self.snapButton = qt.QPushButton()
    self.snapButton.setIcon(pencilIcon)
    self.snapButton.setIconSize(pencilPixmap.rect().size())
    self.snapButton.setCheckable(True)
    self.snapButton.setEnabled(False)
    self.snapButton.setAutoExclusive(True)
    optionsFrame.layout().addWidget(self.snapButton)

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
    self.historyList.setCurrentRow(0)
    self.historyList.setMaximumHeight(redoPixmap.rect().size().height() * 3.3)

    historyGridLayout.addWidget(self.historyList,0,1)

    #
    # View Area
    #
    viewCollapsibleButton = ctk.ctkCollapsibleButton()
    viewCollapsibleButton.text = "View"
    viewCollapsibleButton.collapsed=True
    self.layout.addWidget(viewCollapsibleButton)

    # Layout within the dummy collapsible button
    viewFormLayout = qt.QFormLayout(viewCollapsibleButton)

    #
    # Display
    #


    #
    # Modality
    #
    self.modalityComboBox = qt.QComboBox()
    self.modalityComboBox.addItems(['T1','T2'])
    self.modalityComboBox.setEnabled(False)
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
    viewFormLayout.addRow("&Template Image:", self.templateSlider)

    #
    # Atlas
    #


    self.atlasesComboBox = qt.QComboBox()
    self.atlasesComboBox.addItem('Choose Atlas')
    self.atlasesComboBox.setMaximumWidth(350)
    viewFormLayout.addRow("Import Atlas:", self.atlasesComboBox)


    #
    # Modles Area
    #

    modelsCollapsibleButton = ctk.ctkCollapsibleButton()
    modelsCollapsibleButton.text = "Model Control"
    self.layout.addWidget(modelsCollapsibleButton, 1)

    modelsGridLayout = qt.QGridLayout(modelsCollapsibleButton)    

    self.dataTreeWidget = slicer.qMRMLSubjectHierarchyTreeView(slicer.util.mainWindow())
    self.dataTreeWidget.setMRMLScene(slicer.mrmlScene)
    self.dataTreeWidget.setColumnHidden(self.dataTreeWidget.model().idColumn, True)
    self.dataTreeWidget.setColumnHidden(self.dataTreeWidget.model().transformColumn, True)
    self.dataTreeWidget.nodeTypes = ('vtkMRMLModelNode','vtkMRMLFolderDisplayNode')

    modelsGridLayout.addWidget(self.dataTreeWidget,0,0)

    self.layout.addStretch(0)

    #
    # Save Area
    #

    self.saveButton = qt.QPushButton("Save and Next")
    self.saveButton.setMinimumHeight(30)
    self.saveButton.setStyleSheet("background-color: green")
    self.saveButton.setVisible(False)
    self.layout.addWidget(self.saveButton)

    # connections
    self.warpSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.exit) # deselect effect
    self.warpSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateMRMLFromGUI)
    self.smudgeButton.toggled.connect(self.toogleTools)
    self.snapButton.toggled.connect(self.toogleTools)
    self.smudgeButton.connect('clicked(bool)', self.onSmudgeButton)
    self.snapButton.connect('clicked(bool)', self.onSnapButton)
    
    self.radiusSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.blurrSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.hardnessSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.flattenButton.connect("clicked(bool)", self.onFlattenButton)
    self.saveButton.connect("clicked(bool)", self.onSaveButton)
    self.undoButton.connect("clicked(bool)", self.onUndoButton)
    self.redoButton.connect("clicked(bool)", self.onRedoButton)
    self.historyList.itemSelectionChanged.connect(self.historyItemChanged)

    self.modalityComboBox.connect('currentIndexChanged(int)', self.onModalityChanged)
    self.templateSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.templateSlider.connect('valueChanged(double)', self.updateTemplateView)
    self.atlasesComboBox.connect('currentIndexChanged(int)', self.onAtlasChanged)
    self.dataTreeWidget.doubleClicked.connect(self.onDataTreeDoubleClicked)

    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)    
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.NodeAddedEvent, self.onSceneNodeAdded)    
    interactionNode = slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLInteractionNode')
    self.addObserver(interactionNode, interactionNode.InteractionModeChangedEvent, self.onInteractionModeChanged) 


    # Refresh
    qt.QApplication.processEvents()
    if self.updateMRMLFromArgs(): # was called from command line
      SmudgeModuleLogic().loadSubjectData(anat=False)
      self.onModalityChanged(0)
      self.atlasesComboBox.addItems(ImportAtlas.ImportAtlasLogic().getValidAtlases(self.parameterNode.GetParameter("MNIAtlasPath")))
      self.onAtlasChanged(1,'DISTAL Minimal (Ewert 2017)')
      self.showSingleModule()
      tb = Toolbar.reducedToolbar()
      slicer.util.mainWindow().addToolBar(tb)

    self.updateGuiFromMRML()  
    self.onSceneNodeAdded()
    self.toogleTools()

  

  #
  # Methods
  #
  
  def showSingleModule(self):
    
    singleModule = True

    keepToolbars = []
    slicer.util.setToolbarsVisible(not singleModule, keepToolbars)
    slicer.util.setMenuBarsVisible(not singleModule)
    slicer.util.setApplicationLogoVisible(not singleModule)
    slicer.util.setModuleHelpSectionVisible(not singleModule)
    slicer.util.setModulePanelTitleVisible(not singleModule)
    slicer.util.setDataProbeVisible(not singleModule)
    #slicer.util.setViewControllersVisible(not singleModule)

    if singleModule:
      slicer.util.setPythonConsoleVisible(False)

    self.saveButton.setVisible(singleModule)
    self.modalityComboBox.setEnabled(singleModule)
    self.inputsCollapsibleButton.setVisible(not singleModule)
    self.reloadCollapsibleButton.setVisible(not singleModule)

    slicer.util.mainWindow().setWindowTitle("Name goes here")


  def onSaveButton(self):
    self.exit()
    SmudgeModuleLogic().applyChanges()

    # remove nodes
    slicer.mrmlScene.RemoveNode(slicer.util.getNode(self.parameterNode.GetParameter("affineTransformID")))
    slicer.mrmlScene.RemoveNode(slicer.util.getNode(self.parameterNode.GetParameter("warpID")))
    layoutManager = slicer.app.layoutManager()
    compositeNode = layoutManager.sliceWidget('Red').sliceLogic().GetSliceCompositeNode()
    slicer.mrmlScene.RemoveNode(slicer.util.getNode(compositeNode.GetBackgroundVolumeID()))

    nextSubjectN = int(self.parameterNode.GetParameter("subjectN"))+1
    subjectPaths = self.parameterNode.GetParameter("subjectPaths").split(' ')
    
    if nextSubjectN < len(subjectPaths):
      self.parameterNode.SetParameter("subjectN", str(nextSubjectN))
      self.parameterNode.SetParameter("subjectPath", subjectPaths[nextSubjectN])
      imageNode = SmudgeModuleLogic().loadSubjectData()
      SmudgeModuleLogic().initialize(imageNode)
      slicer.util.setSliceViewerLayers(background=imageNode.GetID())
      self.updateGuiFromMRML()
      self.toogleTools()
    else:
      slicer.util.exit()


  def onModalityChanged(self, index):
    modality = self.modalityComboBox.itemText(index).lower() # get modality
    self.parameterNode.SetParameter("modality",modality)
    # find old node and delete
    layoutManager = slicer.app.layoutManager()
    compositeNode = layoutManager.sliceWidget('Red').sliceLogic().GetSliceCompositeNode()
    if compositeNode.GetBackgroundVolumeID():
      slicer.mrmlScene.RemoveNode(slicer.util.getNode(compositeNode.GetBackgroundVolumeID()))
    # initialize new image and init
    imageNode = SmudgeModuleLogic().loadSubjectData(anat=True, warp=False, affine=False)
    SmudgeModuleLogic().initialize(imageNode)
    slicer.util.setSliceViewerLayers(background=imageNode.GetID())
    # load new temaplate image
    if compositeNode.GetForegroundVolumeID():
      slicer.mrmlScene.RemoveNode(slicer.util.getNode(compositeNode.GetForegroundVolumeID()))
    templateNode = slicer.util.loadVolume(os.path.join(self.parameterNode.GetParameter("MNIPath"), modality + ".nii"), properties={'show':False})
    templateNode.GetDisplayNode().AutoWindowLevelOff()
    templateNode.GetDisplayNode().SetWindow(100)
    templateNode.GetDisplayNode().SetLevel(70)
    slicer.util.setSliceViewerLayers(foreground=templateNode.GetID())


  def updateStructures(self,value):
    modelsInScene = slicer.mrmlScene.GetNodesByClass('vtkMRMLModelNode')
    for i in range(modelsInScene.GetNumberOfItems()):
      modelNode = modelsInScene.GetItemAsObject(i)
      try:
        modelNode.GetDisplayNode().SetSliceIntersectionOpacity(value)
      except:
        pass

  def updateTemplateView(self,value):
    slicer.util.setSliceViewerLayers(foregroundOpacity=value)

  def updateMRMLFromArgs(self): 
    args = sys.argv
    if (len(sys.argv) > 2) and os.path.isfile(os.path.join(sys.argv[1],'lead.m')):
      subjectPaths = ' '.join(sys.argv[2:])
      subjectPath = subjectPaths.split(' ')[0]
      MNIPath = os.path.join(sys.argv[1],'templates','space','MNI_ICBM_2009b_NLIN_ASYM')
      MNIAtlasPath = os.path.join(MNIPath,'atlases')
      if sys.platform == "darwin":
        ext = "maci64"
      elif sys.platform.startswith('win'):
        ext = 'exe'
      else:
        ext = 'glnxa64'
      antsApplyTransformsPath = os.path.join(sys.argv[1],'ext_libs','ANTs','antsApplyTransforms.' + ext)
      # set parameter node
      self.parameterNode.SetParameter("subjectPaths", subjectPaths)
      self.parameterNode.SetParameter("subjectN", "0")
      self.parameterNode.SetParameter("subjectPath", subjectPath)
      self.parameterNode.SetParameter("MNIPath", MNIPath)
      self.parameterNode.SetParameter("MNIAtlasPath", MNIAtlasPath)
      self.parameterNode.SetParameter("antsApplyTransformsPath", antsApplyTransformsPath)
      return True
    else:
      return False

  def updateGuiFromMRML(self,caller=None,event=None):
    # radius - blur - hardness
    radius = float(self.parameterNode.GetParameter("radius"))
    self.radiusSlider.setValue( radius )
    if radius < self.radiusSlider.minimum or radius > self.radiusSlider.maximum:
      self.updateMRMLFromGUI()
    self.blurrSlider.setValue(float(self.parameterNode.GetParameter("blurr")))
    self.hardnessSlider.setValue(float(self.parameterNode.GetParameter("hardness")))
    # get warp node and set selector and buttons
    warpID = self.parameterNode.GetParameter("warpID")
    warpNode = slicer.util.getNode(warpID) if warpID != "" else None
    warpNumberOfComponents = TransformsUtil.TransformsUtilLogic().getNumberOfLayers(warpNode)
    self.warpSelector.setCurrentNode(warpNode)
    self.smudgeButton.enabled = bool(warpNode)
    self.snapButton.enabled = bool(warpNode)
    self.historyList.enabled = bool(warpNode) 
    # undo redo button
    self.undoButton.setEnabled(warpNumberOfComponents > 1) 
    self.redoButton.setEnabled(self.parameterNode.GetParameter("redoTransformID") != "") 
    self.flattenButton.setEnabled(warpNumberOfComponents > 2)
    # update warp history. add number of warp component plus redo transform (if available)
    self.historyList.clear()
    self.historyList.addItems(['Component ' + str(i) for i in range(warpNumberOfComponents + int(self.redoButton.enabled))])
    self.historyList.setCurrentRow(warpNumberOfComponents - 1)
    # structure outline - temaplate
    self.updateStructures(float(self.parameterNode.GetParameter("structureOutline")))
    self.templateSlider.value = float(self.parameterNode.GetParameter("templateOpacity"))
    # subject text
    subjectN = int(self.parameterNode.GetParameter("subjectN"))
    subjectPaths = self.parameterNode.GetParameter("subjectPaths").split(' ')
    self.saveButton.text = 'Save and Exit' if subjectN == len(subjectPaths)-1 else 'Save and Next'


  def updateMRMLFromGUI(self):
    self.parameterNode.SetParameter("radius", str(self.radiusSlider.value) )
    self.parameterNode.SetParameter("blurr", str(self.blurrSlider.value) )
    self.parameterNode.SetParameter("hardness", str(self.hardnessSlider.value) )
    self.parameterNode.SetParameter("warpID", self.warpSelector.currentNode().GetID() if self.warpSelector.currentNode() else "")
    self.parameterNode.SetParameter("templateOpacity", str(self.templateSlider.value))

  def historyItemChanged(self):
    warpID = self.parameterNode.GetParameter("warpID")
    if warpID != "":
      warpNode = slicer.util.getNode(warpID)
      self.historyList.setCurrentRow(TransformsUtil.TransformsUtilLogic().getNumberOfLayers(warpNode) - 1) # keep same value


  def iterItems(self, item):
    item.setFlags(33)
    for row in range(item.rowCount()):
      self.iterItems(item.child(row))

  def onSceneNodeAdded(self,caller=None,event=None):
    sceneItem = self.dataTreeWidget.model().item(0,0)
    self.iterItems(sceneItem)

  def onAtlasChanged(self, index, atlasName = None):
    if not atlasName:
      atlasName = self.atlasesComboBox.itemText(index)
    if index != 0:
      atlasPath = os.path.join(self.parameterNode.GetParameter("MNIAtlasPath"), atlasName)
      folderID = ImportAtlas.ImportAtlasLogic().run(atlasPath)


  def onDataTreeDoubleClicked(self, i):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    node = shNode.GetItemDataNode(self.dataTreeWidget.currentItem())
    if isinstance(node, slicer.vtkMRMLModelNode):
      # center
      pd = node.GetPolyData()
      center = vtk.vtkCenterOfMass()
      center.SetInputData(pd)
      center.Update()
      segCenter = center.GetCenter()
      # create markups node, add center as fiducial and jump and center slices
      markupsNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
      markupsNode.GetDisplayNode().SetVisibility(False)
      markupsNode.AddFiducialFromArray(np.array(segCenter),'')
      markupsLogic = slicer.modules.markups.logic()
      markupsLogic.JumpSlicesToNthPointInMarkup(markupsNode.GetID(),0,True)
      slicer.mrmlScene.RemoveNode(markupsNode)

  def onUndoButton(self):
    SmudgeModuleLogic().removeRedoTransform()
    redoTransformID = TransformsUtil.TransformsUtilLogic().removeLastLayer(slicer.util.getNode(self.parameterNode.GetParameter("warpID")))
    self.parameterNode.SetParameter("redoTransformID", redoTransformID)

  def onRedoButton(self):
    # see if smudging. aux transform causes problems here
    smudgeEnabled = self.smudgeButton.checked
    if smudgeEnabled:
      SmudgeModuleLogic().effectOff()
    # get nodes
    warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))
    redoTransformNode = slicer.util.getNode(self.parameterNode.GetParameter("redoTransformID"))
    # apply
    warpNode.SetAndObserveTransformNodeID(redoTransformNode.GetID())
    warpNode.HardenTransform()
    # delete redo transform
    slicer.mrmlScene.RemoveNode(redoTransformNode)
    self.parameterNode.SetParameter("redoTransformID","")
    # restore state
    if smudgeEnabled:
      SmudgeModuleLogic().effectOn('Smudge')

  def onFlattenButton(self):
    SmudgeModuleLogic().removeRedoTransform()
    warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))
    TransformsUtil.TransformsUtilLogic().flattenTransform(warpNode, False, [193,229,193],[-96.0, -132.0, -78.0],[1.0, 1.0, 1.0])
    self.updateGuiFromMRML() # update history

  def toogleTools(self):
    SmudgeModuleLogic().effectOff() # deactivate previous effect (if any)
    self.radiusSlider.setEnabled(False)
    self.blurrSlider.setEnabled(False)
    self.hardnessSlider.setEnabled(False)
    interactionNode = slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLInteractionNode')
    self.removeObserver(interactionNode, interactionNode.InteractionModeChangedEvent, self.onInteractionModeChanged) # so that changing mode doesnt affect
    if not self.noneButton.checked:
      interactionNode.SetCurrentInteractionMode(2)
      self.addObserver(interactionNode, interactionNode.InteractionModeChangedEvent, self.onInteractionModeChanged) 
    else:
      SmudgeModuleLogic().effectOn('None')

  def onSmudgeButton(self, buttonDown):
    if buttonDown:      
      self.radiusSlider.setEnabled(True)
      self.blurrSlider.setEnabled(True)
      self.hardnessSlider.setEnabled(True)
      SmudgeModuleLogic().effectOn('Smudge')

  def onSnapButton(self, buttonDown):
    if buttonDown:
      self.blurrSlider.setEnabled(True)
      SmudgeModuleLogic().effectOn('Snap')
      
  def exit(self):
    self.noneButton.setChecked(True)
    SmudgeModuleLogic().effectOff()

  def onInteractionModeChanged(self, caller, event):
    self.exit()

  def cleanup(self):
    self.exit()

  def enter(self):
    self.toogleTools()
      
  def onSceneStartClose(self, caller, event):
    pass



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


  def createParameterNode(self):
    node = ScriptedLoadableModuleLogic.createParameterNode(self)
    node.SetParameter("warpID", "")
    node.SetParameter("redoTransformID", "")
    node.SetParameter("affineTransformID", "")
    node.SetParameter("imageID", "")
    node.SetParameter("templateID", "")
    node.SetParameter("radius", "25")
    node.SetParameter("blurr", "60")
    node.SetParameter("hardness", "100")
    node.SetParameter("modality", "t1")
    node.SetParameter("subjectPath", "")
    node.SetParameter("subjectN", "0")
    node.SetParameter("MNIPath", ".")
    node.SetParameter("MNIAtlasPath", ".")
    node.SetParameter("antsApplyTransformsPath", "")
    node.SetParameter("warpModified", "1")
    node.SetParameter("structureOutline", "1")
    node.SetParameter("templateOpacity", "0")
    return node

  def removeRedoTransform(self):
    parameterNode = self.getParameterNode()
    redoTransformID = parameterNode.GetParameter("redoTransformID")
    if redoTransformID != "":
      slicer.mrmlScene.RemoveNode(slicer.util.getNode(redoTransformID))
      parameterNode.SetParameter("redoTransformID","")


  def loadSubjectData(self, anat=True, warp=True, affine=True):
    parameterNode = self.getParameterNode()
    subjectPath = parameterNode.GetParameter("subjectPath")

    if warp: # load warp
      if ImportSubject.ImportSubjectLogic().ish5Transform(subjectPath):
        ImportSubject.ImportSubjectLogic().updateTranform(subjectPath, parameterNode.GetParameter("antsApplyTransformsPath"))

      warpNode = ImportSubject.ImportSubjectLogic().importFile(subjectPath, 'glanatComposite.nii.gz')
      parameterNode.SetParameter("warpID", warpNode.GetID())

    if affine: # load affine
      affineNode = ImportSubject.ImportSubjectLogic().importFile(subjectPath, 'glanat0GenericAffine_backup.mat')
      if not affineNode:
        affineNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLinearTransformNode')
      affineNode.Inverse()
      affineNode.CreateDefaultDisplayNodes()
      parameterNode.SetParameter("affineTransformID", affineNode.GetID())
    if anat: # load image
      imageNode = ImportSubject.ImportSubjectLogic().importFile(subjectPath, 'anat_' + parameterNode.GetParameter("modality") + '.nii')
      if not imageNode:
        imageNode = ImportSubject.ImportSubjectLogic().importFile(subjectPath, 'anat_t1.nii')
      return imageNode


  def initialize(self, imageNode):
    parameterNode = self.getParameterNode()
    # apply affine transform to image
    affineNode = slicer.util.getNode(parameterNode.GetParameter("affineTransformID"))
    imageNode.ApplyTransform(affineNode.GetTransformFromParent())
    # set to image 
    imageNode.SetAndObserveTransformNodeID(parameterNode.GetParameter("affineTransformID"))
    # set transform to affine
    affineNode.SetAndObserveTransformNodeID(parameterNode.GetParameter("warpID"))


  def effectOn(self, effectName):
    
    if effectName == 'Smudge':
      auxTransformNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLGridTransformNode')
      TransformsUtil.TransformsUtilLogic().emptyTransfrom(auxTransformNode)

    parameterNode = self.getParameterNode()
    for color in ['Red','Green','Yellow']:
      sliceWidget = slicer.app.layoutManager().sliceWidget(color)
      if effectName == 'None':
        WarpEffect.NoneEffect(parameterNode, sliceWidget)
      elif effectName == 'Smudge':
        WarpEffect.SmudgeEffectTool(parameterNode, sliceWidget, auxTransformNode)
      elif effectName == 'Snap':
        WarpEffect.SnapEffectTool(parameterNode, sliceWidget)
  
  def effectOff(self):
    WarpEffect.WarpEffectTool.empty()

  def applyChanges(self):

    parameterNode = self.getParameterNode()
    warpNode = slicer.util.getNode(parameterNode.GetParameter("warpID"))
    subjectPath = parameterNode.GetParameter("subjectPath")

    if TransformsUtil.TransformsUtilLogic().getNumberOfLayers(warpNode) == 1:
      msgBox = qt.QMessageBox()
      msgBox.setText('No modifications in warp')
      msgBox.setInformativeText('Save subject as approved?')
      msgBox.setStandardButtons(qt.QMessageBox().Save | qt.QMessageBox().Discard)
      ret = msgBox.exec_()
      if ret == qt.QMessageBox().Save:
        FunctionsUtil.saveApprovedData(subjectPath)
      return
    
    TransformsUtil.TransformsUtilLogic().flattenTransform(warpNode, True, [], [], [])

    slicer.util.saveNode(warpNode, os.path.join(subjectPath,'glanatComposite.nii.gz'))

    warpNode.Inverse()

    # load the volume again to set as reference
    imageNode = slicer.util.loadVolume(os.path.join(subjectPath,'anat_t1.nii'))

    outNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
    transformsLogic = slicer.modules.transforms.logic()
    transformsLogic.ConvertToGridTransform(warpNode, imageNode, outNode)

    slicer.util.saveNode(outNode, os.path.join(subjectPath,'glanatInverseComposite.nii.gz'))

    # delete
    slicer.mrmlScene.RemoveNode(outNode)
    slicer.mrmlScene.RemoveNode(imageNode)







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




    





