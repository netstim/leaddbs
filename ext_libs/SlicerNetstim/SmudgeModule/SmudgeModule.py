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
import uuid

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


    blurPixmap = qt.QPixmap(self.resourcePath(os.path.join('Icons','blurIcon.png')))
    blurIcon = qt.QIcon(blurPixmap)
    self.blurButton = qt.QPushButton()
    self.blurButton.setIcon(blurIcon)
    self.blurButton.setIconSize(blurPixmap.rect().size())
    self.blurButton.setCheckable(True)
    self.blurButton.setEnabled(False)
    self.blurButton.setAutoExclusive(True)
    optionsFrame.layout().addWidget(self.blurButton)


    toolsFormLayout.addRow(optionsFrame)


    #
    # radius value
    #
    self.radiusSlider = ctk.ctkSliderWidget()
    self.radiusSlider.singleStep = 0.1
    self.radiusSlider.minimum = 5
    self.radiusSlider.maximum = float(self.parameterNode.GetParameter("maxRadius"))
    self.radiusSlider.decimals = 1
    self.radiusSlider.value = float(self.parameterNode.GetParameter("radius"))
    toolsFormLayout.addRow("Radius (mm):", self.radiusSlider)

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
    # force
    #
    self.forceSlider = ctk.ctkSliderWidget()
    self.forceSlider.singleStep = 1
    self.forceSlider.minimum = 0
    self.forceSlider.maximum = 100
    self.forceSlider.decimals = 0
    self.forceSlider.value = float(self.parameterNode.GetParameter("force"))
    toolsFormLayout.addRow("Force (%):", self.forceSlider)

    #
    # sigma
    #
    self.sigmaSlider = ctk.ctkSliderWidget()
    self.sigmaSlider.singleStep = 1
    self.sigmaSlider.minimum = 0
    self.sigmaSlider.maximum = 30
    self.sigmaSlider.decimals = 0
    self.sigmaSlider.value = float(self.parameterNode.GetParameter("sigma"))
    toolsFormLayout.addRow("Sigma (mm):", self.sigmaSlider)


    #
    # Advanced Area
    #
    advancedCollapsibleButton = ctk.ctkCollapsibleButton()
    advancedCollapsibleButton.text = "Advanced"
    advancedCollapsibleButton.collapsed=True
    self.layout.addWidget(advancedCollapsibleButton)

    # Layout within the dummy collapsible button
    advancedFormLayout = qt.QFormLayout(advancedCollapsibleButton)

    #
    # expand edge value
    #
    self.expandEdgeSlider = ctk.ctkSliderWidget()
    self.expandEdgeSlider.singleStep = 1
    self.expandEdgeSlider.minimum = 0
    self.expandEdgeSlider.maximum = 100
    self.expandEdgeSlider.decimals = 0
    self.expandEdgeSlider.value = float(self.parameterNode.GetParameter("expandEdge"))
    advancedFormLayout.addRow("Expand Edge (mm):", self.expandEdgeSlider)


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

    self.firstComponentCheckBox = qt.QCheckBox('Component 0')
    self.firstComponentCheckBox.toolTip = 'If checked, Component 0 will be included when flattening'

    undoredoFrame.layout().addWidget(self.undoButton)
    undoredoFrame.layout().addWidget(self.redoButton)
    undoredoFrame.layout().addWidget(self.flattenButton)
    undoredoFrame.layout().addWidget(self.firstComponentCheckBox)
    
    historyGridLayout.addWidget(undoredoFrame,0,0)

    # History Stack

    self.historyList = qt.QListWidget()
    self.historyList.setCurrentRow(0)
    self.historyList.setMaximumHeight(redoPixmap.rect().size().height() * 3.3 + self.firstComponentCheckBox.height)

    historyGridLayout.addWidget(self.historyList,0,1)



    #
    # Modles Area
    #

    modelsCollapsibleButton = ctk.ctkCollapsibleButton()
    modelsCollapsibleButton.text = "Data Control"
    self.layout.addWidget(modelsCollapsibleButton, 1)

    modelsGridLayout = qt.QGridLayout(modelsCollapsibleButton)    

    # Type selector
    self.sceneRadioButton = qt.QRadioButton('Scene')
    self.atlasesRadioButton = qt.QRadioButton('Atlases')
    self.drawingsRadioButton = qt.QRadioButton('Drawings')
    self.fiducialsRadioButton = qt.QRadioButton('Fiducials')
    self.atlasesRadioButton.setChecked(True)

    # delete selected
    self.deleteTreeElement = qt.QPushButton('Delete')

    # Tree widget
    self.dataTreeWidget = slicer.qMRMLSubjectHierarchyTreeView(slicer.util.mainWindow())
    self.dataTreeWidget.setMRMLScene(slicer.mrmlScene)
    self.dataTreeWidget.setColumnHidden(self.dataTreeWidget.model().idColumn, True)
    self.dataTreeWidget.setColumnHidden(self.dataTreeWidget.model().transformColumn, True)
    self.dataTreeWidget.setColumnHidden(self.dataTreeWidget.model().descriptionColumn, True)
    self.dataTreeWidget.contextMenuEnabled = False

    modelsGridLayout.addWidget(self.sceneRadioButton,0,0)
    modelsGridLayout.addWidget(self.atlasesRadioButton,0,1)
    modelsGridLayout.addWidget(self.drawingsRadioButton,0,2)
    modelsGridLayout.addWidget(self.fiducialsRadioButton,0,3)
    modelsGridLayout.addWidget(self.deleteTreeElement,0,4)
    modelsGridLayout.addWidget(self.dataTreeWidget,1,0,1,5)

    self.layout.addStretch(0)


    # connections
    self.warpSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.exit) # deselect effect
    self.warpSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateMRMLFromGUI)
    self.smudgeButton.toggled.connect(self.toogleTools)
    self.snapButton.toggled.connect(self.toogleTools)
    self.blurButton.toggled.connect(self.toogleTools)
    self.smudgeButton.connect('clicked(bool)', self.onSmudgeButton)
    self.snapButton.connect('clicked(bool)', self.onSnapButton)
    self.blurButton.connect('clicked(bool)', self.onBlurButton)
    
    self.radiusSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.hardnessSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.forceSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.sigmaSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.expandEdgeSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.flattenButton.connect("clicked(bool)", self.onFlattenButton)
    self.undoButton.connect("clicked(bool)", self.onUndoButton)
    self.redoButton.connect("clicked(bool)", self.onRedoButton)
    self.historyList.itemSelectionChanged.connect(self.historyItemChanged)

    
    self.sceneRadioButton.connect("toggled(bool)", self.dataTreeTypeToggle)
    self.drawingsRadioButton.connect("toggled(bool)", self.dataTreeTypeToggle)
    self.atlasesRadioButton.connect("toggled(bool)", self.dataTreeTypeToggle)
    self.fiducialsRadioButton.connect("toggled(bool)", self.dataTreeTypeToggle)
    self.deleteTreeElement.connect("clicked(bool)", lambda i: self.dataTreeWidget.deleteSelectedItems())
    self.dataTreeWidget.doubleClicked.connect(self.onDataTreeDoubleClicked)

    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)    
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.NodeAddedEvent, self.onSceneNodeAdded)    
    interactionNode = slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLInteractionNode')
    self.addObserver(interactionNode, interactionNode.InteractionModeChangedEvent, self.onInteractionModeChanged) 


    # Refresh
    qt.QApplication.processEvents()
    if self.updateMRMLFromArgs(): # was called from command line
      self.showSingleModule()
      tb = Toolbar.reducedToolbar()
      slicer.util.mainWindow().addToolBar(tb)

    self.updateGuiFromMRML()  
    self.onSceneNodeAdded()
    self.toogleTools()
    self.dataTreeTypeToggle(1)

  

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

    if singleModule:
      slicer.util.setPythonConsoleVisible(False)

    self.inputsCollapsibleButton.setVisible(not singleModule)
    if self.developerMode:
      self.reloadCollapsibleButton.setVisible(not singleModule)

    self.sceneRadioButton.setEnabled(False)

    slicer.util.mainWindow().setWindowTitle("Name goes here")


  def updateMRMLFromArgs(self): 
    args = sys.argv
    print(args)
    if (len(sys.argv) > 2) and os.path.isfile(os.path.join(sys.argv[1],'lead.m')):
      subjectPaths = self.parameterNode.GetParameter("separator").join(sys.argv[2:])
      subjectPath = subjectPaths.split(self.parameterNode.GetParameter("separator"))[0]
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
    # radius - blur - force
    radius = float(self.parameterNode.GetParameter("radius"))
    self.radiusSlider.setValue( radius )
    if radius < self.radiusSlider.minimum or radius > self.radiusSlider.maximum:
      self.updateMRMLFromGUI()
    self.hardnessSlider.setValue(float(self.parameterNode.GetParameter("hardness")))
    self.forceSlider.setValue(float(self.parameterNode.GetParameter("force")))
    self.sigmaSlider.setValue(float(self.parameterNode.GetParameter("sigma")))
    # get warp node and set selector and buttons
    warpID = self.parameterNode.GetParameter("warpID")
    warpNode = slicer.util.getNode(warpID) if warpID != "" else None
    warpNumberOfComponents = TransformsUtil.TransformsUtilLogic().getNumberOfLayers(warpNode)
    self.warpSelector.setCurrentNode(warpNode)
    self.smudgeButton.enabled = bool(warpNode)
    self.snapButton.enabled = bool(warpNode)
    self.blurButton.enabled = bool(warpNumberOfComponents == 1)
    self.historyList.enabled = bool(warpNode) 
    # undo redo button
    self.undoButton.setEnabled(warpNumberOfComponents > 1) 
    self.redoButton.setEnabled(self.parameterNode.GetParameter("redoTransformID") != "") 
    self.flattenButton.setEnabled(warpNumberOfComponents > 1)
    # update warp history. add number of warp component plus redo transform (if available)
    self.historyList.clear()
    self.historyList.addItems(['Component ' + str(i) for i in range(warpNumberOfComponents + int(self.redoButton.enabled))])
    self.historyList.setCurrentRow(warpNumberOfComponents - 1)
    # resolution change
    if float(self.parameterNode.GetParameter("resolution")) != TransformsUtil.TransformsUtilLogic().getGridDefinition(warpNode)[2][0]:
      self.exit()
    # get subject hierarchy node
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    # if subject is changed
    if bool(int(self.parameterNode.GetParameter("subjectChanged"))):
      self.exit()
      shNode.RemoveItem(int(self.parameterNode.GetParameter("drawingsRootItem")))
      self.parameterNode.SetParameter("drawingsRootItem","0") # reset drawings
      self.parameterNode.SetParameter("subjectChanged","0")
    # drawings hierarchy
    if not bool(int(self.parameterNode.GetParameter("drawingsRootItem"))):
      folderID = shNode.CreateFolderItem(shNode.GetSceneItemID(), 'AnchorDrawings')
      self.parameterNode.SetParameter("drawingsRootItem", str(folderID))
    self.dataTreeTypeToggle(1) # update

  def updateMRMLFromGUI(self):
    self.parameterNode.SetParameter("radius", str(self.radiusSlider.value) )
    self.parameterNode.SetParameter("hardness", str(self.hardnessSlider.value) )
    self.parameterNode.SetParameter("force", str(self.forceSlider.value) )
    self.parameterNode.SetParameter("sigma", str(self.sigmaSlider.value) )
    self.parameterNode.SetParameter("expandEdge", str(self.expandEdgeSlider.value) )
    self.parameterNode.SetParameter("warpID", self.warpSelector.currentNode().GetID() if self.warpSelector.currentNode() else "")

  def historyItemChanged(self):
    warpID = self.parameterNode.GetParameter("warpID")
    if warpID != "":
      warpNode = slicer.util.getNode(warpID)
      self.historyList.setCurrentRow(TransformsUtil.TransformsUtilLogic().getNumberOfLayers(warpNode) - 1) # keep same value

  def setItemChildrenFlags(self, item):
    item.setFlags(qt.Qt.ItemIsSelectable + qt.Qt.ItemIsEnabled)
    for row in range(item.rowCount()):
      self.setItemChildrenFlags(item.child(row))      

  def onSceneNodeAdded(self,caller=None,event=None):
    sceneItem = self.dataTreeWidget.model().item(0,0)
    self.setItemChildrenFlags(sceneItem)

  def onDataTreeDoubleClicked(self, i):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    node = shNode.GetItemDataNode(self.dataTreeWidget.currentItem())
    # get center position of model/drawing
    if isinstance(node, slicer.vtkMRMLModelNode):
      pd = node.GetPolyData()
      center = vtk.vtkCenterOfMass()
      center.SetInputData(pd)
      center.Update()
      centerList = center.GetCenter()
    elif isinstance(node, slicer.vtkMRMLMarkupsCurveNode):
      centerList = [0] * 3
      node.GetNthControlPointPosition(round(node.GetNumberOfControlPoints()/2),centerList)
    else:
      return
    SmudgeModuleLogic().centerPosition(centerList)
    

  def dataTreeTypeToggle(self, b):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    if self.sceneRadioButton.isChecked():
      self.dataTreeWidget.nodeTypes = ('vtkMRMLModelNode','vtkMRMLFolderDisplayNode')
      self.dataTreeWidget.attributeNameFilter = ('')
    elif self.atlasesRadioButton.isChecked():
      self.dataTreeWidget.nodeTypes = ('vtkMRMLModelNode','vtkMRMLFolderDisplayNode')
      self.dataTreeWidget.attributeNameFilter = ('atlas')
    elif self.drawingsRadioButton.isChecked():
      self.dataTreeWidget.nodeTypes = ('vtkMRMLMarkupsCurveNode','vtkMRMLFolderDisplayNode')
      self.dataTreeWidget.attributeNameFilter = ('drawing')
    elif self.fiducialsRadioButton.isChecked():
      self.dataTreeWidget.attributeNameFilter = ('fiducial')
    # reset settings
    self.dataTreeWidget.expandToDepth(0)


  def onUndoButton(self):
    SmudgeModuleLogic().removeRedoTransform()
    redoTransformID = TransformsUtil.TransformsUtilLogic().removeLastLayer(slicer.util.getNode(self.parameterNode.GetParameter("warpID")))
    self.parameterNode.SetParameter("redoTransformID", redoTransformID)
    # apply to anchor points
    redoTransformNode = slicer.util.getNode(redoTransformID)
    SmudgeModuleLogic().applyChangesToAnchorPoints(redoTransformNode.GetTransformFromParent())

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
    SmudgeModuleLogic().applyChangesToAnchorPoints(redoTransformNode.GetTransformToParent())
    # delete redo transform
    slicer.mrmlScene.RemoveNode(redoTransformNode)
    self.parameterNode.SetParameter("redoTransformID","")
    # restore state
    if smudgeEnabled:
      SmudgeModuleLogic().effectOn('Smudge')

  def onFlattenButton(self):
    SmudgeModuleLogic().removeRedoTransform()
    warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))
    includeFirstLayer = self.firstComponentCheckBox.checked
    TransformsUtil.TransformsUtilLogic().flattenTransform(warpNode, includeFirstLayer)
    self.updateGuiFromMRML() # update history


  def warpNodeModified(self, caller=None, event=None):
    # update gui
    self.updateGuiFromMRML()
    self.parameterNode.SetParameter("warpModified", str(int(self.parameterNode.GetParameter("warpModified"))+1) )


  def toogleTools(self):
    SmudgeModuleLogic().effectOff() # deactivate previous effect (if any)
    self.radiusSlider.setEnabled(False)
    self.hardnessSlider.setEnabled(False)
    self.forceSlider.setEnabled(False)
    self.sigmaSlider.setEnabled(False)
    self.expandEdgeSlider.setEnabled(True)
    interactionNode = slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLInteractionNode')
    self.removeObserver(interactionNode, interactionNode.InteractionModeChangedEvent, self.onInteractionModeChanged) # so that changing mode doesnt affect
    if not self.noneButton.checked:
      interactionNode.SetCurrentInteractionMode(2)
      self.addObserver(interactionNode, interactionNode.InteractionModeChangedEvent, self.onInteractionModeChanged) 
      warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))
      eventTag = warpNode.AddObserver(slicer.vtkMRMLGridTransformNode.TransformModifiedEvent, self.warpNodeModified)
    else:
      SmudgeModuleLogic().effectOn('None')

  def onSmudgeButton(self, buttonDown):
    if buttonDown:      
      self.radiusSlider.setEnabled(True)
      self.hardnessSlider.setEnabled(True)
      self.forceSlider.setEnabled(True)
      self.sigmaSlider.setEnabled(True)
      self.expandEdgeSlider.setEnabled(False)
      SmudgeModuleLogic().effectOn('Smudge')

  def onSnapButton(self, buttonDown):
    if buttonDown:
      self.radiusSlider.setEnabled(True)
      self.sigmaSlider.setEnabled(True)
      SmudgeModuleLogic().effectOn('Snap')

  def onBlurButton(self, buttonDown):
    if buttonDown:
      self.radiusSlider.setEnabled(True)
      self.hardnessSlider.setEnabled(True)
      self.forceSlider.setEnabled(True)
      self.sigmaSlider.setEnabled(True)
      SmudgeModuleLogic().effectOn('Blur')


  def exit(self):
    self.noneButton.setChecked(True)
    self.toogleTools()
    SmudgeModuleLogic().removeRedoTransform()
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
    node.SetParameter("radius", "25")
    node.SetParameter("maxRadius", "50")
    node.SetParameter("hardness", "40")
    node.SetParameter("force", "100")
    node.SetParameter("sigma", "2")
    node.SetParameter("expandEdge", "0")
    node.SetParameter("warpModified","0")
    node.SetParameter("drawingsRootItem","0")
    node.SetParameter("lastOperation","")
    # lead dbs specific
    node.SetParameter("affineTransformID", "")
    node.SetParameter("templateID", "")
    node.SetParameter("modality", "t1")
    node.SetParameter("subjectPath", "")
    node.SetParameter("subjectN", "0")
    node.SetParameter("separator",uuid.uuid4().hex)
    node.SetParameter("MNIPath", ".")
    node.SetParameter("MNIAtlasPath", ".")
    node.SetParameter("antsApplyTransformsPath", "")
    node.SetParameter("subjectChanged","0")
    node.SetParameter("resolution","1")
    return node

  def removeRedoTransform(self):
    parameterNode = self.getParameterNode()
    redoTransformID = parameterNode.GetParameter("redoTransformID")
    if redoTransformID != "":
      slicer.mrmlScene.RemoveNode(slicer.util.getNode(redoTransformID))
      parameterNode.SetParameter("redoTransformID","")

  def removeLastDrawing(self):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    drawingsRootItem = int(self.getParameterNode().GetParameter("drawingsRootItem"))
    ids = vtk.vtkIdList()
    shNode.GetItemChildren(drawingsRootItem,ids,False)
    shNode.RemoveItem(ids.GetId(ids.GetNumberOfIds()-1))

  def effectOn(self, effectName):
    
    if effectName in ['Smudge']:
      size, origin, spacing = self.getExpandedGrid()
      auxTransformNode = TransformsUtil.TransformsUtilLogic().emptyGridTransfrom(size, origin, spacing)

    for color in ['Red','Green','Yellow']:
      sliceWidget = slicer.app.layoutManager().sliceWidget(color)
      if effectName == 'None':
        WarpEffect.NoneEffect(sliceWidget)
      elif effectName == 'Smudge':
        WarpEffect.SmudgeEffectTool(sliceWidget, auxTransformNode)
      elif effectName == 'Blur':
        WarpEffect.BlurEffectTool(sliceWidget)
      elif effectName == 'Snap':
        WarpEffect.SnapEffectTool(sliceWidget)
  

  def effectOff(self):
    WarpEffect.WarpEffectTool.empty()


  def applyChangesToAnchorPoints(self, transform):
    # get drawings
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    drawingsRootItem = int(self.getParameterNode().GetParameter("drawingsRootItem"))
    childrenIDs = vtk.vtkIdList()
    shNode.GetItemChildren(drawingsRootItem, childrenIDs, True)
    # apply transform to all source points
    for i in range(childrenIDs.GetNumberOfIds()):
      dataNode = shNode.GetItemDataNode(childrenIDs.GetId(i))
      if dataNode and dataNode.GetName() == 'source':
        dataNode.ApplyTransform(transform)


  def getExpandedGrid(self):
    # create aux transform with same grid as warp
    parameterNode = self.getParameterNode()
    warpNode = slicer.util.getNode(parameterNode.GetParameter("warpID"))
    size,origin,spacing = TransformsUtil.TransformsUtilLogic().getGridDefinition(warpNode)
    # expand aux transform to deal with borders
    expandEdge = float(parameterNode.GetParameter("expandEdge"))
    origin = [o - expandEdge for o in origin]
    size = [int(round(s+expandEdge*2/spacing[0])) for s in size]

    return size, origin, spacing

  def centerPosition(self, centerList):
    # create markups node, add center as fiducial and jump and center slices
    markupsNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    markupsNode.GetDisplayNode().SetVisibility(False)
    markupsNode.AddFiducialFromArray(np.array(centerList),'')
    markupsLogic = slicer.modules.markups.logic()
    markupsLogic.JumpSlicesToNthPointInMarkup(markupsNode.GetID(),0,True)
    slicer.mrmlScene.RemoveNode(markupsNode)




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




    





