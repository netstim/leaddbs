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


    # connections
    self.warpSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.exit) # deselect effect
    self.warpSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateMRMLFromGUI)
    self.smudgeButton.toggled.connect(self.toogleTools)
    self.snapButton.toggled.connect(self.toogleTools)
    self.smudgeButton.connect('clicked(bool)', self.onSmudgeButton)
    self.snapButton.connect('clicked(bool)', self.onSnapButton)
    
    self.radiusSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.hardnessSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.forceSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.flattenButton.connect("clicked(bool)", self.onFlattenButton)
    self.undoButton.connect("clicked(bool)", self.onUndoButton)
    self.redoButton.connect("clicked(bool)", self.onRedoButton)
    self.historyList.itemSelectionChanged.connect(self.historyItemChanged)

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

    slicer.util.mainWindow().setWindowTitle("Name goes here")


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
    # radius - blur - force
    radius = float(self.parameterNode.GetParameter("radius"))
    self.radiusSlider.setValue( radius )
    if radius < self.radiusSlider.minimum or radius > self.radiusSlider.maximum:
      self.updateMRMLFromGUI()
    self.hardnessSlider.setValue(float(self.parameterNode.GetParameter("hardness")))
    self.forceSlider.setValue(float(self.parameterNode.GetParameter("force")))
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
    # if subject is changed
    if bool(int(self.parameterNode.GetParameter("subjectChanged"))):
      self.exit()
      self.toogleTools()
      self.parameterNode.SetParameter("subjectChanged","0")


  def updateMRMLFromGUI(self):
    self.parameterNode.SetParameter("radius", str(self.radiusSlider.value) )
    self.parameterNode.SetParameter("hardness", str(self.hardnessSlider.value) )
    self.parameterNode.SetParameter("force", str(self.forceSlider.value) )
    self.parameterNode.SetParameter("warpID", self.warpSelector.currentNode().GetID() if self.warpSelector.currentNode() else "")

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
    self.hardnessSlider.setEnabled(False)
    self.forceSlider.setEnabled(False)
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
      self.hardnessSlider.setEnabled(True)
      self.forceSlider.setEnabled(True)
      SmudgeModuleLogic().effectOn('Smudge')

  def onSnapButton(self, buttonDown):
    if buttonDown:
      self.radiusSlider.setEnabled(True)
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
    node.SetParameter("hardness", "60")
    node.SetParameter("force", "100")
    node.SetParameter("modality", "t1")
    node.SetParameter("subjectPath", "")
    node.SetParameter("subjectN", "0")
    node.SetParameter("MNIPath", ".")
    node.SetParameter("MNIAtlasPath", ".")
    node.SetParameter("antsApplyTransformsPath", "")
    node.SetParameter("warpModified", "1")
    node.SetParameter("subjectChanged","0")
    return node

  def removeRedoTransform(self):
    parameterNode = self.getParameterNode()
    redoTransformID = parameterNode.GetParameter("redoTransformID")
    if redoTransformID != "":
      slicer.mrmlScene.RemoveNode(slicer.util.getNode(redoTransformID))
      parameterNode.SetParameter("redoTransformID","")

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
    imageNode = slicer.util.loadVolume(os.path.join(subjectPath,'anat_' + parameterNode.GetParameter("modality") + '.nii'))

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




    





