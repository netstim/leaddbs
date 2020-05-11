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
from PythonQt import BoolResult

# netstim helpers
from Helpers import WarpEffect, FunctionsUtil, Toolbar, WarpEffectParameters, treeView

# netstim modules
import TransformsUtil
import ImportAtlas


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
    # I / O
    #
    self.inputsCollapsibleButton = ctk.ctkCollapsibleButton()
    self.inputsCollapsibleButton.text = "I / O"
    self.layout.addWidget(self.inputsCollapsibleButton)

    # Layout within the dummy collapsible button
    inputsFormLayout = qt.QFormLayout(self.inputsCollapsibleButton)


    #
    # clean up checkbox
    #
    self.cleanUpOnNodeChange = qt.QCheckBox('')
    self.cleanUpOnNodeChange.setChecked(True)
    inputsFormLayout.addRow("Clean Up On Master Change: ", self.cleanUpOnNodeChange)


    #
    # master node selector
    #
    self.masterNodeSelector = slicer.qMRMLNodeComboBox()
    self.masterNodeSelector.nodeTypes = ["vtkMRMLGridTransformNode", "vtkMRMLScalarVolumeNode"]
    self.masterNodeSelector.selectNodeUponCreation = False
    self.masterNodeSelector.addEnabled = False
    self.masterNodeSelector.removeEnabled = False
    self.masterNodeSelector.noneEnabled = True
    self.masterNodeSelector.showHidden = False
    self.masterNodeSelector.showChildNodeTypes = False
    self.masterNodeSelector.setMRMLScene( slicer.mrmlScene )
    self.masterNodeSelector.setToolTip( "Pick the warp to refine or volume to apply warp to." )
    inputsFormLayout.addRow("Master Warp or Volume: ", self.masterNodeSelector)

    #
    # output
    #
    self.hardenOutputPushButton = qt.QPushButton('Harden Current Output Warp on Master')
    self.hardenOutputPushButton.setEnabled(False)
    inputsFormLayout.addRow("Output: ", self.hardenOutputPushButton)

    #
    # Tools Area
    #
    toolsCollapsibleButton = ctk.ctkCollapsibleButton()
    toolsCollapsibleButton.text = "Tools"
    self.layout.addWidget(toolsCollapsibleButton)

    # Layout within the dummy collapsible button
    toolsFormLayout = qt.QFormLayout(toolsCollapsibleButton)

    toolsFrame = qt.QFrame()
    toolsFrame.setLayout(qt.QHBoxLayout())
    toolsFormLayout.addRow(toolsFrame)

    warpEffects = [WarpEffectParameters.NoneEffectParameters(), 
                  WarpEffectParameters.SmudgeEffectParameters(), 
                  WarpEffectParameters.DrawEffectParameters(),
                  WarpEffectParameters.SmoothEffectParameters()]

    for warpEffectParametersWidget in warpEffects:
      toolsFrame.layout().addWidget(warpEffectParametersWidget.effectButton)
      toolsFormLayout.addRow(warpEffectParametersWidget.parametersFrame)


    #
    # Undo Redo
    #   

    undoRedoGroupBox = qt.QGroupBox('Edit')
    undoRedoGroupBox.setLayout(qt.QHBoxLayout())
    toolsFormLayout.addRow(undoRedoGroupBox)

    undoAllButton =   {'text':'Undo All',  'icon':'UndoAll',   'toolTip':'Undo all user modifications. Fixed points won\'t be deleted.'}
    undoButton =      {'text':'Undo',      'icon':'Undo',      'toolTip':'Undo last operation. In case it was a drawing, corresponding fixed points will be deleted.'}
    redoButton =      {'text':'Redo',      'icon':'Redo',      'toolTip':'Redo'}

    # dont use QToolButton in order to use QPushButton's pressed and release signals
    buttonStyleSheet = "QPushButton { \
                          background-image: url(%s); \
                          font-size: 10px; \
                          text-align: bottom; \
                          border-radius: 3px; \
                          border-style: solid; \
                          border-color: rgb(182, 182, 182); \
                          border-width: 1px; } \
                        QPushButton:disabled { \
                          background-image: url(%s); }\
                        QPushButton:pressed { \
                          background-color: rgb(232, 232, 232); }"

    for b in [undoAllButton, undoButton, redoButton]:
      buttonIconPath = self.resourcePath(os.path.join('Icons', b['icon'] + '%s.png'))
      buttonPixmap = qt.QPixmap(buttonIconPath %'')
      button = qt.QPushButton(b['text'])
      button.setStyleSheet(buttonStyleSheet % (buttonIconPath %'', buttonIconPath %'_disabled'))
      button.setFixedSize(buttonPixmap.rect().size())
      button.setEnabled(False)
      button.setToolTip(b['toolTip'])
      undoRedoGroupBox.layout().addWidget(button)
      b['widget'] = button

    self.undoAllButton = undoAllButton['widget']
    self.undoButton = undoButton['widget']
    self.redoButton = redoButton['widget']


    #
    # Modles Area
    #

    modelsCollapsibleButton = ctk.ctkCollapsibleButton()
    modelsCollapsibleButton.text = "Data Control"
    self.layout.addWidget(modelsCollapsibleButton, 1)

    modelsFormLayout = qt.QGridLayout(modelsCollapsibleButton)    
    modelsFormLayout.addWidget(treeView.WarpDriveTreeView(),0,0)

    self.layout.addStretch(0)

    # connections
    self.hardenOutputPushButton.connect("clicked(bool)", self.onHardenOutputPushButton) 

    self.masterNodeSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.exit) # deselect effect
    self.masterNodeSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onMasterNodeSelectionChanged)

    self.undoAllButton.connect("clicked(bool)", self.onUndoAllButton)
    self.undoButton.connect("clicked(bool)", self.onUndoButton)
    self.redoButton.connect("clicked(bool)", self.onRedoButton)

    for button in [self.undoAllButton, self.undoButton, self.redoButton]:
      button.connect("pressed()", self.onEditButtonPressed)
      button.connect("released()", self.onEditButtonReleased)
    
    for effect in warpEffects:
      effect.addEditButtonListeners(self)

    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)    

    # Refresh
    qt.QApplication.processEvents()

    # check dependencies
    if self.checkExtensionInstall(extensionName = 'SlicerRT'):
      return

    # Lead-DBS call
    if self.updateMRMLFromArgs(): # was called from command line
      self.showSingleModule()
      slicer.util.mainWindow().addToolBar(Toolbar.reducedToolbar())

    self.updateGuiFromMRML()  

  

  #
  # Methods
  #
  
  def checkExtensionInstall(self, extensionName):
    em = slicer.app.extensionsManagerModel()
    if not em.isExtensionInstalled(extensionName):
      extensionMetaData = em.retrieveExtensionMetadataByName(extensionName)
      url = os.path.join(em.serverUrl().toString(), 'download', 'item', extensionMetaData['item_id'])
      extensionPackageFilename = os.path.join(slicer.app.temporaryPath, extensionMetaData['md5'])
      slicer.util.downloadFile(url, extensionPackageFilename)
      em.installExtension(extensionPackageFilename)
      qt.QMessageBox.information(qt.QWidget(), '', 'Slicer will install %s and quit.\nPlease restart.' % extensionName)
      slicer.util.exit()
      return True

  def showSingleModule(self):
    
    singleModule = True

    # toolbars
    slicer.util.setToolbarsVisible(not singleModule, [])

    # customize view
    viewToolBar = slicer.util.mainWindow().findChild('QToolBar', 'ViewToolBar')
    viewToolBar.setVisible(1)
    layoutMenu = viewToolBar.widgetForAction(viewToolBar.actions()[0]).menu()
    for action in layoutMenu.actions():
      if action.text not in ['Four-Up', 'Tabbed slice']:
        layoutMenu.removeAction(action)

    # customize mouse mode
    mouseModeToolBar = slicer.util.mainWindow().findChild('QToolBar', 'MouseModeToolBar')
    mouseModeToolBar.setVisible(1)
    mouseModeToolBar.removeAction(mouseModeToolBar.actions()[2]) # remove place markups

    # viewers
    viewersToolBar = slicer.util.mainWindow().findChild('QToolBar', 'ViewersToolBar')
    viewersToolBar.setVisible(1)

    # slicer window
    slicer.util.setMenuBarsVisible(not singleModule)
    slicer.util.setApplicationLogoVisible(not singleModule)
    slicer.util.setModuleHelpSectionVisible(not singleModule)
    slicer.util.setModulePanelTitleVisible(not singleModule)
    slicer.util.setDataProbeVisible(not singleModule)
    slicer.util.setPythonConsoleVisible(not singleModule)

    # inputs area
    self.inputsCollapsibleButton.setVisible(not singleModule)

    # reload area
    if self.developerMode:
      self.reloadCollapsibleButton.setVisible(not singleModule)

    # slice controllers
    for color in ["Red","Green","Yellow"]:
      sliceController = slicer.app.layoutManager().sliceWidget(color).sliceController()
      sliceController.pinButton().hide()
      sliceController.viewLabel().hide()

    # data probe
    for i in range(slicer.mrmlScene.GetNumberOfNodesByClass("vtkMRMLScriptedModuleNode")):
      n  = slicer.mrmlScene.GetNthNodeByClass( i, "vtkMRMLScriptedModuleNode" )
      if n.GetModuleName() == "DataProbe":
        n.SetParameter('sliceViewAnnotationsEnabled','0')

    # set name
    slicer.util.mainWindow().setWindowTitle("Warp Drive")


  def updateMRMLFromArgs(self): 
    args = sys.argv
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


  def updateGuiFromMRML(self,caller=None,event=None):
    # get warp node and set selector and buttons
    warpID = self.parameterNode.GetParameter("warpID")
    warpNode = slicer.util.getNode(warpID) if warpID != "" else None
    warpNumberOfComponents = TransformsUtil.TransformsUtilLogic().getNumberOfLayers(warpNode)
    # undo redo button
    self.undoButton.setEnabled(warpNumberOfComponents > 1 and self.parameterNode.GetParameter("redoTransformID") == "" and self.parameterNode.GetParameter("lastOperation") != "UndoAll") 
    self.redoButton.setEnabled(self.parameterNode.GetParameter("redoTransformID") != "") 
    self.undoAllButton.setEnabled(warpNumberOfComponents > 1)
    # resolution change
    if float(self.parameterNode.GetParameter("resolution")) != TransformsUtil.TransformsUtilLogic().getGridDefinition(warpNode)[2][0]:
      self.exit()
    # get subject hierarchy node
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    # if subject is changed
    if bool(int(self.parameterNode.GetParameter("subjectChanged"))):
      self.exit()
      self.parameterNode.SetParameter("subjectChanged","0")


  def onMasterNodeSelectionChanged(self):
    # clean up
    if self.cleanUpOnNodeChange.checked:
      SmudgeModuleLogic().cleanUp()
    # get node
    currentNode = self.masterNodeSelector.currentNode()
    if currentNode:
      # get size origin spacing
      if isinstance(currentNode, slicer.vtkMRMLGridTransformNode):
        size,origin,spacing = TransformsUtil.TransformsUtilLogic().getGridDefinition(currentNode)
      elif isinstance(currentNode, slicer.vtkMRMLScalarVolumeNode):
        size = currentNode.GetImageData().GetDimensions()
        origin = currentNode.GetOrigin()
        spacing = currentNode.GetSpacing()
      # create warp
      warpNode = TransformsUtil.TransformsUtilLogic().emptyGridTransform(size,origin,spacing)
      warpNode.CreateDefaultDisplayNodes()
      warpNode.GetDisplayNode().SetVisibility2D(True)
      warpNode.SetDescription('Current')
      warpNode.SetName(slicer.mrmlScene.GenerateUniqueName('Initial'))
      shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
      shNode.SetItemAttribute(shNode.GetItemByDataNode(warpNode), 'savedWarp', '1')
      # set
      currentNode.SetAndObserveTransformNodeID(warpNode.GetID())
      self.parameterNode.SetParameter("resolution", str(spacing[0]))
      self.parameterNode.SetParameter("warpID", warpNode.GetID())
    else:
      self.parameterNode.SetParameter("warpID", "")
    # enable / disable harden option
    self.hardenOutputPushButton.enabled = self.parameterNode.GetParameter("warpID") != ""

  def onHardenOutputPushButton(self, b):
    self.masterNodeSelector.currentNode().HardenTransform()
    self.masterNodeSelector.setCurrentNode(None)

  def onEditButtonPressed(self):
    qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
    qt.QApplication.processEvents()

  def onEditButtonReleased(self):
    qt.QApplication.setOverrideCursor(qt.QCursor(qt.Qt.ArrowCursor))

  def onUndoAllButton(self):
    self.parameterNode.SetParameter("lastOperation","UndoAll")
    warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))
    if TransformsUtil.TransformsUtilLogic().getNumberOfLayers(warpNode) > 2:
      TransformsUtil.TransformsUtilLogic().flattenTransform(warpNode, False)
    self.onUndoButton()

  def onUndoButton(self):
    # remove redo nodes
    SmudgeModuleLogic().removeRedoNodes()
    # apply and save redo transform
    redoTransformID = TransformsUtil.TransformsUtilLogic().removeLastLayer(slicer.util.getNode(self.parameterNode.GetParameter("warpID")))
    self.parameterNode.SetParameter("redoTransformID", redoTransformID)
    # disable last drawing if was a drawing operation
    if self.parameterNode.GetParameter("lastOperation") == 'Draw':
      SmudgeModuleLogic().disableLastDrawing()

  def onRedoButton(self):
    # get nodes
    warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))
    redoTransformNode = slicer.util.getNode(self.parameterNode.GetParameter("redoTransformID"))
    # apply
    warpNode.SetAndObserveTransformNodeID(redoTransformNode.GetID())
    warpNode.HardenTransform()
    # delete redo transform
    slicer.mrmlScene.RemoveNode(redoTransformNode)
    self.parameterNode.SetParameter("redoTransformID","")
    # re enable drawing
    if self.parameterNode.GetParameter("lastOperation") == 'Draw':
      SmudgeModuleLogic().enableLastDrawing()


  def exit(self):
    WarpEffectParameters.NoneEffectParameters.activateNoneEffect()
    SmudgeModuleLogic().removeRedoNodes()

  def cleanup(self):
    self.exit()

  def enter(self):
    WarpEffectParameters.NoneEffectParameters.activateNoneEffect()
      
  def onSceneStartClose(self, caller, event):
    self.parameterNode.SetParameter("warpID", "")



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
    node.SetParameter("lastDrawingID", "-1")
    node.SetParameter("warpModified","0")
    node.SetParameter("lastOperation","")
    node.SetParameter("currentEfect","None")
    # smudge 
    node.SetParameter("SmudgeRadius", "25")
    node.SetParameter("SmudgeHardness", "40")
    node.SetParameter("SmudgeForce", "100")
    node.SetParameter("SmudgePostSmoothing", "0")
    node.SetParameter("SmudgeSigma", "10")
    node.SetParameter("expandGrid", "0")
    node.SetParameter("maxRadius", "50")
    node.SetParameter("gridBoundsROIID", "")
    # draw
    node.SetParameter("DrawSpread", "15")
    node.SetParameter("DrawSampleDistance", "2")
    node.SetParameter("DrawStiffness", "0.1")
    # Smooth
    node.SetParameter("SmoothRadius", "25")
    node.SetParameter("SmoothHardness", "50")
    node.SetParameter("SmoothSigma", "5")
    node.SetParameter("SmoothUseRadius", "1")
    # lead dbs specific
    node.SetParameter("glanatCompositeID", "")
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

  def removeRedoNodes(self):
    parameterNode = self.getParameterNode()
    # redo transform
    redoTransformID = parameterNode.GetParameter("redoTransformID")
    if redoTransformID != "":
      slicer.mrmlScene.RemoveNode(slicer.util.getNode(redoTransformID))
      parameterNode.SetParameter("redoTransformID","")
    # last drawing
    lastDrawingID = int(parameterNode.GetParameter("lastDrawingID"))
    if lastDrawingID != -1:
      shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
      shNode.RemoveItem(lastDrawingID)
      parameterNode.SetParameter("lastDrawingID","-1")

  def disableLastDrawing(self):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    nMarkups = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLMarkupsFiducialNode')
    for i in range(nMarkups-1,-1,-1):
      markupNode = slicer.mrmlScene.GetNthNodeByClass(i, 'vtkMRMLMarkupsFiducialNode')
      if 'drawing' in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(markupNode)) and markupNode.GetNumberOfControlPoints()>1:
        lastDrawingID = shNode.GetItemByDataNode(markupNode)
        shNode.SetItemParent(lastDrawingID, shNode.GetSceneItemID())
        shNode.SetItemDisplayVisibility(lastDrawingID, 0)
        shNode.SetItemAttribute(lastDrawingID, 'drawing', '0')
        self.getParameterNode().SetParameter("lastDrawingID", str(lastDrawingID))
        return

  def enableLastDrawing(self):
    parameterNode = self.getParameterNode()
    lastDrawingID = int(parameterNode.GetParameter("lastDrawingID"))
    if lastDrawingID != -1:
      shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
      shNode.SetItemDisplayVisibility(lastDrawingID, 1)
      shNode.SetItemAttribute(lastDrawingID, 'drawing', '1')
      parameterNode.SetParameter("lastDrawingID", "-1")

  def cleanUp(self):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    # delete warps
    transformNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLTransformNode')
    transformNodes.UnRegister(slicer.mrmlScene)
    for i in range(transformNodes.GetNumberOfItems()):
      transformNode = transformNodes.GetItemAsObject(i)
      if 'savedWarp' in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(transformNode)):
        slicer.mrmlScene.RemoveNode(transformNode)
    self.getParameterNode().SetParameter("warpID","")

    # delete fiducials
    markupsNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLMarkupsFiducialNode')
    markupsNodes.UnRegister(slicer.mrmlScene)
    for i in range(markupsNodes.GetNumberOfItems()):
      markupNode = markupsNodes.GetItemAsObject(i)
      if 'drawing' in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(markupNode)):
        slicer.mrmlScene.RemoveNode(markupNode)



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




    





