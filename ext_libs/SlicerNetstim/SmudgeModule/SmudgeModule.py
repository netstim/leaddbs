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


    self.noneButton = warpEffects[0].effectButton

    #
    # History Area
    #
    editCollapsibleButton = ctk.ctkCollapsibleButton()
    editCollapsibleButton.text = "Edit"
    self.layout.addWidget(editCollapsibleButton)

    editFormLayout = qt.QFormLayout(editCollapsibleButton)  

    #
    # Undo Redo
    #   

    undoredoFrame = qt.QFrame()
    undoredoFrame.setLayout(qt.QHBoxLayout())

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
      undoredoFrame.layout().addWidget(button)
      b['widget'] = button

    self.undoAllButton = undoAllButton['widget']
    self.undoButton = undoButton['widget']
    self.redoButton = redoButton['widget']
    editFormLayout.addRow(undoredoFrame)


    #
    # Modles Area
    #

    modelsCollapsibleButton = ctk.ctkCollapsibleButton()
    modelsCollapsibleButton.text = "Data Control"
    self.layout.addWidget(modelsCollapsibleButton, 1)

    modelsFormLayout = qt.QGridLayout(modelsCollapsibleButton)    
    tv = treeView.WarpDriveTreeView()
    modelsFormLayout.addWidget(tv,0,0)

    self.layout.addStretch(0)

    # connections
    self.warpSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.exit) # deselect effect
    self.warpSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onWarpSelectionChanged)

    self.undoAllButton.connect("clicked(bool)", self.onUndoAllButton)
    self.undoButton.connect("clicked(bool)", self.onUndoButton)
    self.redoButton.connect("clicked(bool)", self.onRedoButton)
    
    for effect in warpEffects:
      effect.addEditButtonListeners(self)

    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)    

    # Refresh
    qt.QApplication.processEvents()
    if self.updateMRMLFromArgs(): # was called from command line
      self.showSingleModule()
      tb = Toolbar.reducedToolbar()
      slicer.util.mainWindow().addToolBar(tb)
      tv.radioButtons[0].setVisible(False)
    else:
      tv.radioButtons[1].setVisible(False)

    self.updateGuiFromMRML()  
    self.noneButton.setEnabled(True)

  

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

    slicer.util.setPythonConsoleVisible(not singleModule)

    self.inputsCollapsibleButton.setVisible(not singleModule)
    if self.developerMode:
      self.reloadCollapsibleButton.setVisible(not singleModule)

    for color in ["Red","Green","Yellow"]:
      sliceController = slicer.app.layoutManager().sliceWidget(color).sliceController()
      sliceController.pinButton().hide()
      sliceController.viewLabel().hide()

    # data probe
    for i in range(slicer.mrmlScene.GetNumberOfNodesByClass("vtkMRMLScriptedModuleNode")):
      n  = slicer.mrmlScene.GetNthNodeByClass( i, "vtkMRMLScriptedModuleNode" )
      if n.GetModuleName() == "DataProbe":
        n.SetParameter('sliceViewAnnotationsEnabled','0')

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
      shNode.RemoveItem(int(self.parameterNode.GetParameter("drawingsRootItem")))
      self.parameterNode.SetParameter("drawingsRootItem","0") # reset drawings
      self.parameterNode.SetParameter("subjectChanged","0")
    # drawings hierarchy
    if not bool(int(self.parameterNode.GetParameter("drawingsRootItem"))):
      folderID = shNode.CreateFolderItem(shNode.GetSceneItemID(), 'Fixed Points')
      shNode.SetItemAttribute(folderID, 'drawing', '1')
      self.parameterNode.SetParameter("drawingsRootItem", str(folderID))


  def onWarpSelectionChanged(self):
    self.parameterNode.SetParameter("warpID", self.warpSelector.currentNode().GetID() if self.warpSelector.currentNode() else "")


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
    self.noneButton.setChecked(True)
    SmudgeModuleLogic().removeRedoNodes()


  def cleanup(self):
    self.exit()

  def enter(self):
    self.noneButton.setChecked(True)
      
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
    node.SetParameter("lastDrawingID", "-1")
    node.SetParameter("warpModified","0")
    node.SetParameter("drawingsRootItem","0")
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
    node.SetParameter("DrawSpread", "25")
    # Smooth
    node.SetParameter("SmoothRadius", "25")
    node.SetParameter("SmoothHardness", "40")
    node.SetParameter("SmoothSigma", "10")
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
    parameterNode = self.getParameterNode()
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    drawingsRootItem = int(parameterNode.GetParameter("drawingsRootItem"))
    ids = vtk.vtkIdList()
    shNode.GetItemChildren(drawingsRootItem,ids,False)
    for j in range(ids.GetNumberOfIds()-1,-1,-1):
      if isinstance(shNode.GetItemDataNode(ids.GetId(j)), slicer.vtkMRMLMarkupsCurveNode):
        lastDrawingID = ids.GetId(j)
        shNode.SetItemParent(lastDrawingID, shNode.GetSceneItemID())
        shNode.SetItemDisplayVisibility(lastDrawingID, 0)
        shNode.SetItemAttribute(lastDrawingID, 'drawing', '0')
        parameterNode.SetParameter("lastDrawingID", str(lastDrawingID))
        return

  def enableLastDrawing(self):
    parameterNode = self.getParameterNode()
    lastDrawingID = int(parameterNode.GetParameter("lastDrawingID"))
    if lastDrawingID != -1:
      shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
      drawingsRootItem = int(parameterNode.GetParameter("drawingsRootItem"))
      shNode.SetItemParent(lastDrawingID, drawingsRootItem)
      shNode.SetItemDisplayVisibility(lastDrawingID, 1)
      shNode.SetItemAttribute(lastDrawingID, 'drawing', '1')
      parameterNode.SetParameter("lastDrawingID", "-1")
      # add and delete aux fiducial to refresh view #TODO: why doesnt it work
      auxFid = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
      shNode.SetItemAttribute(shNode.GetItemByDataNode(auxFid), 'drawing', '1')
      slicer.mrmlScene.RemoveNode(auxFid)



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




    





