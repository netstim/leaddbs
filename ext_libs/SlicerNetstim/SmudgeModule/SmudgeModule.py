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
from Helpers import WarpEffect, FunctionsUtil, Toolbar, WarpEffectParameters

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
                  WarpEffectParameters.BlurEffectParameters()]

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
    # Undo Redo Flatten
    #   

    undoredoFrame = qt.QFrame()
    undoredoFrame.setLayout(qt.QHBoxLayout())

    buttonNames = ['Undo All', 'Undo', 'Redo', 'Overwrite']
    buttons = []

    for name in buttonNames:
      buttonIconPath = self.resourcePath(os.path.join('Icons', name + '.png'))
      buttonPixmap = qt.QPixmap(buttonIconPath)
      button = qt.QPushButton(name)
      button.setStyleSheet("QPushButton { background-image: url(" + buttonIconPath + "); font-size: 10px; text-align: bottom; border-radius: 3px; border-style: solid; border-color: rgb(182, 182, 182); border-width: 1px; } QPushButton:disabled { background-image: url(" + self.resourcePath(os.path.join('Icons', name + '_d.png')) + "); }")
      button.setFixedSize(buttonPixmap.rect().size())
      button.setEnabled(False)
      button.setToolTip('')
      undoredoFrame.layout().addWidget(button)
      buttons.append(button)

    self.undoAllButton = buttons[0]
    self.undoButton = buttons[1]
    self.redoButton = buttons[2]
    self.overwriteButton = buttons[3]
    editFormLayout.addRow(undoredoFrame)


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
    self.drawingsRadioButton = qt.QRadioButton('Fixed Points')
    self.atlasesRadioButton.setChecked(True)

    # add
    addPixmap = qt.QPixmap(self.resourcePath(os.path.join('Icons','Add.png')))
    addIcon = qt.QIcon(addPixmap)
    self.addTreeElement = qt.QPushButton()
    self.addTreeElement.setIcon(addIcon)
    self.addTreeElement.setIconSize(addPixmap.rect().size())
    # delete
    deletePixmap = qt.QPixmap(self.resourcePath(os.path.join('Icons','Delete.png')))
    deleteIcon = qt.QIcon(deletePixmap)
    self.deleteTreeElement = qt.QPushButton()
    self.deleteTreeElement.setIcon(deleteIcon)
    self.deleteTreeElement.setIconSize(deletePixmap.rect().size())
    # rename
    renamePixmap = qt.QPixmap(self.resourcePath(os.path.join('Icons','Rename.png')))
    renameIcon = qt.QIcon(renamePixmap)
    self.renameTreeElement = qt.QPushButton()
    self.renameTreeElement.setIcon(renameIcon)
    self.renameTreeElement.setIconSize(renamePixmap.rect().size())

    modelEditFrame = qt.QFrame()
    modelEditFrame.setLayout(qt.QHBoxLayout())
    modelEditFrame.layout().addWidget(self.addTreeElement)
    modelEditFrame.layout().addWidget(self.deleteTreeElement)
    modelEditFrame.layout().addWidget(self.renameTreeElement)

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
    modelsGridLayout.addWidget(modelEditFrame,0,3)
    modelsGridLayout.addWidget(self.dataTreeWidget,1,0,1,4)

    self.layout.addStretch(0)


    # connections
    self.warpSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.exit) # deselect effect
    self.warpSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onWarpSelectionChanged)

    self.overwriteButton.connect("clicked(bool)", self.onOverwriteButton)
    self.undoAllButton.connect("clicked(bool)", self.onUndoAllButton)
    self.undoButton.connect("clicked(bool)", self.onUndoButton)
    self.redoButton.connect("clicked(bool)", self.onRedoButton)
    
    for effect in warpEffects:
      effect.addEditButtonListeners(self)

    self.sceneRadioButton.connect("toggled(bool)", self.dataTreeTypeToggle)
    self.drawingsRadioButton.connect("toggled(bool)", self.dataTreeTypeToggle)
    self.atlasesRadioButton.connect("toggled(bool)", self.dataTreeTypeToggle)
    self.deleteTreeElement.connect("clicked(bool)", lambda i: self.dataTreeWidget.deleteSelectedItems())
    self.addTreeElement.connect("clicked(bool)", self.onAddTreeElement)
    self.renameTreeElement.connect("clicked(bool)", self.onRenameTreeElement)
    self.dataTreeWidget.doubleClicked.connect(self.onDataTreeDoubleClicked)
    self.dataTreeWidget.currentItemModified.connect(self.onDataTreeItemModified)

    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)    
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.NodeAddedEvent, self.onSceneNodeAdded)    
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.NodeAboutToBeAddedEvent, self.onSceneNodeAboutToBeAdded)    

    # Refresh
    qt.QApplication.processEvents()
    if self.updateMRMLFromArgs(): # was called from command line
      self.showSingleModule()
      tb = Toolbar.reducedToolbar()
      slicer.util.mainWindow().addToolBar(tb)

    self.updateGuiFromMRML()  
    self.onSceneNodeAdded()
    self.dataTreeTypeToggle(1)
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
    # get warp node and set selector and buttons
    warpID = self.parameterNode.GetParameter("warpID")
    warpNode = slicer.util.getNode(warpID) if warpID != "" else None
    warpNumberOfComponents = TransformsUtil.TransformsUtilLogic().getNumberOfLayers(warpNode)
    # undo redo button
    self.undoButton.setEnabled(warpNumberOfComponents > 1 and self.parameterNode.GetParameter("redoTransformID") == "" and self.parameterNode.GetParameter("lastOperation") != "undoall") 
    self.redoButton.setEnabled(self.parameterNode.GetParameter("redoTransformID") != "") 
    self.overwriteButton.setEnabled(warpNumberOfComponents > 1)
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
      self.parameterNode.SetParameter("drawingsRootItem", str(folderID))


  def onWarpSelectionChanged(self):
    self.parameterNode.SetParameter("warpID", self.warpSelector.currentNode().GetID() if self.warpSelector.currentNode() else "")

  def setItemChildrenFlags(self, item, flags):
    item.setFlags(flags)
    for row in range(item.rowCount()):
      self.setItemChildrenFlags(item.child(row), flags)      

  def onSceneNodeAdded(self,caller=None,event=None):
    sceneItem = self.dataTreeWidget.model().item(0,0)
    self.setItemChildrenFlags(sceneItem, qt.Qt.ItemIsSelectable + qt.Qt.ItemIsEnabled)

  def onDataTreeDoubleClicked(self, i):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    node = shNode.GetItemDataNode(self.dataTreeWidget.currentItem())
    centerList = [0] * 3
    # get center position of model/drawing
    if isinstance(node, slicer.vtkMRMLModelNode):
      pd = node.GetPolyData()
      center = vtk.vtkCenterOfMass()
      center.SetInputData(pd)
      center.Update()
      centerList = center.GetCenter()
    elif isinstance(node, slicer.vtkMRMLMarkupsCurveNode):
      node.GetNthControlPointPosition(round(node.GetNumberOfControlPoints()/2),centerList)
    elif isinstance(node, slicer.vtkMRMLMarkupsFiducialNode):
      node.GetNthFiducialPosition(0,centerList)
    else:
      return
    SmudgeModuleLogic().centerPosition(centerList)
    
  def onSceneNodeAboutToBeAdded(self,caller=None,event=None):
    # when adding fixed points while one of them is selected the new one is not set in the correct parent folder
    # this is overdoing, but fixes the problem
    self.dataTreeWidget.setCurrentItem(int(self.parameterNode.GetParameter("drawingsRootItem")))

  def dataTreeTypeToggle(self, b):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    if self.sceneRadioButton.isChecked():
      self.dataTreeWidget.nodeTypes = ('vtkMRMLModelNode','vtkMRMLFolderDisplayNode')
      self.dataTreeWidget.attributeNameFilter = ('')
      self.dataTreeWidget.attributeValueFilter = ('')
    elif self.atlasesRadioButton.isChecked():
      self.dataTreeWidget.nodeTypes = ('vtkMRMLModelNode','vtkMRMLFolderDisplayNode')
      self.dataTreeWidget.attributeNameFilter = ('atlas')
      self.dataTreeWidget.attributeValueFilter = ('')
    elif self.drawingsRadioButton.isChecked():
      self.dataTreeWidget.nodeTypes = ('vtkMRMLMarkupsCurveNode','vtkMRMLMarkupsFiducialNode','vtkMRMLFolderDisplayNode')
      self.dataTreeWidget.attributeNameFilter = ('drawing')
      self.dataTreeWidget.attributeValueFilter = ('1')
    # reset settings
    self.dataTreeWidget.expandToDepth(0)

  def onAddTreeElement(self):
    if self.atlasesRadioButton.isChecked() and self.parameterNode.GetParameter("MNIAtlasPath") != ".":
      items = ImportAtlas.ImportAtlasLogic().getValidAtlases(self.parameterNode.GetParameter("MNIAtlasPath"))
      result = BoolResult()
      atlasName = qt.QInputDialog.getItem(qt.QWidget(),'Select Atlas','',items,0,0,result)
      if result:
        ImportAtlas.ImportAtlasLogic().run(os.path.join(self.parameterNode.GetParameter("MNIAtlasPath"), atlasName))
    elif self.drawingsRadioButton.isChecked():
      # interaction node
      interactionNode = slicer.app.applicationLogic().GetInteractionNode()
      selectionNode = slicer.app.applicationLogic().GetSelectionNode()
      selectionNode.SetReferenceActivePlaceNodeClassName("vtkMRMLMarkupsFiducialNode")
      # create aux marpus node
      fiducialNode = slicer.vtkMRMLMarkupsFiducialNode()
      slicer.mrmlScene.AddNode(fiducialNode)
      fiducialNode.CreateDefaultDisplayNodes() 
      fiducialNode.GetDisplayNode().SetGlyphScale(2)
      fiducialNode.SetLocked(1)
      fiducialNode.SetName('Fixed Point')
      # add to subject hierarchy
      shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
      shNode.SetItemAttribute(shNode.GetItemByDataNode(fiducialNode), 'drawing', '1')
      shNode.SetItemParent(shNode.GetItemByDataNode(fiducialNode), int(self.parameterNode.GetParameter("drawingsRootItem")))
      # activate placement
      selectionNode.SetActivePlaceNodeID(fiducialNode.GetID())
      interactionNode.SetCurrentInteractionMode(interactionNode.Place)
      
  def onRenameTreeElement(self):
    sceneItem = self.dataTreeWidget.model().item(0,0)
    self.setItemChildrenFlags(sceneItem, qt.Qt.ItemIsSelectable + qt.Qt.ItemIsEnabled + qt.Qt.ItemIsEditable)
    self.dataTreeWidget.edit(self.dataTreeWidget.currentIndex())

  def onDataTreeItemModified(self):
    sceneItem = self.dataTreeWidget.model().item(0,0)
    self.setItemChildrenFlags(sceneItem, qt.Qt.ItemIsSelectable + qt.Qt.ItemIsEnabled)
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    node = shNode.GetItemDataNode(self.dataTreeWidget.currentItem())
    if isinstance(node, slicer.vtkMRMLMarkupsFiducialNode):
      node.SetNthControlPointLabel(0,node.GetName())


  def onUndoAllButton(self):
    self.parameterNode.SetParameter("lastOperation","undoall")
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
    if self.parameterNode.GetParameter("lastOperation") == 'snap':
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
    if self.parameterNode.GetParameter("lastOperation") == 'snap':
      SmudgeModuleLogic().enableLastDrawing()


  def onOverwriteButton(self):
    SmudgeModuleLogic().removeRedoNodes()
    warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))
    TransformsUtil.TransformsUtilLogic().flattenTransform(warpNode, True)
    self.parameterNode.SetParameter("lastOperation","flatten")
    self.updateGuiFromMRML() # update history


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
    node.SetParameter("SmudgeSigma", "0")
    node.SetParameter("expandEdge", "0")
    node.SetParameter("maxRadius", "50")
    # draw
    node.SetParameter("DrawSpread", "25")
    # blur
    node.SetParameter("BlurRadius", "25")
    node.SetParameter("BlurHardness", "40")
    node.SetParameter("BlurSigma", "10")
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




    





