import os
import unittest
import logging
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin

import numpy as np

from WarpDriveLib.Tools import NoneTool, SmudgeTool, DrawTool, PointToPointTool
from WarpDriveLib.Helpers import GridNodeHelper, LeadDBSCall
from WarpDriveLib.Widgets import Tables, Toolbar, Buttons

#
# WarpDrive
#

class WarpDrive(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "WarpDrive" 
    self.parent.categories = ["Netstim"]
    self.parent.dependencies = ["FiducialRegistrationVariableRBF"]
    self.parent.contributors = ["Simon Oxenford (Netstim Berlin)"]
    self.parent.helpText = """
This module provides tools to manually fix misalignments after non linear registration
"""
    self.parent.helpText += self.getDefaultModuleDocumentationLink()  # TODO: verify that the default URL is correct or change it to the actual documentation
    self.parent.acknowledgementText = "" 
    slicer.app.connect("startupCompleted()", setUpSliceNames)

def setUpSliceNames():
  if slicer.app.mainApplicationName == 'SlicerForLeadDBS':
    for color,name in zip(['Red','Green','Yellow'],['Axial','Coronal','Sagittal']):
      sliceWidget = slicer.app.layoutManager().sliceWidget(color)
      if not sliceWidget:
        continue
      sliceWidget.mrmlSliceNode().SetName(name)

#
# WarpDriveWidget
#

class WarpDriveWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent=None):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.__init__(self, parent)
    VTKObservationMixin.__init__(self)  # needed for parameter node observation
    self.logic = None
    self._parameterNode = None
    self._updatingGUIFromParameterNode = False

  def setup(self):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.setup(self)

    # Load widget from .ui file (created by Qt Designer)
    uiWidget = slicer.util.loadUI(self.resourcePath('UI/WarpDrive.ui'))
    self.layout.addWidget(uiWidget)
    self.ui = slicer.util.childWidgetVariables(uiWidget)

    # Add tools buttons
    toolsLayout = qt.QHBoxLayout(self.ui.toolsFrame)

    toolWidgets = [NoneTool.NoneToolWidget(),
                   SmudgeTool.SmudgeToolWidget(),
                   DrawTool.DrawToolWidget(),
                   PointToPointTool.PointToPointToolWidget()]

    for toolWidget in toolWidgets:
      toolsLayout.addWidget(toolWidget.effectButton)

    self.ui.drawModeMenu = toolWidgets[2].effectButton.menu()

    # Add Tree View
    correctionsLayout = qt.QVBoxLayout(self.ui.correctionsFrame)
    correctionsLayout.addWidget(Tables.WarpDriveCorrectionsManager())

    self.atlasesTable = Tables.AtlasesTable()
    atlasesLayout = qt.QVBoxLayout(self.ui.atlasesFrame)
    atlasesLayout.addWidget(self.atlasesTable)
    self.ui.tabWidget.currentChanged.connect(lambda i,a=self.atlasesTable: a.updateTable())

    # add cli progress bar
    self.ui.landwarpWidget = slicer.modules.fiducialregistrationvariablerbf.createNewWidgetRepresentation()
    self.ui.calculateFrame.layout().addWidget(self.ui.landwarpWidget.children()[3], 2, 0, 1, 2) # progress bar

    # Set scene in MRML widgets. Make sure that in Qt designer
    # "mrmlSceneChanged(vtkMRMLScene*)" signal in is connected to each MRML widget's.
    # "setMRMLScene(vtkMRMLScene*)" slot.
    uiWidget.setMRMLScene(slicer.mrmlScene)

    # Create a new parameterNode
    # This parameterNode stores all user choices in parameter values, node selections, etc.
    # so that when the scene is saved and reloaded, these settings are restored.
    self.logic = WarpDriveLogic()

    # Connections
    self.ui.calculateButton.connect('clicked(bool)', self.onCalculateButton)
    self.ui.spacingSameAsInputCheckBox.toggled.connect(lambda b: self.ui.spacingSpinBox.setEnabled(not b))

    # These connections ensure that whenever user changes some settings on the GUI, that is saved in the MRML scene
    # (in the selected parameter node).
    self.ui.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    self.ui.sourceFiducialsComboBox.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    self.ui.targetFiducialsComboBox.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    self.ui.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onOutputNodeChanged)
    self.ui.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    self.ui.radiusSlider.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
    self.ui.spacingSameAsInputCheckBox.connect("toggled(bool)", self.updateParameterNodeFromGUI)
    self.ui.spacingSpinBox.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
    self.ui.stiffnessSpinBox.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
    self.ui.drawModeMenu.triggered.connect(self.updateParameterNodeFromGUI)
    
    # MRML Scene
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)
    # self.addObserver(slicer.mrmlScene, slicer.mrmlScene.NodeAddedEvent, correctionsTree.updateTree)

    # Initial GUI update
    self.updateGUIFromParameterNode()

  def initializeCustomUI(self):
    
    # toolbars
    slicer.util.setToolbarsVisible(False, [])

    # customize view
    viewToolBar = slicer.util.mainWindow().findChild('QToolBar', 'ViewToolBar')
    viewToolBar.setVisible(1)
    layoutMenu = viewToolBar.widgetForAction(viewToolBar.actions()[0]).menu()
    for action in layoutMenu.actions():
      if action.text not in ['Four-Up', 'Tabbed slice']:
        layoutMenu.removeAction(action)

    # viewers
    viewersToolBar = slicer.util.mainWindow().findChild('QToolBar', 'ViewersToolBar')
    viewersToolBar.setVisible(1)

    # slicer window
    slicer.util.setMenuBarsVisible(False)
    slicer.util.setApplicationLogoVisible(False)
    slicer.util.setModuleHelpSectionVisible(False)
    slicer.util.setModulePanelTitleVisible(False)
    slicer.util.setDataProbeVisible(False)
    slicer.util.setPythonConsoleVisible(False)

    # inputs area
    self.ui.IOCollapsibleButton.setVisible(False)

    # data probe
    for i in range(slicer.mrmlScene.GetNumberOfNodesByClass("vtkMRMLScriptedModuleNode")):
      n  = slicer.mrmlScene.GetNthNodeByClass( i, "vtkMRMLScriptedModuleNode" )
      if n.GetModuleName() == "DataProbe":
        n.SetParameter('sliceViewAnnotationsEnabled','0')

    # set name
    if int(WarpDriveLogic().getParameterNode().GetParameter('SegmentMode')):
      slicer.util.mainWindow().setWindowTitle("Warp Drive - Segment")
    else:
      slicer.util.mainWindow().setWindowTitle("Warp Drive")
    slicer.util.mainWindow().showMaximized()
    qt.QApplication.processEvents()

    # Set linked slice views  in all existing slice composite nodes and in the default node
    sliceCompositeNodes = slicer.util.getNodesByClass("vtkMRMLSliceCompositeNode")
    defaultSliceCompositeNode = slicer.mrmlScene.GetDefaultNodeByClass("vtkMRMLSliceCompositeNode")
    if not defaultSliceCompositeNode:
      defaultSliceCompositeNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLSliceCompositeNode")
      defaultSliceCompositeNode.UnRegister(None)  # CreateNodeByClass is factory method, need to unregister the result to prevent memory leaks
      slicer.mrmlScene.AddDefaultNode(defaultSliceCompositeNode)
    sliceCompositeNodes.append(defaultSliceCompositeNode)
    for sliceCompositeNode in sliceCompositeNodes:
      sliceCompositeNode.SetLinkedControl(True)

    # start-up view
    for name in ['Axial','Coronal','Sagittal']:
      sliceWidget = slicer.app.layoutManager().sliceWidget(name)
      if not sliceWidget:
        continue
      fov = sliceWidget.mrmlSliceNode().GetFieldOfView()
      sliceWidget.mrmlSliceNode().SetFieldOfView(fov[0]/4,fov[1]/4,fov[2])
      if name == 'Axial':
        sliceWidget.mrmlSliceNode().SetXYZOrigin(0,-12,0)
        sliceWidget.sliceLogic().SetSliceOffset(-8.4)

    slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutTabbedSliceView)

    # Add custom ToolBar
    slicer.util.mainWindow().addToolBar(Toolbar.reducedToolbar())
    # Add add segmentation button
    button = Buttons.addSegmentationButton()
    self.atlasesTable.buttonsFrame.layout().insertWidget(1,button,1)
    button.clicked.connect(self.atlasesTable.updateTable)
    # Put all toolbars in same row
    for tb in slicer.util.mainWindow().findChildren('QToolBar'):
      slicer.util.mainWindow().removeToolBarBreak(tb)
    # Customize mouse mode
    mouseModeToolBar = slicer.util.mainWindow().findChild('QToolBar', 'MouseModeToolBar')
    mouseModeToolBar.setVisible(1)
    for a in mouseModeToolBar.actions():
      if a.text in ["Fiducial", "Toggle Markups Toolbar"]:
        mouseModeToolBar.removeAction(a)

    qt.QApplication.processEvents()

  def cleanup(self):
    """
    Called when the application closes and the module widget is destroyed.
    """
    self.cleanTools()
    self.removeObservers()

  def exit(self):
    """
    Called each time the user opens a different module.
    """
    # Do not react to parameter node changes (GUI will be updated when the user enters into the module)
    self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)
    self.cleanTools()

  def enter(self):
    """
    Called each time the user opens this module.
    """
    # Make sure parameter node exists and observed
    self.initializeParameterNode()

  def onSceneStartClose(self, caller=None, event=None):
    """
    Called just before the scene is closed.
    """
    # Parameter node will be reset, do not use it anymore
    self.setParameterNode(None)
    self.cleanTools()
        
  def cleanTools(self):
    # uncheck tools and cleanup
    for child in self.ui.toolsFrame.children():
      if isinstance(child, qt.QToolButton):
        child.setAutoExclusive(False)
        child.setChecked(False)
        child.setAutoExclusive(True)

  def onSceneEndClose(self, caller, event):
    """
    Called just after the scene is closed.
    """
    # If this module is shown while the scene is closed then recreate a new parameter node immediately
    if self.parent.isEntered:
      self.initializeParameterNode()

  def initializeParameterNode(self):
    """
    Ensure parameter node exists and observed.
    """
    # Parameter node stores all user choices in parameter values, node selections, etc.
    # so that when the scene is saved and reloaded, these settings are restored.

    self.setParameterNode(self.logic.getParameterNode())

  def setParameterNode(self, inputParameterNode):
    """
    Adds observers to the selected parameter node. Observation is needed because when the
    parameter node is changed then the GUI must be updated immediately.
    """

    if inputParameterNode:
      self.logic.setDefaultParameters(inputParameterNode)

    # Unobserve previously selected parameter node and add an observer to the newly selected.
    # Changes of parameter node are observed so that whenever parameters are changed by a script or any other module
    # those are reflected immediately in the GUI.
    if self._parameterNode is not None:
      self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)
    if inputParameterNode is not None:
      self.addObserver(inputParameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)
    self._parameterNode = inputParameterNode

    # Initial GUI update
    self.updateGUIFromParameterNode()

  def customUIWasInitialized(self):
    return len(list(filter(lambda x: isinstance(x,qt.QToolBar) and x.windowTitle=='LeadDBS', slicer.util.mainWindow().children())))

  def updateGUIFromParameterNode(self, caller=None, event=None):
    """
    This method is called whenever parameter node is changed.
    The module GUI is updated to show the current state of the parameter node.
    """

    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return

    # Make sure GUI changes do not call updateParameterNodeFromGUI (it could cause infinite loop)
    self._updatingGUIFromParameterNode = True

    # Set up toolbar for LeadDBS call
    if self._parameterNode.GetParameter("LeadSubjects") and not self.customUIWasInitialized():
      self.initializeCustomUI()

    # Update each widget from parameter node
    # Need to temporarily block signals to prevent infinite recursion (MRML node update triggers
    # GUI update, which triggers MRML node update, which triggers GUI update, ...)

    self.ui.inputSelector.setCurrentNode(self._parameterNode.GetNodeReference("InputNode"))
    self.ui.sourceFiducialsComboBox.setCurrentNode(self._parameterNode.GetNodeReference("SourceFiducial"))
    self.ui.targetFiducialsComboBox.setCurrentNode(self._parameterNode.GetNodeReference("TargetFiducial"))
    self.ui.outputSelector.setCurrentNode(self._parameterNode.GetNodeReference("OutputGridTransform"))

    radius = float(self._parameterNode.GetParameter("Radius"))
    self.ui.radiusSlider.value = radius
    if radius < self.ui.radiusSlider.minimum or radius > self.ui.radiusSlider.maximum:
      self.updateParameterNodeFromGUI()

    self.ui.spacingSpinBox.value = float(self._parameterNode.GetParameter("Spacing"))
    self.ui.stiffnessSpinBox.value = float(self._parameterNode.GetParameter("Stiffness"))

    self.ui.outputSelector.enabled = self._parameterNode.GetNodeReference("InputNode")
    self.ui.toolsCollapsibleButton.enabled = self._parameterNode.GetNodeReference("InputNode") and self._parameterNode.GetNodeReference("OutputGridTransform")
    self.ui.tabWidget.enabled = self._parameterNode.GetNodeReference("InputNode") and self._parameterNode.GetNodeReference("OutputGridTransform")
    self.ui.outputCollapsibleButton.enabled = self._parameterNode.GetNodeReference("InputNode") and self._parameterNode.GetNodeReference("OutputGridTransform")
    self.ui.calculateButton.enabled = self._parameterNode.GetNodeReference("InputNode") and self._parameterNode.GetNodeReference("OutputGridTransform")

    next(filter(lambda a: a.text == self._parameterNode.GetParameter("DrawMode"), self.ui.drawModeMenu.actions())).setChecked(True)

    # calculate warp
    if self._parameterNode.GetParameter("Update") == "true" and self._parameterNode.GetParameter("Running") == "false" and self.ui.autoUpdateCheckBox.checked:
      self.ui.calculateButton.animateClick()
    
    # set update to false
    self._parameterNode.SetParameter("Update", "false")

    # All the GUI updates are done
    self._updatingGUIFromParameterNode = False

  def updateParameterNodeFromGUI(self, caller=None, event=None):
    """
    This method is called when the user makes any change in the GUI.
    The changes are saved into the parameter node (so that they are restored when the scene is saved and loaded).
    """

    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return

    currentInputNode = self.ui.inputSelector.currentNode()
    if isinstance(currentInputNode, slicer.vtkMRMLTransformNode):
      if not (isinstance(currentInputNode.GetTransformFromParent(), slicer.vtkOrientedGridTransform) or isinstance(currentInputNode.GetTransformToParent(), slicer.vtkOrientedGridTransform)):
        qt.QMessageBox().warning(qt.QWidget(), "", "Select a Transform Node containing a Grid Transform")
        return

    wasModified = self._parameterNode.StartModify()  # Modify all properties in a single batch

    self._parameterNode.SetNodeReferenceID("InputNode", self.ui.inputSelector.currentNodeID)
    self._parameterNode.SetNodeReferenceID("SourceFiducial", self.ui.sourceFiducialsComboBox.currentNodeID)
    self._parameterNode.SetNodeReferenceID("TargetFiducial", self.ui.targetFiducialsComboBox.currentNodeID)
    self._parameterNode.SetNodeReferenceID("OutputGridTransform", self.ui.outputSelector.currentNodeID)
    self._parameterNode.SetParameter("DrawMode", next(filter(lambda a: a.checked, self.ui.drawModeMenu.actions())).text)
    self._parameterNode.SetParameter("Radius", "%.2f" % self.ui.radiusSlider.value)
    self._parameterNode.SetParameter("Stiffness", str(self.ui.stiffnessSpinBox.value))
    # spacing
    if self.ui.spacingSameAsInputCheckBox.checked:
      size,origin,spacing,directionMatrix = GridNodeHelper.getGridDefinition(currentInputNode)
    else:
      spacing = [self.ui.spacingSpinBox.value]
    self._parameterNode.SetParameter("Spacing", str(spacing[0])) 

    self._parameterNode.EndModify(wasModified)


  def onOutputNodeChanged(self):
    # unset if output is the same as input
    currentNodeID = self.ui.outputSelector.currentNodeID
    if currentNodeID == self.ui.inputSelector.currentNodeID:
      wasBlocked = self.ui.outputSelector.blockSignals(True)
      self.ui.outputSelector.currentNodeID = None
      self.ui.outputSelector.blockSignals(wasBlocked)
    # observe
    if self.ui.inputSelector.currentNodeID:
      self.ui.inputSelector.currentNode().SetAndObserveTransformNodeID(self.ui.outputSelector.currentNodeID)


  def onCalculateButton(self):
    """
    Run processing when user clicks "Apply" button.
    """
    # cursor
    qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
    qt.QApplication.processEvents()

    # create nodes
    sourceFiducial = self._parameterNode.GetNodeReference("SourceFiducial")
    targetFiducial = self._parameterNode.GetNodeReference("TargetFiducial")
    # reference
    size,origin,spacing,directionMatrix = GridNodeHelper.getGridDefinition(self._parameterNode.GetNodeReference("InputNode"))
    userSpacing = np.ones(3) * float(self._parameterNode.GetParameter("Spacing"))
    size = size * (spacing / userSpacing)
    auxVolumeNode = GridNodeHelper.emptyVolume(size.astype(int), origin, userSpacing, directionMatrix)
    # output
    outputNode = self._parameterNode.GetNodeReference("OutputGridTransform")
    # params
    RBFRadius = []
    for i in range(targetFiducial.GetNumberOfControlPoints()):
      if targetFiducial.GetNthControlPointSelected(i):
        RBFRadius.append(targetFiducial.GetNthControlPointDescription(i))
    RBFRadius = ",".join(RBFRadius)
    stiffness = float(self._parameterNode.GetParameter("Stiffness"))

    # save current state if leadDBS call in case of error
    if self._parameterNode.GetParameter("subjectPath") != '':
      LeadDBSCall.saveSourceTarget(self._parameterNode.GetParameter("subjectPath"), sourceFiducial, targetFiducial)

    # preview
    visualizationNodes = self.logic.previewWarp(sourceFiducial, targetFiducial)
    qt.QApplication.processEvents()

    self._parameterNode.SetParameter("Running", "true")
    cliNode = self.logic.run(auxVolumeNode, outputNode, sourceFiducial, targetFiducial, RBFRadius, stiffness)

    if cliNode is not None:
      # set up for UI
      self.ui.landwarpWidget.setCurrentCommandLineModuleNode(cliNode)
      # add observer
      cliNode.AddObserver(slicer.vtkMRMLCommandLineModuleNode.StatusModifiedEvent, \
        lambda c,e,o=outputNode,v=visualizationNodes,a=auxVolumeNode: self.onStatusModifiedEvent(c,o,v,a))
    else:
      self.onStatusModifiedEvent(None,outputNode,visualizationNodes,auxVolumeNode)
    
  
  def onStatusModifiedEvent(self, caller, outputNode, visualizationNodes, auxVolumeNode):
    
    if isinstance(caller, slicer.vtkMRMLCommandLineModuleNode):
      if caller.GetStatusString() == 'Completed':
        # delete cli Node
        qt.QTimer.singleShot(1000, lambda: slicer.mrmlScene.RemoveNode(caller))
      else:
        return

    self._parameterNode.GetNodeReference("InputNode").SetAndObserveTransformNodeID(outputNode.GetID())
    self._parameterNode.GetNodeReference("InputNode").Modified()

    # remove aux
    for node in visualizationNodes:
      slicer.mrmlScene.RemoveNode(node)
    slicer.mrmlScene.RemoveNode(auxVolumeNode)

    qt.QApplication.setOverrideCursor(qt.Qt.ArrowCursor)

    self._parameterNode.SetParameter("Running", "false")



#
# WarpDriveLogic
#

class WarpDriveLogic(ScriptedLoadableModuleLogic):
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
    if slicer.util.settingsValue('Developer/DeveloperMode', False, converter=slicer.util.toBool):
      import WarpDriveLib
      import importlib
      import glob
      warpDrivePath = os.path.split(__file__)[0]
      G = glob.glob(os.path.join(warpDrivePath, 'WarpDriveLib','**','*.py'))
      for g in G:
        relativePath = os.path.relpath(g, warpDrivePath) # relative path
        relativePath = os.path.splitext(relativePath)[0] # get rid of .py
        moduleParts = relativePath.split(os.path.sep) # separate
        importlib.import_module('.'.join(moduleParts)) # import module
        module = WarpDriveLib
        for i in range(1,len(moduleParts)): # iterate over parts in order to load subpkgs
          module = getattr(module, moduleParts[i])
        importlib.reload(module) # reload
    

  def setDefaultParameters(self, parameterNode):
    """
    Initialize parameter node with default settings.
    """
    if not parameterNode.GetParameter("Radius"):
      parameterNode.SetParameter("Radius", "15.0")
    if not parameterNode.GetParameter("Spacing"):
      parameterNode.SetParameter("Spacing", "2.0")
    if not parameterNode.GetParameter("RBFRadius"):
      parameterNode.SetParameter("RBFRadius", "30")
    if not parameterNode.GetParameter("Stiffness"):
      parameterNode.SetParameter("Stiffness", "0.1")
    if not parameterNode.GetParameter("DrawMode"):
      parameterNode.SetParameter("DrawMode", 'To Nearest Model')
    if not parameterNode.GetParameter("Running"):
      parameterNode.SetParameter("Running", "false")

  def run(self, referenceVolume, outputNode, sourceFiducial, targetFiducial, RBFRadius, stiffness):

    # run landmark registration if points available
    if RBFRadius != "":
      cliNode = self.computeWarp(referenceVolume, outputNode, sourceFiducial, targetFiducial, RBFRadius, stiffness)
    else:
      size, origin, spacing, directionMatrix = GridNodeHelper.getGridDefinition(referenceVolume)
      GridNodeHelper.emptyGridTransform(size, origin, spacing, directionMatrix, outputNode)
      return
    return cliNode

  def computeWarp(self, referenceVolume, outputNode, sourceFiducial, targetFiducial, RBFRadius, stiffness):

    # Compute the warp with FiducialRegistrationVariableRBF
    cliParams = {
      "referenceVolume" : referenceVolume.GetID(),
      "fixedFiducials" : targetFiducial.GetID(),
      "movingFiducials" : sourceFiducial.GetID(),
      "outputDisplacementField" : outputNode.GetID(),
      "RBFRadius" : RBFRadius,
      "stiffness" : stiffness,
      } 

    cliNode = slicer.cli.run(slicer.modules.fiducialregistrationvariablerbf, None, cliParams, wait_for_completion=False, update_display=False)

    return cliNode

  def previewWarp(self, source, target):
    if isinstance(source, slicer.vtkMRMLMarkupsFiducialNode) and isinstance(target, slicer.vtkMRMLMarkupsFiducialNode):
      sourcePoints = vtk.vtkPoints()
      targetPoints = vtk.vtkPoints()
      for i in range(target.GetNumberOfControlPoints()):
        if target.GetNthControlPointSelected(i):
          sourcePoints.InsertNextPoint(source.GetNthControlPointPosition(i))
          targetPoints.InsertNextPoint(target.GetNthControlPointPosition(i))
    else:
      sourcePoints = source
      targetPoints = target
    sourceDisplayFiducial = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    sourceDisplayFiducial.GetDisplayNode().SetVisibility(0)
    sourceDisplayFiducial.SetControlPointPositionsWorld(sourcePoints)
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
    transformNode.GetDisplayNode().SetAndObserveGlyphPointsNode(sourceDisplayFiducial)
    transformNode.GetDisplayNode().SetVisibility2D(1)
    return transformNode, sourceDisplayFiducial


#
# WarpDriveTest
#

class WarpDriveTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)
    import WarpDrive
    WarpDrive.WarpDriveLogic().getParameterNode().SetParameter('LeadSubjects','')


  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_WarpDrive1()

  def test_WarpDrive1(self):
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

    json_txt = '[{"id":"leadDCM","warpdrive_path":"/Users/simon/Documents/leadDS/derivatives/leaddbs/sub-leadDCM/warpdrive","normlog_file":"/Users/simon/Documents/leadDS/derivatives/leaddbs/sub-leadDCM/normalization/log/sub-leadDCM_desc-normmethod.json","anat_files":{"iso_T1w":"/Users/simon/Documents/leadDS/derivatives/leaddbs/sub-leadDCM/coregistration/anat/sub-leadDCM_ses-preop_space-anchorNative_desc-preproc_acq-iso_T1w.nii"},"forward_transform":"/Users/simon/Documents/leadDS/derivatives/leaddbs/sub-leadDCM/normalization/transformations/sub-leadDCM_from-anchorNative_to-MNI152NLin2009bAsym_desc-ants.nii.gz","inverse_transform":"/Users/simon/Documents/leadDS/derivatives/leaddbs/sub-leadDCM/normalization/transformations/sub-leadDCM_from-MNI152NLin2009bAsym_to-anchorNative_desc-ants.nii.gz"},{"id":"leadDCM","warpdrive_path":"/Users/simon/Documents/leadDS/derivatives/leaddbs/sub-leadDCM/warpdrive","normlog_file":"/Users/simon/Documents/leadDS/derivatives/leaddbs/sub-leadDCM/normalization/log/sub-leadDCM_desc-normmethod.json","anat_files":{"iso_T1w":"/Users/simon/Documents/leadDS/derivatives/leaddbs/sub-leadDCM/coregistration/anat/sub-leadDCM_ses-preop_space-anchorNative_desc-preproc_acq-iso_T1w.nii"},"forward_transform":"/Users/simon/Documents/leadDS/derivatives/leaddbs/sub-leadDCM/normalization/transformations/sub-leadDCM_from-anchorNative_to-MNI152NLin2009bAsym_desc-ants.nii.gz","inverse_transform":"/Users/simon/Documents/leadDS/derivatives/leaddbs/sub-leadDCM/normalization/transformations/sub-leadDCM_from-MNI152NLin2009bAsym_to-anchorNative_desc-ants.nii.gz"}]'
    # json_txt = '{"id":"leadDCM","warpdrive_path":"/Users/simon/Documents/leadDS/derivatives/leaddbs/sub-leadDCM/warpdrive","normlog_file":"/Users/simon/Documents/leadDS/derivatives/leaddbs/sub-leadDCM/normalization/log/sub-leadDCM_desc-normmethod.json","anat_files":{"iso_T1w":"/Users/simon/Documents/leadDS/derivatives/leaddbs/sub-leadDCM/coregistration/anat/sub-leadDCM_ses-preop_space-anchorNative_desc-preproc_acq-iso_T1w.nii"},"forward_transform":"/Users/simon/Documents/leadDS/derivatives/leaddbs/sub-leadDCM/normalization/transformations/sub-leadDCM_from-anchorNative_to-MNI152NLin2009bAsym_desc-ants.nii.gz","inverse_transform":"/Users/simon/Documents/leadDS/derivatives/leaddbs/sub-leadDCM/normalization/transformations/sub-leadDCM_from-MNI152NLin2009bAsym_to-anchorNative_desc-ants.nii.gz"}'
    
    # import json
    # json_txt = json.dumps(json.load(open("C:\\Users\\simon\\Desktop\\.warpdrive_tmp.json")))
    # json_txt = json.dumps(json.load(open("C:\\Users\\simon\\Desktop\\.warpdrive_tmp2.json")))

    import WarpDrive
    parameterNode = WarpDrive.WarpDriveLogic().getParameterNode()
    wasModified = parameterNode.StartModify()
    parameterNode.SetParameter('CurrentSubject','')
    parameterNode.SetParameter('LeadSubjects',json_txt)
    parameterNode.SetParameter('MNIPath','/Users/simon/repo/leaddbs/templates/space/MNI152NLin2009bAsym/')
    parameterNode.EndModify(wasModified)

    self.delayDisplay('Test passed')
