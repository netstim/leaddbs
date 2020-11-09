import os
import unittest
import logging
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin

import numpy as np
import importlib
import glob

from WarpDriveLib.Tools import NoneTool, SmudgeTool, DrawTool, PointToPointTool
from WarpDriveLib.Helpers import GridNodeHelper, WarpDriveUtil, LeadDBSCall
from WarpDriveLib.Widgets import TreeView, Toolbar

#
# WarpDrive
#

class WarpDrive(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "WarpDrive"  # TODO: make this more human readable by adding spaces
    self.parent.categories = ["Netstim"]  # TODO: set categories (folders where the module shows up in the module selector)
    self.parent.dependencies = []  # TODO: add here list of module names that this module requires
    self.parent.contributors = ["John Doe (AnyWare Corp.)"]  # TODO: replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
This is an example of scripted loadable module bundled in an extension.
"""  # TODO: update with short description of the module
    self.parent.helpText += self.getDefaultModuleDocumentationLink()  # TODO: verify that the default URL is correct or change it to the actual documentation
    self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc., Andras Lasso, PerkLab,
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
"""  # TODO: replace with organization, grant and thanks.

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

    # Add Tree View
    dataControlTree = TreeView.WarpDriveTreeView()
    dataControlLayout = qt.QVBoxLayout(self.ui.dataControlFrame)
    dataControlLayout.addWidget(dataControlTree)

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
    self.ui.autoRBFRadiusCheckBox.toggled.connect(lambda b: self.ui.RBFRadiusSpinBox.setEnabled(not b))

    # These connections ensure that whenever user changes some settings on the GUI, that is saved in the MRML scene
    # (in the selected parameter node).
    self.ui.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    self.ui.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onOutputNodeChanged)
    self.ui.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    self.ui.spreadSlider.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
    self.ui.spacingSameAsInputCheckBox.connect("toggled(bool)", self.updateParameterNodeFromGUI)
    self.ui.spacingSpinBox.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
    self.ui.autoRBFRadiusCheckBox.connect("toggled(bool)", self.updateParameterNodeFromGUI)
    self.ui.RBFRadiusSpinBox.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
    self.ui.stiffnessSpinBox.connect("valueChanged(double)", self.updateParameterNodeFromGUI)

    # MRML Scene
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.NodeAddedEvent, dataControlTree.updateTree)

    # check dependencies
    if slicer.app.mainApplicationName != 'SlicerCustom':
      if LeadDBSCall.checkExtensionInstall(extensionName = 'SlicerRT'):
        return
      if LeadDBSCall.checkExtensionInstall(extensionName = 'MarkupsToModel'):
        return

    # Initial GUI update
    self.updateGUIFromParameterNode()



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
    list(map(lambda a: mouseModeToolBar.removeAction(a), filter(lambda a: a.text=="Place Fiducial", mouseModeToolBar.actions())))

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
    self.ui.IOCollapsibleButton.setVisible(not singleModule)

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
    # Do not react to parameter node changes (GUI wlil be updated when the user enters into the module)
    self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)
    self.cleanTools()

  def enter(self):
    """
    Called each time the user opens this module.
    """
    # Make sure parameter node exists and observed
    self.initializeParameterNode()
    # Lead-DBS call
    if LeadDBSCall.updateParameterNodeFromArgs(self._parameterNode): # was called from command line
      self.showSingleModule()
      slicer.util.mainWindow().addToolBar(Toolbar.reducedToolbar())

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

    # Unobserve previusly selected parameter node and add an observer to the newly selected.
    # Changes of parameter node are observed so that whenever parameters are changed by a script or any other module
    # those are reflected immediately in the GUI.
    if self._parameterNode is not None:
      self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)
    if inputParameterNode is not None:
      self.addObserver(inputParameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)
    self._parameterNode = inputParameterNode

    # Initial GUI update
    self.updateGUIFromParameterNode()

  def updateGUIFromParameterNode(self, caller=None, event=None):
    """
    This method is called whenever parameter node is changed.
    The module GUI is updated to show the current state of the parameter node.
    """

    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return

    # Make sure GUI changes do not call updateParameterNodeFromGUI (it could cause infinite loop)
    self._updatingGUIFromParameterNode = True

    # Update each widget from parameter node
    # Need to temporarily block signals to prevent infinite recursion (MRML node update triggers
    # GUI update, which triggers MRML node update, which triggers GUI update, ...)

    self.ui.inputSelector.setCurrentNode(self._parameterNode.GetNodeReference("InputNode"))
    self.ui.outputSelector.setCurrentNode(self._parameterNode.GetNodeReference("OutputGridTransform"))

    spread = float(self._parameterNode.GetParameter("Spread"))
    self.ui.spreadSlider.value = spread
    if spread < self.ui.spreadSlider.minimum or spread > self.ui.spreadSlider.maximum:
      self.updateParameterNodeFromGUI()

    self.ui.spacingSpinBox.value = float(self._parameterNode.GetParameter("Spacing"))
    self.ui.RBFRadiusSpinBox.value = float(self._parameterNode.GetParameter("RBFRadius"))
    self.ui.stiffnessSpinBox.value = float(self._parameterNode.GetParameter("Stiffness"))

    self.ui.outputSelector.enabled = self._parameterNode.GetNodeReference("InputNode")
    self.ui.toolsCollapsibleButton.enabled = self._parameterNode.GetNodeReference("InputNode") and self._parameterNode.GetNodeReference("OutputGridTransform")
    self.ui.outputCollapsibleButton.enabled = self._parameterNode.GetNodeReference("InputNode") and self._parameterNode.GetNodeReference("OutputGridTransform")
    self.ui.calculateButton.enabled = self._parameterNode.GetNodeReference("InputNode") and self._parameterNode.GetNodeReference("OutputGridTransform")

    # calculate warp
    if self._parameterNode.GetParameter("Update") == "true":
      self.ui.calculateButton.animateClick()
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

    wasModified = self._parameterNode.StartModify()  # Modify all properties in a single batch

    self._parameterNode.SetNodeReferenceID("InputNode", self.ui.inputSelector.currentNodeID)
    self._parameterNode.SetNodeReferenceID("OutputGridTransform", self.ui.outputSelector.currentNodeID)
    self._parameterNode.SetParameter("Spread", str(self.ui.spreadSlider.value))
    self._parameterNode.SetParameter("Stiffness", str(self.ui.stiffnessSpinBox.value))
    # spacing
    if self.ui.spacingSameAsInputCheckBox.checked:
      size,origin,spacing = GridNodeHelper.getGridDefinition(self.ui.inputSelector.currentNode())
    else:
      spacing = [self.ui.spacingSpinBox.value]
    self._parameterNode.SetParameter("Spacing", str(spacing[0]))
    # # RBF radius
    if self.ui.autoRBFRadiusCheckBox.checked:
      radius = WarpDriveUtil.getMaxSpread()
    else:
      radius = self.ui.RBFRadiusSpinBox.value
    self._parameterNode.SetParameter("RBFRadius", str(radius))    

    self._parameterNode.EndModify(wasModified)


  def onOutputNodeChanged(self):
    # unset if output is the same as input
    currentNodeID = self.ui.outputSelector.currentNodeID
    if currentNodeID == self.ui.inputSelector.currentNodeID:
      wasBlocked = self.ui.outputSelector.blockSignals(True)
      self.ui.outputSelector.currentNodeID = None
      self.ui.outputSelector.blockSignals(wasBlocked)
    # obvserve
    if self.ui.inputSelector.currentNodeID:
      self.ui.inputSelector.currentNode().SetAndObserveTransformNodeID(self.ui.outputSelector.currentNodeID)


  def onCalculateButton(self):
    """
    Run processing when user clicks "Apply" button.
    """
    # cursor
    qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
    qt.QApplication.processEvents()

    # update (sets the RBF value)
    self.updateParameterNodeFromGUI()
    # get source and target from points
    sourcePoints = WarpDriveUtil.getPointsFromAttribute('source')
    targetPoints = WarpDriveUtil.getPointsFromAttribute('target')
    fixedPoints = WarpDriveUtil.getPointsFromAttribute('fixed')
    # add fixed to source and target
    sourcePoints.InsertPoints(sourcePoints.GetNumberOfPoints(), fixedPoints.GetNumberOfPoints(), 0, fixedPoints)
    targetPoints.InsertPoints(targetPoints.GetNumberOfPoints(), fixedPoints.GetNumberOfPoints(), 0, fixedPoints)
    # create nodes
    sourceFiducial = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    sourceFiducial.SetControlPointPositionsWorld(sourcePoints)
    sourceFiducial.GetDisplayNode().SetTextScale(0)
    targetFiducial = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    targetFiducial.SetControlPointPositionsWorld(targetPoints)
    targetFiducial.GetDisplayNode().SetVisibility(0)
    # reference
    size,origin,spacing = GridNodeHelper.getGridDefinition(self._parameterNode.GetNodeReference("InputNode"))
    userSpacing = [float(self._parameterNode.GetParameter("Spacing"))] * 3
    size = np.array(size) * (np.array(spacing) / np.array(userSpacing))
    auxVolumeNode = GridNodeHelper.emptyVolume([int(s) for s in size], origin, userSpacing)
    # output
    outputNode = self._parameterNode.GetNodeReference("OutputGridTransform")
    if outputNode.GetDisplayNode():
      visibility = outputNode.GetDisplayNode().GetVisibility()
      outputNode.GetDisplayNode().SetVisibility(0)
    else:
      visibility = None
    # params
    RBFRadius = float(self._parameterNode.GetParameter("Spread"))
    stiffness = float(self._parameterNode.GetParameter("Stiffness"))
    # mask
    maskVolume = WarpDriveUtil.getMaskVolume(auxVolumeNode)

    # save current state if leadDBS call in case of error
    if self._parameterNode.GetParameter("subjectPath") != '':
      LeadDBSCall.saveCurrentScene(self._parameterNode.GetParameter("subjectPath"))

    # unset current warp
    self._parameterNode.GetNodeReference("InputNode").SetAndObserveTransformNodeID(None)

    # preview
    visualizationNode = WarpDriveUtil.previewWarp(sourceFiducial, targetFiducial, outputNode)

    # run
    self.logic.run(auxVolumeNode, outputNode, sourceFiducial, targetFiducial, RBFRadius, stiffness, maskVolume)

    # set new warp
    self._parameterNode.GetNodeReference("InputNode").SetAndObserveTransformNodeID(outputNode.GetID())
    # set visibility
    if visibility:
      outputNode.GetDisplayNode().SetVisibility(visibility)

    # remove aux
    slicer.mrmlScene.RemoveNode(maskVolume) 
    slicer.mrmlScene.RemoveNode(visualizationNode)
    slicer.mrmlScene.RemoveNode(sourceFiducial)
    slicer.mrmlScene.RemoveNode(targetFiducial)
    slicer.mrmlScene.RemoveNode(auxVolumeNode)

    # cursor
    qt.QApplication.setOverrideCursor(qt.Qt.ArrowCursor)
    
    


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
    if not parameterNode.GetParameter("Spread"):
      parameterNode.SetParameter("Spread", "15.0")
    if not parameterNode.GetParameter("Spacing"):
      parameterNode.SetParameter("Spacing", "2.0")
    if not parameterNode.GetParameter("RBFRadius"):
      parameterNode.SetParameter("RBFRadius", "30")
    if not parameterNode.GetParameter("Stiffness"):
      parameterNode.SetParameter("Stiffness", "0.1")

  def run(self, referenceVolume, outputNode, sourceFiducial, targetFiducial, RBFRadius, stiffness, maskVolume):

    # run landmark registration if points available
    if sourceFiducial.GetNumberOfControlPoints():
      self.computeWarp(referenceVolume, outputNode, sourceFiducial, targetFiducial, RBFRadius, stiffness)
    else:
      GridNodeHelper.emptyGridTransform(referenceVolume.GetImageData().GetDimensions(), referenceVolume.GetOrigin(), referenceVolume.GetSpacing(), outputNode)
      return

    # get arrays
    transformArray = slicer.util.array(outputNode.GetID())
    maskArray = slicer.util.array(maskVolume.GetID())
    # mask
    transformArray[:] = np.stack([transformArray[:,:,:,i] * maskArray for i in range(3)], 3).squeeze()

    # modified
    outputNode.Modified()
    referenceVolume.Modified()


  def computeWarp(self, referenceVolume, outputNode, sourceFiducial, targetFiducial, RBFRadius, stiffness):
    """
    Run the processing algorithm.
    Can be used without GUI widget.
    :param referenceVolume: Used to set grid definition
    :param outputNode: output warp. will be observed by input node
    :param sourceFiducial: source fiducials
    :param targetFiducial: target fiducials
    :param spread: used for RBF radius
    """

    if not referenceVolume or not outputNode:
      raise ValueError("Input or output is invalid")

    logging.info('Processing started')

    # Compute the warp with plastimatch landwarp
    cliParams = {
      "plmslc_landwarp_fixed_volume" : referenceVolume.GetID(),
      "plmslc_landwarp_moving_volume" : referenceVolume.GetID(),
      "plmslc_landwarp_fixed_fiducials" : targetFiducial.GetID(),
      "plmslc_landwarp_moving_fiducials" : sourceFiducial.GetID(),
      "plmslc_landwarp_output_vf" : outputNode.GetID(),
      "plmslc_landwarp_rbf_type" : "gauss",
      "plmslc_landwarp_rbf_radius" : RBFRadius,
      "plmslc_landwarp_stiffness" : stiffness,
      } 

    cliNode = slicer.cli.run(slicer.modules.plastimatch_slicer_landwarp, None, cliParams, wait_for_completion=True, update_display=False)

    logging.info('Processing completed')



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

    # Get/create input data

    import SampleData
    InputNode = SampleData.downloadFromURL(
      nodeNames='MRHead',
      fileNames='MR-Head.nrrd',
      uris='https://github.com/Slicer/SlicerTestingData/releases/download/MD5/39b01631b7b38232a220007230624c8e',
      checksums='MD5:39b01631b7b38232a220007230624c8e')[0]
    self.delayDisplay('Finished with download and loading')

    inputScalarRange = InputNode.GetImageData().GetScalarRange()
    self.assertEqual(inputScalarRange[0], 0)
    self.assertEqual(inputScalarRange[1], 279)

    OutputGridTransform = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
    threshold = 50

    # Test the module logic

    logic = WarpDriveLogic()

    # Test algorithm with non-inverted threshold
    logic.run(InputNode, OutputGridTransform, threshold, True)
    outputScalarRange = OutputGridTransform.GetImageData().GetScalarRange()
    self.assertEqual(outputScalarRange[0], inputScalarRange[0])
    self.assertEqual(outputScalarRange[1], threshold)

    # Test algorithm with inverted threshold
    logic.run(InputNode, OutputGridTransform, threshold, False)
    outputScalarRange = OutputGridTransform.GetImageData().GetScalarRange()
    self.assertEqual(outputScalarRange[0], inputScalarRange[0])
    self.assertEqual(outputScalarRange[1], inputScalarRange[1])

    self.delayDisplay('Test passed')
