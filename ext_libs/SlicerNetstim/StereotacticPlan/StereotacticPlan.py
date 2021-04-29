import os
import unittest
import logging
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin

import numpy as np
import StereotacticPlanLib.util

#
# StereotacticPlan
#

class StereotacticPlan(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "StereotacticPlan"  # TODO: make this more human readable by adding spaces
    self.parent.categories = ["Netstim"]  # TODO: set categories (folders where the module shows up in the module selector)
    self.parent.dependencies = []  # TODO: add here list of module names that this module requires
    self.parent.contributors = ["Simon Oxenford (Charite Berlin.)"]  # TODO: replace with "Firstname Lastname (Organization)"
    # TODO: update with short description of the module and a link to online module documentation
    self.parent.helpText = """
This is an example of scripted loadable module bundled in an extension.
See more information in <a href="https://github.com/organization/projectname#StereotacticPlan">module documentation</a>.
"""
    # TODO: replace with organization, grant and thanks
    self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc., Andras Lasso, PerkLab,
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
"""

    # Additional initialization step after application startup is complete
    # slicer.app.connect("startupCompleted()", registerSampleData)

#
# Register sample data sets in Sample Data module
#

def registerSampleData():
  """
  Add data sets to Sample Data module.
  """
  # It is always recommended to provide sample data for users to make it easy to try the module,
  # but if no sample data is available then this method (and associated startupCompeted signal connection) can be removed.

  import SampleData
  iconsPath = os.path.join(os.path.dirname(__file__), 'Resources/Icons')

  # To ensure that the source code repository remains small (can be downloaded and installed quickly)
  # it is recommended to store data sets that are larger than a few MB in a Github release.

  # StereotacticPlan1
  SampleData.SampleDataLogic.registerCustomSampleDataSource(
    # Category and sample name displayed in Sample Data module
    category='StereotacticPlan',
    sampleName='StereotacticPlan1',
    # Thumbnail should have size of approximately 260x280 pixels and stored in Resources/Icons folder.
    # It can be created by Screen Capture module, "Capture all views" option enabled, "Number of images" set to "Single".
    thumbnailFileName=os.path.join(iconsPath, 'StereotacticPlan1.png'),
    # Download URL and target file name
    uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95",
    fileNames='StereotacticPlan1.nrrd',
    # Checksum to ensure file integrity. Can be computed by this command:
    #  import hashlib; print(hashlib.sha256(open(filename, "rb").read()).hexdigest())
    checksums = 'SHA256:998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95',
    # This node name will be used when the data set is loaded
    nodeNames='StereotacticPlan1'
  )

  # StereotacticPlan2
  SampleData.SampleDataLogic.registerCustomSampleDataSource(
    # Category and sample name displayed in Sample Data module
    category='StereotacticPlan',
    sampleName='StereotacticPlan2',
    thumbnailFileName=os.path.join(iconsPath, 'StereotacticPlan2.png'),
    # Download URL and target file name
    uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97",
    fileNames='StereotacticPlan2.nrrd',
    checksums = 'SHA256:1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97',
    # This node name will be used when the data set is loaded
    nodeNames='StereotacticPlan2'
  )

#
# StereotacticPlanWidget
#

class StereotacticPlanWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
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

    # Load widget from .ui file (created by Qt Designer).
    # Additional widgets can be instantiated manually and added to self.layout.
    uiWidget = slicer.util.loadUI(self.resourcePath('UI/StereotacticPlan.ui'))
    self.layout.addWidget(uiWidget)
    self.ui = slicer.util.childWidgetVariables(uiWidget)

    # Set scene in MRML widgets. Make sure that in Qt designer the top-level qMRMLWidget's
    # "mrmlSceneChanged(vtkMRMLScene*)" signal in is connected to each MRML widget's.
    # "setMRMLScene(vtkMRMLScene*)" slot.
    uiWidget.setMRMLScene(slicer.mrmlScene)

    # Create logic class. Logic implements all computations that should be possible to run
    # in batch mode, without a graphical user interface.
    self.logic = StereotacticPlanLogic()

    # Connections

    # These connections ensure that we update parameter node when scene is closed
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)

    # Preview line CHeckbox
    self.ui.previewLineCheckBox.connect('toggled(bool)', self.onPreviewLineToggled)
    self.ui.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updatePreviewLineTransform)

    # These connections ensure that whenever user changes some settings on the GUI, that is saved in the MRML scene
    # (in the selected parameter node).
    self.ui.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    self.ui.arcAngleSliderWidget.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
    self.ui.ringAngleSliderWidget.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
    self.ui.headringCoordinatesWidget.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
    self.ui.mountingComboBox.currentIndexChanged.connect(self.updateParameterNodeFromGUI)
    self.ui.autoUpdateCheckBox.connect('toggled(bool)', self.updateParameterNodeFromGUI)

    # Buttons
    self.ui.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.ui.importFromPDFPushButton.connect('clicked(bool)', self.onImportFromPDFPushButton)
    self.ui.setDefaultResliceDriverPushButton.connect('clicked(bool)', self.onSetDefaultResliceDriver)

    # Make sure parameter node is initialized (needed for module reload)
    self.initializeParameterNode()

  def cleanup(self):
    """
    Called when the application closes and the module widget is destroyed.
    """
    self.removeObservers()

  def enter(self):
    """
    Called each time the user opens this module.
    """
    # Make sure parameter node exists and observed
    self.initializeParameterNode()

  def exit(self):
    """
    Called each time the user opens a different module.
    """
    # Remove Preview node
    self.ui.previewLineCheckBox.setChecked(False)
    # Do not react to parameter node changes (GUI wlil be updated when the user enters into the module)
    self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

  def onSceneStartClose(self, caller, event):
    """
    Called just before the scene is closed.
    """
    # Parameter node will be reset, do not use it anymore
    self.setParameterNode(None)

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

    # Select default input nodes if nothing is selected yet to save a few clicks for the user
    if not self._parameterNode.GetNodeReference("InputVolume"):
      firstVolumeNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLScalarVolumeNode")
      if firstVolumeNode:
        self._parameterNode.SetNodeReferenceID("InputVolume", firstVolumeNode.GetID())

  def setParameterNode(self, inputParameterNode):
    """
    Set and observe parameter node.
    Observation is needed because when the parameter node is changed then the GUI must be updated immediately.
    """

    if inputParameterNode:
      self.logic.setDefaultParameters(inputParameterNode)

    # Unobserve previously selected parameter node and add an observer to the newly selected.
    # Changes of parameter node are observed so that whenever parameters are changed by a script or any other module
    # those are reflected immediately in the GUI.
    if self._parameterNode is not None:
      self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)
    self._parameterNode = inputParameterNode
    if self._parameterNode is not None:
      self.addObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

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

    # Update node selectors and sliders
    self.ui.outputSelector.setCurrentNode(self._parameterNode.GetNodeReference("OutputTransform"))
    self.ui.arcAngleSliderWidget.value = float(self._parameterNode.GetParameter("ArcAngle"))
    self.ui.ringAngleSliderWidget.value = float(self._parameterNode.GetParameter("RingAngle"))
    self.ui.headringCoordinatesWidget.coordinates = self._parameterNode.GetParameter("HeadringCoordinates")
    self.ui.mountingComboBox.setCurrentText(self._parameterNode.GetParameter("Mounting"))
    self.ui.autoUpdateCheckBox.setChecked(int(self._parameterNode.GetParameter("AutoUpdate")))

    # Update buttons states and tooltips
    self.ui.applyButton.enabled = bool(self._parameterNode.GetNodeReference("OutputTransform")) and not int(self._parameterNode.GetParameter("AutoUpdate"))
    self.ui.setDefaultResliceDriverPushButton.enabled = bool(self._parameterNode.GetNodeReference("OutputTransform"))

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

    self._parameterNode.SetNodeReferenceID("OutputTransform", self.ui.outputSelector.currentNodeID)
    self._parameterNode.SetParameter("ArcAngle", str(self.ui.arcAngleSliderWidget.value))
    self._parameterNode.SetParameter("RingAngle", str(self.ui.ringAngleSliderWidget.value))
    self._parameterNode.SetParameter("HeadringCoordinates", self.ui.headringCoordinatesWidget.coordinates)
    self._parameterNode.SetParameter("Mounting", self.ui.mountingComboBox.currentText)
    self._parameterNode.SetParameter("AutoUpdate", str(int(self.ui.autoUpdateCheckBox.checked)))

    self._parameterNode.EndModify(wasModified)

    # Run in case of auto update
    if int(self._parameterNode.GetParameter("AutoUpdate")):
      self.onApplyButton()

  def onImportFromPDFPushButton(self):
    # select PDF
    filePath = qt.QFileDialog.getOpenFileName(qt.QWidget(), 'Select Planning PDF', '', '*.pdf')
    if filePath == '':
      return
    # get planning
    planningDictionary = StereotacticPlanLib.util.exctractPlanningFromPDF(filePath)
    if not planningDictionary:
      return
    # set values
    wasModified = self._parameterNode.StartModify()  # Modify all properties in a single batch
    self._parameterNode.SetParameter("ArcAngle", planningDictionary["Arc Angle"])
    self._parameterNode.SetParameter("RingAngle", planningDictionary["Ring Angle"])
    self._parameterNode.SetParameter("HeadringCoordinates", planningDictionary["Headring Coordinates"])
    self._parameterNode.SetParameter("Mounting", planningDictionary["Mounting"])
    self._parameterNode.EndModify(wasModified)

  def onSetDefaultResliceDriver(self):
    StereotacticPlanLib.util.setDefaultResliceDriver(self._parameterNode.GetNodeReferenceID("OutputTransform"))
      
  def onPreviewLineToggled(self, state):
    if state:
      # add node and default points
      markupsLineNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsLineNode')
      self._parameterNode.SetNodeReferenceID("PreviewLine", markupsLineNode.GetID())
      markupsLineNode.AddControlPointWorld(vtk.vtkVector3d(0,0,0), 'Target')
      markupsLineNode.AddControlPointWorld(vtk.vtkVector3d(0,0,80), '2')
      self.updatePreviewLineTransform(self._parameterNode.GetNodeReference("OutputTransform"))
    elif self._parameterNode.GetNodeReferenceID("PreviewLine"):
      # remove node
      slicer.mrmlScene.RemoveNode(self._parameterNode.GetNodeReference("PreviewLine"))

  def updatePreviewLineTransform(self, node):
    if self._parameterNode.GetNodeReferenceID("PreviewLine"):
      self._parameterNode.GetNodeReference("PreviewLine").SetAndObserveTransformNodeID(node.GetID() if node else None)


  def onApplyButton(self):
    """
    Run processing when user clicks "Apply" button.
    """
    try:

      # Compute output
      self.logic.process(self.ui.outputSelector.currentNode(),
                         self.ui.arcAngleSliderWidget.value,
                         self.ui.ringAngleSliderWidget.value,
                         np.array(self.ui.headringCoordinatesWidget.coordinates.split(','), dtype=float),
                         self.ui.mountingComboBox.currentText)

    except Exception as e:
      slicer.util.errorDisplay("Failed to compute transform: "+str(e))
      import traceback
      traceback.print_exc()


#
# StereotacticPlanLogic
#

class StereotacticPlanLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self):
    """
    Called when the logic class is instantiated. Can be used for initializing member variables.
    """
    ScriptedLoadableModuleLogic.__init__(self)

  def setDefaultParameters(self, parameterNode):
    """
    Initialize parameter node with default settings.
    """
    if not parameterNode.GetParameter("ArcAngle"):
      parameterNode.SetParameter("ArcAngle", "90.0")
    if not parameterNode.GetParameter("RingAngle"):
      parameterNode.SetParameter("RingAngle", "90.0")
    if not parameterNode.GetParameter("HeadringCoordinates"):
      parameterNode.SetParameter("HeadringCoordinates", "100.0,100.0,100.0")
    if not parameterNode.GetParameter("Mounting"):
      parameterNode.SetParameter("Mounting", "lateral-left")
    if not parameterNode.GetParameter("AutoUpdate"):
      parameterNode.SetParameter("AutoUpdate", "0")

  def process(self, outputTransform, arcAngle, ringAngle, headringCoordinates, mounting):
    """
    Run the processing algorithm.
    Can be used without GUI widget.
    """

    if not outputTransform:
      raise ValueError("output transform is invalid")

    # Headring coordinates to Slicer world (matching center)
    headringToRAS = np.array([[ -1,  0,  0,  100],
                              [  0,  1,  0, -100],
                              [  0,  0, -1,  100],
                              [  0,  0,  0,    1]])
    
    headringCoordinatesRAS = np.dot(headringToRAS, np.append(headringCoordinates, 1))[:3]

    # Get ring and arc directions
    if mounting == 'lateral-right':
      initDirection = [0, 1, 0]
      ringDirection = [1, 0, 0]
      arcDirection =  [0, -np.sin(np.deg2rad(ringAngle)), np.cos(np.deg2rad(ringAngle))]
    elif mounting == 'lateral-left':
      initDirection = [0, -1, 0]
      ringDirection = [-1, 0, 0]
      arcDirection  = [0, np.sin(np.deg2rad(ringAngle)), np.cos(np.deg2rad(ringAngle))]
    elif mounting == 'sagittal-anterior':
      initDirection = [-1, 0, 0]
      ringDirection = [0, 1, 0]
      arcDirection  = [np.sin(np.deg2rad(ringAngle)), 0, np.cos(np.deg2rad(ringAngle))]
    elif mounting == 'sagittal-posterior':
      initDirection = [1, 0, 0]
      ringDirection = [0, -1, 0]
      arcDirection  = [-np.sin(np.deg2rad(ringAngle)), 0, np.cos(np.deg2rad(ringAngle))]

    # Create vtk Transform
    vtkTransform = vtk.vtkTransform()
    vtkTransform.Translate(headringCoordinatesRAS)
    vtkTransform.RotateWXYZ(arcAngle, arcDirection[0], arcDirection[1], arcDirection[2])
    vtkTransform.RotateWXYZ(ringAngle, ringDirection[0], ringDirection[1], ringDirection[2])
    vtkTransform.RotateWXYZ(90, initDirection[0], initDirection[1], initDirection[2])

    # Set to node
    outputTransform.SetAndObserveTransformToParent(vtkTransform)

#
# StereotacticPlanTest
#

class StereotacticPlanTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear()

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_StereotacticPlan1()

  def test_StereotacticPlan1(self):
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

    # TODO

    # Get/create input data

    # import SampleData
    # registerSampleData()
    # inputVolume = SampleData.downloadSample('StereotacticPlan1')
    # self.delayDisplay('Loaded test data set')

    # inputScalarRange = inputVolume.GetImageData().GetScalarRange()
    # self.assertEqual(inputScalarRange[0], 0)
    # self.assertEqual(inputScalarRange[1], 695)

    # OutputTransform = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
    # threshold = 100

    # # Test the module logic

    # logic = StereotacticPlanLogic()

    # # Test algorithm with non-inverted threshold
    # logic.process(inputVolume, OutputTransform, threshold, True)
    # outputScalarRange = OutputTransform.GetImageData().GetScalarRange()
    # self.assertEqual(outputScalarRange[0], inputScalarRange[0])
    # self.assertEqual(outputScalarRange[1], threshold)

    # # Test algorithm with inverted threshold
    # logic.process(inputVolume, OutputTransform, threshold, False)
    # outputScalarRange = OutputTransform.GetImageData().GetScalarRange()
    # self.assertEqual(outputScalarRange[0], inputScalarRange[0])
    # self.assertEqual(outputScalarRange[1], inputScalarRange[1])

    self.delayDisplay('Test passed')
