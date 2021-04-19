import os
import unittest
import logging
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin

import LeadORLib
import numpy as np
import importlib
import glob
import re

#
# LeadOR
#

class LeadOR(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "LeadOR"  # TODO: make this more human readable by adding spaces
    self.parent.categories = ["Netstim"]  # TODO: set categories (folders where the module shows up in the module selector)
    self.parent.dependencies = []  # TODO: add here list of module names that this module requires
    self.parent.contributors = ["John Doe (AnyWare Corp.)"]  # TODO: replace with "Firstname Lastname (Organization)"
    # TODO: update with short description of the module and a link to online module documentation
    self.parent.helpText = """
This is an example of scripted loadable module bundled in an extension.
See more information in <a href="https://github.com/organization/projectname#LeadOR">module documentation</a>.
"""
    # TODO: replace with organization, grant and thanks
    self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc., Andras Lasso, PerkLab,
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
"""

    # Additional initialization step after application startup is complete
    slicer.app.connect("startupCompleted()", registerSampleData)

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

  # LeadOR1
  SampleData.SampleDataLogic.registerCustomSampleDataSource(
    # Category and sample name displayed in Sample Data module
    category='LeadOR',
    sampleName='LeadOR1',
    # Thumbnail should have size of approximately 260x280 pixels and stored in Resources/Icons folder.
    # It can be created by Screen Capture module, "Capture all views" option enabled, "Number of images" set to "Single".
    thumbnailFileName=os.path.join(iconsPath, 'LeadOR1.png'),
    # Download URL and target file name
    uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95",
    fileNames='LeadOR1.nrrd',
    # Checksum to ensure file integrity. Can be computed by this command:
    #  import hashlib; print(hashlib.sha256(open(filename, "rb").read()).hexdigest())
    checksums = 'SHA256:998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95',
    # This node name will be used when the data set is loaded
    nodeNames='LeadOR1'
  )

  # LeadOR2
  SampleData.SampleDataLogic.registerCustomSampleDataSource(
    # Category and sample name displayed in Sample Data module
    category='LeadOR',
    sampleName='LeadOR2',
    thumbnailFileName=os.path.join(iconsPath, 'LeadOR2.png'),
    # Download URL and target file name
    uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97",
    fileNames='LeadOR2.nrrd',
    checksums = 'SHA256:1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97',
    # This node name will be used when the data set is loaded
    nodeNames='LeadOR2'
  )

#
# LeadORWidget
#

class LeadORWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
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
    uiWidget = slicer.util.loadUI(self.resourcePath('UI/LeadOR.ui'))
    self.layout.addWidget(uiWidget)
    self.ui = slicer.util.childWidgetVariables(uiWidget)

    # AO Channels actions to ToolButton
    for child in self.ui.microElectrodeLayoutFrame.children():
      if isinstance(child, qt.QToolButton):
        # channel
        AOChannelsAction = qt.QAction('AO Channels', self.layout)
        AOChannelsAction.setEnabled(False)
        AOChannelsAction.setMenu(qt.QMenu(child))
        child.addAction(AOChannelsAction)


    # Set scene in MRML widgets. Make sure that in Qt designer the top-level qMRMLWidget's
    # "mrmlSceneChanged(vtkMRMLScene*)" signal in is connected to each MRML widget's.
    # "setMRMLScene(vtkMRMLScene*)" slot.
    uiWidget.setMRMLScene(slicer.mrmlScene)

    # Create logic class. Logic implements all computations that should be possible to run
    # in batch mode, without a graphical user interface.
    self.logic = LeadORLogic()

    # Connections

    # These connections ensure that we update parameter node when scene is closed
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)


    # These connections ensure that whenever user changes some settings on the GUI, that is saved in the MRML scene
    # (in the selected parameter node).
    self.ui.distanceToTargetComboBox.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    self.ui.trajectoryTransformComboBox.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    self.ui.startPushButton.connect("clicked(bool)", self.updateParameterNodeFromGUI)

    # transforms hierarchy
    self.ui.distanceToTargetComboBox.connect("currentNodeChanged(vtkMRMLNode*)", self.setTransformsHierarchy)
    self.ui.trajectoryTransformComboBox.connect("currentNodeChanged(vtkMRMLNode*)", self.setTransformsHierarchy)

    # reset ME when distance to target node change
    self.ui.distanceToTargetComboBox.currentNodeChanged.connect(lambda b,enabledList=[0]*9: self.setMELayout(enabledList))

    # distance to target slider
    self.ui.distanceToTargetComboBox.connect("currentNodeChanged(vtkMRMLNode*)", lambda node: self.ui.distanceToTargetSlider.setMRMLTransformNode(node))
    self.ui.distanceToTargetSlider.connect("valueChanged(double)", lambda value: self.ui.distanceToTargetSlider.applyTransformation(value))

    # micro electrode layouts
    self.ui.MECenterLayoutPushButton.clicked.connect(lambda b,enabledList=[0,0,0,0,1,0,0,0,0]: self.setMELayout(enabledList))
    self.ui.MEPlusLayoutPushButton.clicked.connect(  lambda b,enabledList=[0,1,0,1,1,1,0,1,0]: self.setMELayout(enabledList))
    self.ui.MEXLayoutPushButton.clicked.connect(     lambda b,enabledList=[1,0,1,0,1,0,1,0,1]: self.setMELayout(enabledList))

    # add connection for each micro electro toggle button 
    for child in self.ui.microElectrodeLayoutFrame.children():
      if isinstance(child, qt.QToolButton):
        child.toggled.connect(lambda b,N=int(child.objectName.split('_')[-1]): self.microElectrodeLayoutToggle(b,N))

    # ME visibility
    self.ui.MEModelVisCheckBox.toggled.connect(lambda b: LeadORLib.util.setMEVisibility('Model', b))
    self.ui.MELineVisCheckBox.toggled.connect(lambda b: LeadORLib.util.setMEVisibility('Line', b))
    self.ui.METipVisCheckBox.toggled.connect(lambda b: LeadORLib.util.setMEVisibility('Tip', b))
    self.ui.METraceVisCheckBox.toggled.connect(lambda b: LeadORLib.util.setMEVisibility('Trace',b))


    # start
    self.ui.startPushButton.clicked.connect(self.onStartPushButton)

    # Make sure parameter node is initialized (needed for module reload)
    self.initializeParameterNode()

    # AlphaOmega parameter node (Singleton)
    try:
      AOParameterNode = slicer.modules.alphaomega.logic().getParameterNode()
      self.addObserver(AOParameterNode, vtk.vtkCommand.ModifiedEvent, self.updateParameterNodeFromAO)
      self.updateParameterNodeFromAO(AOParameterNode, None)
    except:
      pass

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
    # Do not react to parameter node changes (GUI wlil be updated when the user enters into the module)
    self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

  def onSceneStartClose(self, caller, event):
    """
    Called just before the scene is closed.
    """
    # Parameter node will be reset, do not use it anymore
    self.setParameterNode(None)
    # reset ME state
    self.setMELayout([0]*9)

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
    # if not self._parameterNode.GetNodeReference("DistanceToTargetTransform"):
    #   distanceToTargetTransform = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLLinearTransformNode")
    #   if distanceToTargetTransform:
    #     self._parameterNode.SetNodeReferenceID("DistanceToTargetTransform", distanceToTargetTransform.GetID())
    # if not self._parameterNode.GetNodeReference("TrajectoryTransform"):
    #   self._parameterNode.SetNodeReferenceID("TrajectoryTransform", slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLinearTransformNode').GetID())
      

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

    self.ui.distanceToTargetComboBox.setCurrentNode(self._parameterNode.GetNodeReference("DistanceToTargetTransform"))
    self.ui.trajectoryTransformComboBox.setCurrentNode(self._parameterNode.GetNodeReference("TrajectoryTransform"))

    # self.ui.leftMELabel.text = 'Lateral' if self._parameterNode.GetParameter('Side') == 'Left' else 'Medial'
    # self.ui.rightMELabel.text = 'Medial' if self._parameterNode.GetParameter('Side') == 'Left' else 'Lateral'


    # disable ME and trajectory modifications while running
    transformsAvailable = bool(self._parameterNode.GetNodeReference("DistanceToTargetTransform") and self._parameterNode.GetNodeReference("TrajectoryTransform"))
    running = bool(int(self._parameterNode.GetParameter("running")))

    self.ui.layoutToggleFrame.enabled            = not running and transformsAvailable
    self.ui.microElectrodeLayoutFrame.enabled    = not running and transformsAvailable
    self.ui.distanceToTargetComboBox.enabled     = not running
    self.ui.trajectoryTransformComboBox.enabled  = not running
    self.ui.AOSignalsCollapsibleButton.enabled   = transformsAvailable

    # start pb text
    self.ui.startPushButton.text = 'Stop Live Update' if running else 'Start Live Update'

    # enable start
    self.ui.startPushButton.enabled = self._parameterNode.GetParameter("AORunning") == "True"

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

    self._parameterNode.SetNodeReferenceID("DistanceToTargetTransform", self.ui.distanceToTargetComboBox.currentNodeID)
    self._parameterNode.SetNodeReferenceID("TrajectoryTransform", self.ui.trajectoryTransformComboBox.currentNodeID)
    self._parameterNode.SetParameter("running", str(int(self.ui.startPushButton.checked)))

    self._parameterNode.EndModify(wasModified)

  def updateParameterNodeFromAO(self, caller, event):
    """
    This method is called when there are changes in the AlphaOmega module
    The parameters are copied
    """

    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return

    wasModified = self._parameterNode.StartModify()  # Modify all properties in a single batch

    self._parameterNode.SetNodeReferenceID("DistanceToTargetTransform", caller.GetNodeReferenceID("DistanceToTargetTransform"))
    self._parameterNode.SetParameter("AORunning", caller.GetParameter("Running"))
    self._parameterNode.SetParameter("AOSignalSavePath", caller.GetParameter("SignalSavePath"))
    self._parameterNode.SetParameter("AOSamplingFrequency", caller.GetParameter("SamplingFrequency"))
    # channels names
    prevChannelsNames = self._parameterNode.GetParameter("AOChannelsNames")
    self._parameterNode.SetParameter("AOChannelsNames", caller.GetParameter("ChannelsNames"))
    if prevChannelsNames != self._parameterNode.GetParameter("AOChannelsNames"):
      self.updateChannelsNames()
    
    self._parameterNode.EndModify(wasModified)

  def updateChannelsNames(self):
    channelsNames = self._parameterNode.GetParameter("AOChannelsNames").split(",")
    for child in self.ui.microElectrodeLayoutFrame.children():
      if isinstance(child, qt.QToolButton):
        # init
        AOChannelsActionGroup = qt.QActionGroup(child)
        # add none
        noneChannelAction = qt.QAction('None', child)
        noneChannelAction.setCheckable(True)
        noneChannelAction.setChecked(True)
        AOChannelsActionGroup.addAction(noneChannelAction)
        # add for each channel
        for name in channelsNames:
          channelAction = qt.QAction(name, child)
          channelAction.setCheckable(True)
          AOChannelsActionGroup.addAction(channelAction)
        # set menu
        child.actions()[0].menu().clear()
        child.actions()[0].menu().addActions(AOChannelsActionGroup.actions())
        child.actions()[0].menu().triggered(noneChannelAction)

  def setTransformsHierarchy(self):
    if self._parameterNode and self._parameterNode.GetNodeReference("DistanceToTargetTransform") and self._parameterNode.GetNodeReference("TrajectoryTransform"):
      self._parameterNode.GetNodeReference("DistanceToTargetTransform").SetAndObserveTransformNodeID(self._parameterNode.GetNodeReferenceID("TrajectoryTransform"))


  def microElectrodeLayoutToggle(self, enabled, N):
    # enable/disable tool button actions
    toolButton = getattr(self.ui, 'METoolButton_'+str(N))
    list(map(lambda a: a.setEnabled(enabled), toolButton.actions()))
    # create/remove micro electrode nodes
    if enabled:
      LeadORLib.util.initNthMicroElectrode(N, self._parameterNode.GetNodeReferenceID("DistanceToTargetTransform"), toolButton)
    else:
      LeadORLib.util.removeNthMicroElectrode(N)

  def setMELayout(self, enabledList):
    for enabled, N in zip(enabledList, range(len(enabledList))):  
      getattr(self.ui, 'METoolButton_'+str(N)).checked = enabled

  def onStartPushButton(self):
    # get nodes
    samplingFrequencyStr = self._parameterNode.GetParameter("AOSamplingFrequency")
    signalSavePath = self._parameterNode.GetParameter("AOSignalSavePath")
    channelsNames = self._parameterNode.GetParameter("AOChannelsNames").split(',')
    distanceToTargetNode = self._parameterNode.GetNodeReference("DistanceToTargetTransform")
    METraceFiducialIDs = LeadORLib.util.GetMETraceFiducials()

    if len(METraceFiducialIDs) == 0:
      return

    # start/stop
    if self.ui.startPushButton.checked:
      self.logic.run(distanceToTargetNode, METraceFiducialIDs, channelsNames,  samplingFrequencyStr, signalSavePath)
    else:
      self.logic.stop(distanceToTargetNode, METraceFiducialIDs)



#
# LeadORLogic
#

class LeadORLogic(ScriptedLoadableModuleLogic):
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
    if slicer.util.settingsValue('Developer/DeveloperMode', False, converter=slicer.util.toBool):
      import LeadORLib
      import LeadORLib.util
      importlib.reload(LeadORLib.util)

    
  def setDefaultParameters(self, parameterNode):
    """
    Initialize parameter node with default settings.
    """

    if not parameterNode.GetParameter("Side"):
      parameterNode.SetParameter("Side", "Left")
    if not parameterNode.GetParameter("running"):
      parameterNode.SetParameter("running", "0")

  def getVTARadius(self, I, pw=60): 
    # I: amplitude in Ampere
    # pw: pulse width in micro seconds
    # returns radius in meter
    return ((pw/90)**0.3) * np.sqrt(0.72*I/165)

  def run(self, distanceToTargetNode, METraceFiducialIDs, channelsNames, samplingFrequencyStr, signalSavePath):

    # init matlab
    slicer.cli.run(slicer.modules.matlabcommander, None, {'cmd': "cd('"+os.path.join(os.path.split(__file__)[0],'LeadORLib')+"')"}, update_display=False)
    slicer.cli.run(slicer.modules.matlabcommander, None, {'cmd': "warning('off')"}, update_display=False)

    # observe distance to target transform modified
    obs = distanceToTargetNode.AddObserver(slicer.vtkMRMLTransformNode.TransformModifiedEvent, lambda c,e: self.onDistanceToTargetModified(c, METraceFiducialIDs, channelsNames, samplingFrequencyStr, signalSavePath))

    # save observer in distance to target transform node
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    shNode.SetItemAttribute(shNode.GetItemByDataNode(distanceToTargetNode), 'Observer', str(obs))

    # set fiducial description to active
    for ID in METraceFiducialIDs:
      fiducialNode = slicer.util.getNode(ID)
      channelName = shNode.GetItemAttribute(shNode.GetItemByDataNode(fiducialNode), 'AO Channel')
      fiducialNode.SetDescription('Active' if channelName in channelsNames else '')

  def onDistanceToTargetModified(self, distanceToTargetNode, METraceFiducialIDs, channelsNames, samplingFrequencyStr, signalSavePath):
    # update previous dtt
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    prevDTT = shNode.GetItemAttribute(shNode.GetItemByDataNode(distanceToTargetNode), 'PrevDTT')
    shNode.SetItemAttribute(shNode.GetItemByDataNode(distanceToTargetNode), 'PrevDTT', str('{:.3f}'.format(distanceToTargetNode.GetMatrixTransformToParent().GetElement(2,3))))
    if prevDTT  == '': # first run
      return
    # look for file path
    g = glob.glob(os.path.join(signalSavePath, 'AO*'+prevDTT+'*.csv'))
    g.sort(key=os.path.getmtime)
    if g == []:
      return # TODO: print warning
    dttSignalFilePath = g[-1]
    # run cli node
    cliNode = slicer.cli.run(slicer.modules.matlabcommander, None, {'cmd': "getMFRwave_clus('"+dttSignalFilePath+"',"+samplingFrequencyStr+")"}, update_display=False)
    # for each ME trace fiducial, get index that matches channel name and add observer to cli node
    for ID in METraceFiducialIDs:
      fiducialNode = slicer.util.getNode(ID)
      channelName = shNode.GetItemAttribute(shNode.GetItemByDataNode(fiducialNode), 'AO Channel')
      if fiducialNode.GetDescription() == 'Active':
        index = channelsNames.index(channelName)
        cliNode.AddObserver(slicer.vtkMRMLCommandLineModuleNode.StatusModifiedEvent, lambda c,e,n=fiducialNode,i=index: self.updateTrace(c, n, i, prevDTT))


  def updateTrace(self, cliNode, fiducialNode, index, dtt): 
    # update corresponding fiducial description with output from the cli node   
    if cliNode.GetStatusString() == 'Completed':
      reply = cliNode.GetParameterAsString('reply')
      matlabAns = re.search('(?<=ans [=]).*',reply).group(0)
      matlabAnsStr = re.split('\s+',matlabAns)[1:] # empty string is returned
      for i in range(fiducialNode.GetNumberOfControlPoints()-1,-1,-1):
        if fiducialNode.GetNthControlPointLabel(i) == ('D = ' + dtt):
          fiducialNode.SetNthControlPointDescription(i, matlabAnsStr[index])
          break

  def onCLIStatusModified(self, caller, event):
    # run cli node after it finishes or is cancelled
    if caller.GetStatusString() in ['Completed', 'Cancelled']:
      qt.QTimer.singleShot(1000, lambda: slicer.cli.run(slicer.modules.merspike, caller, update_display=False))
      
  def stop(self, distanceToTarget, METraceFiducialIDs):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    ID = shNode.GetItemByDataNode(distanceToTarget)
    # remove obs from transform
    distanceToTarget.RemoveObserver(int(shNode.GetItemAttribute(ID, 'Observer')))
    # stop and remove observers from cli node
    cliNode = slicer.util.getNode(shNode.GetItemAttribute(ID, 'CLI Node'))
    cliNode.RemoveAllObservers()
    cliNode.Cancel()
    # set fiducial description
    list(map(lambda ID: slicer.util.getNode(ID).SetDescription(''), METraceFiducialIDs))


#
# LeadORTest
#

class LeadORTest(ScriptedLoadableModuleTest):
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
    self.test_LeadOR1()

  def test_LeadOR1(self):
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
    registerSampleData()
    inputVolume = SampleData.downloadSample('LeadOR1')
    self.delayDisplay('Loaded test data set')

    inputScalarRange = inputVolume.GetImageData().GetScalarRange()
    self.assertEqual(inputScalarRange[0], 0)
    self.assertEqual(inputScalarRange[1], 695)

    outputVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
    threshold = 100

    # Test the module logic

    logic = LeadORLogic()

    # Test algorithm with non-inverted threshold
    logic.process(inputVolume, outputVolume, threshold, True)
    outputScalarRange = outputVolume.GetImageData().GetScalarRange()
    self.assertEqual(outputScalarRange[0], inputScalarRange[0])
    self.assertEqual(outputScalarRange[1], threshold)

    # Test algorithm with inverted threshold
    logic.process(inputVolume, outputVolume, threshold, False)
    outputScalarRange = outputVolume.GetImageData().GetScalarRange()
    self.assertEqual(outputScalarRange[0], inputScalarRange[0])
    self.assertEqual(outputScalarRange[1], inputScalarRange[1])

    self.delayDisplay('Test passed')
