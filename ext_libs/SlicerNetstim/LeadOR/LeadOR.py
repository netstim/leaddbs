import os
import unittest
import logging
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin

import LeadORLib

#
# LeadOR
#

class LeadOR(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Lead-OR"
    self.parent.categories = ["Netstim"]
    self.parent.dependencies = [] 
    self.parent.contributors = ["Simon Oxenford (Netstim Berlin)"]
    self.parent.helpText = """
This module controls micro electrode settings for deep brain stimulation surgery
"""
    self.parent.acknowledgementText = ""


    # Additional initialization step after application startup is complete
    slicer.app.connect("startupCompleted()", registerSampleData)
    slicer.app.connect("startupCompleted()", addCustomLayout)

#
# Register sample data sets in Sample Data module
#

def registerSampleData():
  """
  Add data sets to Sample Data module.
  """

  import SampleData
  iconsPath = os.path.join(os.path.dirname(__file__), 'Resources/Icons')

  SampleData.SampleDataLogic.registerCustomSampleDataSource(
    category='LeadOR',
    sampleName='STN Planning',
    thumbnailFileName=os.path.join(iconsPath, 'LeadOR1.png'),
    uris="https://github.com/netstim/SlicerNetstim/releases/download/SampleData/Lead-OR_STN.mrb",
    fileNames='Lead-OR_STN.mrb',
    loadFiles=True,
    loadFileType='SceneFile'
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

    if hasattr(slicer,'vtkMRMLFiberBundleNode') and hasattr(slicer.vtkMRMLFiberBundleNode,'GetExtractFromROI'):
      self.ui.stimulationCollapsibleButton.enabled = True
    else:
      self.ui.stimulationCollapsibleButton.enabled = False
      self.ui.stimulationCollapsibleButton.collapsed = True
      self.ui.stimulationCollapsibleButton.setToolTip('Updated SlicerDMRI Extension needed for stimulation module')
      
    # AO Channels actions to ToolButton
    stimulationActionGroup = qt.QActionGroup(self.layout)
    for child in self.ui.microElectrodeLayoutFrame.children():
      if isinstance(child, qt.QToolButton):
        # channel
        AOChannelsAction = qt.QAction('AO Channels', self.layout)
        AOChannelsAction.setEnabled(False)
        AOChannelsAction.setMenu(qt.QMenu(child))
        child.addAction(AOChannelsAction)
        # stim
        stimulationAction = qt.QAction('Stim Source', self.layout)
        stimulationAction.setCheckable(True)
        stimulationAction.setEnabled(True)
        stimulationAction.toggled.connect(self.updateStimulationTransform)
        stimulationActionGroup.addAction(stimulationAction)
        child.addAction(stimulationAction)

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
    self.ui.trajectoryTransformComboBox.connect("currentNodeChanged(vtkMRMLNode*)", self.guessSideFromTransform)

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
    self.ui.MEModelVisCheckBox.toggled.connect(lambda b: self.logic.setMEVisibility('microElectrodeModel', b))
    self.ui.MELineVisCheckBox.toggled.connect(lambda b: self.logic.setMEVisibility('trajectoryLine', b))
    self.ui.METipVisCheckBox.toggled.connect(lambda b: self.logic.setMEVisibility('tipFiducial', b))
    self.ui.METraceVisCheckBox.toggled.connect(lambda b: self.logic.setMEVisibility('traceModel', b))
    self.ui.MEClustersCheckBox.toggled.connect(lambda b: self.logic.setMEVisibility('clusterFiducials', b))

    # Stimulation
    self.ui.stimulationActiveCheckBox.connect('toggled(bool)', self.onStimulationActivate)
    self.ui.stimulationAmplitudeSpinBox.valueChanged.connect(self.updateStimulationRadius)

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

    transformsAvailable = bool(self._parameterNode.GetNodeReference("DistanceToTargetTransform") and self._parameterNode.GetNodeReference("TrajectoryTransform"))

    self.ui.layoutToggleFrame.enabled            = transformsAvailable
    self.ui.microElectrodeLayoutFrame.enabled    = transformsAvailable
    self.ui.stimulationCollapsibleButton.enabled = transformsAvailable and hasattr(slicer,'vtkMRMLFiberBundleNode') and hasattr(slicer.vtkMRMLFiberBundleNode,'GetExtractFromROI')

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
    self._parameterNode.EndModify(wasModified)


  def setTransformsHierarchy(self):
    if self._parameterNode and self._parameterNode.GetNodeReference("DistanceToTargetTransform") and self._parameterNode.GetNodeReference("TrajectoryTransform"):
      self._parameterNode.GetNodeReference("DistanceToTargetTransform").SetAndObserveTransformNodeID(self._parameterNode.GetNodeReferenceID("TrajectoryTransform"))


  def microElectrodeLayoutToggle(self, enabled, N):
    toolButton = getattr(self.ui, 'METoolButton_'+str(N))
    if enabled:
      self.logic.initializeNthTrajectory(N, self._parameterNode.GetNodeReferenceID("DistanceToTargetTransform"))
      self.setToolButtonMenu(toolButton, N)
    else:
      self.logic.removeNthTrajectory(N)
      toolButton.actions()[0].menu().clear()

  def setToolButtonMenu(self, toolButton, N):
    channelsNames = []
    for i in range(slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLAlphaOmegaChannelNode')):
      channelsNames.append(slicer.mrmlScene.GetNthNodeByClass(i,'vtkMRMLAlphaOmegaChannelNode').GetChannelName())
    if not channelsNames:
      return
    # init
    AOChannelsActionGroup = qt.QActionGroup(toolButton)
    # add none
    noneChannelAction = qt.QAction('None', toolButton)
    noneChannelAction.setCheckable(True)
    noneChannelAction.setChecked(True)
    AOChannelsActionGroup.addAction(noneChannelAction)
    # add for each channel
    for name in channelsNames:
      channelAction = qt.QAction(name, toolButton)
      channelAction.setCheckable(True)
      AOChannelsActionGroup.addAction(channelAction)
    # set menu
    toolButton.actions()[0].menu().clear()
    toolButton.actions()[0].menu().addActions(AOChannelsActionGroup.actions())
    toolButton.actions()[0].menu().triggered.connect(lambda action, trajectoryN=N: self.logic.setNthTrajectoryAOChannelNode(trajectoryN, action.text))

  def setMELayout(self, enabledList):
    for enabled, N in zip(enabledList, range(len(enabledList))):  
      getattr(self.ui, 'METoolButton_'+str(N)).checked = enabled

  def guessSideFromTransform(self, transformNode):
    if not transformNode:
      return
    currentPoint = [0.0] * 4
    matrix = vtk.vtkMatrix4x4()
    transformNode.GetMatrixTransformToWorld(matrix)
    matrix.MultiplyPoint([0.0, 0.0, 0.0, 1.0], currentPoint)
    guessRightSide = currentPoint[0] > 0
    self.ui.leftMELabel.text = 'Medial' if guessRightSide else 'Lateral'  
    self.ui.rightMELabel.text = 'Lateral' if guessRightSide else 'Medial' 

  def onStimulationActivate(self, active):
    if active: 
      self.logic.VTASource = LeadORLib.util.VTASource()
      self.updateStimulationRadius(self.ui.stimulationAmplitudeSpinBox.value)
      self.updateStimulationTransform()
    else:
      self.logic.VTASource.cleanup()
      self.logic.VTASource = None
      self.ui.amplitudeRadiusLabel.setText('-')
        
  def updateStimulationTransform(self):
    if not self.logic.VTASource:
      return
    # get current active stim
    N = next(filter(lambda n: getattr(self.ui, 'METoolButton_'+str(n)).actions()[1].checked, range(9)), None)
    if N is None or N not in self.logic.trajectories.keys():
      self.ui.stimulationActiveCheckBox.checked = False
      qt.QMessageBox().warning(qt.QWidget(),'','Set a Stimulation Source')
      return
    # set transform
    self.logic.VTASource.SetAndObserveTransformNodeID(self.logic.trajectories[N].translationTransform.GetID())

  def updateStimulationRadius(self, value):
    if not self.logic.VTASource:
      return
    # set  radius
    radius = self.logic.getVTARadius(value * 1e-3) * 1e3 
    self.logic.VTASource.SetRadius(radius)
    self.ui.amplitudeRadiusLabel.setText('%.2f' % radius)


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
      import importlib
      importlib.reload(LeadORLib.util)

    self.trajectories = {}
    self.VTASource = None

    
  def setDefaultParameters(self, parameterNode):
    """
    Initialize parameter node with default settings.
    """
    pass

  def initializeNthTrajectory(self, N, distanceToTargetTransformID):
    self.trajectories[N] = LeadORLib.util.Trajectory(N, distanceToTargetTransformID)
  
  def removeNthTrajectory(self, N):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    IDs = vtk.vtkIdList()
    shNode.GetItemChildren(self.trajectories[N].folderID, IDs, True)
    for i in range(IDs.GetNumberOfIds()):
      shNode.RemoveItem(IDs.GetId(i))
    shNode.RemoveItem(self.trajectories[N].folderID)
    del self.trajectories[N]

  def setNthTrajectoryAOChannelNode(self, N, channelName):
    if N not in self.trajectories.keys() or channelName == 'None':
      return
    channelNode = self.getAOChannelNodeFromChannelName(channelName)
    self.trajectories[N].setAlphaOmegaChannelNode(channelNode)
    
  def getAOChannelNodeFromChannelName(self, channelName):
    for i in range(slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLAlphaOmegaChannelNode')):
      AONode = slicer.mrmlScene.GetNthNodeByClass(i,'vtkMRMLAlphaOmegaChannelNode')
      if AONode.GetChannelName() == channelName:
        return AONode

  def setMEVisibility(self, modelType, visible):
    for trajectory in self.trajectories.values():
      getattr(trajectory, modelType).GetDisplayNode().SetVisibility3D(visible)

  def getVTARadius(self, I, pw=60): 
    # I: amplitude in Ampere
    # pw: pulse width in micro seconds
    # returns radius in meter
    from numpy import sqrt
    return ((pw/90)**0.3) * sqrt(0.8*I/165) # 0.72



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

#
# Custom Layout
#


def addCustomLayout():

  customLayout = """
  <layout type="horizontal" split="true">
  <item>
    <view class="vtkMRMLViewNode" singletontag="1">
    <property name="viewlabel" action="default">1</property>
    </view>
  </item>
  <item>
    <layout type="vertical">
    <item>
      <layout type="horizontal">
      <item>
        <view class="vtkMRMLSliceNode" singletontag="Red">
        <property name="orientation" action="default">Axial</property>
        <property name="viewlabel" action="default">R</property>
        <property name="viewcolor" action="default">#F34A33</property>
        </view>
      </item>
      <item>
        <layout type="vertical">
        <item>
          <view class="vtkMRMLPlotViewNode" singletontag="PlotView1">
          <property name="viewlabel" action="default">1</property>
          </view>
        </item>
        <item>
          <view class="vtkMRMLPlotViewNode" singletontag="PlotView2">
          <property name="viewlabel" action="default">2</property>
          </view>
        </item>
        <item>
          <view class="vtkMRMLPlotViewNode" singletontag="PlotView3">
          <property name="viewlabel" action="default">3</property>
          </view>
        </item>
        </layout>
      </item>
      </layout>
    </item>
    <item>
      <layout type="horizontal">
      <item>
        <view class="vtkMRMLSliceNode" singletontag="Green">
        <property name="orientation" action="default">Coronal</property>
        <property name="viewlabel" action="default">G</property>
        <property name="viewcolor" action="default">#6EB04B</property>
        </view>
      </item>
      <item>
        <view class="vtkMRMLSliceNode" singletontag="Yellow">
        <property name="orientation" action="default">Sagittal</property>
        <property name="viewlabel" action="default">Y</property>
        <property name="viewcolor" action="default">#EDD54C</property>
        </view>
      </item>
      </layout>
    </item>
    </layout>
  </item>
  </layout>
  """

  customLayoutId = 509

  try:
    mainWindow = slicer.util.mainWindow()
    layoutManager = slicer.app.layoutManager()
    layoutNode = layoutManager.layoutLogic().GetLayoutNode()
  except:
    mainWindow = None

  if mainWindow and not layoutNode.GetLayoutDescription(customLayoutId):
    layoutNode.AddLayoutDescription(customLayoutId, customLayout)                                         
    viewToolBar = mainWindow.findChild("QToolBar", "ViewToolBar")
    layoutMenu = viewToolBar.widgetForAction(viewToolBar.actions()[0]).menu()
    layoutSwitchActionParent = layoutMenu
    layoutSwitchAction = layoutSwitchActionParent.addAction("LeadOR")
    layoutSwitchAction.setData(customLayoutId)
    layoutSwitchAction.setIcon(qt.QIcon(":Icons/Go.png"))
