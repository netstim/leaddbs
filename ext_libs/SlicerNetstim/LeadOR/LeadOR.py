import os
import unittest
import logging
import warnings
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin
import json

import LeadORLib
import LeadORLib.util
from LeadORLib.util import Trajectory, Feature
from LeadORLib.Widgets.tables import FeaturesTable

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

    # IGTLink
    if hasattr(slicer.modules,'openigtlinkif'):
      try:
        w = slicer.modules.openigtlinkif.createNewWidgetRepresentation()
        connectorListFrame = next(filter(lambda c: c.name=="ConnectorListFrame", w.children()))
        # Get ListView and ButtonFrame
        connectorListView = next(filter(lambda c: c.name=="ConnectorListView", connectorListFrame.children()))
        connectorButtonFrame = next(filter(lambda c: c.name=="ConnectorButtonFrame", connectorListFrame.children()))
        # Add "active" checkbox to connector ButtonFrame
        connectorPropertyWidget = next(filter(lambda c: c.name=="ConnectorPropertyWidget", connectorListFrame.children()))
        connectorPropertyFrame = next(filter(lambda c: c.name=="ConnectorPropertyFrame", connectorPropertyWidget.children()))
        connectorStateCheckBox = next(filter(lambda c: c.name=="ConnectorStateCheckBox", connectorPropertyFrame.children()))
        connectorButtonFrame.layout().addWidget(connectorStateCheckBox)
        # Add custom button to go to module
        goToModuleButton = qt.QPushButton('Go To Module')
        goToModuleButton.clicked.connect(lambda: slicer.util.mainWindow().moduleSelector().selectModule('OpenIGTLinkIF'))
        connectorButtonFrame.layout().addWidget(goToModuleButton)
        # Add ListView and ButtonFrame to collapsible button
        self.ui.IGTLinkFrame.setLayout(qt.QVBoxLayout())
        self.ui.IGTLinkFrame.layout().addWidget(connectorListView)
        self.ui.IGTLinkFrame.layout().addWidget(connectorButtonFrame)
        # Adjust view
        for i in range(4):
          connectorListView.header().setSectionResizeMode(i, qt.QHeaderView.Stretch)
      except:
        self.ui.IGTLinkCollapsibleButton.enabled = False
        self.ui.IGTLinkCollapsibleButton.collapsed = True
        self.ui.IGTLinkCollapsibleButton.setToolTip('Unable to set up OpenIGTLinkIF. Use the module directly.')
    else:
      self.ui.IGTLinkCollapsibleButton.enabled = False
      self.ui.IGTLinkCollapsibleButton.collapsed = True
      self.ui.IGTLinkCollapsibleButton.setToolTip('OpenIGTLinkIF Module required')

    # Reslice driver
    volumeResliceDriverPixmap = qt.QPixmap(self.resourcePath('Icons/VolumeResliceDriver.png'))
    self.ui.setDefaultResliceDriverToolButton.setIcon(qt.QIcon(volumeResliceDriverPixmap))
    self.ui.setDefaultResliceDriverToolButton.setIconSize(volumeResliceDriverPixmap.rect().size())
    self.ui.setDefaultResliceDriverToolButton.clicked.connect(self.setDefaultResliceDriver)

    # Sequences
    recordSequenzePixmap = qt.QPixmap(self.resourcePath('Icons/VcrRecord16.png'))
    self.ui.recordSequenceToolButton.setIcon(qt.QIcon(recordSequenzePixmap))
    self.ui.recordSequenceToolButton.setIconSize(recordSequenzePixmap.rect().size())
    self.ui.recordSequenceToolButton.clicked.connect(self.setUpSequenzeRecording)

    # Stim actions to ToolButton
    stimulationActionGroup = qt.QActionGroup(self.layout)
    for child in self.ui.trajectoriesLayoutFrame.children():
      if isinstance(child, qt.QToolButton):
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
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.NodeAddedEvent, self.onNodeAdded)
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.NodeRemovedEvent, self.onNodeRemoved)

    # Planning and distance to target
    self.ui.planningTransformComboBox.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    self.ui.planningTransformComboBox.connect("currentNodeChanged(vtkMRMLNode*)", self.guessSideFromTransform)
    self.ui.planningTransformComboBox.connect("currentNodeChanged(vtkMRMLNode*)", self.setTransformsHierarchy)

    self.ui.distanceToTargetComboBox.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    self.ui.distanceToTargetComboBox.connect("currentNodeChanged(vtkMRMLNode*)", self.setTransformsHierarchy)
    self.ui.distanceToTargetComboBox.currentNodeChanged.connect(lambda b,enabledList=[0]*9: self.setTrajectoryLayout(enabledList))
    self.ui.distanceToTargetComboBox.connect("currentNodeChanged(vtkMRMLNode*)", lambda node: self.ui.distanceToTargetSlider.setMRMLTransformNode(node))
    
    self.ui.distanceToTargetSlider.connect("valueChanged(double)", lambda value: self.ui.distanceToTargetSlider.applyTransformation(value))

    # Trajectories
    for i in range(self.ui.trajectoryVisualizationComboBox.model().rowCount()):
      index = self.ui.trajectoryVisualizationComboBox.model().index(i,0)
      self.ui.trajectoryVisualizationComboBox.setCheckState(index, qt.Qt.Checked)
    self.ui.trajectoryVisualizationComboBox.checkedIndexesChanged.connect(self.trajectoryVisualizationChanged)

    for i in range(9):
      toolButton = getattr(self.ui, 'TrajectoryToolButton_'+str(i))
      toolButton.toggled.connect(lambda b,N=i: self.trajectoryLayoutToggle(b,N))

    self.ui.trajectoryPresetComboBox.currentTextChanged.connect(lambda t: self.setTrajectoryLayoutPreset(t))
    self.ui.unlinkedChannelsListWidget.itemSelectionChanged.connect(self.onUnlinkedChannelsSelectionChanged)
    self.ui.linkChannelsToTrajectoriesPushButton.clicked.connect(self.onLinkChannelsToTrajectoriesPushButton)

    # Features
    self.ui.featuresTableWidget = FeaturesTable(self.ui.featuresTableView, self.updateParameterNodeFromGUI)

    # Stimulation
    if hasattr(slicer,'vtkMRMLFiberBundleNode') and hasattr(slicer.vtkMRMLFiberBundleNode,'GetExtractFromROI'):
      self.ui.stimulationCollapsibleButton.enabled = True
    else:
      self.ui.stimulationCollapsibleButton.enabled = False
      self.ui.stimulationCollapsibleButton.collapsed = True
      self.ui.stimulationCollapsibleButton.setToolTip('Updated SlicerDMRI Extension needed for stimulation module')
      
    self.ui.stimulationActiveCheckBox.connect('toggled(bool)', self.onStimulationActivate)
    self.ui.stimulationAmplitudeSpinBox.valueChanged.connect(self.updateStimulationRadius)

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
    # Do not react to parameter node changes (GUI wlil be updated when the user enters into the module)
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

  def onSceneStartClose(self, caller, event):
    """
    Called just before the scene is closed.
    """
    # Parameter node will be reset, do not use it anymore
    self.setParameterNode(None)
    # reset trajectory state
    self.setTrajectoryLayout([0]*9)

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
      with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)
    self._parameterNode = inputParameterNode
    if self._parameterNode is not None:
      self.addObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

    # Initial GUI update
    self.updateGUIFromParameterNode()

  @vtk.calldata_type(vtk.VTK_OBJECT)
  def onNodeAdded(self, caller, event, calldata):

    # todo send PR to opigtlink so that markups name is set before adding to scene
    qt.QTimer.singleShot(100, lambda cd=calldata: self.onNodeWithNameAdded(cd))

  def onNodeWithNameAdded(self, node):

    if (not node.GetName().startswith("LeadOR")) or self._parameterNode is None or self._updatingGUIFromParameterNode:
      return  
    wasModified = self._parameterNode.StartModify()  # Modify all properties in a single batch
    
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    folderID = self._parameterNode.GetParameter("LeadORIGTLFolderID")
    if folderID == '':
      folderID = shNode.CreateFolderItem(shNode.GetSceneItemID(), "LeadOR-IGTL")
      self._parameterNode.SetParameter("LeadORIGTLFolderID", str(folderID))
      shNode.SetItemExpanded(folderID, False)
    shNode.SetItemParent(shNode.GetItemByDataNode(node), int(folderID))
    
    subname = node.GetName().split(':')[-1]
    if subname == "ChannelsNames":
      self.addObserver(node, node.TextModifiedEvent, lambda c,e,n=node: self.onChannelsNamesModified(n))
      self.onChannelsNamesModified(node)
    elif subname == "DTT":
      self._parameterNode.SetNodeReferenceID("DistanceToTargetTransform", node.GetID())
    elif isinstance(node,slicer.vtkMRMLTextNode):
      newFeature = {'name':subname, 'sourceNodeID':node.GetID(), 'projectTo':'Tube', 'property':'', 'visible':1}
      featuresList = json.loads(self._parameterNode.GetParameter("FeaturesJson"))
      featuresList.append(newFeature)
      self._parameterNode.SetParameter("FeaturesJson", json.dumps(featuresList))
      self.addObserver(node, node.TextModifiedEvent, lambda c,e,n=node: self.onFeatureTextModified(n))

    self._parameterNode.EndModify(wasModified)

  def onChannelsNamesModified(self, channelsNamesNode):
    wasModified = self._parameterNode.StartModify()
    self._parameterNode.SetParameter("UnlinkedChannels", channelsNamesNode.GetText())
    self._parameterNode.EndModify(wasModified)

  def onFeatureTextModified(self, node):
    featuresList = json.loads(self._parameterNode.GetParameter("FeaturesJson"))
    for feature in featuresList:
      if node.GetID() == feature['sourceNodeID']:
        self.logic.setUpFeature(**feature)
        return

  @vtk.calldata_type(vtk.VTK_OBJECT)
  def onNodeRemoved(self, caller, event, calldata):
    if self._parameterNode is not None:
      wasModified = self._parameterNode.StartModify()  # Modify all properties in a single batch
    
      if isinstance(calldata,slicer.vtkMRMLTextNode):
        featuresList = json.loads(self._parameterNode.GetParameter("FeaturesJson"))
        for i,feature in enumerate(featuresList):
          if feature['sourceNodeID'] == calldata.GetID():
            featuresList.pop(i)
            break
        self._parameterNode.SetParameter("FeaturesJson", json.dumps(featuresList))

      self._parameterNode.EndModify(wasModified)
  
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
    self.ui.planningTransformComboBox.setCurrentNode(self._parameterNode.GetNodeReference("TrajectoryTransform"))
    
    transformsAvailable = bool(self._parameterNode.GetNodeReference("DistanceToTargetTransform") and self._parameterNode.GetNodeReference("TrajectoryTransform"))
    self.ui.recordSequenceToolButton.enabled = transformsAvailable
    self.ui.trajectoryPresetComboBox.enabled = transformsAvailable
    self.ui.trajectoriesLayoutFrame.enabled = transformsAvailable
    self.ui.linkChannelsToTrajectoriesPushButton.enabled = transformsAvailable
    self.ui.stimulationCollapsibleButton.enabled = transformsAvailable and hasattr(slicer,'vtkMRMLFiberBundleNode') and hasattr(slicer.vtkMRMLFiberBundleNode,'GetExtractFromROI')
    self.ui.setDefaultResliceDriverToolButton.enabled = transformsAvailable and hasattr(slicer.modules,'volumereslicedriver')
    
    self.ui.setDefaultResliceDriverToolButton.toolTip = 'Set default reslice driver.' if hasattr(slicer.modules,'volumereslicedriver') else 'Install SlicerIGT to use volume reslice driver.'

    trajectories = json.loads(self._parameterNode.GetParameter("TrajectoriesJson"))
    for i,trajectory in enumerate(trajectories):
      toolButton = getattr(self.ui, 'TrajectoryToolButton_'+str(i))
      toolButton.setChecked(trajectory["active"])
      toolButton.setToolTip(trajectory["channelName"])

    linkedChannels = set([trajectory["channelName"] for trajectory in trajectories])
    unlinkedChannels = set(self._parameterNode.GetParameter("UnlinkedChannels").split(","))
    self.ui.unlinkedChannelsListWidget.clear()
    self.ui.unlinkedChannelsListWidget.addItems(list(unlinkedChannels.difference(linkedChannels)))

    featuresList = json.loads(self._parameterNode.GetParameter("FeaturesJson"))
    self.ui.featuresTableWidget.updateNumberOfRows(len(featuresList))
    for i,feature in enumerate(featuresList):
      self.ui.featuresTableWidget.updateNthRowFromFeature(i, feature)

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
    self._parameterNode.SetNodeReferenceID("TrajectoryTransform", self.ui.planningTransformComboBox.currentNodeID)

    trajectories = json.loads(self._parameterNode.GetParameter("TrajectoriesJson"))
    for i,trajectory in enumerate(trajectories):
      toolButton = getattr(self.ui, 'TrajectoryToolButton_'+str(i))
      trajectory["active"] = toolButton.checked
      trajectory["channelName"] = toolButton.toolTip.replace('<p>','').replace('</p>','')
    self._parameterNode.SetParameter("TrajectoriesJson", json.dumps(trajectories))

    featuresList = json.loads(self._parameterNode.GetParameter("FeaturesJson"))
    for i,feature in enumerate(featuresList):
      featureUpdated = self.ui.featuresTableWidget.updateFeatureFromNthRow(feature, i)
      if featureUpdated:
        self.logic.setUpFeature(**feature)
    self._parameterNode.SetParameter("FeaturesJson", json.dumps(featuresList))

    self._parameterNode.EndModify(wasModified)

  def setTransformsHierarchy(self):
    if self._parameterNode and self._parameterNode.GetNodeReference("DistanceToTargetTransform") and self._parameterNode.GetNodeReference("TrajectoryTransform"):
      self._parameterNode.GetNodeReference("DistanceToTargetTransform").SetAndObserveTransformNodeID(self._parameterNode.GetNodeReferenceID("TrajectoryTransform"))

  def onUnlinkedChannelsSelectionChanged(self):
    selection = self.ui.unlinkedChannelsListWidget.selectedItems()
    if len(selection)>1:
      selection[0].setSelected(0)

  def trajectoryLayoutToggle(self, enabled, N):
    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return
    if enabled and len(self.ui.unlinkedChannelsListWidget.selectedItems()):
      linkChannelName = self.ui.unlinkedChannelsListWidget.selectedItems()[0].text()
    else:
      linkChannelName = ''
    trajectories = json.loads(self._parameterNode.GetParameter("TrajectoriesJson"))
    trajectories[N]["active"] = enabled
    trajectories[N]["channelName"] = linkChannelName
    self._parameterNode.SetParameter("TrajectoriesJson", json.dumps(trajectories))
    self.setUpTrajectories()

  def onLinkChannelsToTrajectoriesPushButton(self):
    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return
    trajectories = json.loads(self._parameterNode.GetParameter("TrajectoriesJson"))
    standardChannels = ['Anterio'+self.ui.leftTrajectoryLabel.text, 'Anterior', 'Anterio'+self.ui.rightTrajectoryLabel.text,\
                        self.ui.leftTrajectoryLabel.text, 'Central', self.ui.rightTrajectoryLabel.text,\
                        'Posterio'+self.ui.leftTrajectoryLabel.text, 'Posterior', 'Posterio'+self.ui.rightTrajectoryLabel.text]
    standardChannels = [standardChannel.lower() for standardChannel in standardChannels]
    unlinkedChannels = self._parameterNode.GetParameter("UnlinkedChannels").split(",")
    for unlinkedChannel in unlinkedChannels:
      unlinkedChannelLower = unlinkedChannel.lower()
      idx = standardChannels.index(unlinkedChannelLower) if unlinkedChannelLower in standardChannels else None
      if idx is not None:
        trajectories[idx]["active"] = 1
        trajectories[idx]["channelName"] = unlinkedChannel       
    self._parameterNode.SetParameter("TrajectoriesJson", json.dumps(trajectories))
    self.setUpTrajectories()
    
  def setTrajectoryLayoutPreset(self, text):
    if text == "Cross (x)":
      enabledList=[1,0,1,0,1,0,1,0,1]
    elif text == "Plus (+)":
      enabledList=[0,1,0,1,1,1,0,1,0]
    elif text == "Center (.)":
      enabledList=[0,0,0,0,1,0,0,0,0]
    if text != "Select...":
      self.setTrajectoryLayout(enabledList)
    self.ui.trajectoryPresetComboBox.setCurrentText("Select...")

  def setTrajectoryLayout(self, enabledList):
    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return
    trajectories = json.loads(self._parameterNode.GetParameter("TrajectoriesJson"))
    for N,enabled in enumerate(enabledList):
      trajectories[N]["active"] = enabled
      trajectories[N]["channelName"] = ''
    self._parameterNode.SetParameter("TrajectoriesJson", json.dumps(trajectories))
    self.setUpTrajectories()

  def trajectoryVisualizationChanged(self):
    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return
    trajectories = json.loads(self._parameterNode.GetParameter("TrajectoriesJson"))
    for i in range(self.ui.trajectoryVisualizationComboBox.model().rowCount()):
      index = self.ui.trajectoryVisualizationComboBox.model().index(i,0)
      key = index.data().lower() + 'Visibility'
      value = bool(self.ui.trajectoryVisualizationComboBox.checkState(index))
      for trajectory in trajectories:
        trajectory[key] = value
    self._parameterNode.SetParameter("TrajectoriesJson", json.dumps(trajectories))
    self.setUpTrajectories()

  def setUpTrajectories(self):
    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return
    distanceToTargetTransformID = self.ui.distanceToTargetComboBox.currentNodeID
    trajectories = json.loads(self._parameterNode.GetParameter("TrajectoriesJson"))
    for i,trajectory in enumerate(trajectories):
      self.logic.setUpTrajectory(trajectoryNumber=i, distanceToTargetTransformID=distanceToTargetTransformID, **trajectory)

  def guessSideFromTransform(self, transformNode):
    if not transformNode:
      return
    currentPoint = [0.0] * 4
    matrix = vtk.vtkMatrix4x4()
    transformNode.GetMatrixTransformToWorld(matrix)
    matrix.MultiplyPoint([0.0, 0.0, 0.0, 1.0], currentPoint)
    guessRightSide = currentPoint[0] > 0
    self.ui.leftTrajectoryLabel.text = 'Medial' if guessRightSide else 'Lateral'  
    self.ui.rightTrajectoryLabel.text = 'Lateral' if guessRightSide else 'Medial' 

  def setDefaultResliceDriver(self):
    # Get Reslice Driver Logic
    try:    
        logic = slicer.modules.volumereslicedriver.logic()
    except:
        qt.QMessageBox.warning(qt.QWidget(),'','Reslice Driver Module not Found')
        return
    transformNodeID = self._parameterNode.GetNodeReferenceID("DistanceToTargetTransform")
    # Settings
    redSettings    = {'node':slicer.util.getNode('vtkMRMLSliceNodeRed'),    'mode':6, 'angle':90 , 'flip':True}
    yellowSettings = {'node':slicer.util.getNode('vtkMRMLSliceNodeYellow'), 'mode':5, 'angle':180, 'flip':False}
    greenSettings  = {'node':slicer.util.getNode('vtkMRMLSliceNodeGreen'),  'mode':4, 'angle':180, 'flip':False}
    # Set
    for settings in [redSettings, yellowSettings, greenSettings]:
        logic.SetDriverForSlice(    transformNodeID,    settings['node'])
        logic.SetModeForSlice(      settings['mode'],   settings['node'])
        logic.SetRotationForSlice(  settings['angle'],  settings['node'])
        logic.SetFlipForSlice(      settings['flip'],   settings['node'])

  def setUpSequenzeRecording(self):
    browserNode = self._parameterNode.GetNodeReference("browserNode")
    if browserNode is None:
      browserNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSequenceBrowserNode', 'LeadOR')
      self._parameterNode.SetNodeReferenceID("browserNode", browserNode.GetID())

    if self.ui.recordSequenceToolButton.isChecked():
      sequencesLogic = slicer.modules.sequences.logic()
      featuresList = json.loads(self._parameterNode.GetParameter("FeaturesJson"))
      nodesToSync = [feature['sourceNodeID'] for feature in featuresList]
      nodesToSync.insert(0, self._parameterNode.GetNodeReferenceID("DistanceToTargetTransform"))
      for nodeID in nodesToSync:
        sequenceNode = browserNode.GetSequenceNode(slicer.util.getNode(nodeID))
        if sequenceNode is None:
          sequenceNode = sequencesLogic.AddSynchronizedNode(None, slicer.util.getNode(nodeID), browserNode)
        browserNode.SetRecording(sequenceNode, True)

    browserNode.SetRecordingActive(self.ui.recordSequenceToolButton.isChecked())

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
    N = next(filter(lambda n: getattr(self.ui, 'TrajectoryToolButton_'+str(n)).actions()[0].checked, range(9)), None)
    trajectory = Trajectory.GetNthTrajectory(N)
    if N is None or trajectory is None:
      self.ui.stimulationActiveCheckBox.checked = False
      qt.QMessageBox().warning(qt.QWidget(),'','Set a Stimulation Source')
      return
    # set transform
    self.logic.VTASource.SetAndObserveTransformNodeID(trajectory.translationTransformNodeID)

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
      import glob
      import importlib
      import LeadORLib
      LeadORLibPath = os.path.join(os.path.dirname(__file__), 'LeadORLib')
      G = glob.glob(os.path.join(LeadORLibPath, '**','*.py'))
      for g in G:
        relativePath = os.path.relpath(g, LeadORLibPath) # relative path
        relativePath = os.path.splitext(relativePath)[0] # get rid of .py
        moduleParts = relativePath.split(os.path.sep) # separate
        importlib.import_module('.'.join(['LeadORLib']+moduleParts)) # import module
        module = LeadORLib
        for modulePart in moduleParts: # iterate over parts in order to load subpkgs
          module = getattr(module, modulePart)
        importlib.reload(module) # reload

    self.VTASource = None

    
  def setDefaultParameters(self, parameterNode):
    """
    Initialize parameter node with default settings.
    """
    if not parameterNode.GetParameter("UnlinkedChannels"):
      parameterNode.SetParameter("UnlinkedChannels", "")
    if not parameterNode.GetParameter("FeatureNames"):
      parameterNode.SetParameter("FeatureNames", "")
    if not parameterNode.GetParameter("FeaturesJson"):
      parameterNode.SetParameter("FeaturesJson", json.dumps([]))
    if not parameterNode.GetParameter("TrajectoriesJson"):
      parameterNode.SetParameter("TrajectoriesJson", json.dumps([{"active":0,"channelName":'',"modelVisibility":1,"lineVisibility":1,"tipVisibility":1}]*9))

  def setUpTrajectory(self, trajectoryNumber, distanceToTargetTransformID, active=True, channelName='', modelVisibility=1, lineVisibility=1, tipVisibility=1):
    if active:
      trajectory = Trajectory.InitOrGetNthTrajectory(trajectoryNumber)
      trajectory.setDistanceToTargetTransformID(distanceToTargetTransformID)
      trajectory.setChannelName(channelName)
      trajectory.setModelVisibility(modelVisibility)
      trajectory.setLineVisibility(lineVisibility)
      trajectory.setTipVisibility(tipVisibility)
    else:
      Trajectory.RemoveNthTrajectory(trajectoryNumber)

  def setUpFeature(self, sourceNodeID=None, name='', projectTo='Tube', property='', visible=1):
    if sourceNodeID is None:
      return
    feature = Feature(projectTo)
    feature.addSourceNode(sourceNodeID, property, visible)
    feature.update()

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
    # Close Open-Ephys
    import socket
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    result = sock.connect_ex(('127.0.0.1',37497))
    if result == 0:
      import requests
      url = "http://localhost:37497/api/"
      requests.put(url + "status", json={"mode" : "IDLE"})
      requests.put(url + "processors/103/config", json={"text" : "LOR IGTLDISCONNECT"})
      requests.put(url + "window", json={"command" : "quit"})
    # Clear Slicer
    while slicer.mrmlScene.GetNumberOfNodesByClass("vtkMRMLIGTLConnectorNode"):
      slicer.mrmlScene.GetNthNodeByClass(0,"vtkMRMLIGTLConnectorNode").Stop()
      slicer.mrmlScene.RemoveNode(slicer.mrmlScene.GetNthNodeByClass(0,"vtkMRMLIGTLConnectorNode"))
    slicer.mrmlScene.Clear()

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    # self.test_LeadORWithOpenEphys()
    # self.test_LeadORFeatures()
    self.test_LeadORFeaturesWithNan()
    # self.test_LeadORFeaturesBasic()

  def test_LeadORFeaturesBasic(self):
    self.delayDisplay("Starting the test")

    planningNode =slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLinearTransformNode','Planning')
    dttNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLinearTransformNode','DTT')
    dttNode.SetAndObserveTransformNodeID(planningNode.GetID())

    channelName = "Central"

    featureNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTextNode','feature1')
    featureNode.SetText("RecordingSiteDTT,"+channelName)
    
    logic = LeadORLogic()

    logic.setUpTrajectory(4, dttNode.GetID(), True, channelName, 0, 0, 0)
    logic.setUpFeature(featureNode.GetID(), 'feature1', 'Tube', 'RadiusAndColor')

    self.delayDisplay('Test passed')

  def test_LeadORFeatures(self):
    self.delayDisplay("Starting the test")

    planningNode =slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLinearTransformNode','Planning')
    dttNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLinearTransformNode','DTT')
    dttNode.SetAndObserveTransformNodeID(planningNode.GetID())

    channelName = "Central"

    featureNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTextNode','feature1')
    featureNode.SetText("RecordingSiteDTT,"+channelName+"\n\
                          10.0,80\n\
                          9.0,80\n\
                          8.0,80\n\
                          7.0,80\n\
                          6.0,160\n\
                          5.0,200\n\
                          4.0,140\n\
                          3.0,300\n\
                          2.0,200")
    
    logic = LeadORLogic()

    logic.setUpTrajectory(4, dttNode.GetID(), True, channelName, 0, 0, 0)
    logic.setUpFeature(featureNode.GetID(), 'feature1', 'Tube', 'RadiusAndColor')

    self.delayDisplay('Test passed')

  def test_LeadORFeaturesWithNan(self):
    self.delayDisplay("Starting the test")

    planningNode =slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLinearTransformNode','Planning')
    dttNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLinearTransformNode','DTT')
    dttNode.SetAndObserveTransformNodeID(planningNode.GetID())

    channelName = "Central"

    featureNode1 = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTextNode','feature1')
    featureNode1.SetText("RecordingSiteDTT,"+channelName+"\n\
                          10.0,80\n\
                          8.0,80\n\
                          7.0,80\n\
                          6.0,160\n\
                          5.0,200\n\
                          4.0,300\n\
                          2.0,200")

    featureNode2 = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTextNode','feature2')
    featureNode2.SetText("RecordingSiteDTT,"+channelName+"\n\
                          10.0,80\n\
                          9.0,90\n\
                          8.0,80\n\
                          6.0,300\n\
                          5.0,150\n\
                          4.0,80\n\
                          3.0,80\n\
                          2.0,80")
    
    logic = LeadORLogic()

    logic.setUpTrajectory(4, dttNode.GetID(), True, channelName, 0, 0, 0)
    logic.setUpFeature(featureNode1.GetID(), 'feature1', 'Tube', 'Radius')
    logic.setUpFeature(featureNode2.GetID(), 'feature2', 'Tube', 'Color')

    self.delayDisplay('Test passed')


  def test_LeadORWithOpenEphys(self):

    self.delayDisplay("Starting the test")

    # Get/create input data

    # Currently use local data.
    # This is sensitive data recorded during surgery.
    # TODO: see how to share an example dataset for other users.
    from sys import platform
    if platform == "win32":
      test_dir = "C:\\Users\\simon\\Desktop\\143UA53-test"
    elif platform == "darwin":
      test_dir = "/Users/simon/Desktop/143UA53-test"
    if not os.path.isdir(test_dir):
      return

    slicer.util.loadScene(os.path.join(test_dir, "ORScene.mrb"))

    # Add an IGTL Connector

    n = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLIGTLConnectorNode')
    n.SetServerPort(18944)
    n.Start()

    # Open Open-ephys instance with test config, connect to igtl and start aquisition

    if platform == "win32":
      open_ephys_exe = "C:\\Users\\simon\\repo\\plugin-GUI\\Build\\Release\\open-ephys.exe"
    elif platform == "darwin":
      open_ephys_exe = '/Users/simon/repo/plugin-GUI/Build/Release/open-ephys.app/Contents/MacOS/open-ephys'
    open_ephys_config = os.path.join(test_dir, "LeadORConfig.xml")
    import subprocess
    subprocess.Popen([open_ephys_exe, open_ephys_config])

    import requests
    url = "http://localhost:37497/api/"
    r = None
    while (r == None) or (r.json()['mode'] != 'IDLE'):
      try:
        r = requests.get(url + "status")
        break
      except:
        import time
        time.sleep(0.5)

    r = None
    while (r == None) or (r.json()['info'] != 'Connected!'):
      try:
        r = requests.put(url + "processors/103/config", json={"text" : "LOR IGTLCONNECT 18944"})
        break
      except:
        import time
        time.sleep(0.5)
    
    self.assertEqual(r.json()['info'], 'Connected!')

    r = requests.put(url + "status", json={"mode" : "ACQUIRE"})
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

  if mainWindow and not layoutNode.IsLayoutDescription(customLayoutId):
    layoutNode.AddLayoutDescription(customLayoutId, customLayout)                                         
    viewToolBar = mainWindow.findChild("QToolBar", "ViewToolBar")
    layoutMenu = viewToolBar.widgetForAction(viewToolBar.actions()[0]).menu()
    layoutSwitchActionParent = layoutMenu
    layoutSwitchAction = layoutSwitchActionParent.addAction("LeadOR")
    layoutSwitchAction.setData(customLayoutId)
    layoutSwitchAction.setIcon(qt.QIcon(":Icons/Go.png"))
