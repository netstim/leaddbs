import logging
import os
import glob
import numpy as np
import vtk, qt

import slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin

from StereotacticPlanLib.Widgets.CustomWidgets import CustomCoordinatesWidget, TransformableCoordinatesWidget

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
        self.parent.contributors = ["Simon Oxenford (Charite Berlin)"]  # TODO: replace with "Firstname Lastname (Organization)"
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

        # Custom Widgets        
        self.trajectoryCoordinateWidgets = {}
        for name in ['Entry', 'Target']:
            self.trajectoryCoordinateWidgets[name] =  TransformableCoordinatesWidget(name, self.setTransformableWidgetsState)
            self.trajectoryCoordinateWidgets[name].coordinatesChanged.connect(self.updateParameterNodeFromGUI)
            self.ui.trajectoriesCollapsibleButton.layout().insertRow(2,name + ':', self.trajectoryCoordinateWidgets[name])
        for widget in [self.trajectoryCoordinateWidgets['Entry'], self.ui.rollAngleSliderWidget]:
            widget.setVisible(False)
            self.ui.trajectoriesCollapsibleButton.layout().labelForField(widget).setVisible(False)
            self.ui.trajectoryModeComboBox.currentTextChanged.connect(lambda t,w=widget,target_t='Target Entry Roll': [w.setVisible(t==target_t), self.ui.trajectoriesCollapsibleButton.layout().labelForField(w).setVisible(t==target_t)])
        for widget in [self.ui.ringAngleSliderWidget, self.ui.arcAngleSliderWidget, self.ui.mountingComboBox]:
            self.ui.trajectoryModeComboBox.currentTextChanged.connect(lambda t,w=widget,target_t='Target Mounting Ring Arc': [w.setVisible(t==target_t), self.ui.trajectoriesCollapsibleButton.layout().labelForField(w).setVisible(t==target_t)])
        
        self.referenceToFrameCoordinateWidgets = {}
        for name in ['Reference MS', 'Reference PC', 'Reference AC']:
            self.referenceToFrameCoordinateWidgets[name] =  TransformableCoordinatesWidget(name, self.setTransformableWidgetsState)
            self.referenceToFrameCoordinateWidgets[name].coordinatesChanged.connect(self.updateParameterNodeFromGUI)
            self.ui.referenceToFrameCollapsibleButton.layout().insertRow(1, name + ':', self.referenceToFrameCoordinateWidgets[name])
        for name in ['Frame MS', 'Frame PC', 'Frame AC']:
            self.referenceToFrameCoordinateWidgets[name] =  CustomCoordinatesWidget(name)
            self.referenceToFrameCoordinateWidgets[name].coordinatesChanged.connect(self.updateParameterNodeFromGUI)
            self.referenceToFrameCoordinateWidgets[name].setVisible(False)
            self.ui.referenceToFrameCollapsibleButton.layout().insertRow(5, name + ':', self.referenceToFrameCoordinateWidgets[name])
            self.ui.referenceToFrameCollapsibleButton.layout().labelForField(self.referenceToFrameCoordinateWidgets[name]).setVisible(False)
            self.ui.referenceToFrameModeComboBox.currentTextChanged.connect(lambda t,w=self.referenceToFrameCoordinateWidgets[name]: [w.setVisible(t=='ACPC Register'), self.ui.referenceToFrameCollapsibleButton.layout().labelForField(w).setVisible(t=='ACPC Register')])

        buttonSize = self.trajectoryCoordinateWidgets['Entry'].transformButton.height
        transformReferenceVolumeAction = qt.QAction()
        transformReferenceVolumeAction.setIcon(qt.QIcon(":/Icons/Transforms.png"))
        transformReferenceVolumeAction.setCheckable(True)
        self.transformReferenceVolumeButton = qt.QToolButton()
        self.transformReferenceVolumeButton.setDefaultAction(transformReferenceVolumeAction)
        self.transformReferenceVolumeButton.setToolButtonStyle(qt.Qt.ToolButtonIconOnly)
        self.transformReferenceVolumeButton.setFixedSize(buttonSize, buttonSize)
        self.transformReferenceVolumeButton.toggled.connect(self.updateParameterNodeFromGUI)
        self.transformReferenceVolumeButton.toggled.connect(self.setTransformableWidgetsState)
        self.ui.referenceVolumeLayout.addWidget(self.transformReferenceVolumeButton)

        viewTrajectoryAction = qt.QAction()
        viewTrajectoryAction.setIcon(qt.QIcon(":/Icons/Small/SlicerVisible.png"))
        viewTrajectoryAction.setCheckable(True)
        self.ui.viewTrajectoryToolButton.setDefaultAction(viewTrajectoryAction)
        self.ui.viewTrajectoryToolButton.setFixedSize(buttonSize, buttonSize)
        self.ui.viewTrajectoryToolButton.connect("toggled(bool)", self.onViewTrajectoryToggled)

        resliceDriverAction = qt.QAction()
        resliceDriverAction.setIcon(qt.QIcon(qt.QPixmap(self.resourcePath('Icons/VolumeResliceDriver.png'))))
        resliceDriverAction.setCheckable(True)
        self.ui.resliceDriverToolButton.setDefaultAction(resliceDriverAction)
        self.ui.resliceDriverToolButton.connect("toggled(bool)", self.setDefaultResliceDriver)
        self.ui.resliceDriverToolButton.setFixedSize(buttonSize, buttonSize)

        importFromActionGroup = qt.QActionGroup(self.ui.importFromToolButton)
        importFromOptions = [os.path.basename(g).replace('.py','') for g in glob.glob(os.path.join(os.path.dirname(__file__), 'StereotacticPlanLib', 'ImportFrom', 'Import_From_*.py'))]
        for option in importFromOptions:
            importAction = qt.QAction(option.replace('_',' '), self.ui.importFromToolButton)
            importAction.triggered.connect(lambda b,o=option: self.importTrajectoryFrom(o))
            importFromActionGroup.addAction(importAction)
        importFromMenu = qt.QMenu(self.ui.importFromToolButton)
        importFromMenu.addActions(importFromActionGroup.actions())
        self.ui.importFromToolButton.setMenu(importFromMenu)
        self.ui.importFromToolButton.setFixedSize(buttonSize, buttonSize)

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

        # These connections ensure that whenever user changes some settings on the GUI, that is saved in the MRML scene
        # (in the selected parameter node).
        self.ui.trajectoryTransformNodeComboBox.connect("currentNodeChanged(vtkMRMLNode*)", lambda n,w=self.ui.resliceDriverToolButton: self.setDefaultResliceDriver(w.checked))
        self.ui.trajectoryTransformNodeComboBox.connect("currentNodeChanged(vtkMRMLNode*)", self.updatePreviewLineTransform)
        self.ui.trajectoryTransformNodeComboBox.connect("currentNodeChanged(vtkMRMLNode*)", self.updateCurrentTrajectory)
        self.ui.referenceToFrameTransformNodeComboBox.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
        self.ui.referenceVolumeNodeComboBox.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
        self.ui.trajectoryModeComboBox.connect("currentTextChanged(QString)",  self.updateParameterNodeFromGUI)
        self.ui.arcAngleSliderWidget.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
        self.ui.ringAngleSliderWidget.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
        self.ui.rollAngleSliderWidget.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
        self.ui.mountingComboBox.currentIndexChanged.connect(self.updateParameterNodeFromGUI)
        self.ui.referenceToFrameModeComboBox.currentIndexChanged.connect(self.updateParameterNodeFromGUI)

        # Buttons
        self.ui.calculateReferenceToFramePushButton.connect('clicked(bool)', self.onCalculateReferenceToFrame)

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
        self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

    def onSceneStartClose(self, caller, event):
        """
        Called just before the scene is closed.
        """
        # Parameter node will be reset, do not use it anymore
        self.setTransformableWidgetsState(False)
        if self._parameterNode:
            self._parameterNode.SetNodeReferenceID("ReferenceToFrameTransform","")
        self.setParameterNode(None)

    def onSceneEndClose(self, caller, event):
        """
        Called just after the scene is closed.
        """
        # If this module is shown while the scene is closed then recreate a new parameter node immediately
        if self.parent.isEntered:
            self.initializeParameterNode()
        if hasattr(self,'referenceToFrameCoordinateWidgets'):
            for widget in self.referenceToFrameCoordinateWidgets.values():
                widget.setUpMarkupsNode()
            for widget in self.trajectoryCoordinateWidgets.values():
                widget.setUpMarkupsNode()

    def initializeParameterNode(self):
        """
        Ensure parameter node exists and observed.
        """
        # Parameter node stores all user choices in parameter values, node selections, etc.
        # so that when the scene is saved and reloaded, these settings are restored.

        self.setParameterNode(self.logic.getParameterNode())

        # Select default input nodes if nothing is selected yet to save a few clicks for the user
        # if not self._parameterNode.GetNodeReference("InputVolume"):
        #     firstVolumeNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLScalarVolumeNode")
        #     if firstVolumeNode:
        #         self._parameterNode.SetNodeReferenceID("InputVolume", firstVolumeNode.GetID())

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

    def createNewTrajectory(self, nodeID):
        node = slicer.util.getNode(nodeID)
        node.SetAttribute('NetstimStereotacticPlan', '1')
        node.SetAttribute('Entry', '0,0,0;RAS;0')
        node.SetAttribute('Target', '0,0,0;RAS;0')
        node.SetAttribute('Mounting', 'lateral-right')
        node.SetAttribute('Ring', '90')
        node.SetAttribute('Arc', '90')
        node.SetAttribute('Roll', '0')
        node.SetAttribute('Mode', 'Target Mounting Ring Arc')

    def getTrajectoryNodesIDsInScene(self):
        trajectoryNodesIDs = []
        for i in range(slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLLinearTransformNode')):
            node = slicer.mrmlScene.GetNthNodeByClass(i,'vtkMRMLLinearTransformNode')
            if 'NetstimStereotacticPlan' in node.GetAttributeNames():
                trajectoryNodesIDs.append(node.GetID())
        return trajectoryNodesIDs

    def importTrajectoryFrom(self, importerName):
        # Get importer module
        import StereotacticPlanLib.ImportFrom
        import importlib
        importlib.import_module('.'.join(['StereotacticPlanLib', 'ImportFrom', importerName]))
        importerModule = getattr(StereotacticPlanLib.ImportFrom, importerName)
        importInFrameSpace = self.transformReferenceVolumeButton.checked
        dialog = importerModule.ImporterDialog()
        dialog.run(importInFrameSpace)
        if dialog.selectedFile is None:
            return
        wasModified = self._parameterNode.StartModify()
        importer = importerModule.Importer(dialog.selectedFile)
        if dialog.importACPCCoordinates:
            importer.setACPCCoordinatesToParameterNode(self._parameterNode)
        if dialog.computeReferenceToFrame:
            node = importer.getReferenceToFrameTransform()
            self._parameterNode.SetNodeReferenceID("ReferenceToFrameTransform", node.GetID())
        if dialog.importDICOM:
            node = importer.getReferenceVolumeFromDICOM(dialog.DICOMDir)
            self._parameterNode.SetNodeReferenceID("ReferenceVolume", node.GetID())
        loadedNodeIDs = importer.getTrajectoryTransforms(importInFrameSpace)
        self._parameterNode.SetNodeReferenceID("CurrentTrajectoryTransform", loadedNodeIDs[0])
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

        currentTrajectoryAvailable = self._parameterNode.GetNodeReference("CurrentTrajectoryTransform") is not None

        self.ui.trajectoryModeComboBox.setEnabled(currentTrajectoryAvailable)
        self.ui.mountingComboBox.setEnabled(currentTrajectoryAvailable)
        self.ui.ringAngleSliderWidget.setEnabled(currentTrajectoryAvailable)
        self.ui.arcAngleSliderWidget.setEnabled(currentTrajectoryAvailable)
        self.ui.rollAngleSliderWidget.setEnabled(currentTrajectoryAvailable)

        if currentTrajectoryAvailable:
            node = self._parameterNode.GetNodeReference("CurrentTrajectoryTransform")
            if node.GetID() not in self.getTrajectoryNodesIDsInScene():
                self.createNewTrajectory(node.GetID())
            self.ui.trajectoryModeComboBox.setCurrentText(node.GetAttribute('Mode'))
            self.ui.mountingComboBox.setCurrentText(node.GetAttribute('Mounting'))
            self.ui.ringAngleSliderWidget.value = float(node.GetAttribute('Ring'))
            self.ui.arcAngleSliderWidget.value = float(node.GetAttribute('Arc'))
            self.ui.rollAngleSliderWidget.value = float(node.GetAttribute('Roll'))
            self.updateCoordinatesWidgetFromTrajectory(node)
        else:
            self.ui.referenceToFrameTransformNodeComboBox.setCurrentNode(None)

        for widget in self.trajectoryCoordinateWidgets.values():
            widget.setTransformNodeID(self._parameterNode.GetNodeReferenceID("ReferenceToFrameTransform"))
            if not currentTrajectoryAvailable:
                widget.reset()
            widget.setEnabled(currentTrajectoryAvailable)

        for name, widget in self.referenceToFrameCoordinateWidgets.items():
            if isinstance(widget, TransformableCoordinatesWidget):
                widget.setTransformNodeID(self._parameterNode.GetNodeReferenceID("ReferenceToFrameTransform"))
            coords, system = self._parameterNode.GetParameter(name).split(';')
            widget.setSystem(system)
            widget.coordinates = coords

        self.ui.viewTrajectoryToolButton.setEnabled(currentTrajectoryAvailable)
        self.ui.resliceDriverToolButton.setEnabled(currentTrajectoryAvailable)
        self.ui.calculateReferenceToFramePushButton.setEnabled(bool(self._parameterNode.GetNodeReference("ReferenceToFrameTransform")))

        self.ui.referenceToFrameTransformNodeComboBox.setCurrentNode(self._parameterNode.GetNodeReference("ReferenceToFrameTransform"))
        self.ui.referenceVolumeNodeComboBox.setCurrentNode(self._parameterNode.GetNodeReference("ReferenceVolume"))
        self.ui.trajectoryTransformNodeComboBox.setCurrentNode(self._parameterNode.GetNodeReference("CurrentTrajectoryTransform"))

        self.transformReferenceVolumeButton.setEnabled(self._parameterNode.GetNodeReference("ReferenceToFrameTransform") and self._parameterNode.GetNodeReference("ReferenceVolume"))
        self.ui.referenceToFrameModeComboBox.currentText = self._parameterNode.GetParameter("ReferenceToFrameMode")

        # All the GUI updates are done
        self._updatingGUIFromParameterNode = False

    def updateCurrentTrajectory(self, caller=None, event=None):
        if self._parameterNode is None or self._updatingGUIFromParameterNode:
            return
        self._parameterNode.SetNodeReferenceID("CurrentTrajectoryTransform", self.ui.trajectoryTransformNodeComboBox.currentNodeID)
        if self._parameterNode.GetNodeReference("CurrentTrajectoryTransform") is not None:
            self.onCalculateTrajectory()

    def updateParameterNodeFromGUI(self, caller=None, event=None):
        """
        This method is called when the user makes any change in the GUI.
        The changes are saved into the parameter node (so that they are restored when the scene is saved and loaded).
        """

        if self._parameterNode is None or self._updatingGUIFromParameterNode:
            return

        wasModified = self._parameterNode.StartModify()  # Modify all properties in a single batch

        if self._parameterNode.GetNodeReference("CurrentTrajectoryTransform") is not None:
            node = self._parameterNode.GetNodeReference("CurrentTrajectoryTransform")
            node.SetAttribute('Mode', self.ui.trajectoryModeComboBox.currentText)
            node.SetAttribute('Mounting', self.ui.mountingComboBox.currentText)
            node.SetAttribute('Ring', str(self.ui.ringAngleSliderWidget.value))
            node.SetAttribute('Arc', str(self.ui.arcAngleSliderWidget.value))
            node.SetAttribute('Roll', str(self.ui.rollAngleSliderWidget.value))
            self.updateTrajectoryFromCoordinatesWidget(node)

        for name, widget in self.referenceToFrameCoordinateWidgets.items():
            self._parameterNode.SetParameter(name, '%s;%s' % (widget.coordinates, widget.getSystem()) )

        self._parameterNode.SetNodeReferenceID("ReferenceToFrameTransform", self.ui.referenceToFrameTransformNodeComboBox.currentNodeID)
        self._parameterNode.SetNodeReferenceID("ReferenceVolume", self.ui.referenceVolumeNodeComboBox.currentNodeID)

        if self.ui.referenceVolumeNodeComboBox.currentNodeID != "":
            slicer.util.getNode(self.ui.referenceVolumeNodeComboBox.currentNodeID).SetAndObserveTransformNodeID(self.ui.referenceToFrameTransformNodeComboBox.currentNodeID if self.transformReferenceVolumeButton.checked else None)

        self._parameterNode.SetParameter("ReferenceToFrameMode", self.ui.referenceToFrameModeComboBox.currentText)

        if self._parameterNode.GetNodeReference("CurrentTrajectoryTransform") is not None:
            self.onCalculateTrajectory()

        self._parameterNode.EndModify(wasModified)


    def updateCoordinatesWidgetFromTrajectory(self, node):
        for name, widget in self.trajectoryCoordinateWidgets.items():
            coords, system, inFrameSpace = node.GetAttribute(name).split(';')
            widget.setSystem(system)
            widget.coordinates = coords

    def updateTrajectoryFromCoordinatesWidget(self, node):
        for name, widget in self.trajectoryCoordinateWidgets.items():
             node.SetAttribute(name, '%s;%s;%d' % (widget.coordinates, widget.getSystem(), widget.transformButton.checked))

    def setTransformableWidgetsState(self, state):           
        for widget in self.trajectoryCoordinateWidgets.values():
            widget.transformButton.setChecked(state)
        for widget in self.referenceToFrameCoordinateWidgets.values():
            if isinstance(widget, TransformableCoordinatesWidget):
                widget.transformButton.setChecked(state)
        self.transformReferenceVolumeButton.setChecked(state)
        # Also update coords from other trajectories
        trajectoriesNodesIDs = self.getTrajectoryNodesIDsInScene()
        for nodeID in trajectoriesNodesIDs:
            if nodeID == self.ui.trajectoryTransformNodeComboBox.currentNodeID:
                continue
            for name in ['Entry', 'Target']:
                coords, system, inFrameSpace = slicer.util.getNode(nodeID).GetAttribute(name).split(';')
                coords = np.fromstring(coords, dtype=float, sep=',')
                coords = coords if system=='RAS' else self.logic.transformCoordsFromXYZToRAS(coords)
                inFrameSpace = bool(int(inFrameSpace))
                if inFrameSpace and not state:
                    coords = self._parameterNode.GetNodeReference("ReferenceToFrameTransform").GetTransformFromParent().TransformFloatPoint(coords)
                elif not inFrameSpace and state:
                    coords = self._parameterNode.GetNodeReference("ReferenceToFrameTransform").GetTransformToParent().TransformFloatPoint(coords)
                coords = coords if system=='RAS' else self.logic.transformCoordsFromRASToXYZ(coords)
                slicer.util.getNode(nodeID).SetAttribute(name, '%s;%s;%d' % (','.join([str(x) for x in coords]), system, state))

    def onCalculateReferenceToFrame(self):
        if self.ui.referenceToFrameModeComboBox.currentText == 'ACPC Align':
            self.logic.runACPCAlignment(self.ui.referenceToFrameTransformNodeComboBox.currentNode(),
                                        self.referenceToFrameCoordinateWidgets['Reference AC'].getNumpyCoordinates(system='RAS'),
                                        self.referenceToFrameCoordinateWidgets['Reference PC'].getNumpyCoordinates(system='RAS'),
                                        self.referenceToFrameCoordinateWidgets['Reference MS'].getNumpyCoordinates(system='RAS'))

        elif self.ui.referenceToFrameModeComboBox.currentText == 'ACPC Register':
            sourceCoordinates = [self.referenceToFrameCoordinateWidgets[name].getNumpyCoordinates(system='RAS') for name in ['Reference MS', 'Reference PC', 'Reference AC']]
            targetCoordinates = [self.referenceToFrameCoordinateWidgets[name].getNumpyCoordinates(system='RAS') for name in ['Frame MS', 'Frame PC', 'Frame AC']]
            self.logic.runFiducialRegistration(self.ui.referenceToFrameTransformNodeComboBox.currentNode(),
                                                sourceCoordinates,
                                                targetCoordinates)


    def onCalculateTrajectory(self):
        try:
            if self.ui.trajectoryModeComboBox.currentText == 'Target Mounting Ring Arc':
                self.logic.computeTrajectoryFromTargetMountingRingArc(self.ui.trajectoryTransformNodeComboBox.currentNode(),
                                self.trajectoryCoordinateWidgets['Target'].getNumpyCoordinates(system='RAS'),
                                self.ui.mountingComboBox.currentText,
                                self.ui.ringAngleSliderWidget.value,
                                self.ui.arcAngleSliderWidget.value)
            else:
                self.logic.computeTrajectoryFromTargetEntryRoll(self.ui.trajectoryTransformNodeComboBox.currentNode(),
                                self.trajectoryCoordinateWidgets['Target'].getNumpyCoordinates(system='RAS'),
                                self.trajectoryCoordinateWidgets['Entry'].getNumpyCoordinates(system='RAS'),
                                self.ui.rollAngleSliderWidget.value)
        except Exception as e:
            slicer.util.errorDisplay("Failed to compute transform: "+str(e))
            import traceback
            traceback.print_exc()


    def onViewTrajectoryToggled(self, state):
        if state:
            # add node and default points
            markupsLineNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsLineNode')
            self._parameterNode.SetNodeReferenceID("PreviewLine", markupsLineNode.GetID())
            markupsLineNode.AddControlPointWorld(vtk.vtkVector3d(0,0,0), 'Target')
            markupsLineNode.AddControlPointWorld(vtk.vtkVector3d(0,0,80), '2')
            self.updatePreviewLineTransform(self.ui.trajectoryTransformNodeComboBox.currentNode())
        elif self._parameterNode.GetNodeReferenceID("PreviewLine"):
            # remove node
            slicer.mrmlScene.RemoveNode(self._parameterNode.GetNodeReference("PreviewLine"))

    def updatePreviewLineTransform(self, node):
        if self._parameterNode and self._parameterNode.GetNodeReferenceID("PreviewLine"):
            self._parameterNode.GetNodeReference("PreviewLine").SetAndObserveTransformNodeID(node.GetID() if node else None)

    def setDefaultResliceDriver(self, state):
        if state:
            # Get Reslice Driver Logic
            try:    
                logic = slicer.modules.volumereslicedriver.logic()
            except:
                qt.QMessageBox.warning(qt.QWidget(),'','Reslice Driver Module not Found')
                return
            transformNodeID = self.ui.trajectoryTransformNodeComboBox.currentNodeID
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
        else:
            sliceNodes = slicer.util.getNodesByClass("vtkMRMLSliceNode")
            for sliceNode in sliceNodes:
                if sliceNode.GetName() == 'Red':
                    sliceNode.SetOrientationToAxial()
                elif sliceNode.GetName() == 'Green':
                    sliceNode.SetOrientationToCoronal()
                elif sliceNode.GetName() == 'Yellow':
                    sliceNode.SetOrientationToSagittal()


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

    if slicer.util.settingsValue('Developer/DeveloperMode', False, converter=slicer.util.toBool):
      import glob
      import importlib
      import StereotacticPlanLib
      StereotacticPlanLibPath = os.path.join(os.path.dirname(__file__), 'StereotacticPlanLib')
      G = glob.glob(os.path.join(StereotacticPlanLibPath, '**','*.py'))
      for g in G:
        relativePath = os.path.relpath(g, StereotacticPlanLibPath) # relative path
        relativePath = os.path.splitext(relativePath)[0] # get rid of .py
        moduleParts = relativePath.split(os.path.sep) # separate
        importlib.import_module('.'.join(['StereotacticPlanLib']+moduleParts)) # import module
        module = StereotacticPlanLib
        for modulePart in moduleParts: # iterate over parts in order to load subpkgs
          module = getattr(module, modulePart)
        importlib.reload(module) # reload


    def setDefaultParameters(self, parameterNode):
        """
        Initialize parameter node with default settings.
        """
        for name in ["Reference AC", "Reference PC", "Reference MS", "Frame AC", "Frame PC", "Frame MS"]:
            if not parameterNode.GetParameter(name):
                parameterNode.SetParameter(name, "0,0,0;RAS")
        if not parameterNode.GetParameter("TrajectoryIndex"):
            parameterNode.SetParameter("TrajectoryIndex", "-1")
        if not parameterNode.GetParameter("TransformableWidgetsChecked"):
            parameterNode.SetParameter("TransformableWidgetsChecked", "0")
        if not parameterNode.GetParameter("ReferenceToFrameMode"):
            parameterNode.SetParameter("ReferenceToFrameMode","ACPC Align")
        if not parameterNode.GetNodeReference("CurrentTrajectoryTransform"):
            parameterNode.SetNodeReferenceID("CurrentTrajectoryTransform","")
        if not parameterNode.GetNodeReference("ReferenceToFrameTransform"):
            parameterNode.SetNodeReferenceID("ReferenceToFrameTransform","")

    def transformCoordsFromXYZToRAS(self, coords):
        return np.dot(self.getFrameXYZToRASTransform(), np.append(coords, 1))[:3]

    def transformCoordsFromRASToXYZ(self, coords):
        return np.dot(np.linalg.inv(self.getFrameXYZToRASTransform()), np.append(coords, 1))[:3]

    def getFrameXYZToRASTransform(self):
        # Headring coordinates to Slicer world (matching center)
        frameToRAS = np.array([[ -1,  0,  0,  100],
                                [  0,  1,  0, -100],
                                [  0,  0, -1,  100],
                                [  0,  0,  0,    1]])
        return frameToRAS                       

    def runFiducialRegistration(self, outputTransform, sourceCoords, targetCoords):
        auxSourceNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
        auxSourceNode.GetDisplayNode().SetVisibility(False)
        for coord in sourceCoords:
            auxSourceNode.AddFiducialFromArray(coord)

        auxTargetNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
        auxTargetNode.GetDisplayNode().SetVisibility(False)
        for coord in targetCoords:
            auxTargetNode.AddFiducialFromArray(coord)
            
        parameters = {}
        parameters['fixedLandmarks']  = auxTargetNode.GetID()
        parameters['movingLandmarks'] = auxSourceNode.GetID()
        parameters['saveTransform']   = outputTransform.GetID()
        parameters['transformType']   = 'Rigid'
        slicer.cli.run(slicer.modules.fiducialregistration, None, parameters, wait_for_completion=True, update_display=False)

        slicer.mrmlScene.RemoveNode(auxTargetNode)
        slicer.mrmlScene.RemoveNode(auxSourceNode)


    def runACPCAlignment(self, outputTransform, AC, PC, MS):

        auxLineNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsLineNode')
        auxLineNode.GetDisplayNode().SetVisibility(False)
        auxLineNode.AddControlPointWorld(AC, 'AC')
        auxLineNode.AddControlPointWorld(PC, 'PC')

        ACPCMSNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
        ACPCMSNode.GetDisplayNode().SetVisibility(False)
        ACPCMSNode.AddFiducialFromArray(AC)
        ACPCMSNode.AddFiducialFromArray(PC)
        ACPCMSNode.AddFiducialFromArray(MS)

        parameters = {}
        parameters['ACPC']  = auxLineNode.GetID()
        parameters['Midline'] = ACPCMSNode.GetID()
        parameters['centerVolume'] = True 
        parameters['OutputTransform'] = outputTransform.GetID()
        slicer.cli.run(slicer.modules.acpctransform, None, parameters, wait_for_completion=True, update_display=False)

        slicer.mrmlScene.RemoveNode(auxLineNode)
        slicer.mrmlScene.RemoveNode(ACPCMSNode)

    def computeTrajectoryFromTargetMountingRingArc(self, outputTransform, frameTargetCoordinates, mounting, ringAngle, arcAngle):
        """
        Run the processing algorithm.
        Can be used without GUI widget.
        """

        if not outputTransform:
            raise ValueError("output transform is invalid")

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
        vtkTransform.Translate(frameTargetCoordinates)
        vtkTransform.RotateWXYZ(arcAngle, arcDirection[0], arcDirection[1], arcDirection[2])
        vtkTransform.RotateWXYZ(ringAngle, ringDirection[0], ringDirection[1], ringDirection[2])
        vtkTransform.RotateWXYZ(90, initDirection[0], initDirection[1], initDirection[2])

        # Set to node
        outputTransform.SetAndObserveTransformToParent(vtkTransform)

    def computeTrajectoryFromTargetEntryRoll(self, outputTransform, frameTargetCoordinates, frameEntryCoordinates, rollAngle):

        entryTargetDirection = frameEntryCoordinates - frameTargetCoordinates
        vtk.vtkMath().Normalize(entryTargetDirection)
        superiorInferiorDirection = np.array([0,0,1])

        ang_rad = np.arccos(vtk.vtkMath().Dot(entryTargetDirection, superiorInferiorDirection))
        ang_deg = np.rad2deg(ang_rad)

        cross = np.zeros(3)
        vtk.vtkMath().Cross(entryTargetDirection, superiorInferiorDirection, cross)

        if vtk.vtkMath().Dot(cross,superiorInferiorDirection) >= 0:
            ang_deg = -1 * ang_deg
        
        vtkTransform = vtk.vtkTransform()
        vtkTransform.Translate(frameTargetCoordinates)
        vtkTransform.RotateWXYZ(rollAngle, entryTargetDirection[0], entryTargetDirection[1], entryTargetDirection[2])
        vtkTransform.RotateWXYZ(ang_deg, cross[0], cross[1], cross[2])

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
        self.test_BrainlabImport()
    
    def test_BrainlabImport(self):
        import StereotacticPlanLib.ImportFrom.Brainlab as bl
        logic = StereotacticPlanLogic()
        parameterNode = logic.getParameterNode()
        reportFilePath = "C:\\Users\\simon\\Desktop\\test\\StereotaxyReport.pdf"
        bl.setParameterNodeFromDevice(parameterNode, reportFilePath)
        volumeFilePath = "C:\\Users\\simon\\Desktop\\test\\anat_t1.nii"
        volumeNode = slicer.util.loadVolume(volumeFilePath)
        parameterNode.SetNodeReferenceID("ReferenceVolume", volumeNode.GetID())

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

        # Get/create input data

        import SampleData
        registerSampleData()
        inputVolume = SampleData.downloadSample('StereotacticPlan1')
        self.delayDisplay('Loaded test data set')

        inputScalarRange = inputVolume.GetImageData().GetScalarRange()
        self.assertEqual(inputScalarRange[0], 0)
        self.assertEqual(inputScalarRange[1], 695)

        outputVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
        threshold = 100

        # Test the module logic

        logic = StereotacticPlanLogic()

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
