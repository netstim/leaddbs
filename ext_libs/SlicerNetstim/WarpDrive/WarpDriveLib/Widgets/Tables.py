import os
import qt, vtk, slicer
from PythonQt import BoolResult
from slicer.util import VTKObservationMixin
import WarpDrive, ImportAtlas
import numpy as np

class TextEditDelegate(qt.QItemDelegate):
  def __init__(self, parent, renameControlPointsFunction):
    qt.QItemDelegate.__init__(self, parent)
    self.renameControlPointsFunction = renameControlPointsFunction

  def createEditor(self, parent, option, index):
    lineEdit = qt.QLineEdit(parent)
    return lineEdit

  def setEditorData(self, editor, index):
    editor.blockSignals(True)
    editor.text = index.model().data(index) if index.model().data(index) else ""
    editor.blockSignals(False)
  
  def setModelData(self, editor, model, index):
    previousName = index.model().data(index)
    newName = editor.text
    model.setData(index, newName)
    self.renameControlPointsFunction(previousName, newName)

class SpinBoxDelegate(qt.QItemDelegate):
  def __init__(self, parent, updateRadiusFunction):
    qt.QItemDelegate.__init__(self, parent)
    self.updateRadiusFunction = updateRadiusFunction

  def createEditor(self, parent, option, index):
    spinBox = qt.QDoubleSpinBox(parent)
    spinBox.setSingleStep(1)
    spinBox.maximum = 50
    spinBox.minimum = 5
    return spinBox

  def setEditorData(self, editor, index):
    editor.blockSignals(True)
    editor.value = float(index.model().data(index)) if index.model().data(index) else 30
    editor.blockSignals(False)

  def setModelData(self, editor, model, index):
    newValue = editor.value
    model.setData(index, newValue)
    controlPointName = index.model().data(index.siblingAtColumn(1))
    self.updateRadiusFunction(controlPointName, newValue)

class firstColumnCheckableModel(qt.QStandardItemModel):
  def __init__(self , *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.updateSelectedFuntion = None

  def flags(self, index):
    baseFlags = qt.Qt.ItemIsEnabled | qt.Qt.ItemIsSelectable
    if index.column() == 0:
      return baseFlags | qt.Qt.ItemIsUserCheckable
    else:
      return baseFlags | qt.Qt.ItemIsEditable

  def setData(self , *args, **kwargs):
    index = args[0] if args else None
    if isinstance(index, qt.QModelIndex) and args[-1] == qt.Qt.CheckStateRole:
      controlPointName = index.model().data(index.siblingAtColumn(1))
      self.updateSelectedFuntion(controlPointName, args[1])
    qt.QStandardItemModel.setData(self , *args, **kwargs)


class baseTable(qt.QWidget):

  def __init__(self):
    super().__init__()

    effectPixmap = qt.QPixmap(os.path.join(os.path.split(WarpDrive.__file__)[0], 'Resources', 'Icons', 'Add.png'))
    effectIcon = qt.QIcon(effectPixmap)
    self.addButton = qt.QToolButton()
    self.addButton.setToolButtonStyle(qt.Qt.ToolButtonTextUnderIcon)
    self.addButton.setIcon(effectIcon)
    self.addButton.setIconSize(effectPixmap.rect().size())
    self.addButton.setEnabled(True)
    self.addButton.setSizePolicy(qt.QSizePolicy.MinimumExpanding,qt.QSizePolicy.Maximum)
    self.addButton.clicked.connect(self.onAddButton)

    effectPixmap = qt.QPixmap(os.path.join(os.path.split(WarpDrive.__file__)[0], 'Resources', 'Icons', 'Delete.png'))
    effectIcon = qt.QIcon(effectPixmap)
    self.removeButton = qt.QToolButton()
    self.removeButton.setToolButtonStyle(qt.Qt.ToolButtonTextUnderIcon)
    self.removeButton.setIcon(effectIcon)
    self.removeButton.setText('Delete')
    self.removeButton.setIconSize(effectPixmap.rect().size())
    self.removeButton.setEnabled(True)
    self.removeButton.setSizePolicy(qt.QSizePolicy.MinimumExpanding,qt.QSizePolicy.Maximum)
    self.removeButton.clicked.connect(self.onRemoveButton)

    self.buttonsFrame = qt.QFrame()
    self.buttonsFrame.setSizePolicy(qt.QSizePolicy.Preferred, qt.QSizePolicy.Minimum)
    self.buttonsFrame.setLayout(qt.QHBoxLayout())
    self.buttonsFrame.layout().addWidget(self.addButton,1)
    self.buttonsFrame.layout().addWidget(self.removeButton,1)

    self.parentLayout = qt.QVBoxLayout(self)
    self.parentLayout.addWidget(self.buttonsFrame)
    self.parentLayout.addWidget(self.view)
  
  def onAddButton(self):
    pass

  def onRemoveButton(self):
    pass

class AtlasesTable(baseTable):

  def __init__(self):

    self.view = slicer.qMRMLSubjectHierarchyTreeView(slicer.util.mainWindow())
    self.view.setMRMLScene(slicer.mrmlScene)
    self.view.contextMenuEnabled = False
    self.view.setEditTriggers(0) # disable double click to edit
    self.view.nodeTypes = ('vtkMRMLModelNode','vtkMRMLFolderDisplayNode')
    self.view.attributeNameFilter = ('atlas')
    for key in ['idColumn', 'transformColumn', 'descriptionColumn']:
      self.view.setColumnHidden(eval('self.view.model().'+key), True)

    self.view.doubleClicked.connect(self.onDoubleClick)

    super().__init__()

    self.addButton.setText('Atlas')
    self.addButton.setToolTip('Add atlas')

    self.buttonsFrame.layout().addStretch(2)

    self.updateTable()

  def updateTable(self):
    currentValue = self.view.attributeNameFilter
    self.view.attributeNameFilter = ('')
    self.view.attributeNameFilter = currentValue

  def onAddButton(self):
    leadDBSPath = slicer.util.settingsValue("NetstimPreferences/leadDBSPath", "", converter=str)
    if leadDBSPath == "":
      qt.QMessageBox().warning(qt.QWidget(), "", "Add Lead-DBS path to Slicer preferences")
      return
    validAtlasesNames = ImportAtlas.ImportAtlasLogic().getValidAtlases()
    if not validAtlasesNames:
      return
    result = BoolResult()
    atlasName = qt.QInputDialog.getItem(qt.QWidget(),'Select Atlas','',validAtlasesNames,0,0,result) 
    if result:
      qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
      qt.QApplication.processEvents()
      try:
        ImportAtlas.ImportAtlasLogic().readAtlas(os.path.join(ImportAtlas.ImportAtlasLogic().getAtlasesPath(), atlasName, 'atlas_index.mat'))    
      finally:
        qt.QApplication.restoreOverrideCursor()
    tb = next(filter(lambda x: isinstance(x,qt.QToolBar) and x.windowTitle=='LeadDBS', slicer.util.mainWindow().children()))
    tb.invertAtlases(WarpDrive.WarpDriveLogic().getParameterNode().GetNodeReference("InputNode"))
    self.updateTable()

  def onRemoveButton(self):
    nodeID = self.view.currentItem()
    if not nodeID:
      return
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    removeIDs = vtk.vtkIdList()
    shNode.GetItemChildren(nodeID, removeIDs, True)
    removeIDs.InsertNextId(nodeID)
    for i in range(removeIDs.GetNumberOfIds()):
      shNode.RemoveItem(removeIDs.GetId(i))

  def onDoubleClick(self):
    nodeID = self.view.currentItem()
    if not nodeID:
      return
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    node = shNode.GetItemDataNode(nodeID)
    if not isinstance(node, slicer.vtkMRMLModelNode):
      return
    pd = node.GetPolyData()
    center = vtk.vtkCenterOfMass()
    center.SetInputData(pd)
    center.Update()
    centerList = center.GetCenter()
    # create markups node, add center as fiducial and jump and center slices
    markupsNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    markupsNode.GetDisplayNode().SetVisibility(False)
    markupsNode.AddFiducialFromArray(np.array(centerList),'')
    markupsLogic = slicer.modules.markups.logic()
    markupsLogic.JumpSlicesToNthPointInMarkup(markupsNode.GetID(),0,True)
    slicer.mrmlScene.RemoveNode(markupsNode)

class WarpDriveCorrectionsTable(baseTable):


  def __init__(self):

    columnNames = ["Include", "Name", "Radius"]
    self.model = firstColumnCheckableModel(1, len(columnNames))
    self.model.updateSelectedFuntion = self.updateSelected
    for i, columnName in enumerate(columnNames):
      self.model.setHeaderData(i, qt.Qt.Horizontal, columnName)

    self.view = qt.QTableView()
    self.view.setEditTriggers(self.view.DoubleClicked)
    self.view.setSelectionMode(self.view.SingleSelection)
    self.view.setSelectionBehavior(self.view.SelectRows)
    self.view.horizontalHeader().setStretchLastSection(True)
    self.view.verticalHeader().setDefaultSectionSize(self.view.verticalHeader().defaultSectionSize*3/4)
    self.view.setHorizontalScrollMode(self.view.ScrollPerPixel)
    self.view.setModel(self.model)
    self.view.clicked.connect(self.onSelectionChanged)

    self.view.setItemDelegateForColumn(1, TextEditDelegate(self.model, self.renameControlPoints))
    self.view.setItemDelegateForColumn(2, SpinBoxDelegate(self.model, self.updateRadius))

    super().__init__()

    effectPixmap = qt.QPixmap(os.path.join(os.path.split(WarpDrive.__file__)[0], 'Resources', 'Icons', 'SlicerVisible.png'))
    effectIcon = qt.QIcon(effectPixmap)
    self.sourceVisibleButton = qt.QToolButton()
    self.sourceVisibleButton.setToolButtonStyle(qt.Qt.ToolButtonTextUnderIcon)
    self.sourceVisibleButton.setIcon(effectIcon)
    self.sourceVisibleButton.setIconSize(effectPixmap.rect().size())
    self.sourceVisibleButton.setText('Source')
    self.sourceVisibleButton.setToolTip('Toggle source fiducials visibility')
    self.sourceVisibleButton.setCheckable(True)
    self.sourceVisibleButton.setChecked(False)
    self.sourceVisibleButton.setEnabled(True)
    self.sourceVisibleButton.setSizePolicy(qt.QSizePolicy.MinimumExpanding,qt.QSizePolicy.Maximum)
    self.sourceVisibleButton.toggled.connect(self.onSourceVisibleToggled)
    
    effectPixmap = qt.QPixmap(os.path.join(os.path.split(WarpDrive.__file__)[0], 'Resources', 'Icons', 'SlicerVisible.png'))
    effectIcon = qt.QIcon(effectPixmap)
    self.targetVisibleButton = qt.QToolButton()
    self.targetVisibleButton.setToolButtonStyle(qt.Qt.ToolButtonTextUnderIcon)
    self.targetVisibleButton.setIcon(effectIcon)
    self.targetVisibleButton.setIconSize(effectPixmap.rect().size())
    self.targetVisibleButton.setText('Target')
    self.targetVisibleButton.setToolTip('Toggle target fiducials visibility')    
    self.targetVisibleButton.setCheckable(True)
    self.targetVisibleButton.setChecked(False)
    self.targetVisibleButton.setEnabled(True)
    self.targetVisibleButton.setSizePolicy(qt.QSizePolicy.MinimumExpanding,qt.QSizePolicy.Maximum)
    self.targetVisibleButton.toggled.connect(self.onTargetVisibleToggled)

    effectPixmap = qt.QPixmap(os.path.join(os.path.split(WarpDrive.__file__)[0], 'Resources', 'Icons', 'SlicerUndo.png'))
    effectIcon = qt.QIcon(effectPixmap)
    self.undoButton = qt.QToolButton()
    self.undoButton.setToolButtonStyle(qt.Qt.ToolButtonTextUnderIcon)
    self.undoButton.setIcon(effectIcon)
    self.undoButton.setIconSize(effectPixmap.rect().size())
    self.undoButton.setText('Undo')
    self.undoButton.setToolTip('Undo last correction')    
    self.undoButton.setCheckable(False)
    self.undoButton.setChecked(False)
    self.undoButton.setEnabled(True)
    self.undoButton.setSizePolicy(qt.QSizePolicy.MinimumExpanding,qt.QSizePolicy.Maximum)
    self.undoButton.clicked.connect(self.onUndoClicked)

    self.buttonsFrame.layout().addWidget(self.sourceVisibleButton,1)
    self.buttonsFrame.layout().addWidget(self.targetVisibleButton,1)
    self.buttonsFrame.layout().addWidget(self.undoButton,1)

    self.addButton.setText('Fixed point')
    self.addButton.setToolTip('Add fixed point')

    self.previewSelectedCheckBox = qt.QCheckBox('Preview selected correction')
    self.previewSelectedCheckBox.setChecked(True)
    self.parentLayout.addWidget(self.previewSelectedCheckBox)
  
  def onSourceVisibleToggled(self):
    pass

  def onTargetVisibleToggled(self):
    pass

  def onSelectionChanged(self):
    pass

  def onUndoClicked(self):
    pass

  def clearTable(self):
    while self.model.rowCount() > 0:
      self.model.removeRow(self.model.rowCount()-1)

  def getSelectedRow(self):
    selectedRows = self.view.selectionModel().selectedRows()
    for selectedRow in selectedRows:
      return selectedRow.row() # is a single selection view


class WarpDriveCorrectionsManager(VTKObservationMixin, WarpDriveCorrectionsTable):
  def __init__(self):
    VTKObservationMixin.__init__(self)
    WarpDriveCorrectionsTable.__init__(self)
    self.sourceFiducialNodeID = ""
    self.targetFiducialNodeID = ""
    self._updatingFiducials = False
    self.targetFiducialObservers = []
    self.sourceFiducialObservers = []
    self.parameterNode = WarpDrive.WarpDriveLogic().getParameterNode()
    self.addObserver(self.parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateNodesListeners)

  def onSourceVisibleToggled(self):
    if self.sourceFiducialNodeID != "":
      sourceFiducialNode = slicer.mrmlScene.GetNodeByID(self.sourceFiducialNodeID)
      sourceFiducialNode.GetDisplayNode().SetVisibility(self.sourceVisibleButton.checked)

  def onTargetVisibleToggled(self):
    if self.targetFiducialNodeID != "":
      targetFiducialNode = slicer.mrmlScene.GetNodeByID(self.targetFiducialNodeID)
      targetFiducialNode.GetDisplayNode().SetVisibility(self.targetVisibleButton.checked)

  def onAddButton(self):
    if self.targetFiducialNodeID == "":
      return
    interactionNode = slicer.app.applicationLogic().GetInteractionNode()
    selectionNode = slicer.app.applicationLogic().GetSelectionNode()
    selectionNode.SetReferenceActivePlaceNodeClassName("vtkMRMLMarkupsFiducialNode")
    selectionNode.SetActivePlaceNodeID(self.targetFiducialNodeID)
    interactionNode.SetCurrentInteractionMode(interactionNode.Place)

  def updateNodesListeners(self, caller, event):
    sourceFiducialNode = self.parameterNode.GetNodeReference("SourceFiducial")
    if sourceFiducialNode and self.sourceFiducialNodeID != sourceFiducialNode.GetID():
      previousNode = slicer.mrmlScene.GetNodeByID(self.sourceFiducialNodeID)
      if previousNode and bool(self.sourceFiducialObservers):
        for obs in self.sourceFiducialObservers:
          previousNode.RemoveObserver(obs)
      sourceFiducialNode.GetDisplayNode().SetPointLabelsVisibility(0)
      sourceFiducialNode.GetDisplayNode().SetVisibility(0)
      self.sourceFiducialNodeID = sourceFiducialNode.GetID()
      self.sourceFiducialObservers.clear()
      self.sourceFiducialObservers.append(sourceFiducialNode.AddObserver(slicer.vtkMRMLDisplayableNode.DisplayModifiedEvent, self.updateVisibilityWidget))
    
    targetFiducialNode = self.parameterNode.GetNodeReference("TargetFiducial")
    if targetFiducialNode and self.targetFiducialNodeID != targetFiducialNode.GetID():
      previousNode = slicer.mrmlScene.GetNodeByID(self.targetFiducialNodeID)
      if previousNode and bool(self.targetFiducialObservers):
        for obs in self.targetFiducialObservers:
          previousNode.RemoveObserver(obs)
      targetFiducialNode.GetDisplayNode().SetPointLabelsVisibility(0)
      self.targetFiducialNodeID = targetFiducialNode.GetID()
      self.targetFiducialObservers.clear()
      self.targetFiducialObservers.append(targetFiducialNode.AddObserver(slicer.vtkMRMLDisplayableNode.DisplayModifiedEvent, self.updateVisibilityWidget))
      self.targetFiducialObservers.append(targetFiducialNode.AddObserver(targetFiducialNode.PointAddedEvent, self.targetFiducialModified))
      self.targetFiducialObservers.append(targetFiducialNode.AddObserver(targetFiducialNode.PointRemovedEvent, self.targetFiducialModified))
      self.targetFiducialObservers.append(targetFiducialNode.AddObserver(targetFiducialNode.PointModifiedEvent, self.targetFiducialModified))
      self.targetFiducialObservers.append(targetFiducialNode.AddObserver(targetFiducialNode.PointPositionDefinedEvent, self.onPointPositionDefined))
      self.setUpWidget()
      self.updateVisibilityWidget()

  def updateVisibilityWidget(self, caller=None, event=None):
    if self.targetFiducialNodeID != "":
      targetFiducialNode = slicer.mrmlScene.GetNodeByID(self.targetFiducialNodeID)
      if targetFiducialNode.GetDisplayNode():
        self.targetVisibleButton.checked = targetFiducialNode.GetDisplayNode().GetVisibility()
    if self.sourceFiducialNodeID != "":
      sourceFiducialNode = slicer.mrmlScene.GetNodeByID(self.sourceFiducialNodeID)
      if sourceFiducialNode.GetDisplayNode():
        self.sourceVisibleButton.checked = sourceFiducialNode.GetDisplayNode().GetVisibility()

  def targetFiducialModified(self, caller, event):
    self.setUpWidget()

  def setUpWidget(self):
    if self._updatingFiducials:
      return
    self.clearTable()
    targetFiducialNode = slicer.mrmlScene.GetNodeByID(self.targetFiducialNodeID)
    newCorrection = {"include": False, "name": "", "radius":""}
    for i in range(targetFiducialNode.GetNumberOfControlPoints()):
      if newCorrection["name"] == targetFiducialNode.GetNthControlPointLabel(i):
        continue
      else:
        self.addCorrectionToWidget(newCorrection)
        newCorrection["name"] = targetFiducialNode.GetNthControlPointLabel(i)
        newCorrection["radius"] = targetFiducialNode.GetNthControlPointDescription(i)
        newCorrection["include"] = targetFiducialNode.GetNthControlPointSelected(i)
    self.addCorrectionToWidget(newCorrection)

  def addCorrectionToWidget(self, newCorrection):
    if newCorrection["name"] == "":
      return

    row = self.model.rowCount()
    self.model.insertRow(row)
    for col,val in enumerate(newCorrection.values()):
      index = self.model.index(row, col)
      if col == 0:
        val = qt.Qt.Checked if val else qt.Qt.Unchecked
        role = qt.Qt.CheckStateRole
      else:
        role = qt.Qt.DisplayRole

      self.model.setData(index, val, role)

  def onPointPositionDefined(self, caller, event):
    targetFiducialNode = slicer.mrmlScene.GetNodeByID(self.targetFiducialNodeID)
    lastControlPoint = targetFiducialNode.GetNumberOfControlPoints()-1
    name = targetFiducialNode.GetNthControlPointLabel(lastControlPoint)
    if name.startswith(targetFiducialNode.GetName()):
      sourceFiducialNode = slicer.mrmlScene.GetNodeByID(self.sourceFiducialNodeID)
      sourceFiducialNode.AddControlPoint(vtk.vtkVector3d(targetFiducialNode.GetNthControlPointPosition(lastControlPoint)))
      targetFiducialNode.SetNthControlPointLabel(lastControlPoint, slicer.mrmlScene.GenerateUniqueName('fixed point'))
      targetFiducialNode.SetNthControlPointDescription(lastControlPoint, self.parameterNode.GetParameter("Radius"))
      sourceFiducialNode.SetNthControlPointLabel(lastControlPoint, targetFiducialNode.GetNthControlPointLabel(lastControlPoint))
      sourceFiducialNode.SetNthControlPointDescription(lastControlPoint, self.parameterNode.GetParameter("Radius"))
    self.parameterNode.SetParameter("Update","true")

  def getSelectedCorrectionName(self):
    row = self.getSelectedRow()
    if row is None or self.targetFiducialNodeID=="":
      return
    index = self.model.index(row, 1)
    return self.model.itemData(index)[0]

  def onRemoveButton(self):
    correctionName = self.getSelectedCorrectionName()
    if correctionName is None:
      return
    self.removeCorrectionByName(correctionName)

  def onUndoClicked(self):
    if self.model.rowCount() == 0:
      return
    index = self.model.index(self.model.rowCount()-1, 1)
    correctionName = self.model.itemData(index)[0]
    self.removeCorrectionByName(correctionName)

  def removeCorrectionByName(self, correctionName):
    targetFiducialNode = slicer.mrmlScene.GetNodeByID(self.targetFiducialNodeID)
    sourceFiducialNode = slicer.mrmlScene.GetNodeByID(self.sourceFiducialNodeID)
    for i in range(targetFiducialNode.GetNumberOfControlPoints()-1,-1,-1):
      if targetFiducialNode.GetNthControlPointLabel(i) == correctionName:
        targetFiducialNode.RemoveNthControlPoint(i)
        sourceFiducialNode.RemoveNthControlPoint(i)
    self.parameterNode.SetParameter("Update","true")

  def renameControlPoints(self, previousName, newName):
    if self.targetFiducialNodeID == "":
      return
    self._updatingFiducials = True
    targetFiducialNode = slicer.mrmlScene.GetNodeByID(self.targetFiducialNodeID)
    sourceFiducialNode = slicer.mrmlScene.GetNodeByID(self.sourceFiducialNodeID)
    for i in range(targetFiducialNode.GetNumberOfControlPoints()):
      if targetFiducialNode.GetNthControlPointLabel(i) == previousName:
        targetFiducialNode.SetNthControlPointLabel(i, newName)
        sourceFiducialNode.SetNthControlPointLabel(i, newName)
    self._updatingFiducials = False

  def updateRadius(self, controlPointName, value):
    if self.targetFiducialNodeID == "":
      return
    self._updatingFiducials = True
    targetFiducialNode = slicer.mrmlScene.GetNodeByID(self.targetFiducialNodeID)
    sourceFiducialNode = slicer.mrmlScene.GetNodeByID(self.sourceFiducialNodeID)
    for i in range(targetFiducialNode.GetNumberOfControlPoints()):
      if targetFiducialNode.GetNthControlPointLabel(i) == controlPointName:
        targetFiducialNode.SetNthControlPointDescription(i, "%.01f"%value)
        sourceFiducialNode.SetNthControlPointDescription(i, "%.01f"%value)
        self.parameterNode.SetParameter("Update","true")
    self._updatingFiducials = False

  def updateSelected(self, controlPointName, value):
    if self.targetFiducialNodeID == "":
      return    
    self._updatingFiducials = True
    targetFiducialNode = slicer.mrmlScene.GetNodeByID(self.targetFiducialNodeID)
    sourceFiducialNode = slicer.mrmlScene.GetNodeByID(self.sourceFiducialNodeID)
    for i in range(targetFiducialNode.GetNumberOfControlPoints()):
      if targetFiducialNode.GetNthControlPointLabel(i) == controlPointName:
        targetFiducialNode.SetNthControlPointSelected(i, value)
        sourceFiducialNode.SetNthControlPointSelected(i, value)
        self.parameterNode.SetParameter("Update","true")
    self._updatingFiducials = False

  def onSelectionChanged(self):
    if not self.previewSelectedCheckBox.checked:
      return
    correctionName = self.getSelectedCorrectionName()
    if correctionName is None:
      return
    jumped = False
    sourcePoints = vtk.vtkPoints()
    targetPoints = vtk.vtkPoints()
    targetFiducialNode = slicer.mrmlScene.GetNodeByID(self.targetFiducialNodeID)
    sourceFiducialNode = slicer.mrmlScene.GetNodeByID(self.sourceFiducialNodeID)
    for i in range(targetFiducialNode.GetNumberOfControlPoints()):
      if targetFiducialNode.GetNthControlPointLabel(i) == correctionName:
        if not jumped:
          markupsLogic = slicer.modules.markups.logic()
          markupsLogic.JumpSlicesToNthPointInMarkup(self.targetFiducialNodeID,i,False)
          jumped = True
        sourcePoints.InsertNextPoint(sourceFiducialNode.GetNthControlPointPosition(i))
        targetPoints.InsertNextPoint(targetFiducialNode.GetNthControlPointPosition(i))        
    tmpNodes = WarpDrive.WarpDriveLogic().previewWarp(sourcePoints, targetPoints)
    for n in tmpNodes:
      qt.QTimer.singleShot(1000, lambda node=n: slicer.mrmlScene.RemoveNode(node))
