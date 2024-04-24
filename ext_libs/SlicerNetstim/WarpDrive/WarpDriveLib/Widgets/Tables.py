import os
import qt, vtk, slicer
from PythonQt import BoolResult
from slicer.util import VTKObservationMixin
import WarpDrive, ImportAtlas
from ..Helpers import GridNodeHelper
import numpy as np
import glob
import json

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
    spinBox.minimum = 1
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

class AtlasSegmentationBaseTable(baseTable):

  def __init__(self):

    self.view.setMRMLScene(slicer.mrmlScene)
    self.view.contextMenuEnabled = False
    self.view.setEditTriggers(0) # disable double click to edit
    self.view.nodeTypes = ('vtkMRMLModelNode','vtkMRMLFolderDisplayNode')
    for key in ['idColumn', 'transformColumn', 'descriptionColumn']:
      self.view.setColumnHidden(eval('self.view.model().'+key), True)

    self.view.doubleClicked.connect(self.onDoubleClick)

    super().__init__()

    self.saveButton = qt.QToolButton()
    self.saveButton.setVisible(False)
    self.saveButton.setToolButtonStyle(qt.Qt.ToolButtonTextUnderIcon)
    self.saveButton.setText('Save')
    self.saveButton.setIcon(qt.QIcon(":/Icons/Small/SlicerSave.png"))
    self.saveButton.setIconSize(self.removeButton.iconSize)
    self.saveButton.setSizePolicy(qt.QSizePolicy.MinimumExpanding,qt.QSizePolicy.Maximum)
    self.saveButton.clicked.connect(self.onSaveClicked)

    self.buttonsFrame.layout().addWidget(self.saveButton,1)

  def onSaveClicked(self):
    currentSubjectPath = os.path.dirname(json.loads(WarpDrive.WarpDriveLogic().getParameterNode().GetParameter("CurrentSubject"))['warpdrive_path'])
    outputDir = os.path.join(currentSubjectPath, 'segmentations') if int(WarpDrive.WarpDriveLogic().getParameterNode().GetParameter("InverseMode")) else os.path.join(currentSubjectPath, 'segmentations_mni')
    qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
    qt.QApplication.processEvents()
    self.saveSegmentation(outputDir)
    qt.QApplication.restoreOverrideCursor()
  
  def saveSegmentation(self, segmentationPath):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    modelNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLModelNode')
    modelNodes.UnRegister(slicer.mrmlScene)
    for i in range(modelNodes.GetNumberOfItems()):
      modelNode = modelNodes.GetItemAsObject(i)
      modelItem = shNode.GetItemByDataNode(modelNode)
      if modelNode.GetDisplayNode().GetVisibility() and self.view.attributeNameFilter in shNode.GetItemAttributeNames(modelItem):
        # construct name
        name = modelNode.GetName()
        node = modelNode
        while shNode.GetItemParent(shNode.GetItemByDataNode(node)) != shNode.GetSceneItemID():
          node = shNode.GetItemDataNode(shNode.GetItemParent(shNode.GetItemByDataNode(node)))
          name = os.path.join(node.GetName(), name)
        # save segmentation
        clonedItemID = slicer.modules.subjecthierarchy.logic().CloneSubjectHierarchyItem(shNode, shNode.GetItemByDataNode(modelNode))
        tmpNode = shNode.GetItemDataNode(clonedItemID)
        tmpNode.SetAndObserveTransformNodeID(modelNode.GetTransformNodeID())
        tmpNode.HardenTransform()
        labelNode = self.modelToLabel(modelNode)
        if os.path.dirname(name):
          os.makedirs(os.path.dirname(name), exist_ok=True)
        slicer.util.saveNode(labelNode, os.path.join(segmentationPath, name + '.nii.gz'))      
        slicer.mrmlScene.RemoveNode(labelNode)
        slicer.mrmlScene.RemoveNode(tmpNode)

  def modelToLabel(self, model_node):
    segmentationsLogic = slicer.modules.segmentations.logic()
    segmentNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentationNode', 'segmentation')
    labelNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLabelMapVolumeNode', 'label')
    segmentationsLogic.ImportModelToSegmentationNode(model_node, segmentNode)
    segmentationsLogic.ExportAllSegmentsToLabelmapNode(segmentNode, labelNode)
    slicer.mrmlScene.RemoveNode(segmentNode)
    # pad image
    extend_mm = 5
    dimensions, origin, spacing, directionMatrix = GridNodeHelper.getGridDefinition(labelNode)
    origin = [x - extend_mm for x in origin]
    dimensions = [int(d + (2*extend_mm) / s) for (d,s) in zip(dimensions, spacing)]
    auxVolumeNode = GridNodeHelper.emptyVolume(dimensions, origin, spacing, directionMatrix)
    parameters = {}
    parameters['inputVolume'] = labelNode.GetID()
    parameters['referenceVolume'] =  auxVolumeNode.GetID()
    parameters['outputVolume'] = labelNode.GetID()
    parameters['pixelType'] = 'uchar'
    parameters['interpolationMode'] = 'NearestNeighbor'
    slicer.cli.run(slicer.modules.brainsresample, None, parameters, wait_for_completion=True, update_display=False)
    slicer.mrmlScene.RemoveNode(auxVolumeNode)
    return labelNode

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
    markupsNode.AddControlPoint(np.array(centerList),'')
    markupsLogic = slicer.modules.markups.logic()
    markupsLogic.JumpSlicesToNthPointInMarkup(markupsNode.GetID(),0,True)
    slicer.mrmlScene.RemoveNode(markupsNode)

  def updateTable(self):
    currentValue = self.view.attributeNameFilter
    self.view.attributeNameFilter = ('')
    self.view.attributeNameFilter = currentValue

class SegmentationsTable(AtlasSegmentationBaseTable):

  def __init__(self):

    self.view = slicer.qMRMLSubjectHierarchyTreeView(slicer.util.mainWindow())
    self.view.attributeNameFilter = ('segmentation')

    super().__init__()

    self.addButton.setText('Segmentation')
    self.addButton.setToolTip('Add segmentation')

    self.buttonsFrame.layout().addStretch(2)

    self.updateTable()

  def onAddButton(self):
    currentSubject = json.loads(WarpDrive.WarpDriveLogic().getParameterNode().GetParameter("CurrentSubject"))
    segmentationsPath = os.path.join(os.path.dirname(currentSubject['warpdrive_path']), 'segmentations')
    subFolders = [os.path.basename(x) for x in sorted(glob.glob(os.path.join(segmentationsPath, '*')))]
    result = BoolResult()
    segmentationName = qt.QInputDialog.getItem(qt.QWidget(),'Select Atlas','',subFolders,0,0,result) 
    if result:
      selection = os.path.join(segmentationsPath, segmentationName)
      segmentationFiles = [selection] if os.path.isfile(selection) else glob.glob(os.path.join(selection, '**', '*.nii*'), recursive=True)
      qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
      qt.QApplication.processEvents()
      try:
        self.loadSegmentations(segmentationsPath, segmentationFiles)    
      finally:
        qt.QApplication.restoreOverrideCursor()
    self.updateTable()

  def loadSegmentations(self, segmentationPath, niiFiles):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    for file in niiFiles:
      # load
      volumeNode = slicer.util.loadVolume(file,{'show':False})
      modelNode = self.volumeToModel(volumeNode)
      slicer.mrmlScene.RemoveNode(volumeNode)
      modelNode.GetDisplayNode().SetVisibility2D(1)
      # set in hierarchy
      basePath,fileName = os.path.split(file)
      parentFolders = os.path.relpath(basePath, segmentationPath).split(os.sep)
      if '.' in parentFolders:
        parentFolders.remove('.')
      parentID = shNode.GetSceneItemID()
      while(parentFolders):
        parentID = self.getOrCreateFolderItem(parentID, parentFolders.pop(0))
      shNode.SetItemParent(shNode.GetItemByDataNode(modelNode), parentID)
      shNode.SetItemAttribute(shNode.GetItemByDataNode(modelNode), 'segmentation', '1')
      modelNode.SetName(fileName.split('.')[0])
    WarpDrive.WarpDriveLogic().invertAtlases(WarpDrive.WarpDriveLogic().getParameterNode().GetNodeReference("InputNode"), int(WarpDrive.WarpDriveLogic().getParameterNode().GetParameter("InverseMode")))

  def volumeToModel(self, volumeNode):
    labelNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLabelMapVolumeNode')
    volumeArray = slicer.util.array(volumeNode.GetID())
    volumeArray[:] = volumeArray > 0.5
    volumeNode.Modified()
    slicer.modules.volumes.logic().CreateLabelVolumeFromVolume(slicer.mrmlScene, labelNode, volumeNode)
    segmentationsLogic = slicer.modules.segmentations.logic()
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    segmentNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentationNode', 'segmentation')
    segmentationsLogic.ImportLabelmapToSegmentationNode(labelNode, segmentNode)
    itemID = shNode.CreateFolderItem(shNode.GetSceneItemID(), 'aux')
    segmentationsLogic.ExportAllSegmentsToModels(segmentNode, itemID)
    childID = []
    shNode.GetItemChildren(itemID, childID, False)
    shNode.SetItemParent(childID[0], shNode.GetSceneItemID())
    shNode.RemoveItem(itemID)
    slicer.mrmlScene.RemoveNode(segmentNode)
    slicer.mrmlScene.RemoveNode(labelNode)
    return shNode.GetItemDataNode(childID[0])

  def createFolderDisplayNode(self, folderID, color=[0.66,0.66,0.66], opacity=1.0):
    # from qSlicerSubjectHierarchyFolderPlugin.cxx
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    displayNode = slicer.vtkMRMLFolderDisplayNode()
    displayNode.SetName(shNode.GetItemName(folderID))
    displayNode.SetHideFromEditors(0)
    displayNode.SetAttribute('SubjectHierarchy.Folder', "1")
    displayNode.SetColor(*color)
    displayNode.SetOpacity(opacity)
    shNode.GetScene().AddNode(displayNode)
    shNode.SetItemDataNode(folderID, displayNode)
    shNode.ItemModified(folderID)

  def getOrCreateFolderItem(self, parentID, name):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    folderID = shNode.GetItemChildWithName(parentID, name)
    if folderID == 0:
      folderID = shNode.CreateFolderItem(parentID, name)
      self.createFolderDisplayNode(folderID)
      shNode.SetItemAttribute(folderID, 'segmentation', '1')
    return folderID


class AtlasesTable(AtlasSegmentationBaseTable):

  def __init__(self):

    self.view = slicer.qMRMLSubjectHierarchyTreeView(slicer.util.mainWindow())
    self.view.attributeNameFilter = ('atlas')

    super().__init__()

    self.addButton.setText('Atlas')
    self.addButton.setToolTip('Add atlas')

    self.buttonsFrame.layout().addStretch(2)

    self.updateTable()

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
    WarpDrive.WarpDriveLogic().invertAtlases(WarpDrive.WarpDriveLogic().getParameterNode().GetNodeReference("InputNode"), int(WarpDrive.WarpDriveLogic().getParameterNode().GetParameter("InverseMode")))
    self.updateTable()

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

    self.sourceVisibleAction = qt.QAction(self)
    self.sourceVisibleAction.setText('Source')
    self.sourceVisibleAction.setCheckable(True)
    self.sourceVisibleAction.setChecked(False)
    self.targetVisibleAction = qt.QAction(self)
    self.targetVisibleAction.setText('Target')
    self.targetVisibleAction.setCheckable(True)
    self.targetVisibleAction.setChecked(True)
    self.previewSelectedAction = qt.QAction(self)
    self.previewSelectedAction.setText('Preview selected correction')
    self.previewSelectedAction.setCheckable(True)
    self.previewSelectedAction.setChecked(True)
    visibilityGroup = qt.QActionGroup(self)
    visibilityGroup.setExclusive(False)
    visibilityGroup.connect('triggered(QAction*)', self.visibilityChanged)
    visibilityGroup.addAction(self.sourceVisibleAction)
    visibilityGroup.addAction(self.targetVisibleAction)
    visibilityGroup.addAction(self.previewSelectedAction)

    modifiableCorrectionsAction = qt.QAction(self)
    modifiableCorrectionsAction.setText('Modifiable corrections')
    modifiableCorrectionsAction.setCheckable(True)
    modifiableCorrectionsAction.setChecked(False)
    modifiableCorrectionsAction.setToolTip('When checked, corrections are modified by new ones.')
    modifiableCorrectionsAction.connect('toggled(bool)', self.modiableCorrectionsChanged)

    settingsMenu = qt.QMenu(self)
    visibilityMenu = settingsMenu.addMenu("Visibility")
    visibilityMenu.addActions(visibilityGroup.actions())
    settingsMenu.addAction(modifiableCorrectionsAction)

    settingsAction = qt.QAction(self)
    settingsAction.setIcon(qt.QIcon(":/Icons/Small/SlicerConfigure.png"))
    settingsAction.setText('Settings')
    settingsButton = qt.QToolButton()
    settingsButton.setDefaultAction(settingsAction)
    settingsButton.setMenu(settingsMenu)
    settingsButton.setToolButtonStyle(qt.Qt.ToolButtonTextUnderIcon)
    settingsButton.setPopupMode(qt.QToolButton.InstantPopup)
    settingsButton.setIconSize(effectPixmap.rect().size())
    settingsButton.setSizePolicy(qt.QSizePolicy.MinimumExpanding,qt.QSizePolicy.Maximum)

    self.snapAutoApplyCheckBox = qt.QCheckBox('Auto apply')
    self.snapAutoApplyCheckBox.setToolTip('When checked, Snap will run after every correction is made.')
    self.snapAutoApplyCheckBox.setChecked(False)
    self.snapModeComboBox = qt.QComboBox()
    self.snapModeComboBox.addItems(['Last correction', 'All corrections'])
    self.snapSourceComboBox = slicer.qMRMLNodeComboBox()
    self.snapSourceComboBox.nodeTypes = ['vtkMRMLScalarVolumeNode']
    self.snapSourceComboBox.selectNodeUponCreation = False
    self.snapSourceComboBox.noneEnabled = True
    self.snapSourceComboBox.addEnabled = False
    self.snapSourceComboBox.removeEnabled = False
    self.snapSourceComboBox.setMRMLScene(slicer.mrmlScene)
    self.snapSourceComboBox.setToolTip("Select the source volume for snap")
    self.snapTargetComboBox = slicer.qMRMLNodeComboBox()
    self.snapTargetComboBox.nodeTypes = ['vtkMRMLScalarVolumeNode']
    self.snapTargetComboBox.selectNodeUponCreation = False
    self.snapTargetComboBox.noneEnabled = True
    self.snapTargetComboBox.addEnabled = False
    self.snapTargetComboBox.removeEnabled = False
    self.snapTargetComboBox.setMRMLScene(slicer.mrmlScene)
    self.snapTargetComboBox.setToolTip("Select the target volume for snap")
    menu = qt.QMenu(self)
    for widget,signal in zip([self.snapAutoApplyCheckBox, self.snapModeComboBox, self.snapSourceComboBox, self.snapTargetComboBox], ["stateChanged(int)", "currentIndexChanged(int)", "currentNodeChanged(vtkMRMLNode*)", "currentNodeChanged(vtkMRMLNode*)"]):
      widget.connect(signal, self.updateSnapOptionsFromGUI)
      a = qt.QWidgetAction(self)
      a.setDefaultWidget(widget)
      menu.addAction(a) 
    effectPixmap = qt.QPixmap(os.path.join(os.path.split(WarpDrive.__file__)[0], 'Resources', 'Icons', 'Magnet.png'))
    effectIcon = qt.QIcon(effectPixmap)    
    self.sanpAction = qt.QAction(self)
    self.sanpAction.setIcon(effectIcon)
    self.sanpAction.setText('Snap')
    self.sanpAction.setToolTip('Calls ANTs SyN registering the image with modifications to the MNI template. Then applies the transform to the target fiducials.')
    self.sanpAction.setCheckable(False)
    self.sanpAction.connect("triggered(bool)", self.onSnapTriggered)
    snapToolButton = qt.QToolButton()
    snapToolButton.setDefaultAction(self.sanpAction)
    snapToolButton.setMenu(menu)
    snapToolButton.setToolButtonStyle(qt.Qt.ToolButtonTextUnderIcon)
    snapToolButton.setPopupMode(qt.QToolButton.MenuButtonPopup)
    snapToolButton.setIconSize(effectPixmap.rect().size())
    snapToolButton.setSizePolicy(qt.QSizePolicy.MinimumExpanding,qt.QSizePolicy.Maximum)

    self.buttonsFrame.layout().addWidget(self.undoButton,1)
    self.buttonsFrame.layout().addWidget(settingsButton,1)
    if hasattr(slicer.modules,'antsregistration'):
      self.buttonsFrame.layout().addWidget(snapToolButton,1)

    self.addButton.setText('Fixed point')
    self.addButton.setToolTip('Add fixed point')
  

  def updateSnapOptionsFromGUI(self):
    pass

  def onSnapTriggered(self):
    pass

  def visibilityChanged(self, action):
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

  def modiableCorrectionsChanged(self, checked):
    pass

class WarpDriveCorrectionsManager(WarpDriveCorrectionsTable, VTKObservationMixin):
  def __init__(self):
    WarpDriveCorrectionsTable.__init__(self)
    VTKObservationMixin.__init__(self)
    self.sourceFiducialNodeID = ""
    self.targetFiducialNodeID = ""
    self._updatingFiducials = False
    self.targetFiducialObservers = []
    self.sourceFiducialObservers = []
    self.parameterNode = WarpDrive.WarpDriveLogic().getParameterNode()
    self.addObserver(self.parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateNodesListeners)
    self.addObserver(self.parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromSnapOptions)
    self._updatingSnapGUI = False

  def visibilityChanged(self, action):
    if action.text == 'Source':
      nodeID = self.sourceFiducialNodeID  
    elif action.text == 'Target':
      nodeID = self.targetFiducialNodeID
    else:
      nodeID = None
    if nodeID:
      slicer.util.getNode(nodeID).GetDisplayNode().SetVisibility(action.checked)

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
        self.targetVisibleAction.checked = targetFiducialNode.GetDisplayNode().GetVisibility()
    if self.sourceFiducialNodeID != "":
      sourceFiducialNode = slicer.mrmlScene.GetNodeByID(self.sourceFiducialNodeID)
      if sourceFiducialNode.GetDisplayNode():
        self.sourceVisibleAction.checked = sourceFiducialNode.GetDisplayNode().GetVisibility()

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
    undoSnap = self.undoSnapIfPresent()
    if undoSnap:
      return
    if self.model.rowCount() == 0:
      return
    index = self.model.index(self.model.rowCount()-1, 1)
    correctionName = self.model.itemData(index)[0]
    self.removeCorrectionByName(correctionName)

  def undoSnapIfPresent(self):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    markupsNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLMarkupsFiducialNode')
    markupsNodes.UnRegister(slicer.mrmlScene)
    for i in range(markupsNodes.GetNumberOfItems()):
      backupNode = markupsNodes.GetItemAsObject(i)
      if ('SnapBackUp' in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(backupNode))):
        targetFiducialNode = slicer.mrmlScene.GetNodeByID(self.targetFiducialNodeID)
        targetFiducialNode.Copy(backupNode)
        slicer.mrmlScene.RemoveNode(backupNode)
        targetFiducialNode.CreateDefaultDisplayNodes()
        targetFiducialNode.GetDisplayNode().SetGlyphTypeFromString('Sphere3D')
        targetFiducialNode.GetDisplayNode().SetGlyphScale(1)
        targetFiducialNode.GetDisplayNode().SetPointLabelsVisibility(0)
        for sliceNode in slicer.util.getNodesByClass('vtkMRMLSliceNode'):
          sliceNode.Modified()
        snapOptions = json.loads(self.parameterNode.GetParameter("SnapOptions"))
        snapOptions['SnapRun'] = True
        self.parameterNode.SetParameter("SnapOptions", json.dumps(snapOptions))
        self.parameterNode.SetParameter("Update","true")
        return True

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
    if not self.previewSelectedAction.checked:
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

  def updateSnapOptionsFromGUI(self):
    if self._updatingSnapGUI or not self.parameterNode.GetParameter("SnapOptions"):
      return
    options = json.loads(self.parameterNode.GetParameter("SnapOptions"))
    options['AutoApply'] = self.snapAutoApplyCheckBox.isChecked()
    options['Mode'] = self.snapModeComboBox.currentText
    options['SourceID'] = self.snapSourceComboBox.currentNodeID
    options['TargetID'] = self.snapTargetComboBox.currentNodeID
    self.parameterNode.SetParameter("SnapOptions", json.dumps(options))

  def updateGUIFromSnapOptions(self, caller=None, event=None):
    if not self.parameterNode.GetParameter("SnapOptions"):    
      return
    self._updatingSnapGUI = True
    options = json.loads(self.parameterNode.GetParameter("SnapOptions"))
    self.snapAutoApplyCheckBox.setChecked(options['AutoApply'] if 'AutoApply' in options else False)
    self.snapModeComboBox.setCurrentText(options['Mode'] if 'Mode' in options else 'Last correction')
    sourceImageNodeID = options['SourceID'] if 'SourceID' in options else ''
    if not sourceImageNodeID and self.parameterNode.GetNodeReferenceID("ImageNode"):
      sourceImageNodeID = self.parameterNode.GetNodeReferenceID("ImageNode")
    self.snapSourceComboBox.setCurrentNodeID(sourceImageNodeID)
    targetImageNodeID = options['TargetID'] if 'TargetID' in options else ''
    if not targetImageNodeID and self.parameterNode.GetNodeReferenceID("TemplateNode"):
      targetImageNodeID = self.parameterNode.GetNodeReferenceID("TemplateNode")
    self.snapTargetComboBox.setCurrentNodeID(targetImageNodeID)
    self._updatingSnapGUI = False

  def onSnapTriggered(self):
    snapMode = self.snapModeComboBox.currentText
    sourceImageNodeID = self.snapSourceComboBox.currentNodeID
    targetImageNodeID = self.snapTargetComboBox.currentNodeID
    targetFiducialNodeID = self.targetFiducialNodeID
    WarpDrive.WarpDriveLogic().runSnap(snapMode, sourceImageNodeID, targetImageNodeID, targetFiducialNodeID)

  def modiableCorrectionsChanged(self, checked):
    self.parameterNode.SetParameter("ModifiableCorrections",str(int(checked)))