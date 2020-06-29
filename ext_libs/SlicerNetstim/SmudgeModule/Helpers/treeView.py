import qt, vtk, slicer
from PythonQt import BoolResult
import os
import numpy as np

import SmudgeModule, ImportAtlas, TransformsUtil


class treeViewFilter(object):

  def __init__(self):
    # defaults
    self.parameterNode = SmudgeModule.SmudgeModuleLogic().getParameterNode()
    self.toolTip = ''
    self.filterDictionary = {'nodeTypes': (), 'attributeNameFilter': (''), 'attributeValueFilter': ('')}
    self.columnHidden = {'idColumn': True, 'transformColumn': True, 'descriptionColumn': True}

  def deleteFunction(self, node):
    if not node:
      return
    # get subject hierarchy node ID
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    nodeID = shNode.GetItemByDataNode(node)
    # get children
    removeIDs = vtk.vtkIdList()
    shNode.GetItemChildren(nodeID, removeIDs, True)
    # add selected ID
    removeIDs.InsertNextId(nodeID)
    # remove
    for i in range(removeIDs.GetNumberOfIds()):
      shNode.RemoveItem(removeIDs.GetId(i))

  def addFunction(self):
    pass

  def doubleClickFunction(self, node):
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
    self.centerPosition(centerList)

  def centerPosition(self, centerList):
    # create markups node, add center as fiducial and jump and center slices
    markupsNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    markupsNode.GetDisplayNode().SetVisibility(False)
    markupsNode.AddFiducialFromArray(np.array(centerList),'')
    markupsLogic = slicer.modules.markups.logic()
    markupsLogic.JumpSlicesToNthPointInMarkup(markupsNode.GetID(),0,True)
    slicer.mrmlScene.RemoveNode(markupsNode)

  def renameFunction(self, node):
    if not node:
      return
    name = qt.QInputDialog.getText(qt.QWidget(),'Rename','New name:')
    if name != '':
      node.SetName(name)


class treeViewSceneFilter(treeViewFilter):

  def __init__(self):
    super().__init__()
    self.name = 'Scene'
    self.toolTip = ''
    self.addText = ''

  def deleteFunction(self, node):
    pass


class treeViewAtlasFilter(treeViewFilter):

  def __init__(self):
    super().__init__()
    self.name = 'Atlases'
    self.toolTip = 'Lead-DBS Atlases. Press + to import more atlases'
    self.addText = 'Add Atlas'
    self.filterDictionary['nodeTypes'] = ('vtkMRMLModelNode','vtkMRMLFolderDisplayNode')
    self.filterDictionary['attributeNameFilter'] = ('atlas')

  def addFunction(self):
    if self.parameterNode.GetParameter("MNIAtlasPath") != ".":
      items = ImportAtlas.ImportAtlasLogic().getValidAtlases(self.parameterNode.GetParameter("MNIAtlasPath"))
      result = BoolResult()
      atlasName = qt.QInputDialog.getItem(qt.QWidget(),'Select Atlas','',items,0,0,result)
      if result:
        ImportAtlas.ImportAtlasLogic().run(os.path.join(self.parameterNode.GetParameter("MNIAtlasPath"), atlasName))    


class treeViewDrawingsFilter(treeViewFilter):

  def __init__(self):
    super().__init__()
    self.name = 'Fixed Points'
    self.toolTip = 'Points in this list will remian fixed in following drawings. Press + to add fiducials.'
    self.addText = 'Add Fiducial'
    self.filterDictionary['nodeTypes'] = ('vtkMRMLMarkupsCurveNode','vtkMRMLMarkupsFiducialNode','vtkMRMLFolderDisplayNode')
    self.filterDictionary['attributeNameFilter'] = ('drawing')
    self.filterDictionary['attributeValueFilter'] = ('1')

  def deleteFunction(self, node):
    if not isinstance(node, slicer.vtkMRMLFolderDisplayNode):
      super().deleteFunction(node)

  def addFunction(self):
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
    fiducialNode.SetName(slicer.mrmlScene.GenerateUniqueName('Fixed Point'))
    # add to subject hierarchy
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    shNode.SetItemAttribute(shNode.GetItemByDataNode(fiducialNode), 'drawing', '1')
    # activate placement
    selectionNode.SetActivePlaceNodeID(fiducialNode.GetID())
    interactionNode.SetCurrentInteractionMode(interactionNode.Place)

  def renameFunction(self, node):
    super().renameFunction(node)
    if isinstance(node, slicer.vtkMRMLMarkupsFiducialNode):
      node.SetNthControlPointLabel(0, node.GetName())

class treeViewSavedWarpFilter(treeViewFilter):

  def __init__(self):
    super().__init__()
    self.name = 'Output Warp'
    self.toolTip = 'Saved user modifications. Press + to save current state. Double click to change current warp.'
    self.addText = 'Duplicate Current Warp'
    self.filterDictionary['nodeTypes'] = ('vtkMRMLTransformNode','vtkMRMLGridTransformNode')
    self.filterDictionary['attributeNameFilter'] = ('savedWarp')
    self.columnHidden['descriptionColumn'] = False

  def deleteFunction(self, node):
    if node and node.GetDescription() != 'Current':
      super().deleteFunction(node)

  def addFunction(self):
    if not self.parameterNode.GetNodeReferenceID("warpID"):
      return
    qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
    qt.QApplication.processEvents()
    # get node
    warpNode = self.parameterNode.GetNodeReference("warpID")
    # save visibility
    vis = warpNode.GetDisplayNode().GetVisibility()
    warpNode.GetDisplayNode().SetVisibility(False)
    # clone
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    clonedID = slicer.modules.subjecthierarchy.logic().CloneSubjectHierarchyItem(shNode, shNode.GetItemByDataNode(warpNode))
    shNode.SetItemAttribute(clonedID, 'savedWarp', '1')
    # flat new warp
    newWarpNode = shNode.GetItemDataNode(clonedID)
    newWarpNode.SetName(slicer.mrmlScene.GenerateUniqueName('SavedWarp'))
    size, origin, spacing = TransformsUtil.TransformsUtilLogic().getGridDefinition(newWarpNode)
    outNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
    referenceVolume = TransformsUtil.TransformsUtilLogic().createEmpyVolume(size, origin, spacing)    
    transformsLogic = slicer.modules.transforms.logic()
    transformsLogic.ConvertToGridTransform(newWarpNode, referenceVolume, outNode)
    newWarpNode.SetAndObserveTransformFromParent(outNode.GetTransformFromParent())
    # remove aux
    slicer.mrmlScene.RemoveNode(outNode)
    slicer.mrmlScene.RemoveNode(referenceVolume)
    # restore visibility
    warpNode.GetDisplayNode().SetVisibility(vis)
    # simulate double click to change
    self.doubleClickFunction(newWarpNode)
    qt.QApplication.setOverrideCursor(qt.QCursor(qt.Qt.ArrowCursor))

  def doubleClickFunction(self, node):
    SmudgeModule.SmudgeModuleLogic().removeRedoNodes()
    # reset descriptions
    previousWarpNode = self.parameterNode.GetNodeReference("warpID")
    previousWarpNode.SetDescription('')
    # apply current node
    self.parameterNode.GetNodeReference("glanatCompositeID").SetAndObserveTransformNodeID(node.GetID())
    # apply visibility
    node.GetDisplayNode().SetVisibility(previousWarpNode.GetDisplayNode().GetVisibility())
    previousWarpNode.GetDisplayNode().SetVisibility(False)
    # set description
    node.SetDescription('Current')
    # change parameter node
    self.parameterNode.SetNodeReferenceID("warpID", node.GetID())
    




class WarpDriveTreeView(qt.QWidget):

  def __init__(self):

    super().__init__()

    layout = qt.QGridLayout(self)  

    # set up tree view
    self.treeView = slicer.qMRMLSubjectHierarchyTreeView(slicer.util.mainWindow())
    self.treeView.setMRMLScene(slicer.mrmlScene)
    self.treeView.contextMenuEnabled = False
    self.treeView.setEditTriggers(0) # disable double click to edit

    # add
    addPixmap = qt.QPixmap(os.path.join(os.path.split(__file__)[0] ,'Icons', 'Add.png'))
    addIcon = qt.QIcon(addPixmap)
    self.addButton = qt.QToolButton()
    self.addButton.setToolButtonStyle(qt.Qt.ToolButtonTextBesideIcon)
    self.addButton.setSizePolicy(qt.QSizePolicy.MinimumExpanding,qt.QSizePolicy.Maximum)
    self.addButton.setIcon(addIcon)
    self.addButton.setIconSize(addPixmap.rect().size())
    self.addButton.setToolTip('Add')
    # delete
    deletePixmap = qt.QPixmap(os.path.join(os.path.split(__file__)[0] ,'Icons', 'Delete.png'))
    deleteIcon = qt.QIcon(deletePixmap)
    self.deleteButton = qt.QPushButton()
    self.deleteButton.setIcon(deleteIcon)
    self.deleteButton.setIconSize(deletePixmap.rect().size())
    self.deleteButton.setToolTip('Delete')
    # rename
    renamePixmap = qt.QPixmap(os.path.join(os.path.split(__file__)[0] ,'Icons', 'Rename.png'))
    renameIcon = qt.QIcon(renamePixmap)
    self.renameButton = qt.QPushButton()
    self.renameButton.setIcon(renameIcon)
    self.renameButton.setIconSize(renamePixmap.rect().size())
    self.renameButton.setToolTip('Rename')

    # set up filters
    filters = [treeViewSavedWarpFilter(), treeViewAtlasFilter(), treeViewDrawingsFilter(), treeViewSceneFilter()]
    self.radioButtons = []

    for filt,pos in zip(filters,[[0,0],[0,2],[1,0],[1,2]]):
      filterRadioButton = qt.QRadioButton(filt.name)
      filterRadioButton.setToolTip(filt.toolTip)
      filterRadioButton.clicked.connect(lambda b,f=filt: self.onFilterRadioButtonClicked(f))
      layout.addWidget(filterRadioButton,pos[0],pos[1],1,2)
      self.radioButtons.append(filterRadioButton)


    layout.addWidget(self.addButton,0,4,1,2)
    layout.addWidget(self.deleteButton,1,4,1,1)
    layout.addWidget(self.renameButton,1,5,1,1)
    layout.addWidget(self.treeView,2,0,1,6)

    # when adding fixed points while one of them is selected the new one is not set in the correct parent folder
    # this is overdoing, but fixes the problem
    self.treeView.model().rowsAboutToBeInserted.connect(lambda: self.treeView.setCurrentItem(0))

    # init
    self.radioButtons[3].animateClick()


  def onFilterRadioButtonClicked(self, filt):
    # set scene item. if not crushes
    self.treeView.setCurrentItem(0)
    # filter data tree
    for key,value in filt.filterDictionary.items():
      setattr(self.treeView, key, value)
    # columns hidden
    for key,value in filt.columnHidden.items():
      self.treeView.setColumnHidden(eval('self.treeView.model().'+key), value)
    # reset depth
    self.treeView.expandToDepth(0)
    # set add function
    self.addButton.disconnect('clicked(bool)')
    self.addButton.connect('clicked(bool)', filt.addFunction)
    self.addButton.connect('clicked(bool)', self.updateTree)
    self.addButton.setText(filt.addText)
    # set delete function
    self.deleteButton.disconnect('clicked(bool)')
    self.deleteButton.connect('clicked(bool)', lambda: filt.deleteFunction(self.currentNode()))
    # set rename function
    self.renameButton.disconnect('clicked(bool)')
    self.renameButton.connect('clicked(bool)', lambda: filt.renameFunction(self.currentNode()))
    # set double click function
    self.treeView.doubleClicked.disconnect()
    self.treeView.doubleClicked.connect(lambda: filt.doubleClickFunction(self.currentNode()))

  def currentNode(self):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    if self.treeView.currentItem():
      return shNode.GetItemDataNode(self.treeView.currentItem())

  def updateTree(self):
    # annimate click will update nodes in tree - useful when adding saved warps 
    for button in self.radioButtons:
      if button.checked:
        button.animateClick()
