import qt, vtk, slicer
from PythonQt import BoolResult
import os
import numpy as np

import WarpDrive, ImportAtlas

from ..Helpers import WarpDriveUtil

class treeViewFilter(object):

  def __init__(self):
    # defaults
    self.parameterNode = WarpDrive.WarpDriveLogic().getParameterNode()
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

  def doubleClickFunction(self, node, tree=None):
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


class treeViewSegmentationsFilter(treeViewFilter):

  def __init__(self):
    super().__init__()
    self.name = 'Segments'
    self.toolTip = ''
    self.addText = ''
    self.filterDictionary['attributeNameFilter'] = ('Segment')



class treeViewCorrectionsFiducialsFilter(treeViewFilter):

  def __init__(self):
    super().__init__()
    self.name = 'Corrections'
    self.toolTip = ''
    self.addText = 'Add Fixed Point'
    self.filterDictionary['nodeTypes'] = ('vtkMRMLMarkupsFiducialNode','vtkMRMLFolderDisplayNode','vtkMRMLLabelMapVolumeNode')
    self.filterDictionary['attributeNameFilter'] = ('correction')
    self.columnHidden['descriptionColumn'] = False

  def deleteFunction(self, node):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    if shNode.GetItemName(shNode.GetItemParent(shNode.GetItemByDataNode(node))) == 'Fixed Points' or (isinstance(node,slicer.vtkMRMLFolderDisplayNode) and node.GetName()!='Fixed Points'):
      super().deleteFunction(node)

  def doubleClickFunction(self, node, tree):
    if isinstance(node, slicer.vtkMRMLMarkupsFiducialNode):
      super().doubleClickFunction(node)
    elif isinstance(node, slicer.vtkMRMLFolderDisplayNode):
      node.SetDescription('Disabled' if node.GetDescription() == 'Enabled' else 'Enabled')
    elif isinstance(node, slicer.vtkMRMLLabelMapVolumeNode):
      previousValue = int(node.GetDescription())
      newValue = qt.QInputDialog.getInt(qt.QWidget(),'Modify Value', 'Spread:', previousValue, 5, 50, 5)
      if newValue != previousValue:
        shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
        sourceNode = shNode.GetItemDataNode(shNode.GetItemChildWithName(shNode.GetItemParent(shNode.GetItemByDataNode(node)),'source'))
        newNode = WarpDriveUtil.generateMaskFromPointsSpread(sourceNode, newValue, self.parameterNode.GetNodeReference("InputNode"))
        # replace old by new node
        newNode.SetDescription(str(newValue))
        newNode.SetName('spread')
        shNode.SetItemAttribute(shNode.GetItemByDataNode(newNode), 'correction', '1')
        shNode.SetItemParent(shNode.GetItemByDataNode(newNode), shNode.GetItemParent(shNode.GetItemByDataNode(node)))
        slicer.mrmlScene.RemoveNode(node)

  def renameFunction(self, node):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    if shNode.GetItemName(shNode.GetItemParent(shNode.GetItemByDataNode(node))) == 'Fixed Points' or (isinstance(node,slicer.vtkMRMLFolderDisplayNode) and node.GetName()!='Fixed Points'):
      super().renameFunction(node)

  def addFunction(self):
    # interaction node
    interactionNode = slicer.app.applicationLogic().GetInteractionNode()
    selectionNode = slicer.app.applicationLogic().GetSelectionNode()
    selectionNode.SetReferenceActivePlaceNodeClassName("vtkMRMLMarkupsFiducialNode")
    # create aux marpus node
    fiducialNode = slicer.vtkMRMLMarkupsFiducialNode()
    slicer.mrmlScene.AddNode(fiducialNode)
    fiducialNode.CreateDefaultDisplayNodes() 
    fiducialNode.GetDisplayNode().SetGlyphScale(1)
    fiducialNode.SetLocked(1)
    fiducialNode.SetName(slicer.mrmlScene.GenerateUniqueName('Fixed Point'))
    # add to subject hierarchy
    WarpDriveUtil.addFixedPoint(fiducialNode)
    # activate placement
    selectionNode.SetActivePlaceNodeID(fiducialNode.GetID())
    interactionNode.SetCurrentInteractionMode(interactionNode.Place)

class treeViewAtlasFilter(treeViewFilter):

  def __init__(self):
    super().__init__()
    self.name = 'Atlases'
    self.toolTip = 'Lead-DBS Atlases. Press + to import more atlases'
    self.addText = 'Add Atlas'
    self.filterDictionary['nodeTypes'] = ('vtkMRMLModelNode','vtkMRMLFolderDisplayNode')
    self.filterDictionary['attributeNameFilter'] = ('atlas')

  def addFunction(self):
    if self.parameterNode.GetParameter("MNIAtlasPath") not in [".", ""]:
      directory = self.parameterNode.GetParameter("MNIAtlasPath") 
    else:
      with open(os.path.join(os.path.split(ImportAtlas.__file__)[0],'Resources','previousDirectory.txt'), 'r') as f: 
        directory = f.readlines()[0]
        if directory == ".":
          return
    # load
    items = ImportAtlas.ImportAtlasLogic().getValidAtlases(directory)
    result = BoolResult()
    atlasName = qt.QInputDialog.getItem(qt.QWidget(),'Select Atlas','',items,0,0,result)
    if result:
      ImportAtlas.ImportAtlasLogic().readAtlas(os.path.join(directory, atlasName, 'atlas_index.mat'))    


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
    addPixmap = qt.QPixmap(os.path.join(os.path.split(WarpDrive.__file__)[0], 'Resources', 'Icons', 'Add.png'))
    addIcon = qt.QIcon(addPixmap)
    self.addButton = qt.QToolButton()
    self.addButton.setToolButtonStyle(qt.Qt.ToolButtonTextBesideIcon)
    self.addButton.setSizePolicy(qt.QSizePolicy.MinimumExpanding,qt.QSizePolicy.Maximum)
    self.addButton.setIcon(addIcon)
    self.addButton.setIconSize(addPixmap.rect().size())
    self.addButton.setToolTip('Add')
    # delete
    deletePixmap = qt.QPixmap(os.path.join(os.path.split(WarpDrive.__file__)[0], 'Resources', 'Icons', 'Delete.png'))
    deleteIcon = qt.QIcon(deletePixmap)
    self.deleteButton = qt.QPushButton()
    self.deleteButton.setIcon(deleteIcon)
    self.deleteButton.setIconSize(deletePixmap.rect().size())
    self.deleteButton.setToolTip('Delete')
    # rename
    renamePixmap = qt.QPixmap(os.path.join(os.path.split(WarpDrive.__file__)[0], 'Resources', 'Icons', 'Rename.png'))
    renameIcon = qt.QIcon(renamePixmap)
    self.renameButton = qt.QPushButton()
    self.renameButton.setIcon(renameIcon)
    self.renameButton.setIconSize(renamePixmap.rect().size())
    self.renameButton.setToolTip('Rename')

    # set up filters
    #filters = [treeViewSavedWarpFilter(), treeViewAtlasFilter(), treeViewDrawingsFilter(), treeViewSceneFilter()]
    filters = [treeViewCorrectionsFiducialsFilter(), treeViewAtlasFilter(), treeViewSegmentationsFilter()]
    self.radioButtons = []

    for filt,pos in zip(filters,[[0,0],[0,2],[0,4]]):
      filterRadioButton = qt.QRadioButton(filt.name)
      filterRadioButton.setToolTip(filt.toolTip)
      filterRadioButton.clicked.connect(lambda b,f=filt: self.onFilterRadioButtonClicked(f))
      layout.addWidget(filterRadioButton,pos[0],pos[1],1,2)
      self.radioButtons.append(filterRadioButton)


    layout.addWidget(self.addButton,0,6,1,2)
    layout.addWidget(self.deleteButton,0,8,1,1)
    layout.addWidget(self.renameButton,0,9,1,1)
    layout.addWidget(self.treeView,1,0,1,10)

    # when adding fixed points while one of them is selected the new one is not set in the correct parent folder
    # this is overdoing, but fixes the problem
    self.treeView.model().rowsAboutToBeInserted.connect(lambda: self.treeView.setCurrentItem(0))

    # init
    self.radioButtons[0].animateClick()


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
    self.treeView.doubleClicked.connect(lambda: filt.doubleClickFunction(self.currentNode(), self.treeView))

  def currentNode(self):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    if self.treeView.currentItem():
      return shNode.GetItemDataNode(self.treeView.currentItem())

  def updateTree(self, caller=None, event=None):
    # reset the attribute filter so the tree is updated
    currentValue = self.treeView.attributeNameFilter
    self.treeView.attributeNameFilter = ('')
    self.treeView.attributeNameFilter = currentValue
