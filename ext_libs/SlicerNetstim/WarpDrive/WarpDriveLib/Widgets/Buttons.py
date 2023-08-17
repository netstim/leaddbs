import qt, slicer
import os
import WarpDrive
import json
import glob

def volumeToModel(volumeNode):
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

def createFolderDisplayNode(folderID, color=[0.66,0.66,0.66], opacity=1.0):
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

def getOrCreateFolderItem(parentID, name):
  shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
  folderID = shNode.GetItemChildWithName(parentID, name)
  if folderID == 0:
    folderID = shNode.CreateFolderItem(parentID, name)
    createFolderDisplayNode(folderID)
    shNode.SetItemAttribute(folderID, 'atlas', 'native')
  return folderID

def onAddSegmentation():
  shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
  currentSubject = json.loads(WarpDrive.WarpDriveLogic().getParameterNode().GetParameter("CurrentSubject"))
  segmentationsPath = os.path.join(os.path.dirname(currentSubject['warpdrive_path']), 'segmentations')
  niiFiles = glob.glob(os.path.join(segmentationsPath, '**', '*.nii*'), recursive=True)
  for file in niiFiles:
    # load
    volumeNode = slicer.util.loadVolume(file,{'show':False})
    modelNode = volumeToModel(volumeNode)
    slicer.mrmlScene.RemoveNode(volumeNode)
    modelNode.GetDisplayNode().SetVisibility2D(1)
    # set in hierarchy
    basePath,fileName = os.path.split(file)
    parentFolders = os.path.relpath(basePath, os.path.dirname(currentSubject['warpdrive_path'])).split(os.sep)
    folderID = getOrCreateFolderItem(shNode.GetSceneItemID(), parentFolders.pop(0))
    for parentName in parentFolders:
      folderID = getOrCreateFolderItem(folderID, parentName)
    shNode.SetItemParent(shNode.GetItemByDataNode(modelNode), folderID)
    shNode.SetItemAttribute(shNode.GetItemByDataNode(modelNode), 'atlas', 'native')
    modelNode.SetName(fileName.split('.')[0])
  tb = next(filter(lambda x: isinstance(x,qt.QToolBar) and x.windowTitle=='LeadDBS', slicer.util.mainWindow().children()))
  tb.invertAtlases(WarpDrive.WarpDriveLogic().getParameterNode().GetNodeReference("InputNode"))

class addSegmentationButton(qt.QToolButton):
  def __init__(self):
    super().__init__()
    effectPixmap = qt.QPixmap(os.path.join(os.path.split(WarpDrive.__file__)[0], 'Resources', 'Icons', 'Add.png'))
    effectIcon = qt.QIcon(effectPixmap)
    self.setToolButtonStyle(qt.Qt.ToolButtonTextUnderIcon)
    self.setIcon(effectIcon)
    self.setIconSize(effectPixmap.rect().size())
    self.setEnabled(True)
    self.setSizePolicy(qt.QSizePolicy.MinimumExpanding,qt.QSizePolicy.Maximum)
    self.setText('Segmentation')
    self.setToolTip('Add segmentation. Looks for nifti files in the segmentations subfolder.')
    self.clicked.connect(onAddSegmentation)
