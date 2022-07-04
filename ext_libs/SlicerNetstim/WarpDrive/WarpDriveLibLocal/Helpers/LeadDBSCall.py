import vtk, qt, slicer
import os
import json

from . import GridNodeHelper

def saveApprovedData(normalizationMethodFile):
  with open(normalizationMethodFile, 'r') as f:
    normalizationMethod = json.load(f)
    normalizationMethod['approval'] = 1
  with open(normalizationMethodFile, 'w') as f:
    json.dump(normalizationMethod, f)

def queryUserApproveSubject():
  msgBox = qt.QMessageBox()
  msgBox.setText('No corrections made')
  msgBox.setInformativeText('Save subject as approved?')
  msgBox.setStandardButtons(qt.QMessageBox().Save | qt.QMessageBox().Discard | qt.QMessageBox().Cancel)
  ret = msgBox.exec_()
  return (ret != qt.QMessageBox().Cancel)

def applyChanges(inputNode, imageNode, forwardWarpPath, inverseWarpPath):

  qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
  qt.QApplication.processEvents()
  
  # undo changes to image node
  imageNode.SetAndObserveTransformNodeID(None)

  # FORWARD
  
  size, origin, spacing = GridNodeHelper.getGridDefinition(inputNode)
  # harden changes in input
  inputNode.HardenTransform()
  # to grid transform
  outNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
  referenceVolume = GridNodeHelper.emptyVolume(size, origin, spacing)    
  slicer.modules.transforms.logic().ConvertToGridTransform(inputNode, referenceVolume, outNode)
  # set to input and delete aux
  inputNode.SetAndObserveTransformFromParent(outNode.GetTransformFromParent())
  slicer.mrmlScene.RemoveNode(outNode)
  slicer.mrmlScene.RemoveNode(referenceVolume)
  # save
  slicer.util.saveNode(inputNode, forwardWarpPath)

  # BACKWARD

  inputNode.Inverse()
  # to grid
  outNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
  slicer.modules.transforms.logic().ConvertToGridTransform(inputNode, imageNode, outNode)
  # save
  slicer.util.saveNode(outNode, inverseWarpPath)
  # delete aux node
  slicer.mrmlScene.RemoveNode(outNode)
  
  # back to original
  inputNode.Inverse()
  imageNode.SetAndObserveTransformNodeID(inputNode.GetID())
  
  qt.QApplication.setOverrideCursor(qt.QCursor(qt.Qt.ArrowCursor))


def saveSourceTarget(warpDriveSavePath, sourceNode, targetNode):
  """
  Save source and target in subject directory so will be loaded next time
  """
  if not os.path.isdir(warpDriveSavePath):
    os.mkdir(warpDriveSavePath)  
  slicer.util.saveNode(sourceNode, os.path.join(warpDriveSavePath, 'source.json'))
  slicer.util.saveNode(targetNode, os.path.join(warpDriveSavePath, 'target.json'))

def getAtlasesNamesInScene():
  shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
  folderNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLFolderDisplayNode')
  folderNodes.UnRegister(slicer.mrmlScene)
  names = []
  for i in range(folderNodes.GetNumberOfItems()):
    folderNode = folderNodes.GetItemAsObject(i)
    if 'atlas' in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(folderNode)):
      if  shNode.GetItemParent(shNode.GetItemByDataNode(folderNode)) == shNode.GetSceneItemID():
        names.append(folderNode.GetName())
  return names

def saveSceneInfo(warpDriveSavePath):
  info = {}
  info["atlasNames"] = getAtlasesNamesInScene()
  with open(os.path.join(warpDriveSavePath,'info.json'), 'w') as jsonFile:
    json.dump(info, jsonFile)
