import vtk, qt, slicer
import os, sys, shutil
import uuid
from scipy import io
import numpy as np
from subprocess import call
import json
import glob

from . import WarpDriveUtil, GridNodeHelper

def checkExtensionInstall(extensionName):
  em = slicer.app.extensionsManagerModel()
  if not em.isExtensionInstalled(extensionName):
    extensionMetaData = em.retrieveExtensionMetadataByName(extensionName)
    url = os.path.join(em.serverUrl().toString(), 'download', 'item', extensionMetaData['item_id'])
    extensionPackageFilename = os.path.join(slicer.app.temporaryPath, extensionMetaData['md5'])
    slicer.util.downloadFile(url, extensionPackageFilename)
    em.installExtension(extensionPackageFilename)
    qt.QMessageBox.information(qt.QWidget(), '', 'Slicer will install %s and quit.\nPlease restart.' % extensionName)
    slicer.util.exit()
    return True

def updateParameterNodeFromArgs(parameterNode): 
  if parameterNode.GetParameter("MNIPath") != '':
    return # was already called

  args = sys.argv
  if (len(sys.argv) > 2) and os.path.isfile(os.path.join(sys.argv[1],'lead.m')):
    pathsSeparator = uuid.uuid4().hex
    subjectPaths = pathsSeparator.join(sys.argv[2:])
    subjectPath = subjectPaths.split(pathsSeparator)[0]
    MNIPath = os.path.join(sys.argv[1],'templates','space','MNI152NLin2009bAsym')
    MNIAtlasPath = os.path.join(MNIPath,'atlases')
    if sys.platform == "darwin":
      ext = "maci64"
    elif sys.platform.startswith('win'):
      ext = 'exe'
    else:
      ext = 'glnxa64'
    antsApplyTransformsPath = os.path.join(sys.argv[1],'ext_libs','ANTs','antsApplyTransforms.' + ext)
    # set parameter node
    parameterNode.SetParameter("separator", pathsSeparator)
    parameterNode.SetParameter("subjectPaths", subjectPaths)
    parameterNode.SetParameter("subjectN", "0")
    parameterNode.SetParameter("subjectPath", subjectPath)
    parameterNode.SetParameter("MNIPath", MNIPath)
    parameterNode.SetParameter("MNIAtlasPath", MNIAtlasPath)
    parameterNode.SetParameter("antsApplyTransformsPath", antsApplyTransformsPath)
    parameterNode.SetNodeReferenceID("ImageNode", None)
    parameterNode.SetNodeReferenceID("TemplateNode", None)
    return True

class LeadBIDS():
  def __init__(self, subjectPath):
    self.subjectPath = subjectPath
    self.subjectID = os.path.split(self.subjectPath)[-1]
  
  def getNormalizationMethod(self):
    return os.path.join(self.getNormalizationPath(), 'log', self.subjectID + '_desc-normmethod.json')

  def getCoregImages(self, modality='*'):
    return glob.glob(os.path.join(self.subjectPath, 'coregistration', 'anat', self.subjectID + '*ses-preop_' + modality + '.nii'))

  def getNormalizedImages(self):
    return glob.glob(os.path.join(self.subjectPath, 'normalization', 'anat', self.subjectID + '*ses-preop*.nii'))

  def getANTSForwardWarp(self):
    return os.path.join(self.getNormalizationPath(), 'transformations', self.subjectID + '_from-anchorNative_to-MNI152NLin2009bAsym_desc-ants.nii.gz')

  def getANTSInverseWarp(self):
    return os.path.join(self.getNormalizationPath(), 'transformations', self.subjectID + '_from-MNI152NLin2009bAsym_to-anchorNative_desc-ants.nii.gz')

  def getNormalizationPath(self):
    return os.path.join(self.subjectPath, 'normalization')

  def getWarpDrivePath(self):
    return os.path.join(self.subjectPath, 'warpdrive')

def saveApprovedData(subjectPath):
  normalizationMethodFile = LeadBIDS(subjectPath).getNormalizationMethod()
  with open(normalizationMethodFile, 'a+') as f:
    normalizationMethod = json.load(f)
    normalizationMethod['approval'] = 1
    f.seek(0)
    json.dump(normalizationMethod, f)
    f.truncate()

def loadSubjectTransform(subjectPath, antsApplyTransformsPath):

  transformPath = LeadBIDS(subjectPath).getANTSForwardWarp()

  # update subject warp fields to new lead dbs specification
  if os.path.isfile(transformPath.replace('.nii.gz','.h5')):
    updateTranform(subjectPath, antsApplyTransformsPath)

  # load glanat composite
  glanatCompositeNode = slicer.util.loadTransform(transformPath)
  
  return glanatCompositeNode

def updateTranform(directory, antsApplyTransformsPath):
  bidsSubject = LeadBIDS(directory)
  transforms = [bidsSubject.getANTSForwardWarp(), bidsSubject.getANTSInverseWarp()]
  references = [bidsSubject.getNormalizedImages()[0], bidsSubject.getCoregImages()[0]]
  for transform,reference in zip(transforms,references):
    transformFullPath = os.path.join(directory,transform + '.h5') # in case inverse doesnt exist
    if os.path.isfile(transformFullPath):
      h5transform = transform.replace('.nii.gz','.h5')
      command = antsApplyTransformsPath + " -r " + reference + " -t " + h5transform + " -o [" + transform + ",1] -v 1"
      commandOut = call(command, env=slicer.util.startupEnvironment(), shell=True) # run antsApplyTransforms
      os.remove(h5transform)
  return True
  
def queryUserApproveSubject(subjectPath):
  msgBox = qt.QMessageBox()
  msgBox.setText('No corrections made')
  msgBox.setInformativeText('Save subject as approved?')
  msgBox.setStandardButtons(qt.QMessageBox().Save | qt.QMessageBox().Discard | qt.QMessageBox().Cancel)
  ret = msgBox.exec_()
  if ret == qt.QMessageBox().Cancel:
    return False
  if ret == qt.QMessageBox().Save:
    saveApprovedData(subjectPath)
  return True

def applyChanges(subjectPath, inputNode, imageNode):

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
  slicer.util.saveNode(inputNode, LeadBIDS(subjectPath).getANTSForwardWarp())

  # BACKWARD

  inputNode.Inverse()
  # to grid
  outNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
  slicer.modules.transforms.logic().ConvertToGridTransform(inputNode, imageNode, outNode)
  # save
  slicer.util.saveNode(outNode, LeadBIDS(subjectPath).getANTSInverseWarp())
  # delete aux node
  slicer.mrmlScene.RemoveNode(outNode)
  
  # back to original
  inputNode.Inverse()
  imageNode.SetAndObserveTransformNodeID(inputNode.GetID())
  
  qt.QApplication.setOverrideCursor(qt.QCursor(qt.Qt.ArrowCursor))


def setTargetFiducialsAsFixed():
  shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
  fiducialNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLMarkupsFiducialNode')
  fiducialNodes.UnRegister(slicer.mrmlScene)
  for i in range(fiducialNodes.GetNumberOfItems()):
    fiducialNode = fiducialNodes.GetItemAsObject(i)
    if 'target' in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(fiducialNode)):
      # get parent folder
      parentFolder = shNode.GetItemDataNode(shNode.GetItemParent(shNode.GetItemByDataNode(fiducialNode)))
      parentFolderName = parentFolder.GetName()
      # remove target attribute
      shNode.RemoveItemAttribute(shNode.GetItemByDataNode(fiducialNode), 'target')
      # add as fixed point
      WarpDriveUtil.addFixedPoint(fiducialNode)
      # remove correction
      removeNodeAndChildren(parentFolder)
      # change fixed point name
      fiducialNode.SetName(parentFolderName)

def saveCurrentScene(subjectPath):
  """
  Save corrections and fixed points is subject directory so will be loaded next time
  """
  warpDriveSavePath = LeadBIDS(subjectPath).getWarpDrivePath()
  # delete previous
  if os.path.isdir(warpDriveSavePath):
    shutil.rmtree(warpDriveSavePath)
  # create directories
  os.mkdir(warpDriveSavePath)
  os.mkdir(os.path.join(warpDriveSavePath,'Data'))
  # set scene URL
  slicer.mrmlScene.SetURL(os.path.join(warpDriveSavePath, 'WarpDriveScene.mrml'))
  # save corrections
  shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
  for nodeType, nodeExt in zip(['vtkMRMLMarkupsFiducialNode', 'vtkMRMLLabelMapVolumeNode'], ['.fcsv', '.nrrd']):
    nodes = slicer.mrmlScene.GetNodesByClass(nodeType)
    nodes.UnRegister(slicer.mrmlScene)
    for i in range(nodes.GetNumberOfItems()):
      node = nodes.GetItemAsObject(i)
      if 'correction' in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(node)):
        slicer.util.saveNode(node, os.path.join(warpDriveSavePath, 'Data', uuid.uuid4().hex + nodeExt))
  # save scene
  slicer.mrmlScene.Commit()

def DeleteCorrections():
  shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
  # delete folders
  folderNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLFolderDisplayNode')
  folderNodes.UnRegister(slicer.mrmlScene)
  for i in range(folderNodes.GetNumberOfItems()):
    folderNode = folderNodes.GetItemAsObject(i)
    if 'correction' in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(folderNode)):
      removeNodeAndChildren(folderNode)

def removeNodeAndChildren(node):
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