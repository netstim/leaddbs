import vtk, qt, slicer
import os, sys, shutil
import uuid
from scipy import io
import numpy as np
from subprocess import call

from . import WarpDriveUtil, GridNodeHelper

try:
  import h5py
except:
  slicer.util.pip_install('h5py')
  import h5py

try:
  import hdf5storage
except:
  slicer.util.pip_install('hdf5storage')
  import hdf5storage


def saveApprovedData(subjectPath):
  approvedFile = os.path.join(subjectPath,'ea_coreg_approved.mat')
  matfiledata = {}
  if os.path.isfile(approvedFile):
    try:
      # read file and copy data except for glanat
      with h5py.File(approvedFile,'r') as f:
        for k in f.keys():
          if k != 'glanat':
            keyValue = f[k][()]
            matfiledata[k] = keyValue
      # now add approved glanat
      matfiledata[u'glanat'] = np.array([2])

    except: # use other reader for .mat file
      f = io.loadmat(approvedFile)
      for k in f.keys():
        if k != 'glanat':
          keyValue = f[k]
          matfiledata[k] = keyValue
      matfiledata['glanat'] = np.array([[2]],dtype='uint8')
      io.savemat(approvedFile,matfiledata)
      return

  else:
    matfiledata[u'glanat'] = np.array([2])

  # save
  # for some reason putting subject path into hdf5storage.write doesnt work
  currentDir = os.getcwd()
  os.chdir(subjectPath)
  hdf5storage.write(matfiledata, '.', 'ea_coreg_approved.mat', matlab_compatible=True)
  os.chdir(currentDir)


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
    MNIPath = os.path.join(sys.argv[1],'templates','space','MNI_ICBM_2009b_NLIN_ASYM')
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


def loadSubjectTransform(subjectPath, antsApplyTransformsPath):

  # update subject warp fields to new lead dbs specification
  if os.path.isfile(os.path.join(subjectPath, 'glanatComposite.h5')):
    updateTranform(subjectPath, antsApplyTransformsPath)

  # load glanat composite
  glanatCompositeNode = slicer.util.loadTransform(os.path.join(subjectPath, 'glanatComposite.nii.gz'))
  
  return glanatCompositeNode

def updateTranform(directory, antsApplyTransformsPath):
  for transform,reference in zip(['glanatComposite','glanatInverseComposite'],['glanat','anat_t1']):
    transformFullPath = os.path.join(directory,transform + '.h5') # in case inverse doesnt exist
    if os.path.isfile(transformFullPath):
      command = antsApplyTransformsPath + " -r " + os.path.join(directory,reference + '.nii') + " -t " + transformFullPath + " -o [" + os.path.join(directory,transform + '.nii.gz') + ",1] -v 1"
      commandOut = call(command, env=slicer.util.startupEnvironment(), shell=True) # run antsApplyTransforms
      os.remove(transformFullPath)
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
  slicer.util.saveNode(inputNode, os.path.join(subjectPath,'glanatComposite.nii.gz'))

  # BACKWARD

  inputNode.Inverse()
  # to grid
  outNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
  slicer.modules.transforms.logic().ConvertToGridTransform(inputNode, imageNode, outNode)
  # save
  slicer.util.saveNode(outNode, os.path.join(subjectPath,'glanatInverseComposite.nii.gz'))
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
  warpDriveSavePath = os.path.join(subjectPath,'WarpDrive')
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