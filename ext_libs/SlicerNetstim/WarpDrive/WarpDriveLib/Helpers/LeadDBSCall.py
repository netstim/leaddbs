import vtk, qt, ctk, slicer
import os, sys, platform
import json
import shutil
import glob

def getApprovedData(normalizationMethodFile):
  with open(normalizationMethodFile, 'r') as f:
    normalizationMethod = json.load(f)
    value = normalizationMethod['approval']
  return value

def setApprovedData(normalizationMethodFile, value):
  with open(normalizationMethodFile, 'r') as f:
    normalizationMethod = json.load(f)
    normalizationMethod['approval'] = value
  with open(normalizationMethodFile, 'w') as f:
    json.dump(normalizationMethod, f)

def applyChanges(correctionsTransformNodeID, nativeReferencePath, templateReferencePath, forwardWarpPath, inverseWarpPath, subjectWarpDrivePath, useExternalInstance):

  forwardParams = {
    "inputTransform1File": forwardWarpPath,
    "inputTransform2Node": correctionsTransformNodeID,
    "inputReferenceVolumeFile" : templateReferencePath,
    "outputFileName" : forwardWarpPath
    } 

  inverseParams = {
    "inputTransform1Node": correctionsTransformNodeID,
    "inputTransform2File": inverseWarpPath,
    "inputReferenceVolumeFile" : nativeReferencePath,
    "outputFileName" : inverseWarpPath
    }

  forwardCliNode = slicer.cli.createNode(slicer.modules.compositetogridtransform, forwardParams)
  forwardCliNode.SetName('forwardCompositeToGrid')
  inverseCliNode = slicer.cli.createNode(slicer.modules.compositetogridtransform, inverseParams)
  inverseCliNode.SetName('inverseCompositeToGrid')

  subName = os.path.basename(os.path.dirname(subjectWarpDrivePath))
  tmpScenePath = os.path.join(subjectWarpDrivePath, 'tmpScene')
  python_commands = 'slicer.util.mainWindow().hide();\
                    w = slicer.modules.compositetogridtransform.createNewWidgetRepresentation();\
                    w.show();\
                    w.setWindowTitle(\'WarpDrive\');\
                    w.children()[2].hide();\
                    w.children()[4].hide();\
                    w.children()[5].hide();\
                    w.children()[6].hide();\
                    w.children()[7].hide();\
                    w.setCurrentCommandLineModuleNode(forwardCliNode);\
                    txt = ctk.ctkFittedTextBrowser(w);\
                    txt.setHtml(\'Subject: '+subName+'.<br><br>Saving changes to the normalization transformation files.<br><br>This window is independent of Slicer and Lead-DBS and will close when finished.\');\
                    w.children()[1].insertWidget(0,txt);\
                    w.resize(w.width,w.height/2);\
                    qt.QApplication.processEvents();\
                    forwardCliNode.AddObserver(\'ModifiedEvent\', lambda c,e,w=w,icli=inverseCliNode: [slicer.util.getNode(c.GetParameterAsString(\'inputTransform2Node\')).Inverse(), w.setCurrentCommandLineModuleNode(icli), slicer.cli.run(slicer.modules.compositetogridtransform, icli)] if (c.GetStatus() == c.Completed) else None);\
                    inverseCliNode.AddObserver(\'ModifiedEvent\', lambda c,e,w=w,: [shutil.rmtree(r\''+tmpScenePath+'\') if os.path.isdir(r\''+tmpScenePath+'\') else None, w.close(), slicer.util.exit(0)] if (c.GetStatus() == c.Completed) else None);\
                    slicer.cli.run(slicer.modules.compositetogridtransform, forwardCliNode);'

  if useExternalInstance:

    if os.path.isdir(tmpScenePath):
      shutil.rmtree(tmpScenePath)
    os.mkdir(tmpScenePath)

    tmpScene = slicer.vtkMRMLScene()
    slicer.mrmlScene.CopyDefaultNodesToScene(tmpScene)
    tmpScene.AddNode(forwardCliNode)
    tmpScene.AddNode(inverseCliNode)
    tmpScene.AddNode(slicer.util.getNode(correctionsTransformNodeID))
    tmpScene.SaveSceneToSlicerDataBundleDirectory(tmpScenePath)
    tmpScene.Clear()
    del tmpScene

    tmpScriptPath = os.path.join(subjectWarpDrivePath, 'tmpScript.py')
    with open(tmpScriptPath, 'w') as f:
      f.write('import os, shutil, ctk;\
                loadScene(r\''+os.path.join(tmpScenePath,'tmpScene.mrml')+'\');\
                forwardCliNode = slicer.mrmlScene.GetFirstNodeByName(\'forwardCompositeToGrid\');\
                inverseCliNode = slicer.mrmlScene.GetFirstNodeByName(\'inverseCompositeToGrid\');'\
                + python_commands + 'os.remove(r\''+tmpScriptPath+'\')')

    slicerInstallPath = os.path.dirname(os.path.dirname(sys.executable))
    if platform.system() == 'Darwin':
      slicerPath = os.path.join(slicerInstallPath, 'MacOS', slicer.app.mainApplicationName)
    elif platform.system() == 'Linux':
      slicerPath = os.path.join(slicerInstallPath, slicer.app.mainApplicationName)
    elif platform.system() == 'Windows':
      slicerPath = os.path.join(slicerInstallPath, slicer.app.mainApplicationName + '.exe')
    
    commands = [slicerPath, 
                '--ignore-slicerrc', 
                '--no-splash',
                '--python-script', tmpScriptPath]
    
    if slicer.app.mainApplicationName != 'SlicerForLeadDBS':
      slicerNetstimModule = glob.glob(os.path.join(slicerInstallPath,'**','SlicerNetstim','**','cli-modules'),recursive=True)[0]
      commands += ['--disable-settings', '--additional-module-paths', slicerNetstimModule]

    import subprocess
    subprocess.Popen(commands, env=slicer.util.startupEnvironment())
  
  else:
    exec(python_commands)


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
    folderItem = shNode.GetItemByDataNode(folderNode)
    if 'atlas' in shNode.GetItemAttributeNames(folderItem) and shNode.GetItemAttribute(folderItem,'atlas') == 'template':
      if  shNode.GetItemParent(shNode.GetItemByDataNode(folderNode)) == shNode.GetSceneItemID():
        names.append(folderNode.GetName())
  return names

def saveSceneInfo(warpDriveSavePath, inverseApplied):
  info = {}
  info["atlasNames"] = getAtlasesNamesInScene()
  info["inverseApplied"] = inverseApplied
  with open(os.path.join(warpDriveSavePath,'info.json'), 'w') as jsonFile:
    json.dump(info, jsonFile)

def saveSegmentation(segmentationsDir):
  shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
  modelNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLModelNode')
  modelNodes.UnRegister(slicer.mrmlScene)
  for i in range(modelNodes.GetNumberOfItems()):
    modelNode = modelNodes.GetItemAsObject(i)
    modelItem = shNode.GetItemByDataNode(modelNode)
    if modelNode.GetDisplayNode().GetVisibility() and 'atlas' in shNode.GetItemAttributeNames(modelItem) and shNode.GetItemAttribute(modelItem,'atlas') == 'template':
      # construct name
      name = modelNode.GetName()
      node = modelNode
      while shNode.GetItemParent(shNode.GetItemByDataNode(node)) != shNode.GetSceneItemID():
        node = shNode.GetItemDataNode(shNode.GetItemParent(shNode.GetItemByDataNode(node)))
        name = os.path.join(node.GetName(), name)
      # save segmentation
      modelNode.HardenTransform()
      labelNode = modelToLabel(modelNode)
      os.makedirs(os.path.dirname(name), exist_ok=True)
      slicer.util.saveNode(labelNode, os.path.join(segmentationsDir, name + '.nii.gz'))      
      slicer.mrmlScene.RemoveNode(labelNode)
      slicer.mrmlScene.RemoveNode(modelNode)

def modelToLabel(model_node):
  segmentationsLogic = slicer.modules.segmentations.logic()
  segmentNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentationNode', 'segmentation')
  labelNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLabelMapVolumeNode', 'label')
  segmentationsLogic.ImportModelToSegmentationNode(model_node, segmentNode)
  segmentationsLogic.ExportAllSegmentsToLabelmapNode(segmentNode, labelNode)
  slicer.mrmlScene.RemoveNode(segmentNode)
  return labelNode
