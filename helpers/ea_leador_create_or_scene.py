import os
import slicer
import glob
import ImportAtlas
import StereotacticPlanLib
import StereotacticPlanLib.ImportFrom
import importlib

'''
These variables are defined from Lead-DBS when launching this script:
    subjectAnatFiles
    atlasPath
    subjectDICOMFolder
    MNIToAnchorNativeFile
'''

def getImporterData(folder):
    possibleImporterExt = ['.ros', '.pdf']
    possibleImporterName = ['Import_From_ROSA', 'Import_From_Brainlab']
    for ext,name in zip(possibleImporterExt,possibleImporterName):
        g = glob.glob(os.path.join(folder,'*' + ext))
        if len(g):
            return [name] + g
    print('No planning files available in: ' + folder)
    slicer.util.exit()

subjectLeadORFolder = os.path.dirname(__file__)
subjectFolder = os.path.dirname(subjectLeadORFolder)
subjectID = os.path.basename(subjectFolder)
savePrefix = subjectID + '_desc-'

#
# Load planning information
#

importerData = getImporterData(subjectLeadORFolder)
importerName = importerData[0]
planningFiles = importerData[1:]

importlib.import_module('.'.join(['StereotacticPlanLib', 'ImportFrom', importerName]))
importerModule = getattr(StereotacticPlanLib.ImportFrom, importerName)

referenceLoaded = False

for i,file in enumerate(planningFiles):
    importer = importerModule.Importer(file)
    loadedNodeIDs = importer.getTrajectoryTransforms(importInFrameSpace=True)
    for nodeID in loadedNodeIDs:
        node = slicer.util.getNode(nodeID)
        slicer.util.saveNode(node, os.path.join(subjectLeadORFolder, savePrefix + 'planning' + node.GetName().replace(' ','') + '.txt'))
    if not referenceLoaded:
        try:
            referenceToFrameTransformNode = importer.getReferenceToFrameTransform()
            referenceVolumeNode = importer.getReferenceVolumeFromDICOM(subjectDICOMFolder)
            slicer.util.saveNode(referenceToFrameTransformNode, os.path.join(subjectLeadORFolder, savePrefix + 'referenceToFrame.txt'))
            slicer.util.saveNode(referenceVolumeNode, os.path.join(subjectLeadORFolder, savePrefix + 'reference.nii'))
            referenceLoaded = True
        except:
            referenceLoaded = False

if not referenceLoaded:
    print('Couldn''t load reference for: ' + subjectFolder)
    slicer.util.exit()

#
# Load subject anchorNative volumes and transform to frame space
#

subjectAnatNodes = [slicer.util.loadVolume(file) for file in subjectAnatFiles]

anchorNativeToReferenceTransformNode  = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLinearTransformNode')
parameters = {}
parameters['fixedVolume'] 			  = referenceVolumeNode.GetID()
parameters['movingVolume'] 			  = subjectAnatNodes[0].GetID()
parameters['samplingPercentage'] 	  = 0.02
parameters['linearTransform'] 		  = anchorNativeToReferenceTransformNode.GetID()
parameters['initializeTransformMode'] = 'useMomentsAlign'
parameters['useRigid'] 				  = True
parameters['costMetric'] 			  = 'MMI'
cli = slicer.cli.run(slicer.modules.brainsfit, None, parameters, wait_for_completion=True, update_display=False)
slicer.util.saveNode(anchorNativeToReferenceTransformNode, os.path.join(subjectLeadORFolder, savePrefix + 'anchorNativeToReference.txt'))

for node in subjectAnatNodes:
    node.ApplyTransform(anchorNativeToReferenceTransformNode.GetTransformToParent())
    node.ApplyTransform(referenceToFrameTransformNode.GetTransformToParent())
    node.HardenTransform()

#
# Load atlas and transform to frame space
#

ImportAtlas.ImportAtlasLogic().readAtlas(atlasPath)

MNIToAnchorNativeFileTransformNode = slicer.util.loadTransform(MNIToAnchorNativeFile)

shnode = slicer.mrmlScene.GetSubjectHierarchyNode()
for i in range(slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')):
    modelNode = slicer.mrmlScene.GetNthNodeByClass(i, 'vtkMRMLModelNode')
    if shnode.GetItemAttribute(shnode.GetItemByDataNode(modelNode),'atlas') == '1':
        modelNode.ApplyTransform(MNIToAnchorNativeFileTransformNode.GetTransformToParent())
        modelNode.ApplyTransform(anchorNativeToReferenceTransformNode.GetTransformToParent())
        modelNode.ApplyTransform(referenceToFrameTransformNode.GetTransformToParent())
        modelNode.HardenTransform()

#
# Clean-up and save Scene
#

slicer.mrmlScene.RemoveNode(MNIToAnchorNativeFileTransformNode)
slicer.mrmlScene.RemoveNode(anchorNativeToReferenceTransformNode)
slicer.mrmlScene.RemoveNode(referenceToFrameTransformNode)
slicer.mrmlScene.RemoveNode(referenceVolumeNode)

sceneSaveFilename = os.path.join(subjectLeadORFolder, savePrefix + 'ORScene.mrb')
slicer.util.saveScene(sceneSaveFilename)

slicer.util.exit()
