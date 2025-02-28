import os
import slicer

'''
These variables are defined from Lead-DBS when launching this script:
    nativeToMNITransformFile
    MNIToNativeTransformFile
    MNIReferenceFile
'''

MNIToNativeTransformNode = slicer.util.loadTransform(MNIToNativeTransformFile)
MNIToNativeTransformNode.Inverse()

params = {
    "inputTransform1Node": MNIToNativeTransformNode.GetID(),
    "inputTransform2Node": slicer.mrmlScene.AddNewNodeByClass('vtkMRMLGridTransformNode').GetID(),
    "inputReferenceVolumeFile" : MNIReferenceFile,
    "outputFileName" : nativeToMNITransformFile
} 

slicer.cli.run(slicer.modules.compositetogridtransform, None, params, wait_for_completion=True, update_display=False)

slicer.util.exit()
