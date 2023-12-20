
import slicer
import sys

fwd_transform_file = sys.argv[1]
reference_volume_file = sys.argv[2]
inv_transform_file = sys.argv[3]

in_transform_node = slicer.util.loadTransform(fwd_transform_file)
out_transform_node = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
reference_volume_node = slicer.util.loadVolume(reference_volume_file)
in_transform_node.Inverse()
slicer.modules.transforms.logic().ConvertToGridTransform(in_transform_node,reference_volume_node,out_transform_node)
slicer.util.saveNode(out_transform_node, inv_transform_file)

exit()
