import vtk, slicer
import numpy as np

def getGridDefinition(node):

  if isinstance(node, slicer.vtkMRMLScalarVolumeNode):
    directionMatrix = vtk.vtkMatrix4x4()
    node.GetIJKToRASDirectionMatrix(directionMatrix)
    grid = node.GetImageData()
    size = np.array(grid.GetDimensions())
    origin = np.array(node.GetOrigin())
    spacing = np.array(node.GetSpacing())

  elif isinstance(node, slicer.vtkMRMLTransformNode):
    fp = node.GetTransformFromParent()
    tp = node.GetTransformToParent()
    if isinstance(fp, slicer.vtkOrientedGridTransform) and fp.GetDisplacementGrid():
      grid = fp.GetDisplacementGrid()
      directionMatrix = slicer.vtkOrientedGridTransform().SafeDownCast(fp).GetGridDirectionMatrix()
    elif isinstance(tp, slicer.vtkOrientedGridTransform) and tp.GetDisplacementGrid():
      grid = tp.GetDisplacementGrid()
      directionMatrix = slicer.vtkOrientedGridTransform().SafeDownCast(tp).GetGridDirectionMatrix()
    origin = np.array(grid.GetOrigin())
    size = np.array(grid.GetDimensions())
    spacing = np.array(grid.GetSpacing())

  return size,origin,spacing,directionMatrix

def getTransformRASToIJK(transformNode):
  size,origin,spacing,directionMatrix = getGridDefinition(transformNode)
  m = vtk.vtkMatrix4x4()
  for row in range(3):
    m.SetElement(row,3,origin[row])
    for col in range(3):
      m.SetElement(row, col, spacing[col] * directionMatrix.GetElement(row,col))
  m.Invert()
  return m

def emptyGridTransform(transformSize, transformOrigin, transformSpacing, directionMatrix = None, transformNode = None):

  if not directionMatrix:
    directionMatrix = vtk.vtkMatrix4x4()

  if not transformNode:
    transformNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLGridTransformNode')

  voxelType = vtk.VTK_FLOAT
  fillVoxelValue = 0
  # Create an empty image volume, filled with fillVoxelValue
  imageData = vtk.vtkImageData()
  imageData.SetDimensions(transformSize)
  imageData.AllocateScalars(voxelType, 3)
  imageData.GetPointData().GetScalars().Fill(fillVoxelValue)
  # Create transform
  transform = slicer.vtkOrientedGridTransform()
  transform.SetInterpolationModeToCubic()
  transform.SetDisplacementGridData(imageData)
  transform.SetGridDirectionMatrix(directionMatrix)
  # Create transform node
  transformNode.SetAndObserveTransformFromParent(transform)
  transformNode.GetTransformFromParent().GetDisplacementGrid().SetOrigin(transformOrigin)
  transformNode.GetTransformFromParent().GetDisplacementGrid().SetSpacing(transformSpacing)

  return transformNode


def emptyVolume(imageSize, imageOrigin, imageSpacing, directionMatrix):
  voxelType = vtk.VTK_UNSIGNED_CHAR
  fillVoxelValue = 0
  # Create an empty image volume, filled with fillVoxelValue
  imageData = vtk.vtkImageData()
  imageData.SetDimensions(imageSize)
  imageData.AllocateScalars(voxelType, 1)
  imageData.GetPointData().GetScalars().Fill(fillVoxelValue)
  # Create volume node
  volumeNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
  volumeNode.SetOrigin(imageOrigin)
  volumeNode.SetSpacing(imageSpacing)
  volumeNode.SetIJKToRASDirectionMatrix(directionMatrix)
  volumeNode.SetAndObserveImageData(imageData)
  volumeNode.CreateDefaultDisplayNodes()

  return volumeNode