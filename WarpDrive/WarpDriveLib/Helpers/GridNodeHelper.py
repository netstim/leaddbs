import vtk, slicer


def getGridDefinition(node):

  if isinstance(node, slicer.vtkMRMLScalarVolumeNode):
    grid = node.GetImageData()
    # get origin (defined in volume node)
    origin = node.GetOrigin()

  elif isinstance(node, slicer.vtkMRMLGridTransformNode):
    fp = node.GetTransformFromParent()
    tp = node.GetTransformToParent()
    if isinstance(fp, slicer.vtkOrientedGridTransform) and fp.GetDisplacementGrid():
      grid = fp.GetDisplacementGrid()
    elif isinstance(tp, slicer.vtkOrientedGridTransform) and tp.GetDisplacementGrid():
      grid = tp.GetDisplacementGrid()
    # get origin
    origin = grid.GetOrigin()

  size = grid.GetDimensions()
  spacing = grid.GetSpacing()

  return size,origin,spacing

def getTransformRASToIJK(transformNode):
  size,origin,spacing = getGridDefinition(transformNode)
  IJKToRAS = [ 
              [spacing[0],          0 ,         0 ,  origin[0] ],
              [        0 ,  spacing[1],         0 ,  origin[1] ],
              [        0 ,          0 , spacing[2],  origin[2] ],
              [        0 ,          0 ,         0 ,         1  ],
              ] 
  m = vtk.vtkMatrix4x4()
  for i in range(4):
    for j in range(4):
      m.SetElement(i,j,IJKToRAS[i][j])
  m.Invert()
  return m

def emptyGridTransform(transformSize = [193,229,193], transformOrigin = [-96.0, -132.0, -78.0], transformSpacing = [1.0, 1.0, 1.0], transformNode = None):

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
  # Create transform node
  transformNode.SetAndObserveTransformFromParent(transform)
  transformNode.GetTransformFromParent().GetDisplacementGrid().SetOrigin(transformOrigin)
  transformNode.GetTransformFromParent().GetDisplacementGrid().SetSpacing(transformSpacing)

  return transformNode


def emptyVolume(imageSize, imageOrigin, imageSpacing):
  voxelType = vtk.VTK_UNSIGNED_CHAR
  imageDirections = [[1,0,0], [0,1,0], [0,0,1]]
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
  volumeNode.SetIJKToRASDirections(imageDirections)
  volumeNode.SetAndObserveImageData(imageData)
  volumeNode.CreateDefaultDisplayNodes()

  return volumeNode