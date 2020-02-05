import vtk

def getTransformRASToIJK(transformNode):
  origin = transformNode.GetTransformFromParent().GetDisplacementGrid().GetOrigin()
  spacing = transformNode.GetTransformFromParent().GetDisplacementGrid().GetSpacing()
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