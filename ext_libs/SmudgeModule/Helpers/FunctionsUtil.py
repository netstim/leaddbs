import vtk
import os
import slicer
import SimpleITK as sitk
import sitkUtils
import numpy as np
import glob
import sys

try:
  import h5py
  import hdf5storage
  from scipy import io
except:
  import importlib  
  for package in ['h5py','hdf5storage','scipy']: 
    if sys.version_info[0] < 3: 
      from pip._internal import main as pipmain
      failed = pipmain(['install', package])
    else:
      failed = slicer.util.pip_install(package)
    if not failed:
      importlib.import_module(package)


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


def loadAtlas(atlasPath, modelHierarchyNode=None):

  if modelHierarchyNode:
    modelHierarchyNode.RemoveAllChildrenNodes()
  else:
    modelHierarchyNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelHierarchyNode')
  
  with h5py.File(os.path.join(atlasPath,'atlas_index.mat'),'r') as atlasFile:
    fv = atlasFile['atlases']['fv']

    colors = atlasFile['atlases']['colors'][()]   

    try:
      colormap = atlasFile['atlases']['colormap'][()].transpose()
    except: # colormap not present
      #if colors[-1][0]==64:
      colormap = np.array([[0.2422,0.1504,0.6603],[0.2504,0.1650,0.7076],[0.2578,0.1818,0.7511],[0.2647,0.1978,0.7952],[0.2706,0.2147,0.8364],[0.2751,0.2342,0.8710],[0.2783,0.2559,0.8991],[0.2803,0.2782,0.9221],[0.2813,0.3006,0.9414],[0.2810,0.3228,0.9579],[0.2795,0.3447,0.9717],[0.2760,0.3667,0.9829],[0.2699,0.3892,0.9906],[0.2602,0.4123,0.9952],[0.2440,0.4358,0.9988],[0.2206,0.4603,0.9973],[0.1963,0.4847,0.9892],[0.1834,0.5074,0.9798],[0.1786,0.5289,0.9682],[0.1764,0.5499,0.9520],[0.1687,0.5703,0.9359],[0.1540,0.5902,0.9218],[0.1460,0.6091,0.9079],[0.1380,0.6276,0.8973],[0.1248,0.6459,0.8883],[0.1113,0.6635,0.8763],[0.0952,0.6798,0.8598],[0.0689,0.6948,0.8394],[0.0297,0.7082,0.8163],[0.0036,0.7203,0.7917],[0.0067,0.7312,0.7660],[0.0433,0.7411,0.7394],[0.0964,0.7500,0.7120],[0.1408,0.7584,0.6842],[0.1717,0.7670,0.6554],[0.1938,0.7758,0.6251],[0.2161,0.7843,0.5923],[0.2470,0.7918,0.5567],[0.2906,0.7973,0.5188],[0.3406,0.8008,0.4789],[0.3909,0.8029,0.4354],[0.4456,0.8024,0.3909],[0.5044,0.7993,0.3480],[0.5616,0.7942,0.3045],[0.6174,0.7876,0.2612],[0.6720,0.7793,0.2227],[0.7242,0.7698,0.1910],[0.7738,0.7598,0.1646],[0.8203,0.7498,0.1535],[0.8634,0.7406,0.1596],[0.9035,0.7330,0.1774],[0.9393,0.7288,0.2100],[0.9728,0.7298,0.2394],[0.9956,0.7434,0.2371],[0.9970,0.7659,0.2199],[0.9952,0.7893,0.2028],[0.9892,0.8136,0.1885],[0.9786,0.8386,0.1766],[0.9676,0.8639,0.1643],[0.9610,0.8890,0.1537],[0.9597,0.9135,0.1423],[0.9628,0.9373,0.1265],[0.9691,0.9606,0.1064],[0.9769,0.9839,0.0805]])
      #else:
      #  colormap = np.array([[ 0, 0,0.5156],[ 0, 0,0.5312],[ 0, 0,0.5469],[ 0, 0,0.5625],[ 0, 0,0.5781],[ 0, 0,0.5938],[ 0, 0,0.6094],[ 0, 0,0.6250],[ 0, 0,0.6406],[ 0, 0,0.6562],[ 0, 0,0.6719],[ 0, 0,0.6875],[ 0, 0,0.7031],[ 0, 0,0.7188],[ 0, 0,0.7344],[ 0, 0,0.7500],[ 0, 0,0.7656],[ 0, 0,0.7812],[ 0, 0,0.7969],[ 0, 0,0.8125],[ 0, 0,0.8281],[ 0, 0,0.8438],[ 0, 0,0.8594],[ 0, 0,0.8750],[ 0, 0,0.8906],[ 0, 0,0.9062],[ 0, 0,0.9219],[ 0, 0,0.9375],[ 0, 0,0.9531],[ 0, 0,0.9688],[ 0, 0,0.9844],[ 0, 0,1.0000],[ 0,0.0156,1.0000],[ 0,0.0312,1.0000],[ 0,0.0469,1.0000],[ 0,0.0625,1.0000],[ 0,0.0781,1.0000],[ 0,0.0938,1.0000],[ 0,0.1094,1.0000],[ 0,0.1250,1.0000],[ 0,0.1406,1.0000],[ 0,0.1562,1.0000],[ 0,0.1719,1.0000],[ 0,0.1875,1.0000],[ 0,0.2031,1.0000],[ 0,0.2188,1.0000],[ 0,0.2344,1.0000],[ 0,0.2500,1.0000],[ 0,0.2656,1.0000],[ 0,0.2812,1.0000],[ 0,0.2969,1.0000],[ 0,0.3125,1.0000],[ 0,0.3281,1.0000],[ 0,0.3438,1.0000],[ 0,0.3594,1.0000],[ 0,0.3750,1.0000],[ 0,0.3906,1.0000],[ 0,0.4062,1.0000],[ 0,0.4219,1.0000],[ 0,0.4375,1.0000],[ 0,0.4531,1.0000],[ 0,0.4688,1.0000],[ 0,0.4844,1.0000],[ 0,0.5000,1.0000],[ 0,0.5156,1.0000],[ 0,0.5312,1.0000],[ 0,0.5469,1.0000],[ 0,0.5625,1.0000],[ 0,0.5781,1.0000],[ 0,0.5938,1.0000],[ 0,0.6094,1.0000],[ 0,0.6250,1.0000],[ 0,0.6406,1.0000],[ 0,0.6562,1.0000],[ 0,0.6719,1.0000],[ 0,0.6875,1.0000],[ 0,0.7031,1.0000],[ 0,0.7188,1.0000],[ 0,0.7344,1.0000],[ 0,0.7500,1.0000],[ 0,0.7656,1.0000],[ 0,0.7812,1.0000],[ 0,0.7969,1.0000],[ 0,0.8125,1.0000],[ 0,0.8281,1.0000],[ 0,0.8438,1.0000],[ 0,0.8594,1.0000],[ 0,0.8750,1.0000],[ 0,0.8906,1.0000],[ 0,0.9062,1.0000],[ 0,0.9219,1.0000],[ 0,0.9375,1.0000],[ 0,0.9531,1.0000],[ 0,0.9688,1.0000],[ 0,0.9844,1.0000],[ 0,1.0000,1.0000],[0.0156,1.0000,0.9844],[0.0312,1.0000,0.9688],[0.0469,1.0000,0.9531],[0.0625,1.0000,0.9375],[0.0781,1.0000,0.9219],[0.0938,1.0000,0.9062],[0.1094,1.0000,0.8906],[0.1250,1.0000,0.8750],[0.1406,1.0000,0.8594],[0.1562,1.0000,0.8438],[0.1719,1.0000,0.8281],[0.1875,1.0000,0.8125],[0.2031,1.0000,0.7969],[0.2188,1.0000,0.7812],[0.2344,1.0000,0.7656],[0.2500,1.0000,0.7500],[0.2656,1.0000,0.7344],[0.2812,1.0000,0.7188],[0.2969,1.0000,0.7031],[0.3125,1.0000,0.6875],[0.3281,1.0000,0.6719],[0.3438,1.0000,0.6562],[0.3594,1.0000,0.6406],[0.3750,1.0000,0.6250],[0.3906,1.0000,0.6094],[0.4062,1.0000,0.5938],[0.4219,1.0000,0.5781],[0.4375,1.0000,0.5625],[0.4531,1.0000,0.5469],[0.4688,1.0000,0.5312],[0.4844,1.0000,0.5156],[0.5000,1.0000,0.5000],[0.5156,1.0000,0.4844],[0.5312,1.0000,0.4688],[0.5469,1.0000,0.4531],[0.5625,1.0000,0.4375],[0.5781,1.0000,0.4219],[0.5938,1.0000,0.4062],[0.6094,1.0000,0.3906],[0.6250,1.0000,0.3750],[0.6406,1.0000,0.3594],[0.6562,1.0000,0.3438],[0.6719,1.0000,0.3281],[0.6875,1.0000,0.3125],[0.7031,1.0000,0.2969],[0.7188,1.0000,0.2812],[0.7344,1.0000,0.2656],[0.7500,1.0000,0.2500],[0.7656,1.0000,0.2344],[0.7812,1.0000,0.2188],[0.7969,1.0000,0.2031],[0.8125,1.0000,0.1875],[0.8281,1.0000,0.1719],[0.8438,1.0000,0.1562],[0.8594,1.0000,0.1406],[0.8750,1.0000,0.1250],[0.8906,1.0000,0.1094],[0.9062,1.0000,0.0938],[0.9219,1.0000,0.0781],[0.9375,1.0000,0.0625],[0.9531,1.0000,0.0469],[0.9688,1.0000,0.0312],[0.9844,1.0000,0.0156],[1.0000,1.0000, 0],[1.0000,0.9844, 0],[1.0000,0.9688, 0],[1.0000,0.9531, 0],[1.0000,0.9375, 0],[1.0000,0.9219, 0],[1.0000,0.9062, 0],[1.0000,0.8906, 0],[1.0000,0.8750, 0],[1.0000,0.8594, 0],[1.0000,0.8438, 0],[1.0000,0.8281, 0],[1.0000,0.8125, 0],[1.0000,0.7969, 0],[1.0000,0.7812, 0],[1.0000,0.7656, 0],[1.0000,0.7500, 0],[1.0000,0.7344, 0],[1.0000,0.7188, 0],[1.0000,0.7031, 0],[1.0000,0.6875, 0],[1.0000,0.6719, 0],[1.0000,0.6562, 0],[1.0000,0.6406, 0],[1.0000,0.6250, 0],[1.0000,0.6094, 0],[1.0000,0.5938, 0],[1.0000,0.5781, 0],[1.0000,0.5625, 0],[1.0000,0.5469, 0],[1.0000,0.5312, 0],[1.0000,0.5156, 0],[1.0000,0.5000, 0],[1.0000,0.4844, 0],[1.0000,0.4688, 0],[1.0000,0.4531, 0],[1.0000,0.4375, 0],[1.0000,0.4219, 0],[1.0000,0.4062, 0],[1.0000,0.3906, 0],[1.0000,0.3750, 0],[1.0000,0.3594, 0],[1.0000,0.3438, 0],[1.0000,0.3281, 0],[1.0000,0.3125, 0],[1.0000,0.2969, 0],[1.0000,0.2812, 0],[1.0000,0.2656, 0],[1.0000,0.2500, 0],[1.0000,0.2344, 0],[1.0000,0.2188, 0],[1.0000,0.2031, 0],[1.0000,0.1875, 0],[1.0000,0.1719, 0],[1.0000,0.1562, 0],[1.0000,0.1406, 0],[1.0000,0.1250, 0],[1.0000,0.1094, 0],[1.0000,0.0938, 0],[1.0000,0.0781, 0],[1.0000,0.0625, 0],[1.0000,0.0469, 0],[1.0000,0.0312, 0],[1.0000,0.0156, 0],[1.0000, 0, 0],[0.9844, 0, 0],[0.9688, 0, 0],[0.9531, 0, 0],[0.9375, 0, 0],[0.9219, 0, 0],[0.9062, 0, 0],[0.8906, 0, 0],[0.8750, 0, 0],[0.8594, 0, 0],[0.8438, 0, 0],[0.8281, 0, 0],[0.8125, 0, 0],[0.7969, 0, 0],[0.7812, 0, 0],[0.7656, 0, 0],[0.7500, 0, 0],[0.7344, 0, 0],[0.7188, 0, 0],[0.7031, 0, 0],[0.6875, 0, 0],[0.6719, 0, 0],[0.6562, 0, 0],[0.6406, 0, 0],[0.6250, 0, 0],[0.6094, 0, 0],[0.5938, 0, 0],[0.5781, 0, 0],[0.5625, 0, 0],[0.5469, 0, 0],[0.5312, 0, 0],[0.5156, 0, 0],[0.5000, 0, 0]])

    names = []
    atlases = atlasFile['atlases']
    for column in atlases['names']:
      row_data = []
      for row_number in range(len(column)):     
        if sys.version_info[0] < 3: 
          row_data.append(''.join(map(unichr, atlases[column[row_number]][:])))   
        else:
          row_data.append(''.join(map(chr, atlases[column[row_number]][:]))) 
      names.append(row_data)

    types = atlasFile['atlases']['types'][()]

    for index in range(len(names)):

      structureName = os.path.splitext(os.path.splitext(names[index][0])[0])[0]

      if types[index][0] in [3,4]:
        endName = ['_rh', '_lh']
      else:
        endName = ['']

      for sideIndex,sideName in zip(range(len(endName)),endName):

        ref = fv[sideIndex][index]
        b = atlasFile[ref]

        vertices = b['vertices'][()].transpose()
        faces = b['faces'][()].transpose()

        points = vtk.vtkPoints()
        for v in vertices:
          points.InsertNextPoint(*v)

        triangles = vtk.vtkCellArray()
        for f in faces:
          triangle = vtk.vtkTriangle()
          for j in range(3):
            triangle.GetPointIds().SetId(j,int(f[j])-1) # fix 1based matlab index
          triangles.InsertNextCell(triangle)

        pd = vtk.vtkPolyData()
        pd.SetPoints(points)
        pd.SetPolys(triangles)

        m = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode', structureName + sideName)
        m.SetAndObservePolyData(pd)
        m.CreateDefaultDisplayNodes()
        m.GetDisplayNode().SetColor(*colormap[int(colors[index])-1])

        mh = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelHierarchyNode')
        mh.SetModelNodeID(m.GetID())
        mh.SetParentNodeID(modelHierarchyNode.GetID())

  return modelHierarchyNode



def emptyTransform(transformNode):
  # init
  transformSize = [193,229,193]
  voxelType=vtk.VTK_FLOAT
  transformOrigin = [-96.0, -132.0, -78.0]
  transformSpacing = [1.0, 1.0, 1.0]
  transformPathections = [[1,0,0], [0,1,0], [0,0,1]]
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
  #transformNode.CreateDefaultDisplayNodes()
  transformNode.CreateDefaultStorageNode()  