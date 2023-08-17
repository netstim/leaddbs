import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import glob

import sys
import numpy as np

#
# ImportAtlas
#

class ImportAtlas(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "ImportAtlas" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Netstim"]
    self.parent.dependencies = ["NetstimPreferences"]
    self.parent.contributors = ["Simon Oxenford (Netstim Berlin)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
This module loads Lead-DBS atlases into Slicer. Set the Lead-DBS path in the Settings menu to see available atlases.
"""
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = "" # replace with organization, grant and thanks.

#
# ImportAtlasWidget
#

class ImportAtlasWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)


    #
    # atlas combo box
    #
    self.atlasComboBox = qt.QComboBox()
    parametersFormLayout.addRow("Atlas: ", self.atlasComboBox)

    #
    # Import Button
    #
    self.importButton = qt.QPushButton("Import")
    self.importButton.toolTip = "Import atlas."
    self.importButton.enabled = False
    parametersFormLayout.addRow(self.importButton)

    # connections
    self.importButton.connect('clicked(bool)', self.onImportButton)

    # Add vertical spacer
    self.layout.addStretch(1)


  def cleanup(self):
    pass

  def enter(self):
    validAtlases = ImportAtlasLogic().getValidAtlases()
    self.atlasComboBox.clear()
    self.atlasComboBox.addItems(validAtlases)
    self.importButton.enabled = len(validAtlases) > 0
    if 'DISTAL Minimal (Ewert 2017)' in validAtlases:
      self.atlasComboBox.setCurrentText('DISTAL Minimal (Ewert 2017)')
      self.importButton.setToolTip('')
    else:
      self.importButton.setToolTip('Set leaddbs path in setting menu')


  def onImportButton(self):
    logic = ImportAtlasLogic()
    atlasPath = os.path.join(logic.getAtlasesPath(), self.atlasComboBox.currentText)
    qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
    qt.QApplication.processEvents()
    try:
      logic.readAtlas(os.path.join(atlasPath,'atlas_index.mat'))
    finally:
      qt.QApplication.restoreOverrideCursor() 

#
# ImportAtlasLogic
#

class ImportAtlasLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """


  def getAtlasesPath(self):
    leadDBSPath = slicer.util.settingsValue("NetstimPreferences/leadDBSPath", "", converter=str)
    leadDBSSpace = slicer.util.settingsValue("NetstimPreferences/leadDBSSpace", "", converter=str)
    if leadDBSPath and leadDBSSpace:
      return os.path.join(leadDBSPath, "templates", "space", leadDBSSpace, "atlases")
    for possibleSpace in ["MNI152NLin2009bAsym", "MNI_ICBM_2009b_NLIN_ASYM"]:
      possiblePath = os.path.join(leadDBSPath, "templates", "space", possibleSpace, "atlases")
      if os.path.isdir(possiblePath):
        return possiblePath
    return ""

  def getValidAtlases(self):
    validAtlases = glob.glob(os.path.join(self.getAtlasesPath(), '*', 'atlas_index.mat'))
    if not validAtlases:
      qt.QMessageBox().warning(qt.QWidget(), "", "Invalid Lead-DBS path in preferences.")
      return []
    validAtlases = [os.path.basename(os.path.dirname(a)) for a in validAtlases]
    validAtlases.sort()
    return validAtlases

  def createFolderDisplayNode(self, folderID, color=[0.66,0.66,0.66], opacity=1.0):
    # from qSlicerSubjectHierarchyFolderPlugin.cxx
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    displayNode = slicer.vtkMRMLFolderDisplayNode()
    displayNode.SetName(shNode.GetItemName(folderID))
    displayNode.SetHideFromEditors(0)
    displayNode.SetAttribute('SubjectHierarchy.Folder', "1")
    displayNode.SetColor(*color)
    displayNode.SetOpacity(opacity)
    shNode.GetScene().AddNode(displayNode)
    shNode.SetItemDataNode(folderID, displayNode)
    shNode.ItemModified(folderID)

  def readAtlas(self, atlasPath, atlasName=None):
    """
    Run the actual algorithm
    """

    atlasName = atlasName if atlasName else os.path.basename(os.path.dirname(atlasPath))

    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    folderID = shNode.CreateFolderItem(shNode.GetSceneItemID(), atlasName)
    self.createFolderDisplayNode(folderID, opacity=0.8)
    shNode.SetItemAttribute(folderID, 'atlas', 'template')

    atlas = LeadDBSAtlas(atlasPath)

    for structure in atlas.structures:
      
      if structure.isBilateral():
        subName = ['rh', 'lh']
        sideIndexes = [0, 1]
        subFolderID = shNode.CreateFolderItem(folderID, structure.name)
        self.createFolderDisplayNode(subFolderID, structure.color)
        shNode.SetItemDisplayVisibility(subFolderID, structure.visibility)
        shNode.SetItemExpanded(subFolderID, 0)
        shNode.SetItemAttribute(subFolderID, 'atlas', 'template')
      else:
        subName = [structure.name]
        sideIndexes = [1] if structure.type == 2 else [0]
        subFolderID = folderID

      for sideIndex, sideName in zip(sideIndexes, subName):
        node = structure.getStructureNode(sideIndex)
        node.SetName(sideName)
        shNode.SetItemParent(shNode.GetItemChildWithName(shNode.GetSceneItemID(), sideName), subFolderID)
        shNode.SetItemAttribute(shNode.GetItemByDataNode(node), 'atlas', 'template')

    return folderID

#
# Atlas Structure
#

class LeadDBSAtlasStructure:
  def __init__(self):
    self.atlasPath = None
    self.index = None
    self.name = None
    self.color = None
    self.type = None
    self.visibility = None

  def isBilateral(self):
    return self.type in [3,4]

  def getStructureNode(self, sideIndex):
    pass

  def createNode(self, polyData, nodeType):
    node = slicer.mrmlScene.AddNewNodeByClass(nodeType)
    node.SetAndObservePolyData(polyData)
    node.CreateDefaultDisplayNodes()
    node.GetDisplayNode().SetColor(*self.color)
    node.GetDisplayNode().SetBackfaceCulling(0)
    node.GetDisplayNode().SetVisibility(self.visibility)
    return node

class ModelStructure(LeadDBSAtlasStructure):
  def __init__(self):
      super().__init__()

  def getStructureNode(self, sideIndex):
    faces, vertices = self.getFacesVertices(sideIndex)
    polyData = self.getPolyData(faces, vertices)
    modelNode = self.createNode(polyData)
    return modelNode

  def getFacesVertices(self, sideIndex):
    smooth = slicer.util.settingsValue("NetstimPreferences/useSmoothAtlas", True, converter=slicer.util.toBool)
    import h5py
    with h5py.File(self.atlasPath,'r') as atlasFile:
      roi = atlasFile['atlases']['roi']
      ref = roi[sideIndex][self.index]
      fvKey = 'sfv' if (smooth and 'sfv' in atlasFile[ref].keys()) else 'fv'
      fv = atlasFile[ref][fvKey]
      vertices = fv['vertices'][()].transpose()
      faces = fv['faces'][()].transpose()
    return faces, vertices

  def getPolyData(self, faces, vertices):
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

    # Compute surface normals for better appearance
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(pd)
    normals.ConsistencyOn()
    normals.AutoOrientNormalsOn()
    normals.Update()

    return normals.GetOutput()

  def createNode(self, polyData):
    node = super().createNode(polyData, 'vtkMRMLModelNode')
    node.GetDisplayNode().SetVisibility2D(1)
    return node


class FibersStructure(LeadDBSAtlasStructure):
  def __init__(self):
      super().__init__()
  
  def getStructureNode(self, sideIndex):
    points = self.getPointsIdx(sideIndex)
    polyData = self.getPolyData(points)
    fiberNode = self.createNode(polyData)
    return fiberNode

  def getPointsIdx(self, sideIndex):
    if self.isBilateral():
      subFolder = 'rh' if sideIndex == 0 else 'lh'
      fibersPath = os.path.join(os.path.dirname(self.atlasPath), subFolder, self.name+'.mat')
    else:
      import glob
      fibersPath = glob.glob(os.path.join(os.path.dirname(self.atlasPath), '*', self.name+'.mat'))[0]
    import h5py
    with h5py.File(fibersPath,'r') as fibersFile:
      points = fibersFile['fibers'][()].transpose()
    return points

  def getPolyData(self, points):
    vtkPoints = vtk.vtkPoints()
    lines = vtk.vtkCellArray()
    ptsIndex = 0
    lineIndex = 1
    finish = False
    while not finish:
      line = vtk.vtkPolyLine()
      while lineIndex == points[ptsIndex,3]:
        vtkPoints.InsertNextPoint(points[ptsIndex,0:3])
        line.GetPointIds().InsertNextId(ptsIndex)
        ptsIndex = ptsIndex + 1
        if ptsIndex == points.shape[0]:
          finish = True
          break
      lines.InsertNextCell(line)
      lineIndex = lineIndex + 1
    
    pd = vtk.vtkPolyData()
    pd.SetPoints(vtkPoints)
    pd.SetLines(lines)

    return pd

  def createNode(self, polyData):
    node = super().createNode(polyData, 'vtkMRMLFiberBundleNode')
    node.GetDisplayNode().SetColorModeToScalar()
    return node

class DiscFibersStructure(FibersStructure):
  def __init__(self):
      super().__init__()

  def getStructureNode(self, sideIndex):
    points, scalars = self.getPointsIdx(sideIndex)
    polyData = self.getPolyData(points, scalars)
    fiberNode = self.createNode(polyData)
    return fiberNode

  def getPointsIdx(self, sideIndex):
    import glob
    fibersPath = glob.glob(os.path.join(os.path.dirname(self.atlasPath), '*', self.name+'.mat'))[0]
    import h5py
    fiberCount = 1
    with h5py.File(fibersPath,'r') as fibersFile:
      for i in range(fibersFile['fibcell'].shape[0]):
        ref = fibersFile['fibcell'][i,0]
        for j in range(fibersFile[ref].shape[1]):
          fiber = fibersFile[fibersFile[ref][0,j]][()].transpose()
          fiber = np.concatenate((fiber, fiberCount*np.ones((fiber.shape[0],1))), axis=1)
          scalar = fibersFile[fibersFile['vals'][i,0]][0][j]
          if fiberCount == 1:
            points = fiber
            scalars = scalar * np.ones((fiber.shape[0],1))
          else:
            points = np.concatenate((points, fiber), axis=0)
            scalars = np.concatenate((scalars, scalar * np.ones((fiber.shape[0],1))), axis=0)
          fiberCount = fiberCount + 1
    return points, scalars

  def getPolyData(self, points, scalars):
    pd = super().getPolyData(points)
    vtkValuesArray = vtk.vtkDoubleArray()
    vtkValuesArray.SetName('values')
    for scalar in scalars:
      vtkValuesArray.InsertNextTuple((scalar,))
    pd.GetPointData().AddArray(vtkValuesArray)
    pd.GetPointData().SetScalars(vtkValuesArray)
    return pd

  def createNode(self, polyData):
    node = super().createNode(polyData)
    node.GetDisplayNode().SetVisibility(True)
    node.GetDisplayNode().SetColorModeToScalarData()
    node.GetDisplayNode().SetActiveScalarName('values')
    node.GetDisplayNode().SetAndObserveColorNodeID(slicer.util.getNode('DivergingBlueRed').GetID())
    return node

#
# Lead-DBS Atlas
#

class LeadDBSAtlas:
  def __init__(self, atlasPath):

    try:
      import h5py
    except:
      slicer.util.pip_install('h5py')
      import h5py
      
    with h5py.File(atlasPath,'r') as atlasFile:
      names = self.readNames(atlasFile)
      colors = self.readColors(atlasFile)
      types = self.readTypes(atlasFile)
      showIndex = self.readShowIndex(atlasFile)
      pixdimTypes = self.readPixdimType(atlasFile)

    self.structures = []

    for i, pixdimType in enumerate(pixdimTypes):
      if pixdimType == 'numeric':
        structure = ModelStructure()
      elif pixdimType == 'fibers':
        structure = FibersStructure()
      elif pixdimType == 'discfibers':
        structure = DiscFibersStructure()
      if isinstance(structure, FibersStructure) and not hasattr(slicer.modules,'tractographydisplay'):
        qt.QMessageBox().warning(qt.QWidget(), "Error", "Install SlicerDMRI Extension to load fiber atlases")
        continue
      structure.atlasPath = atlasPath
      structure.index = i
      structure.name = names[i]
      structure.color = colors[i]
      structure.type = types[i]
      structure.visibility = i in showIndex
      self.structures.append(structure)

  def readColors(self, atlasFile):
    colorIndex = atlasFile['atlases']['colors'][()].squeeze()
    if not colorIndex.shape:
      colorIndex = np.array([colorIndex])
    colorIndex[np.isnan(colorIndex)] = 1
    try:
      colormap = atlasFile['atlases']['colormap'][()].transpose()
    except: # colormap not present
      colormap = np.array([[0.2422,0.1504,0.6603],[0.2504,0.1650,0.7076],[0.2578,0.1818,0.7511],[0.2647,0.1978,0.7952],[0.2706,0.2147,0.8364],[0.2751,0.2342,0.8710],[0.2783,0.2559,0.8991],[0.2803,0.2782,0.9221],[0.2813,0.3006,0.9414],[0.2810,0.3228,0.9579],[0.2795,0.3447,0.9717],[0.2760,0.3667,0.9829],[0.2699,0.3892,0.9906],[0.2602,0.4123,0.9952],[0.2440,0.4358,0.9988],[0.2206,0.4603,0.9973],[0.1963,0.4847,0.9892],[0.1834,0.5074,0.9798],[0.1786,0.5289,0.9682],[0.1764,0.5499,0.9520],[0.1687,0.5703,0.9359],[0.1540,0.5902,0.9218],[0.1460,0.6091,0.9079],[0.1380,0.6276,0.8973],[0.1248,0.6459,0.8883],[0.1113,0.6635,0.8763],[0.0952,0.6798,0.8598],[0.0689,0.6948,0.8394],[0.0297,0.7082,0.8163],[0.0036,0.7203,0.7917],[0.0067,0.7312,0.7660],[0.0433,0.7411,0.7394],[0.0964,0.7500,0.7120],[0.1408,0.7584,0.6842],[0.1717,0.7670,0.6554],[0.1938,0.7758,0.6251],[0.2161,0.7843,0.5923],[0.2470,0.7918,0.5567],[0.2906,0.7973,0.5188],[0.3406,0.8008,0.4789],[0.3909,0.8029,0.4354],[0.4456,0.8024,0.3909],[0.5044,0.7993,0.3480],[0.5616,0.7942,0.3045],[0.6174,0.7876,0.2612],[0.6720,0.7793,0.2227],[0.7242,0.7698,0.1910],[0.7738,0.7598,0.1646],[0.8203,0.7498,0.1535],[0.8634,0.7406,0.1596],[0.9035,0.7330,0.1774],[0.9393,0.7288,0.2100],[0.9728,0.7298,0.2394],[0.9956,0.7434,0.2371],[0.9970,0.7659,0.2199],[0.9952,0.7893,0.2028],[0.9892,0.8136,0.1885],[0.9786,0.8386,0.1766],[0.9676,0.8639,0.1643],[0.9610,0.8890,0.1537],[0.9597,0.9135,0.1423],[0.9628,0.9373,0.1265],[0.9691,0.9606,0.1064],[0.9769,0.9839,0.0805]])
    colors = [colormap[int(i)-1] for i in colorIndex]
    return colors

  def readNames(self, atlasFile):
    names = []
    atlases = atlasFile['atlases']
    for column in atlases['names']:
      name = ''.join(map(chr, np.squeeze(atlases[column[0]][:])))
      names.append(name.split('.')[0])
    return names
  
  def readTypes(self, atlasFile):
    types = atlasFile['atlases']['types'][()].squeeze()
    if not types.shape:
      types = np.array([types])
    return types

  def readShowIndex(self, atlasFile):
    import h5py
    try:
      if isinstance(atlasFile['atlases']['presets']['show'][0,0], h5py.h5r.Reference):
        showIndexRef = atlasFile['atlases']['presets']['show'][0,0]
        showIndex = atlasFile[showIndexRef][()].squeeze() - 1 # -1 to fix index base
      else:
        showIndex = atlasFile['atlases']['presets']['show'][()].squeeze() - 1 
    except:
      showIndex = np.array(range(len(atlasFile['atlases']['pixdim'][0])))
    if not showIndex.shape:
      showIndex = np.array([showIndex])
    return showIndex

  def readPixdimType(self, atlasFile):
    pixdimType = []
    for i in range(len(atlasFile['atlases']['pixdim'][0])):
      ref = atlasFile['atlases']['pixdim'][0,i]
      data = atlasFile[ref][()]
      if data.dtype == np.dtype('uint16'):
        pixdimType.append(''.join(map(chr, np.squeeze(data))))
      else:
        pixdimType.append('numeric')
    return pixdimType


#
# File Reader
#

class ImportAtlasFileReader:

  def __init__(self, parent):
    self.parent = parent

  def description(self):
    return 'Lead-DBS atlas'

  def fileType(self):
    return 'LeadDBSAtlas'

  def extensions(self):
    return ['Lead-DBS atlas (*.mat)']

  def canLoadFile(self, filePath):
    # filename must be atlas_index.mat
    filename = os.path.split(filePath)[1]
    return filename == "atlas_index.mat"

  def load(self, properties):
    try:

      # Import data
      filePath = properties['fileName']
      logic = ImportAtlasLogic()
      shFolderItem = logic.readAtlas(filePath)

      # Create list of loaded noed IDs
      loadedNodeIDs = []
      shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
      childrenIdList = vtk.vtkIdList()
      shNode.GetItemChildren(shFolderItem, childrenIdList, True)
      for childIndex in range(childrenIdList.GetNumberOfIds()):
        dataNode = shNode.GetItemDataNode(childrenIdList.GetId(childIndex))
        if dataNode:
          loadedNodeIDs.append(dataNode.GetID())

    except Exception as e:
      logging.error('Failed to load file: '+str(e))
      import traceback
      traceback.print_exc()
      return False

    self.parent.loadedNodes = loadedNodeIDs
    return True

#
# Test
#

class ImportAtlasTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_ImportAtlas1()

  def test_ImportAtlas1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    #import SampleData
    #SampleData.downloadFromURL(
    #  nodeNames='FA',
    #  fileNames='FA.nrrd',
    #  uris='http://slicer.kitware.com/midas3/download?items=5767')
    #self.delayDisplay('Finished with download and loading')
    #
    #volumeNode = slicer.util.getNode(pattern="FA")
    #logic = ImportAtlasLogic()
    #self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
