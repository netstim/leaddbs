import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging

import sys
import numpy as np

try:
  import h5py
except:
  slicer.util.pip_install('h5py')
  import h5py


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
    self.parent.dependencies = []
    self.parent.contributors = ["John Doe (AnyWare Corp.)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
This is an example of scripted loadable module bundled in an extension.
It performs a simple thresholding on the input volume and optionally captures a screenshot.
"""
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

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
    # atlases directory
    #
    self.atlasDirectoryButton = ctk.ctkDirectoryButton()
    parametersFormLayout.addRow("Directory: ", self.atlasDirectoryButton)

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
    self.atlasDirectoryButton.directoryChanged.connect(self.onAtlasDirectoryChanged)
    #self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    #self.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh dir

    modulePath = os.path.split(__file__)[0]
    with open(os.path.join(modulePath,'Resources','previousDirectory.txt'), 'r') as f: 
      directory = f.readlines()[0]
    self.atlasDirectoryButton.directory = directory if os.path.isdir(directory) else '.'

  def cleanup(self):
    pass

  def onAtlasDirectoryChanged(self, directory):
    # remove atlases
    self.atlasComboBox.clear()  

    # add new dirs if valid lead paths
    self.atlasComboBox.addItems(ImportAtlasLogic().getValidAtlases(directory))
    
    # if MNI use DISTAL as default
    parentDir, dirName = os.path.split(directory)
    parentDirName = os.path.split(parentDir)[-1]
    if parentDirName == 'MNI_ICBM_2009b_NLIN_ASYM':
      self.atlasComboBox.setCurrentText('DISTAL Minimal (Ewert 2017)')
    
    # change text so its no so large
    self.atlasDirectoryButton.text = os.path.join('[...]',parentDirName,dirName)

    # enable import button if available atlases
    self.importButton.enabled = self.atlasComboBox.itemText(0) != ''

    # save for future
    modulePath = os.path.split(__file__)[0]
    with open(os.path.join(modulePath,'Resources','previousDirectory.txt'), 'w') as f: 
      f.write(directory) 

  def onImportButton(self):
    logic = ImportAtlasLogic()
    atlasPath = os.path.join(self.atlasDirectoryButton.directory, self.atlasComboBox.currentText)
    logic.run(atlasPath)

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

  def getValidAtlases(self, directory):
    validAtlases = []
    atlases = sorted(os.listdir(directory))
    for atlas in atlases:
      if atlas[0] != '.' and os.path.isfile(os.path.join(directory,atlas,'atlas_index.mat')):
        validAtlases.append(atlas)
    return validAtlases

  def getAtlasNames(self,atlasFile):
    """
    Get names from atlas.mat. This is a bit different from the others
    """
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
    return names

  def createPolyData(self, vertices, faces):
    """
     Generate vtk polydata from vertices and faces
    """
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

    return pd

  def createModel(self, polydata, name, color):
    """
    create model node with polydata, color, name and set default view mode
    """
    modelNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode', name )
    modelNode.SetAndObservePolyData(polydata)
    modelNode.CreateDefaultDisplayNodes()
    modelNode.GetDisplayNode().SetColor(*color)
    modelNode.GetDisplayNode().SetVisibility2D(1)
    modelNode.GetDisplayNode().SetBackfaceCulling(0)
    modelNode.GetDisplayNode().SetOpacity(0.8)
    modelNode.GetDisplayNode().SetVisibility(0) # hide by default
    return modelNode

  def createFolderDisplayNode(self, folderID, color=[0.66,0.66,0.66]):
    # from qSlicerSubjectHierarchyFolderPlugin.cxx
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    displayNode = slicer.vtkMRMLFolderDisplayNode()
    displayNode.SetName(shNode.GetItemName(folderID))
    displayNode.SetHideFromEditors(0)
    displayNode.SetAttribute('SubjectHierarchy.Folder', "1")
    displayNode.SetColor(*color)
    shNode.GetScene().AddNode(displayNode)
    shNode.SetItemDataNode(folderID, displayNode)
    shNode.ItemModified(folderID)

  def run(self, atlasPath):
    """
    Run the actual algorithm
    """
    qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
    qt.QApplication.processEvents()
  
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    folderID = shNode.CreateFolderItem(shNode.GetSceneItemID(), os.path.split(atlasPath)[-1])
    self.createFolderDisplayNode(folderID)
    shNode.SetItemAttribute(folderID, 'atlas', '1')
  
    with h5py.File(os.path.join(atlasPath,'atlas_index.mat'),'r') as atlasFile:
      # get .mat data
      roi = atlasFile['atlases']['roi']
      colors = atlasFile['atlases']['colors'][()]   
      try:
        colormap = atlasFile['atlases']['colormap'][()].transpose()
      except: # colormap not present
        colormap = np.array([[0.2422,0.1504,0.6603],[0.2504,0.1650,0.7076],[0.2578,0.1818,0.7511],[0.2647,0.1978,0.7952],[0.2706,0.2147,0.8364],[0.2751,0.2342,0.8710],[0.2783,0.2559,0.8991],[0.2803,0.2782,0.9221],[0.2813,0.3006,0.9414],[0.2810,0.3228,0.9579],[0.2795,0.3447,0.9717],[0.2760,0.3667,0.9829],[0.2699,0.3892,0.9906],[0.2602,0.4123,0.9952],[0.2440,0.4358,0.9988],[0.2206,0.4603,0.9973],[0.1963,0.4847,0.9892],[0.1834,0.5074,0.9798],[0.1786,0.5289,0.9682],[0.1764,0.5499,0.9520],[0.1687,0.5703,0.9359],[0.1540,0.5902,0.9218],[0.1460,0.6091,0.9079],[0.1380,0.6276,0.8973],[0.1248,0.6459,0.8883],[0.1113,0.6635,0.8763],[0.0952,0.6798,0.8598],[0.0689,0.6948,0.8394],[0.0297,0.7082,0.8163],[0.0036,0.7203,0.7917],[0.0067,0.7312,0.7660],[0.0433,0.7411,0.7394],[0.0964,0.7500,0.7120],[0.1408,0.7584,0.6842],[0.1717,0.7670,0.6554],[0.1938,0.7758,0.6251],[0.2161,0.7843,0.5923],[0.2470,0.7918,0.5567],[0.2906,0.7973,0.5188],[0.3406,0.8008,0.4789],[0.3909,0.8029,0.4354],[0.4456,0.8024,0.3909],[0.5044,0.7993,0.3480],[0.5616,0.7942,0.3045],[0.6174,0.7876,0.2612],[0.6720,0.7793,0.2227],[0.7242,0.7698,0.1910],[0.7738,0.7598,0.1646],[0.8203,0.7498,0.1535],[0.8634,0.7406,0.1596],[0.9035,0.7330,0.1774],[0.9393,0.7288,0.2100],[0.9728,0.7298,0.2394],[0.9956,0.7434,0.2371],[0.9970,0.7659,0.2199],[0.9952,0.7893,0.2028],[0.9892,0.8136,0.1885],[0.9786,0.8386,0.1766],[0.9676,0.8639,0.1643],[0.9610,0.8890,0.1537],[0.9597,0.9135,0.1423],[0.9628,0.9373,0.1265],[0.9691,0.9606,0.1064],[0.9769,0.9839,0.0805]])
      names = self.getAtlasNames(atlasFile)
      types = atlasFile['atlases']['types'][()]
      try:
        showIndexRef = atlasFile['atlases']['presets']['show'][0,0]
        showIndex = atlasFile[showIndexRef][()].squeeze() - 1 # -1 to fix index base
      except:
        showIndex = np.array(range(len(names)))

      for index in range(len(names)): # for each structure     
        
        structureName = os.path.splitext(os.path.splitext(names[index][0])[0])[0]
        structureColor = colormap[int(colors[index])-1]

        if types[index][0] in [3,4]:
          subName = ['rh', 'lh']
          subFolderID = shNode.CreateFolderItem(folderID, structureName)
          self.createFolderDisplayNode(subFolderID, structureColor)
          shNode.SetItemDisplayVisibility(subFolderID, index in showIndex)
          shNode.SetItemExpanded(subFolderID, 0)
          shNode.SetItemAttribute(subFolderID, 'atlas', '1')
        else:
          subName = [structureName]
          subFolderID = folderID

        for sideIndex,sideName in zip(range(len(subName)),subName):
          # get faces and vertices data
          ref = roi[sideIndex][index]
          fv = atlasFile[ref]['fv']
          vertices = fv['vertices'][()].transpose()
          faces = fv['faces'][()].transpose()
          # create polydata
          structurePolyData = self.createPolyData(vertices, faces)          
          # add model node
          modelNode = self.createModel(structurePolyData, sideName, structureColor)
          modelNode.GetDisplayNode().SetVisibility(index in showIndex)
          # add as child to parent
          shNode.SetItemParent(shNode.GetItemChildWithName(shNode.GetSceneItemID(), sideName), subFolderID)
          shNode.SetItemAttribute(shNode.GetItemByDataNode(modelNode), 'atlas', '1')
        
    qt.QApplication.setOverrideCursor(qt.QCursor(qt.Qt.ArrowCursor))
  
    return folderID


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
