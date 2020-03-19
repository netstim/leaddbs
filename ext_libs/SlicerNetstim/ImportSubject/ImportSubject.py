import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging


import sys
from subprocess import call
sys.path.append(os.path.join(os.path.dirname(__file__),'..','TransformsUtil'))
import TransformsUtil

#
# ImportSubject
#

class ImportSubject(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "ImportSubject" # TODO make this more human readable by adding spaces
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
# ImportSubjectWidget
#

class ImportSubjectWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    #
    # Subject Directory Area
    #

    self.subjectDirectoryButton = ctk.ctkDirectoryButton()
    self.subjectDirectoryButton.text = "Select Lead-DBS directory"
    self.layout.addWidget(self.subjectDirectoryButton)


    #
    # Images Area
    #
    imagesCollapsibleButton = ctk.ctkCollapsibleButton()
    imagesCollapsibleButton.text = "Images"
    self.layout.addWidget(imagesCollapsibleButton)

    # Layout within the dummy collapsible button
    imagesFormLayout = qt.QFormLayout(imagesCollapsibleButton)

    #
    # check box select images
    #

    self.imageCheckboxes = []
    parameterNode = ImportSubjectLogic().getParameterNode()
    posibleImages = parameterNode.GetParameter("posibleImages").split(' ')

    for posibleImage in posibleImages:
      imageCheckbox = qt.QCheckBox(posibleImage)
      imageCheckbox.checked = False
      imageCheckbox.visible = False
      imagesFormLayout.addRow(imageCheckbox)
      self.imageCheckboxes.append(imageCheckbox)

    #
    # Transforms Area
    #
    transformsCollapsibleButton = ctk.ctkCollapsibleButton()
    transformsCollapsibleButton.text = "Transforms"
    self.layout.addWidget(transformsCollapsibleButton)

    # Layout within the dummy collapsible button
    transformsFormLayout = qt.QFormLayout(transformsCollapsibleButton)

    # converts transforms
    self.updateTransformButton = qt.QPushButton('Update Transform')
    self.updateTransformButton.visible = False
    transformsFormLayout.addRow(self.updateTransformButton)

    #
    # check box select transforms
    #

    self.transformCheckboxes = []
    posibleTransforms = parameterNode.GetParameter("posibleTransforms").split(' ')
    for posibleTransform in posibleTransforms:
      transformCheckbox = qt.QCheckBox(posibleTransform)
      transformCheckbox.checked = False
      transformCheckbox.visible = False
      transformsFormLayout.addRow(transformCheckbox)
      self.transformCheckboxes.append(transformCheckbox)     


    #
    # Electrode reconstruction
    #
    reconstructionCollapsibleButton = ctk.ctkCollapsibleButton()
    reconstructionCollapsibleButton.text = "Electrode Reconstruction"
    self.layout.addWidget(reconstructionCollapsibleButton)

    # Layout within the dummy collapsible button
    reconstructionFormLayout = qt.QFormLayout(reconstructionCollapsibleButton)

    self.recoMNICheckbox = qt.QCheckBox('MNI')
    self.recoMNICheckbox.checked = False
    self.recoMNICheckbox.visible = False

    #
    # Import Button
    #
    self.importButton = qt.QPushButton("Import")
    self.importButton.toolTip = "Import selected options."
    self.importButton.enabled = True
    self.layout.addWidget(self.importButton)

    # connections
    self.importButton.connect('clicked(bool)', self.onImportButton)
    self.updateTransformButton.connect('clicked(bool)', self.onUpdateTransformButton)
    self.subjectDirectoryButton.directoryChanged.connect(self.onSubjectDirectoryChanged)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSubjectDirectoryChanged('.')

  def cleanup(self):
    pass

  def onSubjectDirectoryChanged(self, directory):
    logic = ImportSubjectLogic()
    availableFiles = logic.getAvailableFiles(directory, "Images") + logic.getAvailableFiles(directory, "Transforms")
    for checkbox in self.imageCheckboxes + self.transformCheckboxes:
      checkbox.checked = False
      checkbox.visible = checkbox.text in availableFiles # enable checkbox if image is available
    # check for old transform version
    self.updateTransformButton.visible = logic.ish5Transform(directory)
    # reconstruction
    self.recoMNICheckbox.visible = os.path.isfile(os.path.join(directory, 'ea_reconstruction.mat'))
    # change subject buton text to subject directory name
    subjectName = os.path.basename(directory) if any([ch.visible for ch in self.imageCheckboxes + self.transformCheckboxes]) else "Select Lead-DBS directory"
    self.subjectDirectoryButton.text = subjectName


  def onImportButton(self):
    logic = ImportSubjectLogic()
    # import images and transforms
    for checkbox in self.imageCheckboxes + self.transformCheckboxes:
      if checkbox.visible and checkbox.checked:
        logic.importFile(self.subjectDirectoryButton.directory, checkbox.text)
    # import electrode reconstruction
    if self.recoMNICheckbox.checked:
      logic.importReconstruction(self.subjectDirectoryButton.directory)

  def onUpdateTransformButton(self):
    logic = ImportSubjectLogic()
    directory = self.subjectDirectoryButton.directory
    logic.updateTranform(directory)
    # update
    self.onSubjectDirectoryChanged(directory)

#
# ImportSubjectLogic
#

class ImportSubjectLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def createParameterNode(self):
    node = ScriptedLoadableModuleLogic.createParameterNode(self)
    node.SetParameter("posibleImages", "anat_t1.nii anat_t2.nii anat_pd.nii postop_ct.nii glanat.nii")
    node.SetParameter("posibleTransforms", "glanat0GenericAffine_backup.mat glanatComposite.nii.gz glanatInverseComposite.nii.gz")
    return node 

  def ish5Transform(self, directory):
    return os.path.isfile(os.path.join(directory, 'glanatComposite.h5'))

  def getAvailableFiles(self, directory, fileType):
    """
    Check for files available in directory. fileType is "Images" os "Transforms"
    """
    parameterNode = self.getParameterNode()
    posibleFiles = parameterNode.GetParameter("posible" + fileType).split(' ')
    availableFiles = [posibleFile for posibleFile in posibleFiles if os.path.isfile(os.path.join(directory, posibleFile))]
    return availableFiles

  def importFile(self, directory, fileName):
    """
    Import image or transform from directory and return node
    """
    if fileName in self.getAvailableFiles(directory, "Images"):
      isImage = True
    elif fileName in self.getAvailableFiles(directory, "Transforms"):
      isImage = False
    else:
      return None

    filePath = os.path.join(directory, fileName)

    if isImage:
      node = slicer.util.loadVolume(filePath, properties={'show':False})
    else:
      node = slicer.util.loadTransform(filePath)
    
    subjectName = os.path.basename(directory)
    fileNameNoExt = os.path.splitext(os.path.splitext(fileName)[0])[0]    
    node.SetName(subjectName + '_' + fileNameNoExt)

    return node

  def importReconstruction(self, directory):
    pass


  def saveAffineComponent(self,transformNode):
    # split transforms and get new nodes
    newNodeNames = TransformsUtil.TransformsUtilLogic().splitAndGetNodeNames(transformNode)
    # save affine
    slicer.util.saveNode(transformNode, os.path.join(os.path.dirname(transformNode.GetStorageNode().GetFileName()),'glanat0GenericAffine_backup.mat'))
    # join again and delete created node
    transformNode.HardenTransform()
    slicer.mrmlScene.RemoveNode(slicer.util.getNode(newNodeNames[0]))
    

  def updateTranform(self, directory, antsApplyTransformsPath=None):
    transformNode = slicer.util.loadTransform(os.path.join(directory,'glanatComposite.h5'))
    self.saveAffineComponent(transformNode)
    slicer.mrmlScene.RemoveNode(transformNode)
    # flatten
    if not antsApplyTransformsPath:
      w = qt.QWidget()
      fd = qt.QFileDialog(w,'AntsApplyTransformsPath')
      if fd.exec():
        antsApplyTransformsPath = fd.selectedFiles()[0]
      else:
        return False
    
    for transform,reference in zip(['glanatComposite','glanatInverseComposite'],['glanat','anat_t1']):
      transformFullPath = os.path.join(directory,transform + '.h5') # in case inverse doesnt exist
      if os.path.isfile(transformFullPath):
        command = antsApplyTransformsPath + " -r " + os.path.join(directory,reference + '.nii') + " -t " + transformFullPath + " -o [" + os.path.join(directory,transform + '.nii.gz') + ",1] -v 1"
        commandOut = call(command, env=slicer.util.startupEnvironment(), shell=True) # run antsApplyTransforms
        os.remove(transformFullPath)
    return True


class ImportSubjectTest(ScriptedLoadableModuleTest):
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
    self.test_ImportSubject1()

  def test_ImportSubject1(self):
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
    import SampleData
    SampleData.downloadFromURL(
      nodeNames='FA',
      fileNames='FA.nrrd',
      uris='http://slicer.kitware.com/midas3/download?items=5767',
      checksums='SHA256:12d17fba4f2e1f1a843f0757366f28c3f3e1a8bb38836f0de2a32bb1cd476560')
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = ImportSubjectLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
