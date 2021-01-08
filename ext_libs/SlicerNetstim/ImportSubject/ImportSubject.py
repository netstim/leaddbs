import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import glob

import sys
from subprocess import call
import numpy as np
import SimpleITK as sitk
import sitkUtils

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
    # select images list
    #

    self.imagesList = qt.QListWidget()
    self.imagesList.setSelectionMode(qt.QAbstractItemView.ExtendedSelection)
    imagesFormLayout.addRow(self.imagesList)

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

    self.transformsList = qt.QListWidget()
    self.transformsList.setSelectionMode(qt.QAbstractItemView.ExtendedSelection)
    imagesFormLayout.addRow(self.transformsList)


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
    self.imagesList.clear()
    self.imagesList.addItems(logic.getAvailableModalities(directory))
    self.transformsList.clear()
    self.transformsList.addItems(logic.getAvailableTransforms(directory))
    # check for old transform version
    self.updateTransformButton.visible = logic.ish5Transform(directory)
    # change subject buton text to subject directory name
    subjectName = os.path.basename(directory) if self.imagesList.count else "Select Lead-DBS directory"
    self.subjectDirectoryButton.text = subjectName


  def onImportButton(self):
    logic = ImportSubjectLogic()
    for i in range(self.imagesList.count):
      if self.imagesList.item(i).isSelected():
        logic.importImage(self.subjectDirectoryButton.directory, self.imagesList.item(i).text())
    for i in range(self.transformsList.count):
      if self.transformsList.item(i).isSelected():
        logic.importTransform(self.subjectDirectoryButton.directory, self.transformsList.item(i).text())


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

  def ish5Transform(self, directory):
    return os.path.isfile(os.path.join(directory, 'glanatComposite.h5'))

  def getAvailableModalities(self, directory):
    modalities = []
    listing = glob.glob(os.path.join(directory,'anat_*.nii'))
    for fileName in listing:
      fileName = os.path.split(fileName)[-1] # remove directory
      fileName = os.path.splitext(fileName)[0] # remove extension
      modality = fileName[5:] # remove 'anat_'
      modalities.append(modality)
    return modalities

  def getAvailableTransforms(self, directory):
    posibleTransforms = ["glanat0GenericAffine_backup.mat", "glanatComposite.nii.gz", "glanatInverseComposite.nii.gz"]
    availableTransforms = [pt for pt in posibleTransforms if os.path.isfile(os.path.join(directory,pt))]
    return availableTransforms


  def createNodeName(self, directory, fileName):
    subjectName = os.path.split(os.path.abspath(directory))[-1]
    fileNameNoExt = os.path.splitext(os.path.splitext(fileName)[0])[0]    
    return subjectName + '_' + fileNameNoExt

  def importImage(self, directory, fileName):
    if os.path.splitext(fileName)[-1] != '.nii':
      fileName = 'anat_' + fileName + '.nii'
    filePath = os.path.join(directory, fileName)
    node = slicer.util.loadVolume(filePath, properties={'show':False})
    node.SetName(self.createNodeName(directory,fileName))
    return node

  def importTransform(self, directory, fileName):
    filePath = os.path.join(directory, fileName)
    if os.path.isfile(filePath):
      node = slicer.util.loadTransform(filePath)
      node.SetName(self.createNodeName(directory,fileName))
      return node
    else:
      return None

  def importReconstruction(self, directory):
    pass


  def updateTranform(self, directory, antsApplyTransformsPath=None):

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

  def runBinaryThresholdImageFilter(self, inputNode, outputNode):
    # run Simple ITK threshold Filter
    inputImage = sitkUtils.PullVolumeFromSlicer(inputNode)
    myFilter = sitk.BinaryThresholdImageFilter()
    myFilter.SetLowerThreshold(0.5)
    outputImage = myFilter.Execute(inputImage)
    sitkUtils.PushVolumeToSlicer(outputImage, outputNode)
    # run Simple ITK fill holes Filter
    inputImage = sitkUtils.PullVolumeFromSlicer(outputNode)
    myFilter = sitk.BinaryFillholeImageFilter()
    outputImage = myFilter.Execute(inputImage)
    sitkUtils.PushVolumeToSlicer(outputImage, outputNode)

  def importSegmentations(self, directory):
    # look for nifty files in the segmentations subdirectory
    # load the files and binarize if necesary
    # add the binary segments to a segmentation node
    # segmentation to model node

    listing = glob.glob(os.path.join(directory,'segmentations','*.nii*'))

    if not listing: 
      return

    # init segmentation node
    segmentationNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentationNode')
    segmentColor = [0] * 4
    currentSegment = 0

    for filename, i in zip(listing,range(len(listing))):

      segmentName = os.path.split(filename)[-1].split('.')[0] # name
      volumeNode = slicer.util.loadVolume(filename, properties={'name': 'tmp'}) # load volume node
      labelMapNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLabelMapVolumeNode') # init label map node

      # volume to labelmap
      if volumeNode.GetImageData().GetScalarType() in [vtk.VTK_FLOAT, vtk.VTK_DOUBLE]:
        self.runBinaryThresholdImageFilter(volumeNode, labelMapNode)
      else:
        slicer.modules.volumes.logic().CreateLabelVolumeFromVolume(slicer.mrmlScene, labelMapNode, volumeNode)
      
      # add to segmentation
      slicer.modules.segmentations.logic().ImportLabelmapToSegmentationNode(labelMapNode, segmentationNode)
      for s in range(currentSegment, segmentationNode.GetSegmentation().GetNumberOfSegments()):
        slicer.util.getNode('Labels').GetColor(s+1, segmentColor) # color
        addedSegment = segmentationNode.GetSegmentation().GetNthSegment(s)
        addedSegment.SetName(slicer.mrmlScene.GenerateUniqueName(segmentName))
        addedSegment.SetColor(segmentColor[:-1])
      
      currentSegment = segmentationNode.GetSegmentation().GetNumberOfSegments()
      # remove nodes
      slicer.mrmlScene.RemoveNode(volumeNode)
      slicer.mrmlScene.RemoveNode(labelMapNode)

    # add data attributes
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    shNode.SetItemAttribute(shNode.GetItemByDataNode(segmentationNode), 'Segment', '1')
    # segmentation children
    IDList = vtk.vtkIdList()
    shNode.GetItemChildren(shNode.GetItemByDataNode(segmentationNode), IDList)
    for i in range(IDList.GetNumberOfIds()):
      shNode.SetItemAttribute(IDList.GetId(i), 'Segment', '1')
      
    return segmentationNode



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
