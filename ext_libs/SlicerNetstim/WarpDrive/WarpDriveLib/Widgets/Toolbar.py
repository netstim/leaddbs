import qt, vtk, slicer
from qt import QToolBar
import os
import itertools
from slicer.util import VTKObservationMixin
import glob
import re

import WarpDrive
import ImportAtlas
from ..Helpers import LeadDBSCall, WarpDriveUtil
from ..Widgets import ToolWidget

class reducedToolbar(QToolBar, VTKObservationMixin):

  def __init__(self):

    QToolBar.__init__(self)
    VTKObservationMixin.__init__(self)

    self.parameterNode = WarpDrive.WarpDriveLogic().getParameterNode()
    self.addObserver(self.parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateToolbarFromParameterNode)
    
    self.setWindowTitle(qt.QObject().tr("LeadDBS"))
    self.name = 'LeadDBS'
  

    #
    # Modality
    #
    self.addWidget(qt.QLabel('Modality:'))
    self.modalityComboBox = qt.QComboBox()
    self.modalityComboBox.addItem('t1')
    self.modalityComboBox.view().pressed.connect(self.onModalityPressed)
    self.addWidget(self.modalityComboBox)

    #
    # B <-> F slider
    #
    self.addSeparator()
    self.addWidget(qt.QLabel('Template:'))
    templateSlider = qt.QSlider(1)
    templateSlider.singleStep = 10
    templateSlider.minimum = 0
    templateSlider.maximum = 100
    templateSlider.value = 0
    templateSlider.setFixedWidth(120)
    templateSlider.connect('valueChanged(int)', lambda value: slicer.util.setSliceViewerLayers(foregroundOpacity = value / 100.0))
    self.addWidget(templateSlider)


    #
    # Space Separator
    #
    self.addSeparator()
    empty = qt.QWidget()
    empty.setSizePolicy(qt.QSizePolicy.Expanding,qt.QSizePolicy.Preferred)
    self.addWidget(empty)

    #
    # Subject
    #

    self.subjectNameLabel = qt.QLabel('Subject: ')    
    self.addWidget(self.subjectNameLabel)

    #
    # Harden Changes
    #
    self.addSeparator()
    self.hardenChangesCheckBox = qt.QCheckBox("Harden Changes")
    self.hardenChangesCheckBox.checked = True
    self.addWidget(self.hardenChangesCheckBox)

    #
    # Save
    #
    self.nextButton = qt.QPushButton("Exit")
    self.nextButton.setFixedWidth(75)
    self.nextButton.setStyleSheet("background-color: green")
    self.addWidget(self.nextButton)
    self.nextButton.connect("clicked(bool)", self.onNextButton)

    #
    # Update
    #

    self.updateModalities(self.parameterNode.GetParameter("subjectPath"))
    self.initSubject()
    self.updateToolbarFromParameterNode()



  def initSubject(self):
    if os.path.isfile(os.path.join(self.parameterNode.GetParameter("subjectPath"),'WarpDrive','WarpDriveScene.mrml')):
      self.initFromScene()
    else:
      self.initFromRaw()
    
    self.onModalityPressed([],self.modalityComboBox.currentText)

    # load default atlas if there is none in scene
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    folderNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLFolderDisplayNode')
    folderNodes.UnRegister(slicer.mrmlScene)
    for i in range(folderNodes.GetNumberOfItems()):
      folderNode = folderNodes.GetItemAsObject(i)
      if 'atlas' in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(folderNode)):
        return
    ImportAtlas.ImportAtlasLogic().readAtlas(os.path.join(self.parameterNode.GetParameter("MNIAtlasPath"), 'DISTAL Minimal (Ewert 2017)', 'atlas_index.mat'))


  def initFromRaw(self):
    # set up transform
    inputNode = LeadDBSCall.loadSubjectTransform(self.parameterNode.GetParameter("subjectPath"), self.parameterNode.GetParameter("antsApplyTransformsPath"))
    outputNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLGridTransformNode')
    inputNode.SetAndObserveTransformNodeID(outputNode.GetID())
    # parameter node
    self.parameterNode.SetNodeReferenceID("InputNode", inputNode.GetID())
    self.parameterNode.SetNodeReferenceID("OutputGridTransform", outputNode.GetID())

  def initFromScene(self):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    subjectPaths = self.parameterNode.GetParameter("subjectPaths")
    subjectPath = self.parameterNode.GetParameter("subjectPath")
    subjectPathSep = self.parameterNode.GetParameter("separator")
    MNIPath = self.parameterNode.GetParameter("MNIPath")
    MNIAtlasPath = self.parameterNode.GetParameter("MNIAtlasPath")
    antsApplyTransformsPath = self.parameterNode.GetParameter("antsApplyTransformsPath")

    # load previous scene
    try:
      slicer.util.loadScene(os.path.join(self.parameterNode.GetParameter("subjectPath"),'warpdrive','WarpDriveScene.mrml'))
    except:
      pass

    # restore parameters
    self.parameterNode.SetParameter("subjectPaths", subjectPaths)
    self.parameterNode.SetParameter("subjectPath", subjectPath)
    self.parameterNode.SetParameter("separator", subjectPathSep)
    self.parameterNode.SetParameter("MNIPath", MNIPath)
    self.parameterNode.SetParameter("MNIAtlasPath", MNIAtlasPath)
    self.parameterNode.SetParameter("antsApplyTransformsPath", antsApplyTransformsPath)

    # set input node
    gridTransformNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLGridTransformNode')
    gridTransformNodes.UnRegister(slicer.mrmlScene)
    for i in range(gridTransformNodes.GetNumberOfItems()):
      gridTransformNode = gridTransformNodes.GetItemAsObject(i)
      if gridTransformNode.GetID() != self.parameterNode.GetNodeReferenceID("InputNode"):
        slicer.mrmlScene.RemoveNode(gridTransformNode)
    # set image and template
    volumeNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLVolumeNode')
    volumeNodes.UnRegister(slicer.mrmlScene)
    for i in range(volumeNodes.GetNumberOfItems()):
      volumeNode = volumeNodes.GetItemAsObject(i)
      if 'correction' not in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(volumeNode)):
        slicer.mrmlScene.RemoveNode(volumeNode)
    # get atlas name and delete folders
    atlasName = None
    folderNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLFolderDisplayNode')
    folderNodes.UnRegister(slicer.mrmlScene)
    for i in range(folderNodes.GetNumberOfItems()):
      folderNode = folderNodes.GetItemAsObject(i)
      if 'atlas' in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(folderNode)):
        if  atlasName is None and (shNode.GetItemParent(shNode.GetItemByDataNode(folderNode)) == shNode.GetSceneItemID()):
          atlasName = folderNode.GetName()
        slicer.mrmlScene.RemoveNode(folderNode)
    # delete models
    modelNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLModelNode')
    modelNodes.UnRegister(slicer.mrmlScene)
    for i in range(modelNodes.GetNumberOfItems()):
      slicer.mrmlScene.RemoveNode(modelNodes.GetItemAsObject(i))
    # delete fiducials
    fiducialNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLMarkupsFiducialNode')
    fiducialNodes.UnRegister(slicer.mrmlScene)
    for i in range(fiducialNodes.GetNumberOfItems()):
      fiducialNode = fiducialNodes.GetItemAsObject(i)
      if 'correction' not in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(fiducialNode)):
        slicer.mrmlScene.RemoveNode(fiducialNode)

    # load atlas
    if atlasName is not None and os.path.isdir(os.path.join(self.parameterNode.GetParameter("MNIAtlasPath"), atlasName)):
      ImportAtlas.ImportAtlasLogic().readAtlas(os.path.join(self.parameterNode.GetParameter("MNIAtlasPath"), atlasName, 'atlas_index.mat'))

    # init output
    outputNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLGridTransformNode')
    self.parameterNode.SetNodeReferenceID("OutputGridTransform", outputNode.GetID())

    # in case one not set
    WarpDrive.WarpDriveLogic().setDefaultParameters(self.parameterNode)


  def onModalityPressed(self, item, modality=None):
    if modality is None:
      modality = self.modalityComboBox.itemText(item.row())
    # find old nodes and delete
    slicer.mrmlScene.RemoveNode(self.parameterNode.GetNodeReference("ImageNode"))
    slicer.mrmlScene.RemoveNode(self.parameterNode.GetNodeReference("TemplateNode"))
    # initialize new image and init
    imageNode = slicer.util.loadVolume(LeadDBSCall.LeadBIDS(self.parameterNode.GetParameter("subjectPath")).getCoregImages(modality)[0], properties={'show':False})
    imageNode.SetAndObserveTransformNodeID(self.parameterNode.GetNodeReferenceID("InputNode"))    
    # change to t1 in case modality not present
    modality = re.search(r"[^w]+", modality)[0].lower()
    templateFile = glob.glob(os.path.join(self.parameterNode.GetParameter("MNIPath"), modality+ ".nii"))
    if templateFile:
      templateFile = templateFile[0]
    else:
      templateFile = os.path.join(self.parameterNode.GetParameter("MNIPath"), "t1.nii")
    templateNode = slicer.util.loadVolume(templateFile, properties={'show':False})
    templateNode.GetDisplayNode().AutoWindowLevelOff()
    templateNode.GetDisplayNode().SetWindow(100)
    templateNode.GetDisplayNode().SetLevel(70)
    # set view
    slicer.util.setSliceViewerLayers(background=imageNode.GetID(), foreground=templateNode.GetID())
    # set parameter
    self.parameterNode.SetParameter("modality", modality)
    self.parameterNode.SetNodeReferenceID("ImageNode", imageNode.GetID())
    self.parameterNode.SetNodeReferenceID("TemplateNode", templateNode.GetID())


  def updateToolbarFromParameterNode(self, caller=None, event=None):
    # subject text
    subjectN = int(self.parameterNode.GetParameter("subjectN"))
    subjectPaths = self.parameterNode.GetParameter("subjectPaths").split(self.parameterNode.GetParameter("separator"))
    self.subjectNameLabel.text = 'Subject: ' + os.path.split(os.path.abspath(self.parameterNode.GetParameter("subjectPath")))[-1]
    self.nextButton.text = 'Exit' if subjectN == len(subjectPaths)-1 else 'Next'
    # modality
    self.modalityComboBox.setCurrentText(self.parameterNode.GetParameter("modality"))


  def onNextButton(self):
    ToolWidget.AbstractToolWidget.cleanEffects()
    subjectPath = self.parameterNode.GetParameter("subjectPath")

    if WarpDriveUtil.getPointsFromAttribute('source').GetNumberOfPoints(): # corrections made
      if self.hardenChangesCheckBox.checked:
        LeadDBSCall.applyChanges(subjectPath, self.parameterNode.GetNodeReference("InputNode"), self.parameterNode.GetNodeReference("ImageNode")) # save changes
        LeadDBSCall.setTargetFiducialsAsFixed() # change target fiducial as fixed
        LeadDBSCall.saveCurrentScene(subjectPath) # save scene
      else:
        LeadDBSCall.saveCurrentScene(subjectPath)
    else:
      if not LeadDBSCall.queryUserApproveSubject(subjectPath):
        return # user canceled 

    # clean up

    # remove nodes
    slicer.mrmlScene.RemoveNode(self.parameterNode.GetNodeReference("InputNode"))
    slicer.mrmlScene.RemoveNode(self.parameterNode.GetNodeReference("ImageNode"))
    slicer.mrmlScene.RemoveNode(self.parameterNode.GetNodeReference("OutputGridTransform"))
    
    LeadDBSCall.DeleteCorrections()

    # move to next subject

    nextSubjectN = int(self.parameterNode.GetParameter("subjectN"))+1
    subjectPaths = self.parameterNode.GetParameter("subjectPaths").split(self.parameterNode.GetParameter("separator"))
  
    if nextSubjectN < len(subjectPaths):
      self.updateModalities(subjectPaths[nextSubjectN])
      self.parameterNode.SetParameter("subjectN", str(nextSubjectN))
      self.parameterNode.SetParameter("subjectPath", subjectPaths[nextSubjectN])
      self.initSubject()
      self.updateToolbarFromParameterNode()
    else:
      slicer.util.exit()

    # this sets parameter in case new scene loaded
    slicer.util.moduleSelector().selectModule('Data')
    slicer.util.moduleSelector().selectModule('WarpDrive')

  def updateModalities(self, subjectPath):
    currentModality = self.modalityComboBox.currentText
    subjectModalities = self.getAvailableModalities(subjectPath)
    if currentModality not in subjectModalities:
      self.parameterNode.SetParameter("modality", subjectModalities[0])
    self.modalityComboBox.clear()
    self.modalityComboBox.addItems(subjectModalities)

  def getAvailableModalities(self, directory):
    modalities = []
    fileNames = LeadDBSCall.LeadBIDS(directory).getCoregImages()
    for fileName in fileNames:
      fileName = os.path.split(fileName)[-1] # remove directory
      fileName = os.path.splitext(fileName)[0] # remove extension
      modality = re.search(r"(?<=ses-preop_)\w+", fileName)[0]
      modalities.append(modality)
    return modalities