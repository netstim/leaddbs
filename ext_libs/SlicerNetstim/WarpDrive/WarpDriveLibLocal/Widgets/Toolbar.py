import qt, vtk, slicer
from qt import QToolBar
import os
import json
from slicer.util import VTKObservationMixin
import glob
import re

import WarpDrive
import ImportAtlas
from ..Helpers import LeadDBSCall
from ..Widgets import ToolWidget

class reducedToolbar(QToolBar, VTKObservationMixin):

  def __init__(self):

    QToolBar.__init__(self)
    VTKObservationMixin.__init__(self)

    self.parameterNode = WarpDrive.WarpDriveLogic().getParameterNode()
    self.parameterNode.SetParameter("modality", 'ax_T2w')
    self.addObserver(self.parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateToolbarFromParameterNode)
    
    self.setWindowTitle(qt.QObject().tr("LeadDBS"))
    self.name = 'LeadDBS'
  
    #
    # Space Separator
    #
    empty = qt.QWidget()
    empty.setSizePolicy(qt.QSizePolicy.Expanding,qt.QSizePolicy.Preferred)
    self.addWidget(empty)

    #
    # Modality
    #
    self.addWidget(qt.QLabel('Modality:'))
    self.modalityComboBox = qt.QComboBox()
    self.modalityComboBox.addItem('ax_T2w')
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
    # Subject
    #
    self.addSeparator()
    self.subjectNameLabel = qt.QLabel('Subject: ')    
    self.addWidget(self.subjectNameLabel)

    #
    # Harden Changes
    #
    self.addSeparator()
    self.hardenChangesCheckBox = qt.QCheckBox("Harden Changes")
    self.hardenChangesCheckBox.checked = True
    self.hardenChangesCheckBox.toolTip = 'When checked, the changes will be written to the subject transform files. If not, the transforms will not be modified and source and target points will be saved for the next time warpdrive is opened.'
    self.addWidget(self.hardenChangesCheckBox)

    #
    # Save
    #
    self.nextButton = qt.QPushButton("Next")
    self.nextButton.setFixedWidth(75)
    self.nextButton.setStyleSheet("background-color: green")
    self.addWidget(self.nextButton)
    self.nextButton.connect("clicked(bool)", self.nextSubject)

    #
    # Update
    #

    self.nextSubject()


  def nextSubject(self):
    print("Going to next subject")
    if self.parameterNode.GetParameter("CurrentSubject"):
      keep_same_subject = self.finalizeCurrentSubject()
      if keep_same_subject:
        return
    leadSubjects  = json.loads(self.parameterNode.GetParameter("LeadSubjects"))
    if not leadSubjects:
      print("No more subjects. Terminate Slicer")
      slicer.util.exit()
    if isinstance(leadSubjects, dict):
      leadSubjects = [leadSubjects]
    self.parameterNode.SetParameter("CurrentSubject", json.dumps(leadSubjects.pop(0)))
    self.parameterNode.SetParameter("LeadSubjects", json.dumps(leadSubjects))
    self.initializeCurrentSubject()

  def initializeCurrentSubject(self):
    currentSubject = json.loads(self.parameterNode.GetParameter("CurrentSubject"))
    print("Initialize subject: %s" % currentSubject["id"])

    inputNode = slicer.util.loadTransform(currentSubject["forward_transform"])
    outputNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLGridTransformNode')
    inputNode.SetAndObserveTransformNodeID(outputNode.GetID())

    if os.path.isfile(os.path.join(currentSubject["warpdrive_path"],'target.json')):
      print("Loading previious session")
      targetFiducial = slicer.util.loadMarkups(os.path.join(currentSubject["warpdrive_path"],'target.json'))
      sourceFiducial = slicer.util.loadMarkups(os.path.join(currentSubject["warpdrive_path"],'source.json'))
    else:
      targetFiducial = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
      targetFiducial.GetDisplayNode().SetGlyphTypeFromString('Sphere3D')
      targetFiducial.GetDisplayNode().SetGlyphScale(1)
      sourceFiducial = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
      sourceFiducial.GetDisplayNode().SetGlyphTypeFromString('Sphere3D')
      sourceFiducial.GetDisplayNode().SetGlyphScale(1)

    # parameter node
    self.parameterNode.SetNodeReferenceID("InputNode", inputNode.GetID())
    self.parameterNode.SetNodeReferenceID("OutputGridTransform", outputNode.GetID())
    self.parameterNode.SetNodeReferenceID("SourceFiducial", sourceFiducial.GetID())
    self.parameterNode.SetNodeReferenceID("TargetFiducial", targetFiducial.GetID())

    self.updateModalities()
    self.onModalityPressed([], self.parameterNode.GetParameter("modality"))

    self.setUpAtlases()
    print("Finish loading subject %s" % currentSubject["id"])



  def finalizeCurrentSubject(self):
    ToolWidget.AbstractToolWidget.cleanEffects()
    currentSubject = json.loads(self.parameterNode.GetParameter("CurrentSubject"))
    sourceFiducial = self.parameterNode.GetNodeReference("SourceFiducial")
    targetFiducial = self.parameterNode.GetNodeReference("TargetFiducial")

    if sourceFiducial.GetNumberOfControlPoints(): # corrections made
      if self.hardenChangesCheckBox.checked:
        LeadDBSCall.applyChanges(self.parameterNode.GetNodeReference("InputNode"), self.parameterNode.GetNodeReference("ImageNode"), currentSubject["forward_transform"], currentSubject["inverse_transform"])
        sourceFiducial.Copy(targetFiducial) # set all as fixed points
      LeadDBSCall.saveSourceTarget(currentSubject["warpdrive_path"], sourceFiducial, targetFiducial)
      LeadDBSCall.saveSceneInfo(currentSubject["warpdrive_path"])
    else:
      if LeadDBSCall.queryUserApproveSubject():
        LeadDBSCall.saveApprovedData(currentSubject["normlog_file"])
      else:
        return 1 # user canceled 

    # clean up
    slicer.mrmlScene.RemoveNode(sourceFiducial)
    slicer.mrmlScene.RemoveNode(targetFiducial)
    slicer.mrmlScene.RemoveNode(self.parameterNode.GetNodeReference("InputNode"))
    slicer.mrmlScene.RemoveNode(self.parameterNode.GetNodeReference("ImageNode"))
    slicer.mrmlScene.RemoveNode(self.parameterNode.GetNodeReference("OutputGridTransform"))
    
  def setUpAtlases(self):
    print("Set up atlases")
    currentSubject = json.loads(self.parameterNode.GetParameter("CurrentSubject"))
    jsonFileName = os.path.join(currentSubject["warpdrive_path"],'info.json')
    if os.path.isfile(jsonFileName):
      with open(jsonFileName, 'r') as jsonFile:
        info = json.load(jsonFile)
    else:
      info = {"atlasNames": []}
    atlasNames = info["atlasNames"] if info["atlasNames"] != [] else ['DISTAL Minimal (Ewert 2017)']
    # load atlas if not already in scene
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    folderNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLFolderDisplayNode')
    folderNodes.UnRegister(slicer.mrmlScene)
    for i in range(folderNodes.GetNumberOfItems()):
      folderNode = folderNodes.GetItemAsObject(i)
      if ('atlas' in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(folderNode))) and (folderNode.GetName() in atlasNames):
        atlasNames.pop(atlasNames.index(folderNode.GetName()))
    for name in atlasNames:
      print("Loading atlas %s" % name)
      try:
        ImportAtlas.ImportAtlasLogic().readAtlas(os.path.join(ImportAtlas.ImportAtlasLogic().getAtlasesPath(), name, 'atlas_index.mat'))
      except:
        print("Could not load atlas %s" % name)

  def onModalityPressed(self, item, modality=None):
    if modality is None:
      modality = self.modalityComboBox.itemText(item.row())
    print("Loading %s modality" % modality)
    # find old nodes and delete
    slicer.mrmlScene.RemoveNode(self.parameterNode.GetNodeReference("ImageNode"))
    slicer.mrmlScene.RemoveNode(self.parameterNode.GetNodeReference("TemplateNode"))
    # initialize new image and init
    currentSubject = json.loads(self.parameterNode.GetParameter("CurrentSubject"))
    imageNode = slicer.util.loadVolume(currentSubject["anat_files"][modality], properties={'show':False})
    imageNode.SetAndObserveTransformNodeID(self.parameterNode.GetNodeReferenceID("InputNode"))    
    # change to t1 in case modality not present
    mni_modality = re.findall(r'(?<=T)\d', modality) + ['1']
    mni_modality = mni_modality[0]
    templateFile = glob.glob(os.path.join(self.parameterNode.GetParameter("MNIPath"), "t" + mni_modality + ".nii"))
    templateFile = templateFile[0] if templateFile else os.path.join(self.parameterNode.GetParameter("MNIPath"), "t1.nii")
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
    self.subjectNameLabel.text = 'Subject: ' + json.loads(self.parameterNode.GetParameter("CurrentSubject"))["id"]
    self.subjectNameLabel.toolTip = os.path.dirname(json.loads(self.parameterNode.GetParameter("CurrentSubject"))["warpdrive_path"])
    self.nextButton.text = 'Next' if len(json.loads(self.parameterNode.GetParameter("LeadSubjects"))) else 'Exit'
    self.modalityComboBox.setCurrentText(self.parameterNode.GetParameter("modality"))      


  def updateModalities(self):
    print("Update modalities")
    currentSubject = json.loads(self.parameterNode.GetParameter("CurrentSubject"))
    currentModality = self.modalityComboBox.currentText
    subjectModalities = list(currentSubject["anat_files"].keys())
    self.modalityComboBox.clear()
    self.modalityComboBox.addItems(subjectModalities)
    if currentModality not in subjectModalities:
      self.parameterNode.SetParameter("modality", subjectModalities[0])