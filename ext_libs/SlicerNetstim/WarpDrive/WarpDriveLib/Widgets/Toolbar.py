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
    # Inverse
    #
    self.inverseAction = qt.QAction(self)
    self.inverseAction.setIcon(qt.QIcon(":/Icons/Small/SlicerCheckForUpdates.png"))
    self.inverseAction.setText('Inverse')
    self.inverseAction.setToolTip('Switch between native and template views. Default is with template fixed, and deform the native image. In inverse mode, the native image is fixed, and the template is deformed.')
    self.inverseAction.setCheckable(True)
    self.inverseAction.setChecked(int(self.parameterNode.GetParameter("SegmentMode")))
    self.inverseAction.connect("triggered(bool)", self.onInverseTriggered)

    inverseToolButton = qt.QToolButton()
    inverseToolButton.setDefaultAction(self.inverseAction)
    inverseToolButton.setToolButtonStyle(qt.Qt.ToolButtonTextBesideIcon)
    inverseToolButton.setVisible(not int(self.parameterNode.GetParameter("SegmentMode")))

    #
    # Subject
    #
    self.subjectAction = qt.QAction(self)
    self.subjectAction.setIcon(qt.QIcon(":/Icons/Patient.png"))
    self.subjectAction.setEnabled(False)

    subjectButton = qt.QToolButton()
    subjectButton.setDefaultAction(self.subjectAction)
    subjectButton.setStyleSheet("color: black")
    subjectButton.setToolButtonStyle(qt.Qt.ToolButtonTextBesideIcon)

    #
    # Modality
    #
    a = qt.QAction(self)
    a.setText('ax_T2w')
    a.setCheckable(True)
    a.setChecked(True)

    self.modalitiesGroup = qt.QActionGroup(self)
    self.modalitiesGroup.setExclusive(True)
    self.modalitiesGroup.addAction(a)
    self.modalitiesGroup.connect('triggered(QAction*)', self.modalityChanged)

    self.modalitiesMenu = qt.QMenu(self)
    self.modalitiesMenu.addActions(self.modalitiesGroup.actions())

    modalityAction = qt.QAction(self)
    modalityAction.setIcon(qt.QIcon(":/Icons/Small/SlicerVisible.png"))
    modalityAction.setText('Modality')

    modalityButton = qt.QToolButton()
    modalityButton.setDefaultAction(modalityAction)
    modalityButton.setMenu(self.modalitiesMenu)
    modalityButton.setToolButtonStyle(qt.Qt.ToolButtonTextBesideIcon)
    modalityButton.setPopupMode(qt.QToolButton.InstantPopup)

    #
    # B <-> F slider
    #
    templateSlider = qt.QSlider(1)
    templateSlider.singleStep = 10
    templateSlider.minimum = 0
    templateSlider.maximum = 100
    templateSlider.value = 0

    action = qt.QWidgetAction(self)
    action.setDefaultWidget(templateSlider)

    menu = qt.QMenu(self)
    menu.addAction(action)

    templateAction = qt.QAction(self)
    templateAction.setIcon(qt.QIcon(":/Icons/Small/SlicerVisible.png"))
    templateAction.setText('Template')
    templateAction.setToolTip('Toggle between subject and template image. Long-press to set opacity value.')
    templateAction.setCheckable(True)
    templateAction.connect("triggered(bool)", lambda b: templateSlider.setValue(int(b)*100))
    templateSlider.connect('valueChanged(int)', lambda value: [templateAction.setChecked(value>0), slicer.util.setSliceViewerLayers(foregroundOpacity = value / 100.0)])

    templateButton = qt.QToolButton()
    templateButton.setDefaultAction(templateAction)
    templateButton.setMenu(menu)
    templateButton.setToolButtonStyle(qt.Qt.ToolButtonTextBesideIcon)
    templateButton.setPopupMode(qt.QToolButton.DelayedPopup)

    #
    # Save
    #
    self.hardenChangesAction = qt.QAction(self)
    self.hardenChangesAction.setText('Harden Changes')
    self.hardenChangesAction.setToolTip('When checked, the transformation files will be overwritten. If not, only the corrections are saved and will be loaded next time WarpDrive is used.')
    self.hardenChangesAction.setCheckable(True)
    self.hardenChangesAction.setChecked(True)

    self.saveApprovedAction = qt.QAction(self)
    self.saveApprovedAction.setText('Approve Normalization')
    self.saveApprovedAction.setToolTip('Lead-DBS will only open non-approved normalization files (or when re-touching).')
    self.saveApprovedAction.setCheckable(True)
    self.saveApprovedAction.connect("triggered(bool)", self.setSubjecApproved)

    menu = qt.QMenu(self)
    menu.addAction(self.hardenChangesAction)
    menu.addSeparator()
    menu.addAction(self.saveApprovedAction)

    self.nextAction = qt.QAction(self)
    self.nextAction.setIcon(qt.QIcon(":/Icons/Small/SlicerSave.png"))
    self.nextAction.setText('Save and Exit')
    self.nextAction.connect("triggered(bool)", self.nextSubject)

    nextButton = qt.QToolButton()
    nextButton.setFixedWidth(120)
    nextButton.setMenu(menu)
    nextButton.setDefaultAction(self.nextAction)
    nextButton.setPopupMode(qt.QToolButton.MenuButtonPopup)
    nextButton.setToolButtonStyle(qt.Qt.ToolButtonTextBesideIcon)

    #
    # Set up toolbar
    #
    empty = qt.QWidget()
    empty.setSizePolicy(qt.QSizePolicy.Expanding,qt.QSizePolicy.Preferred)
    self.addWidget(empty)
    self.addSeparator()
    self.addWidget(inverseToolButton)
    self.addSeparator()
    self.addWidget(subjectButton)
    self.addSeparator()
    self.addWidget(modalityButton)
    self.addSeparator()
    self.addWidget(templateButton)
    self.addSeparator()
    self.addWidget(nextButton)

    #
    # Update
    #
    self.updateToolbarFromParameterNode()

  def updateToolbarFromParameterNode(self, caller=None, event=None):
    currentSubjectJson = self.parameterNode.GetParameter("CurrentSubject")
    if currentSubjectJson:
      currentSubject = json.loads(currentSubjectJson)
      self.subjectAction.text = "%s (%s/%s)" % (currentSubject["id"], self.parameterNode.GetParameter("CurrentSubjectNumber"), self.parameterNode.GetParameter("TotalNumberOfSubjects"))
      self.subjectAction.toolTip = os.path.dirname(currentSubject["warpdrive_path"])
      self.nextAction.text = 'Save and Next' if len(json.loads(self.parameterNode.GetParameter("LeadSubjects"))) else 'Save and Exit'
      next(filter(lambda a: a.text==self.parameterNode.GetParameter("modality"), self.modalitiesGroup.actions())).setChecked(True)
    elif self.parameterNode.GetParameter("LeadSubjects"):
      self.nextSubject()

  def nextSubject(self):
    print("Going to next subject")
    wasModified = self.parameterNode.StartModify()
    # Get remaining subjects
    leadSubjects  = json.loads(self.parameterNode.GetParameter("LeadSubjects"))
    if isinstance(leadSubjects, dict):
      leadSubjects = [leadSubjects]
    numberOfSubjects = len(leadSubjects)
    oneOrMoreRemainingSubjects = numberOfSubjects >= 1
    if not self.parameterNode.GetParameter("TotalNumberOfSubjects"):
      self.parameterNode.SetParameter("TotalNumberOfSubjects", str(numberOfSubjects))
    # Save current subject
    slicerWillExitAfterSave = False
    if self.parameterNode.GetParameter("CurrentSubject"):
      slicerWillExitAfterSave = self.saveCurrentSubject(saveInExternalInstance = oneOrMoreRemainingSubjects)
    # Load next subject
    if oneOrMoreRemainingSubjects:
      self.cleanUpNodes()
      self.parameterNode.SetParameter("CurrentSubject", json.dumps(leadSubjects.pop(0)))
      self.parameterNode.SetParameter("LeadSubjects", json.dumps(leadSubjects))
      self.parameterNode.SetParameter("CurrentSubjectNumber", str(int(self.parameterNode.GetParameter("TotalNumberOfSubjects")) - len(leadSubjects)))
      self.initializeCurrentSubject()
    elif (not slicerWillExitAfterSave):
      slicer.util.exit(0)
    self.parameterNode.EndModify(wasModified)

  def cleanUpNodes(self):
    for param in ["SourceFiducial", "TargetFiducial", "InputNode", "ImageNode", "OutputGridTransform"]:
      if self.parameterNode.GetNodeReference(param):
        slicer.mrmlScene.RemoveNode(self.parameterNode.GetNodeReference(param))

  def initializeCurrentSubject(self):
    currentSubject = json.loads(self.parameterNode.GetParameter("CurrentSubject"))
    print("Initialize subject: %s" % currentSubject["id"])
    jsonFileName = os.path.join(currentSubject["warpdrive_path"],'info.json')
    if os.path.isfile(jsonFileName) and not int(self.parameterNode.GetParameter("SegmentMode")):
      with open(jsonFileName, 'r') as jsonFile:
        subjectInfo = json.load(jsonFile)
    else:
      subjectInfo = None

    useInverse = self.inverseAction.checked
    inputNode = slicer.util.loadTransform(currentSubject["inverse_transform" if useInverse else "forward_transform"])
    outputNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLGridTransformNode')
    inputNode.SetAndObserveTransformNodeID(outputNode.GetID())

    if os.path.isfile(os.path.join(currentSubject["warpdrive_path"],'target.json')) and not int(self.parameterNode.GetParameter("SegmentMode")):
      print("Loading previous session")
      targetFiducial = slicer.util.loadMarkups(os.path.join(currentSubject["warpdrive_path"],'target.json'))
      sourceFiducial = slicer.util.loadMarkups(os.path.join(currentSubject["warpdrive_path"],'source.json'))
      if subjectInfo and "inverseApplied" in subjectInfo.keys() and subjectInfo["inverseApplied"]:
        if not useInverse:
          self.invertSourceTargetNodes(sourceFiducial, targetFiducial, inputNode)
      elif useInverse:
        self.invertSourceTargetNodes(sourceFiducial, targetFiducial, inputNode)
    else:
      targetFiducial = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
      targetFiducial.GetDisplayNode().SetGlyphTypeFromString('Sphere3D')
      targetFiducial.GetDisplayNode().SetGlyphScale(1)
      sourceFiducial = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
      sourceFiducial.GetDisplayNode().SetGlyphTypeFromString('Sphere3D')
      sourceFiducial.GetDisplayNode().SetGlyphScale(1)

    # parameter node
    wasModified = self.parameterNode.StartModify()
    self.parameterNode.SetNodeReferenceID("InputNode", inputNode.GetID())

    self.saveApprovedAction.setChecked(LeadDBSCall.getApprovedData(currentSubject["normlog_file"]))
    self.updateModalitiesToolButton()
    self.updateModalitiesImages(self.parameterNode.GetParameter("modality"))
    self.setUpAtlases(subjectInfo)

    self.parameterNode.SetNodeReferenceID("OutputGridTransform", outputNode.GetID())
    self.parameterNode.SetNodeReferenceID("SourceFiducial", sourceFiducial.GetID())
    self.parameterNode.SetNodeReferenceID("TargetFiducial", targetFiducial.GetID())
    self.parameterNode.EndModify(wasModified)
    print("Finish loading subject %s" % currentSubject["id"])

  def saveCurrentSubject(self, saveInExternalInstance):
    slicerWillExit = False
    slicer.util.setSliceViewerLayers(background=None, foreground=None)
    ToolWidget.AbstractToolWidget.cleanEffects()
    currentSubject = json.loads(self.parameterNode.GetParameter("CurrentSubject"))
    sourceFiducial = self.parameterNode.GetNodeReference("SourceFiducial")
    targetFiducial = self.parameterNode.GetNodeReference("TargetFiducial")
    if int(self.parameterNode.GetParameter("SegmentMode")):
      LeadDBSCall.saveSegmentation(os.path.join(os.path.dirname(currentSubject["warpdrive_path"]), 'segmentations'))
    else:
      if sourceFiducial.GetNumberOfControlPoints(): # corrections made
        if self.hardenChangesAction.checked:
          sourceFiducial.Copy(targetFiducial) # set all as fixed points
        LeadDBSCall.saveSourceTarget(currentSubject["warpdrive_path"], sourceFiducial, targetFiducial)
        LeadDBSCall.saveSceneInfo(currentSubject["warpdrive_path"], self.inverseAction.checked)
        if self.hardenChangesAction.checked:
          slicerWillExit = True
          if self.inverseAction.checked:
            self.parameterNode.GetNodeReference("OutputGridTransform").Inverse()
          LeadDBSCall.applyChanges(self.parameterNode.GetNodeReferenceID("OutputGridTransform"), 
                                  self.parameterNode.GetNodeReference("ImageNode").GetStorageNode().GetFileName(), 
                                  os.path.join(self.parameterNode.GetParameter("MNIPath"), "t1.nii"),
                                  currentSubject["forward_transform"], 
                                  currentSubject["inverse_transform"], 
                                  currentSubject["warpdrive_path"],
                                  saveInExternalInstance)
    return slicerWillExit

  def setUpAtlases(self, info):
    print("Set up atlases")
    if info is not None:
      pass
    elif self.parameterNode.GetParameter("LeadAtlas"):
      info = {"atlasNames": [self.parameterNode.GetParameter("LeadAtlas")]}
    else:
      info = {"atlasNames": []}
    atlasNames = info["atlasNames"] if info["atlasNames"] != [] else ['DISTAL Nano (Ewert 2017)']
    # load atlas if not already in scene
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    folderNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLFolderDisplayNode')
    folderNodes.UnRegister(slicer.mrmlScene)
    for i in range(folderNodes.GetNumberOfItems()):
      folderNode = folderNodes.GetItemAsObject(i)
      folderItem = shNode.GetItemByDataNode(folderNode)
      if ('atlas' in shNode.GetItemAttributeNames(folderItem) and shNode.GetItemAttribute(folderItem,'atlas') == 'template') and (folderNode.GetName() in atlasNames):
        atlasNames.pop(atlasNames.index(folderNode.GetName()))
    for name in atlasNames:
      print("Loading atlas %s" % name)
      try:
        ImportAtlas.ImportAtlasLogic().readAtlas(os.path.join(ImportAtlas.ImportAtlasLogic().getAtlasesPath(), name, 'atlas_index.mat'))
      except:
        print("Could not load atlas %s" % name)
    if self.inverseAction.checked:
      self.invertAtlases(self.parameterNode.GetNodeReference("InputNode"), 1)

  def updateModalitiesToolButton(self):
    print("Update modalities")
    currentSubject = json.loads(self.parameterNode.GetParameter("CurrentSubject"))
    currentModality = self.modalitiesGroup.checkedAction().text
    subjectModalities = list(currentSubject["anat_files"].keys())
    for action in self.modalitiesGroup.actions():
      self.modalitiesGroup.removeAction(action)
      self.modalitiesMenu.removeAction(action)
    for subjectModality in subjectModalities:
      a = qt.QAction(self)
      a.setText(subjectModality)
      a.setCheckable(True)
      a.setChecked(subjectModality==currentModality)
      self.modalitiesGroup.addAction(a)
    self.modalitiesMenu.addActions(self.modalitiesGroup.actions())
    if not self.modalitiesGroup.checkedAction():
      a.setChecked(True)
      self.parameterNode.SetParameter("modality", a.text)

  def modalityChanged(self, action):
    self.updateModalitiesImages(action.text)

  def updateModalitiesImages(self, modality=None):
    if modality is None:
      modality = self.modalitiesGroup.checkedAction().text
    print("Loading %s modality" % modality)
    # find old nodes and delete
    slicer.mrmlScene.RemoveNode(self.parameterNode.GetNodeReference("ImageNode"))
    slicer.mrmlScene.RemoveNode(self.parameterNode.GetNodeReference("TemplateNode"))
    # initialize new image and init
    currentSubject = json.loads(self.parameterNode.GetParameter("CurrentSubject"))
    imageNode = slicer.util.loadVolume(currentSubject["anat_files"][modality], properties={'show':self.inverseAction.checked})  
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
    for sliceNode in slicer.util.getNodesByClass('vtkMRMLSliceNode'):
      sliceNode.Modified()
    # set parameter
    self.parameterNode.SetParameter("modality", modality)
    self.parameterNode.SetNodeReferenceID("ImageNode", imageNode.GetID())
    self.parameterNode.SetNodeReferenceID("TemplateNode", templateNode.GetID())
    if self.inverseAction.checked:
      templateNode.SetAndObserveTransformNodeID(self.parameterNode.GetNodeReferenceID("InputNode"))
    else:
      imageNode.SetAndObserveTransformNodeID(self.parameterNode.GetNodeReferenceID("InputNode"))

  def onInverseTriggered(self, useInverse):
    wasModified = self.parameterNode.StartModify()
    outputTransformID = self.parameterNode.GetNodeReferenceID("OutputGridTransform")
    slicer.mrmlScene.RemoveNode(self.parameterNode.GetNodeReference("InputNode"))
    # load
    currentSubject = json.loads(self.parameterNode.GetParameter("CurrentSubject"))
    inputNode = slicer.util.loadTransform(currentSubject["inverse_transform" if useInverse else "forward_transform"])
    # transform corrections
    sourceFiducialNode = self.parameterNode.GetNodeReference("SourceFiducial")
    targetFiducialNode = self.parameterNode.GetNodeReference("TargetFiducial")
    self.invertSourceTargetNodes(sourceFiducialNode, targetFiducialNode, inputNode)
    # transform atlas
    self.invertAtlases(inputNode, useInverse)
    # set output
    inputNode.SetAndObserveTransformNodeID(outputTransformID)
    # set parameter
    self.parameterNode.SetNodeReferenceID("InputNode", inputNode.GetID())
    self.parameterNode.SetNodeReferenceID("OutputGridTransform", outputTransformID)
    self.updateModalitiesImages()
    # set views
    if not useInverse:
      sliceNodes = slicer.util.getNodesByClass("vtkMRMLSliceNode")
      for sliceNode in sliceNodes:
          if sliceNode.GetName() in ['Red', 'Axial']:
              sliceNode.SetOrientationToAxial()
          elif sliceNode.GetName() in ['Green', 'Coronal']:
              sliceNode.SetOrientationToCoronal()
          elif sliceNode.GetName() in ['Yellow', 'Sagittal']:
              sliceNode.SetOrientationToSagittal()
    # update
    self.parameterNode.SetParameter("Update", "true")
    self.parameterNode.EndModify(wasModified)

  def invertAtlases(self, inputNode, useInverse=None):
    useInverse = useInverse if useInverse is not None else self.inverseAction.checked
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    nModels = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
    for i in range(nModels):
      modelNode = slicer.mrmlScene.GetNthNodeByClass(i, 'vtkMRMLModelNode')
      modelItem = shNode.GetItemByDataNode(modelNode)
      if 'atlas' in shNode.GetItemAttributeNames(modelItem):
        if shNode.GetItemAttribute(modelItem,'atlas') == 'template':
          modelNode.SetAndObserveTransformNodeID(inputNode.GetID() if useInverse else None)
        elif shNode.GetItemAttribute(modelItem,'atlas') == 'native':
          modelNode.SetAndObserveTransformNodeID(None if useInverse else inputNode.GetID())

  def invertSourceTargetNodes(self, sourceFiducialNode, targetFiducialNode, transformNode):
    sourcePoints = vtk.vtkPoints()
    sourceFiducialNode.SetAndObserveTransformNodeID(transformNode.GetID())
    sourceFiducialNode.HardenTransform()
    sourceFiducialNode.GetControlPointPositionsWorld(sourcePoints)
    targetPoints = vtk.vtkPoints()
    targetFiducialNode.SetAndObserveTransformNodeID(transformNode.GetID()) 
    targetFiducialNode.HardenTransform()
    targetFiducialNode.GetControlPointPositionsWorld(targetPoints)
    # update points
    sourceFiducialNode.SetControlPointPositionsWorld(targetPoints)
    targetFiducialNode.SetControlPointPositionsWorld(sourcePoints)

  def setSubjecApproved(self, value):
    LeadDBSCall.setApprovedData(json.loads(self.parameterNode.GetParameter("CurrentSubject"))["normlog_file"], int(value))
