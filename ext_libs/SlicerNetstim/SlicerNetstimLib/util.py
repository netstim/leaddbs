import os
import glob
import sys
import json
import re
import numpy as np
import slicer
from datetime import datetime
from DICOMLib import DICOMUtils
import pydicom as dicom

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "StereotacticPlan"))
from StereotacticPlanLib.util import StereotaxyReport

#------------------------------------------------------------------------------------------
class LeadDBSSubject():

  def __init__(self, subjectPath, leaddbsPath=None):
    if not os.path.exists(subjectPath):
      raise RuntimeError(subjectPath + ' Invalid input subject path')
    self.subjectPath = subjectPath
    self.subject = os.path.basename(subjectPath)
    self.leadORPath = os.path.join(self.subjectPath,'leador')
    self.stereotaxyReports = self.getStereotaxyReports()
    self.leaddbsPath = leaddbsPath

  def getStereotaxyReports(self, patientID=None, maxOutput=2):
    out = []
    for possiblePath in self.getPossibleStereotaxyReportPaths(rootPath=self.leadORPath):
      stereotaxyReport = StereotaxyReport(possiblePath)
      if  (patientID is None) or stereotaxyReport.hasPatientID(patientID):
        out.append(stereotaxyReport)
      if len(out) == maxOutput:
        break
    return out

  def getPossibleStereotaxyReportPaths(self, rootPath):
    return glob.glob(os.path.join(rootPath, 'StereotaxyReport*.pdf'))

  def createORScene(self):
    anatVolumeNode = self.getAnatVolumeNode()
    anatToframeTransformNode =  slicer.util.loadTransform(os.path.join(self.leadORPath, 'anatToFrame.txt'))
    glanatInverseTransformNode = self.getInverseNormalizationNode()

    import ImportAtlas
    atlasPath = os.path.join(self.leaddbsPath,'templates','space','MNI152NLin2009bAsym','atlases','DISTAL Nano (Ewert 2017)','atlas_index.mat')
    ImportAtlas.ImportAtlasLogic().readAtlas(atlasPath)

    anatVolumeNode.ApplyTransform(anatToframeTransformNode.GetTransformToParent())
    anatVolumeNode.HardenTransform()

    shnode = slicer.mrmlScene.GetSubjectHierarchyNode()

    for i in range(slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')):
      modelNode = slicer.mrmlScene.GetNthNodeByClass(i, 'vtkMRMLModelNode')
      if shnode.GetItemAttribute(shnode.GetItemByDataNode(modelNode),'atlas') == '1':
        modelNode.ApplyTransform(glanatInverseTransformNode.GetTransformToParent())
        modelNode.ApplyTransform(anatToframeTransformNode.GetTransformToParent())
        modelNode.HardenTransform()

    slicer.mrmlScene.RemoveNode(glanatInverseTransformNode)
    slicer.mrmlScene.RemoveNode(anatToframeTransformNode)

    import StereotacticPlan
    for stereotaxyReport in self.stereotaxyReports:
      name = 'Planning ' + stereotaxyReport.getTrajectoryInformation()['Name']
      outputTransform = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLinearTransformNode', name)
      planningDictionary = stereotaxyReport.getArcSettings()
      StereotacticPlan.StereotacticPlanLogic().process(outputTransform, 
        float(planningDictionary['Arc Angle']), 
        float(planningDictionary['Ring Angle']), 
        np.array(planningDictionary['Headring Coordinates'].split(','), dtype=float), 
        planningDictionary['Mounting'])

    sceneSaveFilename = os.path.join(self.leadORPath, 'ORScene.mrb')
    slicer.util.saveScene(sceneSaveFilename)

  def createAnatToFrameTransform(self):
    # frame fiducials
    frameFidNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode','FrameFid')
    frameFidNode.AddFiducialFromArray(self.stereotaxyReports[0].getCoordinates('AC', 'Headring', 'RAS'), 'frameAC')
    frameFidNode.AddFiducialFromArray(self.stereotaxyReports[0].getCoordinates('PC', 'Headring', 'RAS'), 'framePC')
    frameFidNode.AddFiducialFromArray(self.stereotaxyReports[0].getCoordinates('MS', 'Headring', 'RAS'), 'frameMS')
    # anat fiducials
    anatFidNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode','AnatFid')
    anatFidNode.AddFiducialFromArray(self.stereotaxyReports[0].getCoordinates('AC', 'DICOM', 'RAS'), 'anatAC')
    anatFidNode.AddFiducialFromArray(self.stereotaxyReports[0].getCoordinates('PC', 'DICOM', 'RAS'), 'anatPC')
    anatFidNode.AddFiducialFromArray(self.stereotaxyReports[0].getCoordinates('MS', 'DICOM', 'RAS'), 'anatMS')
    anatFidNode.ApplyTransform(self.getRawAnatToAnatTransform().GetTransformToParent())
    # fiducial registration
    anatToFrameTransformNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLinearTransformNode')
    parameters = {}
    parameters['fixedLandmarks']  = frameFidNode.GetID()
    parameters['movingLandmarks'] = anatFidNode.GetID()
    parameters['saveTransform']   = anatToFrameTransformNode.GetID()
    parameters['transformType']   = 'Rigid'
    cli = slicer.cli.run(slicer.modules.fiducialregistration, None, parameters, wait_for_completion=True, update_display=False)
    slicer.util.saveNode(anatToFrameTransformNode, os.path.join(self.leadORPath,'anatToFrame.txt'))

  def getRawAnatToAnatTransform(self):
    dcmInfo = self.stereotaxyReports[0].getDICOMInformation()
    modality = self.getModalityFromSeriesDescription(dcmInfo['SeriesDescription'])
    # define nodes
    rawAnatVolumeNode = self.getImageFromDICOMInformation(dcmInfo)
    anatVolumeNode 	  = self.getAnatVolumeNode(modality)
    outTransformNode  = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLinearTransformNode')
    # registration
    parameters = {}
    parameters['fixedVolume'] 			      = anatVolumeNode.GetID()
    parameters['movingVolume'] 			      = rawAnatVolumeNode.GetID()
    parameters['samplingPercentage'] 	    = 0.02
    parameters['linearTransform'] 		    = outTransformNode.GetID()
    parameters['initializeTransformMode'] = 'useMomentsAlign'
    parameters['useRigid'] 				        = True
    parameters['costMetric'] 			        = 'MMI'
    cli = slicer.cli.run(slicer.modules.brainsfit, None, parameters, wait_for_completion=True, update_display=False)
    slicer.util.saveNode(outTransformNode, os.path.join(self.leadORPath, rawAnatVolumeNode.GetName().split(' ')[-1] + 'ToAnat.txt'))
    return outTransformNode

  def getInverseNormalizationNode(self):
    g = glob.glob(os.path.join(self.subjectPath, 'normalization', 'transformations', '*from-MNI152NLin2009bAsym*-ants.nii.gz'))
    if g:
      return slicer.util.loadTransform(g[0])
    raise RuntimeError('Unable to find inverse normalization : ' + self.subjectPath)

  def getAnatVolumeNode(self, modality=None):
    modalities = [modality] if modality is not None else ['T2','T1']
    for m in modalities:
      g = glob.glob(os.path.join(self.subjectPath, 'coregistration', 'anat', '*' + m.upper() + 'w.nii'))
      if g:
        return slicer.util.loadVolume(g[0])
    raise RuntimeError('Unable to find anat volume : ' + self.subjectPath)

  def getModalityFromSeriesDescription(self, seriesDescription):
    for g in glob.glob(os.path.join(self.subjectPath, 'coregistration', 'anat', '*preop*w.nii')):
      modality = re.search(r'\w+(?=w.nii)', g).group(0)
      if modality.lower() in seriesDescription.lower():
        return modality
    return 't1'

  def getImageFromDICOMInformation(self, dcmInfo):
    loadedNodeIDs = []
    with DICOMUtils.TemporaryDICOMDatabase() as database:
      DICOMUtils.importDicom(os.path.join(self.subjectPath,'..','..','..','sourcedata',self.subject), database)
      series = SlicerDICOMDatabase().geSeriestMatchingDescriptionAndDateTime(dcmInfo['SeriesDescription'], dcmInfo['AcquisitionDateTime'])
      loadedNodeIDs.extend(DICOMUtils.loadSeriesByUID([series]))

    for nodeID in loadedNodeIDs[::-1]:
      volumeNode = slicer.util.getNode(nodeID)
      if re.search('.*' + dcmInfo['SeriesDescription'] + '.*', volumeNode.GetName()):
        return volumeNode

    raise RuntimeError('Unable to find image in DICOM: ' + self.subjectPath)

#------------------------------------------------------------------------------------------
class SlicerDICOMDatabase():

  def __init__(self):
      self.db = slicer.dicomDatabase

  def geSeriestMatchingDescriptionAndDateTime(self, descriptionIn, dateTimeIn):
    for patient in self.db.patients():
      for study in self.db.studiesForPatient(patient):
        for series in self.db.seriesForStudy(study):
          try:
            seriesDescription = self.getSeriesAcquisitionInformationFromTag(series, 'SeriesDescription')
            acquisitionDateTime = self.getSeriesAcquisitionInformationFromTag(series, 'AcquisitionDateTime')
            if not acquisitionDateTime:
              acquisitionDate = self.getSeriesAcquisitionInformationFromTag(series, 'AcquisitionDate')
              acquisitionTime = self.getSeriesAcquisitionInformationFromTag(series, 'AcquisitionTime')
              acquisitionDateTime = acquisitionDate + acquisitionTime
            dateTime = self.DICOMDateTimeStringToDateTime(acquisitionDateTime)
          except:
            continue
          descriptionMatch = self.seriesDescriptionMatch(seriesDescription, descriptionIn)
          dateTimeMatch = self.seriesDateTimeMatch(dateTime, dateTimeIn)
          if descriptionMatch and dateTimeMatch:
            return series

  def getSeriesAcquisitionInformationFromTag(self, series, tag):
    fileList = self.db.filesForSeries(series)
    tagStr = str(dicom.tag.Tag(tag))[1:-1].replace(' ','')
    return self.db.fileValue(fileList[0],tagStr)

  @staticmethod
  def DICOMDateTimeStringToDateTime(dateTimeString):
    return datetime.strptime(dateTimeString, '%Y%m%d%H%M%S.%f')

  @staticmethod
  def seriesDescriptionMatch(seriesDescriptionA, seriesDescriptionB):
    return seriesDescriptionA == seriesDescriptionB
  
  @staticmethod
  def seriesDateTimeMatch(seriesDateTimeA, seriesDateTimeB):
    timeDelta = seriesDateTimeA - seriesDateTimeB
    return abs(timeDelta.total_seconds()) < 60 # allow 1 min because of different resolutions


#------------------------------------------------------------------------------------------
class NORASubject(LeadDBSSubject):

  def __init__(self, subjectPath):
    self.nodeInfo = self.getNodeInfo(subjectPath)
    leaddbsPath = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'leaddbs')
    LeadDBSSubject.__init__(self, subjectPath, leaddbsPath)

  def getNodeInfo(self, subjectPath):
    with open(os.path.abspath(os.path.join(subjectPath,'..','nodeinfo.json'))) as jsonFile:
      jsonLoad = json.load(jsonFile)
    return jsonLoad

  def getStereotaxyReports(self, maxOutput=2):
    return super().getStereotaxyReports(patientID=self.nodeInfo['PatientID'], maxOutput=maxOutput)

  def getPossibleStereotaxyReportPaths(self, rootPath):
    paths = super().getPossibleStereotaxyReportPaths(rootPath)
    pathsAppend = super().getPossibleStereotaxyReportPaths('/mnt/NCBrainlab')
    pathsAppend.sort(key=os.path.getmtime, reverse=True) # First is newer
    return paths + pathsAppend

  def getImageFromDICOMInformation(self, dcmInfo):
    for seriesKey, seriesMetaData in self.nodeInfo['Series'].items():
      descriptionMatch = SlicerDICOMDatabase.seriesDescriptionMatch(seriesMetaData['SeriesDescription'], dcmInfo['SeriesDescription'])
      aquisitionDateTime = SlicerDICOMDatabase.DICOMDateTimeStringToDateTime(seriesMetaData['AcquisitionDate'] + seriesMetaData['AcquisitionTime'])
      dateTimeMatch = SlicerDICOMDatabase.seriesDateTimeMatch(aquisitionDateTime, dcmInfo['AcquisitionDateTime'])
      if descriptionMatch and dateTimeMatch:
        return slicer.util.loadVolume(os.path.join(self.subjectPath, '..', seriesKey + '.nii'))
    raise RuntimeError('Could not find matching DICOM for: ' + self.subjectPath)
