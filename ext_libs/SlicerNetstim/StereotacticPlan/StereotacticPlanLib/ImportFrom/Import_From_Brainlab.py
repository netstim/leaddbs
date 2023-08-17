import qt, slicer
import numpy as np
import re
import json
from datetime import datetime
from DICOMLib import DICOMUtils
import pydicom as dicom
import StereotacticPlan

from .importerBase import ImporterDialogBase


class ImporterDialog(ImporterDialogBase):
    def __init__(self):
        ImporterDialogBase.__init__(self)
        self.importerName = 'Brainlab'
        self.fileSelectTitle = 'Select Planning PDF'
        self.fileSelectExt = 'pdf'        

class Importer():
    def __init__(self, filePath):
        self.logic =  StereotacticPlan.StereotacticPlanLogic()
        self.stereotaxyReport = StereotaxyReport(filePath)
        self.planningDictionary = self.stereotaxyReport.getArcSettings()

    def setACPCCoordinatesToParameterNode(self, parameterNode):
        parameterNode.SetParameter("Frame AC", self.stereotaxyReport.getCoordinates('AC', 'Headring') + ';XYZ')
        parameterNode.SetParameter("Frame PC", self.stereotaxyReport.getCoordinates('PC', 'Headring') + ';XYZ')
        parameterNode.SetParameter("Frame MS", self.stereotaxyReport.getCoordinates('MS', 'Headring') + ';XYZ')
        parameterNode.SetParameter("Reference AC", self.stereotaxyReport.getCoordinates('AC', 'DICOM') + ';RAS')
        parameterNode.SetParameter("Reference PC", self.stereotaxyReport.getCoordinates('PC', 'DICOM') + ';RAS')
        parameterNode.SetParameter("Reference MS", self.stereotaxyReport.getCoordinates('MS', 'DICOM') + ';RAS')
        parameterNode.SetParameter("ReferenceToFrameMode", "ACPC Register")

    def getReferenceToFrameTransform(self):
        referenceToFrameNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLinearTransformNode", "Reference To Frame")
        sourceCoordinates = [np.fromstring(self.stereotaxyReport.getCoordinates(x, 'DICOM'), dtype=float, sep=',') for x in ['AC', 'PC', 'MS']]
        targetCoordinates = [self.logic.transformCoordsFromXYZToRAS(np.fromstring(self.stereotaxyReport.getCoordinates(x, 'Headring'), dtype=float, sep=',')) for x in ['AC', 'PC', 'MS']]
        self.logic.runFiducialRegistration(referenceToFrameNode, sourceCoordinates, targetCoordinates)
        return referenceToFrameNode

    def getTrajectoryTransforms(self, importInFrameSpace):
        # Set up node
        trajectoryTransform = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLinearTransformNode", "Brainlab Trajectory " + self.stereotaxyReport.getTrajectoryInformation()['Name'])
        trajectoryTransform.SetAttribute('NetstimStereotacticPlan', '1')
        trajectoryTransform.SetAttribute('Mode', 'Target Mounting Ring Arc')
        if importInFrameSpace:
            trajectoryTransform.SetAttribute('Entry', self.stereotaxyReport.getCoordinates('Entry', 'Headring') + ';XYZ;1')
            trajectoryTransform.SetAttribute('Target', self.stereotaxyReport.getCoordinates('Target', 'Headring') + ';XYZ;1')
        else:
            trajectoryTransform.SetAttribute('Entry', self.stereotaxyReport.getCoordinates('Entry', 'DICOM') + ';RAS;0')
            trajectoryTransform.SetAttribute('Target', self.stereotaxyReport.getCoordinates('Target', 'DICOM') + ';RAS;0')
        trajectoryTransform.SetAttribute('Mounting', self.planningDictionary["Mounting"])
        trajectoryTransform.SetAttribute('Ring', self.planningDictionary["Ring Angle"])
        trajectoryTransform.SetAttribute('Arc', self.planningDictionary["Arc Angle"])
        trajectoryTransform.SetAttribute('Roll', '0')
        # Compute
        if importInFrameSpace:
            targetCoordinatesForComputation = self.logic.transformCoordsFromXYZToRAS(np.fromstring(self.stereotaxyReport.getCoordinates('Target', 'Headring'), dtype=float, sep=','))
        else:
            targetCoordinatesForComputation = np.fromstring(self.stereotaxyReport.getCoordinates('Target', 'DICOM'), dtype=float, sep=',') 
        self.logic.computeTrajectoryFromTargetMountingRingArc(trajectoryTransform,
                                                                targetCoordinatesForComputation,
                                                                self.planningDictionary["Mounting"],
                                                                float(self.planningDictionary["Ring Angle"]),
                                                                float(self.planningDictionary["Arc Angle"]))
        return [trajectoryTransform.GetID()]
    
    def getReferenceVolumeFromDICOM(self, DICOMDir):
        DICOMinfo = self.stereotaxyReport.getDICOMInformation()
        loadedNodeIDs = []
        with DICOMUtils.TemporaryDICOMDatabase() as database:
            DICOMUtils.importDicom(DICOMDir, database)
            series = SlicerDICOMDatabase().getSeriestMatchingDescriptionAndDateTime(DICOMinfo['SeriesDescription'], DICOMinfo['AcquisitionDateTime'])
            loadedNodeIDs.extend(DICOMUtils.loadSeriesByUID([series]))

        for nodeID in loadedNodeIDs[::-1]:
            volumeNode = slicer.util.getNode(nodeID)
            if re.search('.*' + DICOMinfo['SeriesDescription'] + '.*', volumeNode.GetName()):
                return volumeNode

        raise RuntimeError('Unable to find image in DICOM')


class SlicerDICOMDatabase():

    def __init__(self):
        self.db = slicer.dicomDatabase

    def getSeriestMatchingDescriptionAndDateTime(self, descriptionIn, dateTimeIn):
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

class StereotaxyReport():

    def __init__(self, PDFPath):
        try:
            import pdfplumber
        except:
            slicer.util.pip_install('pdfplumber')
        import pdfplumber    
        
        self.pdf = pdfplumber.open(PDFPath)
        self.pdfWidth = float(self.pdf.pages[0].width)
        self.pdfHeight = float(self.pdf.pages[0].height)


    def hasPatientID(self, patientID):
        return patientID == self.getPatientInformation()['Patient ID']

    def hasSide(self, side):
        return side in self.getTrajectoryInformation()['Name']

    def getTrajectoryInformation(self):
        cropRegion = (self.pdfWidth/2, 130, self.pdfWidth, 240)
        tableSettings = {
            "vertical_strategy": "text",
            "horizontal_strategy": "lines",
            "intersection_y_tolerance": 20,    
            "keep_blank_chars": True,
            }
        outList = self.pdf.pages[0].crop(cropRegion).extract_table(tableSettings)
        outDict = {r[0]:r[1] for r in outList}
        return outDict

    def getPatientInformation(self):
        cropRegion = (0, 130, self.pdfWidth/2, 240)
        tableSettings = {
            "vertical_strategy": "text",
            "horizontal_strategy": "lines",
            "intersection_y_tolerance": 20,    
            "keep_blank_chars": True,
            }
        outList = self.pdf.pages[0].crop(cropRegion).extract_table(tableSettings)
        outDict = {r[0]:r[1] for r in outList}
        return outDict

    def getArcSettings(self):
        cropRegion = (0, 419, self.pdfWidth, 480)
        tableSettings = {
            "vertical_strategy": "text",
            "horizontal_strategy": "text",
            "min_words_vertical": 0,
            "keep_blank_chars": True,
            }
        outList = self.pdf.pages[0].crop(cropRegion).extract_table(tableSettings)
        outList = [[outList[0][i], outList[1][i]] for i in range(len(outList[0]))] # Transpose
        outDict = {r[0]:r[1].split(' ')[0] for r in outList} # Remove units
        outDict["Headring Coordinates"] = ",".join([outDict[C] for C in ["X","Y","Z"]]) # Join X Y Z
        return outDict

    def getCoordinates(self, queryPoint, queryCoordinateSystem):
        # define crop bounding box and transform to RAS
        if queryCoordinateSystem == 'Headring':
            if queryPoint in ['Entry', 'Target']:
                PDFPage = 0
                cropBoundingBox = (0, 350 , self.pdfWidth, 395)
            elif queryPoint in ['AC', 'PC', 'MS']:
                PDFPage = 1
                cropBoundingBox = (0, self.pdfHeight * 0.57 , self.pdfWidth/2, self.pdfHeight * 0.85)
        elif queryCoordinateSystem == 'DICOM':
            PDFPage = 1
            cropBoundingBox = (self.pdfWidth/2, self.pdfHeight * 0.57 , self.pdfWidth, self.pdfHeight * 0.85)
        else:
            raise RuntimeError('Invalid queryCoordinateSystem: ' + queryCoordinateSystem)
        # extract text
        PDFText = self.pdf.pages[PDFPage].crop(cropBoundingBox).extract_text()
        # extract coords
        queryPoint = queryPoint + ' Point' if queryPoint in ['AC','PC','MS'] else queryPoint
        m = re.search('(?<=' + queryPoint + ')' + r'\s+[-]?\d+[.]\d+ mm' * 3, PDFText)
        xyz_str = m.group(0).split('mm')
        xyz_flt = [float(x) for x in xyz_str[:-1]]
        # transform
        if queryCoordinateSystem == 'DICOM':
            toRAS = np.array([[ -1,  0,  0,  0],
                                [  0, -1,  0,  0],
                                [  0,  0,  1,  0],
                                [  0,  0,  0,  1]])
            xyz_flt = np.dot(toRAS, np.append(xyz_flt, 1))[:3]
        return ','.join([str(x) for x in xyz_flt])

    def getDICOMInformation(self):
        hStart = self.findHeightContainingText(1, self.pdfHeight * 0.5, "DICOM Coordinates") + 15
        hEnd = self.findHeightContainingText(1, self.pdfHeight * 0.61, "X Y Z") - 5
        cropRegion = (self.pdfWidth/2, hStart , self.pdfWidth, hEnd)
        tableSettings = {
            "vertical_strategy": "text",
            "horizontal_strategy": "lines",
            "min_words_vertical": 1,
            "keep_blank_chars": True,
            "intersection_y_tolerance":15,
            "edge_min_length":15,
            "explicit_horizontal_lines":[hEnd],
            "explicit_vertical_lines":[570]
            }
        outList = self.pdf.pages[1].crop(cropRegion).extract_table(tableSettings)
        outDict = {r[0]:r[1].replace('\n','') for r in outList}
        outDict['SeriesDescription'] = outDict['Image Set']
        outDict['AcquisitionDateTime'] = datetime.strptime(outDict['Scanned'], '%m/%d/%Y, %I:%M %p')
        return outDict

    def findHeightContainingText(self, pageNumber, heightStart, matchText):
        t = None
        maxHeight = heightStart
        while not t or t.find(matchText)==-1:
            maxHeight = maxHeight+1
            t = self.pdf.pages[pageNumber].crop((0, heightStart , self.pdfWidth, maxHeight)).extract_text()
        return maxHeight

