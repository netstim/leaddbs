import re
import slicer
from DICOMLib import DICOMUtils
import numpy as np
import StereotacticPlan

from .importerBase import ImporterDialogBase

class ImporterDialog(ImporterDialogBase):
    def __init__(self):
        ImporterDialogBase.__init__(self)
        self.importerName = 'ROSA'
        self.fileSelectTitle = 'Select ROSA file'
        self.fileSelectExt = 'ros'
class Importer():
    def __init__(self, filePath):
        self.logic = StereotacticPlan.StereotacticPlanLogic()
        self.manager = ROSAManager(filePath)

    def setACPCCoordinatesToParameterNode(self, parameterNode):
        parameterNode.SetParameter("Reference AC", self.manager.getCoordinates('AC') + ';RAS')
        parameterNode.SetParameter("Reference PC", self.manager.getCoordinates('PC') + ';RAS')
        parameterNode.SetParameter("Reference MS", self.manager.getCoordinates('IH') + ';RAS')
        parameterNode.SetParameter("ReferenceToFrameMode", "ACPC Align")

    def getReferenceToFrameTransform(self):
        referenceToFrameNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLinearTransformNode", "Reference To Frame")
        self.logic.runACPCAlignment(referenceToFrameNode, 
                                np.fromstring(self.manager.getCoordinates('AC'), dtype=float, sep=','),
                                np.fromstring(self.manager.getCoordinates('PC'), dtype=float, sep=','),
                                np.fromstring(self.manager.getCoordinates('IH'), dtype=float, sep=','))
        return referenceToFrameNode

    def getTrajectoryTransforms(self, importInFrameSpace):
        loadedNodeIDs = []
        for rosa_trajectory in self.manager.getTrajectoriesList():
            # Set up node
            trajectoryTransform = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLinearTransformNode", "ROSA Trajectory " + rosa_trajectory['name'])
            trajectoryTransform.SetAttribute('NetstimStereotacticPlan', '1')
            trajectoryTransform.SetAttribute('Mode', 'Target Entry Roll')
            trajectoryTransform.SetAttribute('Entry', rosa_trajectory['entry'] + ';RAS;0')
            trajectoryTransform.SetAttribute('Target', rosa_trajectory['target'] + ';RAS;0')
            trajectoryTransform.SetAttribute('Mounting', 'lateral-left')
            trajectoryTransform.SetAttribute('Arc', '90')
            trajectoryTransform.SetAttribute('Ring', '90')
            trajectoryTransform.SetAttribute('Roll', '0')
            loadedNodeIDs.append(trajectoryTransform.GetID())
            # Compute
            targetCoordinatesForComputation = np.fromstring(rosa_trajectory['target'], dtype=float, sep=',') 
            entryCoordinatesForComputation = np.fromstring(rosa_trajectory['entry'], dtype=float, sep=',') 
            if importInFrameSpace:
                transformNode = self.getReferenceToFrameTransform()
                transform = slicer.util.array(transformNode.GetID())
                targetCoordinatesForComputation = np.dot(transform, np.append(targetCoordinatesForComputation, 1))[:3]
                entryCoordinatesForComputation = np.dot(transform, np.append(entryCoordinatesForComputation, 1))[:3]
                slicer.mrmlScene.RemoveNode(transformNode)
            self.logic.computeTrajectoryFromTargetEntryRoll(trajectoryTransform,
                                targetCoordinatesForComputation,
                                entryCoordinatesForComputation,
                                0)
        return loadedNodeIDs
    
    def getReferenceVolumeFromDICOM(self, DICOMDir):
        first_series_uid = self.manager.getFirstSeriesUID()
        with DICOMUtils.TemporaryDICOMDatabase() as database:
            DICOMUtils.importDicom(DICOMDir, database)
            rosa_reference_node_ID = DICOMUtils.loadSeriesByUID([first_series_uid])[0]
        slicer.modules.volumes.logic().CenterVolume(slicer.util.getNode(rosa_reference_node_ID))
        if slicer.app.layoutManager():
            slicer.app.layoutManager().resetSliceViews()
        return slicer.util.getNode(rosa_reference_node_ID)
class ROSAManager:
    def __init__(self, ros_file_path):
        with open(ros_file_path, 'r') as f:
            self.ros_txt = f.read()
      
    def getFirstSeriesUID(self):
        return re.search(r"(?<=\[SERIE_UID\]\n).*", self.ros_txt).group()

    def getTrajectoriesList(self):
        pattern = r"(?P<name>\w+) (?P<type>\d) (?P<color>\d+) (?P<entry_point_defined>\d) (?P<entry>-?\d+\.\d+ -?\d+\.\d+ -?\d+\.\d+) (?P<target_point_defined>\d) (?P<target>-?\d+\.\d+ -?\d+\.\d+ -?\d+\.\d+) (?P<instrument_length>\d+\.\d+) (?P<instrument_diameter>\d+\.\d+)\n"
        trajectories = [m.groupdict() for m in re.finditer(pattern, self.ros_txt)]
        for trajectory in trajectories:
          for pos in ['entry', 'target']:
            trajectory[pos] = np.array(list(map(float, trajectory[pos].split(' ')))) # str to array
            trajectory[pos] = trajectory[pos] * np.array([-1, -1, 1]) # LPS to RAS
            trajectory[pos] = ','.join([str(x) for x in trajectory[pos]])
        return trajectories

    def getCoordinates(self, queryPoint):
        pattern = r"(?<=\[ACPC\]).*" + queryPoint + r" \d -?\d+\.\d+ -?\d+\.\d+ -?\d+\.\d+"
        m = re.search(pattern, self.ros_txt, re.DOTALL)
        if not m:
          raise RuntimeError('Unable to find: %s' % queryPoint)
        coords_str = m.group().split(' ')[-3:]
        coords_lps = np.array(list(map(float, coords_str)))
        coords_ras = coords_lps * np.array([-1, -1, 1])
        return ','.join([str(x) for x in coords_ras])
