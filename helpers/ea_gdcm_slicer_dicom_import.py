from DICOMLib import DICOMUtils
import os

loadedNodeIDs = []
with DICOMUtils.TemporaryDICOMDatabase() as db: 
  DICOMUtils.importDicom(dicomDataDir, db)
  for patientUID in db.patients():
    loadedNodeIDs.extend(DICOMUtils.loadPatientByUID(patientUID))
loadedNodes = [getNode(nodeID) for nodeID in loadedNodeIDs if getNodes(nodeID)]
for node in loadedNodes: 
  saveNode(node, os.path.join(outputDir, node.GetName() + '.nii.gz'))
exit()
   