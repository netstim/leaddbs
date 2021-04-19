import slicer, vtk
import numpy as np
import re


#
# Micro Electrode
#

def setMEVisibility(METype, enabled):
  # determine MRML node type
  if METype in ['Model', 'Trace', 'Line']:
    nodeType = 'vtkMRMLModelNode'
  else: # tip
    nodeType = 'vtkMRMLMarkupsFiducialNode'

  shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
  nNodes = slicer.mrmlScene.GetNumberOfNodesByClass(nodeType)
  for i in range(nNodes): # iterate
    node = slicer.mrmlScene.GetNthNodeByClass(i, nodeType)
    if METype in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(node)):
      node.GetDisplayNode().SetVisibility3D(enabled)
  

def GetMETraceFiducials():
  IDs = []
  shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
  nNodes = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLMarkupsFiducialNode')
  for i in range(nNodes): # iterate
    node = slicer.mrmlScene.GetNthNodeByClass(i, 'vtkMRMLMarkupsFiducialNode')
    if 'Trace' in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(node)):
      IDs.append(node.GetID())
  return IDs

def removeNthMicroElectrode(N):
  # remove nodes associated with the micro electrode
  shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
  IDs = vtk.vtkIdList()
  shNode.GetItemChildren(shNode.GetSceneItemID(), IDs, True)
  for i in range(IDs.GetNumberOfIds()):
    if shNode.GetItemAttribute(IDs.GetId(i), 'Micro Electrode Number') == str(N):
      shNode.RemoveItem(IDs.GetId(i))

def initNthMicroElectrode(N, distanceToTargetTransformID, toolButton=None):

  # create folder to store ME nodes
  shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
  folderID = shNode.CreateFolderItem(shNode.GetSceneItemID(), "Micro Electrode " + str(N))

  # translation transform from central ME
  translationTransform = createTranslationTransform(N)
  translationTransform.SetName("Translation Transform - ME: " + str(N))
  translationTransform.SetAndObserveTransformNodeID(distanceToTargetTransformID)
  shNode.SetItemParent(shNode.GetItemByDataNode(translationTransform), folderID)

  # ME 3D model
  MEModel = createMEModel()
  MEModel.SetName("ME Model - ME: " + str(N))
  MEModel.SetAndObserveTransformNodeID(translationTransform.GetID())
  shNode.SetItemParent(shNode.GetItemByDataNode(MEModel), folderID)
  shNode.SetItemAttribute(shNode.GetItemByDataNode(MEModel), 'Model', '1')

  # ME line
  MELine = createTrajectoryLine()
  MELine.SetName("Trajectory Line - ME: " + str(N))
  MELine.SetAndObserveTransformNodeID(translationTransform.GetID())
  shNode.SetItemParent(shNode.GetItemByDataNode(MELine), folderID)
  shNode.SetItemAttribute(shNode.GetItemByDataNode(MELine), 'Line', '1')

  # ME tip
  METipFiducial = createTipFiducial()
  METipFiducial.SetName("Tip Fiducial - ME: " + str(N))
  METipFiducial.SetAndObserveTransformNodeID(translationTransform.GetID())
  shNode.SetItemParent(shNode.GetItemByDataNode(METipFiducial), folderID)
  shNode.SetItemAttribute(shNode.GetItemByDataNode(METipFiducial), 'Tip', '1')

  # ME trace fiducials
  traceFiducials = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
  traceFiducials.SetName("Trace Fiducial - ME: " + str(N))
  traceFiducials.GetDisplayNode().SetVisibility(0)
  shNode.SetItemParent(shNode.GetItemByDataNode(traceFiducials), folderID)
  shNode.SetItemAttribute(shNode.GetItemByDataNode(traceFiducials), 'Trace', '1')

  # ME trace model
  traceModel = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode')
  traceModel.SetName("Trace Model - ME: " + str(N))
  traceModel.CreateDefaultDisplayNodes()
  traceModel.GetDisplayNode().SetAndObserveColorNodeID(slicer.util.getNode('Plasma').GetID())
  traceModel.GetDisplayNode().ScalarVisibilityOn()
  shNode.SetItemParent(shNode.GetItemByDataNode(traceModel), folderID)
  shNode.SetItemAttribute(shNode.GetItemByDataNode(traceModel), 'Trace', '1')

  # observers
  
  # add fiducial every time the transform moves
  translationTransform.AddObserver(slicer.vtkMRMLTransformNode.TransformModifiedEvent, lambda c,e: onTransformModified(c, traceFiducials))
  # update trace model every time the trace fiducials are modified 
  traceFiducials.AddObserver(slicer.vtkMRMLMarkupsNode.PointAddedEvent,    lambda c,e: updateModelFromFiducial(c, traceModel))
  traceFiducials.AddObserver(slicer.vtkMRMLMarkupsNode.PointModifiedEvent, lambda c,e: updateModelFromFiducial(c, traceModel))
  traceFiducials.AddObserver(slicer.vtkMRMLMarkupsNode.PointRemovedEvent,  lambda c,e: updateModelFromFiducial(c, traceModel))
  # update trace fiducial attribute with the AO channel name
  if toolButton:
    toolButton.actions()[0].menu().triggered.connect(lambda action: shNode.SetItemAttribute(shNode.GetItemByDataNode(traceFiducials), 'AO Channel', action.text))

  # set item attribute
  IDs = vtk.vtkIdList()
  shNode.GetItemChildren(folderID, IDs, True)
  IDs = [IDs.GetId(i) for i in range(IDs.GetNumberOfIds())]
  for ID in IDs + [folderID]:
    shNode.SetItemAttribute(ID, 'Micro Electrode Number', str(N))


def onTransformModified(caller, fiducial):
  # check if was set to active
  if fiducial.GetDescription() != 'Active':
    return
  # get current point
  currentPoint = [0.0] * 4
  matrix = vtk.vtkMatrix4x4()
  caller.GetMatrixTransformToWorld(matrix)
  matrix.MultiplyPoint([0.0, 0.0, 0.0, 1.0], currentPoint)
  # add fiducial in current point with distance to target as name
  distanceToTarget = slicer.mrmlScene.GetNodeByID(caller.GetTransformNodeID())
  fiducial.AddFiducialFromArray(currentPoint[:3], 'D = {:.3f}'.format(distanceToTarget.GetMatrixTransformToParent().GetElement(2,3)))

def createTranslationTransform(N):
  index = np.array(np.unravel_index(N, (3,3))) - 1  # from 0:8 to (-1,-1) : (1,1)
  # get translation ammount
  np.seterr(divide='ignore')
  hypotenuse = 2
  alpha = np.arctan(abs(np.divide(index[0],index[1]))) if any(index) else 0
  # create matrix
  m = vtk.vtkMatrix4x4()
  m.SetElement(0, 3, -index[0]*np.sin(alpha)*hypotenuse)
  m.SetElement(1, 3, -index[1]*np.cos(alpha)*hypotenuse)
  # add node
  transformNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLinearTransformNode')
  transformNode.SetMatrixTransformToParent(m)
  return transformNode

def createMEModel():
  # create cylinder (macro)
  cylinder = vtk.vtkCylinderSource()
  cylinder.SetRadius(0.35)
  cylinder.SetHeight(80)
  cylinder.SetCenter(0, cylinder.GetHeight()/2 + 3, 0)
  cylinder.SetResolution(12)
  cylinder.Update()
  # cone line (micro)
  cone = vtk.vtkConeSource()
  cone.SetRadius(0.2)
  cone.SetHeight(3)
  cone.SetDirection(0,-1,0)
  cone.SetCenter(0, 1.5, 0)
  cone.SetResolution(12)
  cone.Update()
  # append cone and cylinder
  input1 = vtk.vtkPolyData()
  input2 = vtk.vtkPolyData()
  input1.ShallowCopy(cylinder.GetOutput())
  input2.ShallowCopy(cone.GetOutput())
  appendFilter = vtk.vtkAppendPolyData()
  appendFilter.AddInputData(input1)
  appendFilter.AddInputData(input2)
  appendFilter.Update()
  cleanFilter = vtk.vtkCleanPolyData()
  cleanFilter.SetInputConnection(appendFilter.GetOutputPort())
  cleanFilter.Update()
  # add model node
  modelsLogic = slicer.modules.models.logic()
  model = modelsLogic.AddModel(cleanFilter.GetOutput())
  model.CreateDefaultSequenceDisplayNodes()
  model.CreateDefaultDisplayNodes()
  model.GetDisplayNode().SetColor(0,1,1)
  model.GetDisplayNode().SetVisibility2D(1)
  model.GetDisplayNode().SetVisibility3D(1)
  setDefaultOrientationToModel(model)
  return model

def createTrajectoryLine():
  ls = vtk.vtkLineSource()
  ls.SetPoint1(0, -10, 0)
  ls.SetPoint2(0,  80, 0)
  ls.Update()
  model = slicer.modules.models.logic().AddModel(ls.GetOutput())
  model.CreateDefaultSequenceDisplayNodes()
  model.CreateDefaultDisplayNodes()
  model.SetDisplayVisibility(1)
  model.GetDisplayNode().SetColor(0.9,0.9,0.9)
  setDefaultOrientationToModel(model)
  return model

def createTipFiducial():
  fiducialNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
  fiducialNode.GetDisplayNode().SetTextScale(0)
  fiducialNode.GetDisplayNode().SetVisibility3D(1)
  fiducialNode.AddControlPointWorld(vtk.vtkVector3d([0, 0, 0]))
  fiducialNode.SetLocked(True)
  return fiducialNode

def setDefaultOrientationToModel(model):
  # put in I-S axis
  vtkTransform = vtk.vtkTransform()
  vtkTransform.RotateWXYZ(90, 1, 0, 0)
  model.ApplyTransform(vtkTransform)


def updateModelFromFiducial(fiducialNode, modelNode):
  # set start modify
  wasModified = fiducialNode.StartModify()
  # get fiducial points to vtkPoints and description to vtkDoubleArray
  pts = vtk.vtkPoints()
  valuesArray = vtk.vtkDoubleArray()
  valuesArray.SetName('values')
  pos = [0.0] * 3
  i = 0
  while i < fiducialNode.GetNumberOfControlPoints():
    if fiducialNode.GetNthControlPointDescription(i) == 'nan':
      fiducialNode.RemoveNthControlPoint(i)
      i -= 1
    elif fiducialNode.GetNthControlPointDescription(i) == '':
      pass
    else:
      fiducialNode.GetNthFiducialPosition(i,pos)
      pts.InsertNextPoint(pos)
      valuesArray.InsertNextTuple((min(float(fiducialNode.GetNthControlPointDescription(i)),50.0),))
    # increase counter
    i += 1
  # end modify
  fiducialNode.EndModify(wasModified)
  # line source
  polyLineSource = vtk.vtkPolyLineSource()
  polyLineSource.SetPoints(pts)
  polyLineSource.Update()
  # poly data
  polyData = polyLineSource.GetOutput()
  polyData.GetPointData().AddArray(valuesArray)
  polyData.GetPointData().SetScalars(valuesArray)
  # run tube filter
  tubeFilter = vtk.vtkTubeFilter()
  tubeFilter.SetInputData(polyData)
  tubeFilter.SetVaryRadiusToVaryRadiusByScalar()
  tubeFilter.SetRadius(0.1)
  tubeFilter.SetRadiusFactor(10 * (valuesArray.GetValueRange()[-1] / 50))
  tubeFilter.SetRadiusFactor(10)
  tubeFilter.SetNumberOfSides(16)
  tubeFilter.CappingOn()
  tubeFilter.Update()
  # update
  modelNode.SetAndObservePolyData(tubeFilter.GetOutput())
  modelNode.GetDisplayNode().SetActiveScalarName('values')
  modelNode.GetDisplayNode().SetAutoScalarRange(False)
  modelNode.GetDisplayNode().SetScalarRange(0.0,50.0)
  modelNode.Modified()
  return 