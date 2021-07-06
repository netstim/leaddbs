import os
import slicer, vtk
import numpy as np

from slicer.util import VTKObservationMixin


class Trajectory(VTKObservationMixin):
  
  def __init__(self, N, distanceToTargetTransformID):
    VTKObservationMixin.__init__(self)

    self.trajectoryNumber = N
    self.alphaOmegaChannelNode = None

    # create folder to store ME nodes
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    self.folderID = shNode.CreateFolderItem(shNode.GetSceneItemID(), "Micro Electrode " + str(N))

    # translation transform from central ME
    self.translationTransform = self.createTranslationTransform()
    self.translationTransform.SetName("Translation Transform - ME: " + str(N))
    self.translationTransform.SetAndObserveTransformNodeID(distanceToTargetTransformID)
    shNode.SetItemParent(shNode.GetItemByDataNode(self.translationTransform), self.folderID)

    # 3D model
    self.microElectrodeModel = self.createMEModel()
    self.microElectrodeModel.SetName("ME Model - ME: " + str(N))
    self.microElectrodeModel.SetAndObserveTransformNodeID(self.translationTransform.GetID())
    shNode.SetItemParent(shNode.GetItemByDataNode(self.microElectrodeModel), self.folderID)

    # line
    self.trajectoryLine = self.createTrajectoryLine()
    self.trajectoryLine.SetName("Trajectory Line - ME: " + str(N))
    self.trajectoryLine.SetAndObserveTransformNodeID(self.translationTransform.GetID())
    shNode.SetItemParent(shNode.GetItemByDataNode(self.trajectoryLine), self.folderID)

    # tip
    self.tipFiducial = self.createTipFiducial()
    self.tipFiducial.SetName("Tip Fiducial - ME: " + str(N))
    self.tipFiducial.SetAndObserveTransformNodeID(self.translationTransform.GetID())
    shNode.SetItemParent(shNode.GetItemByDataNode(self.tipFiducial), self.folderID)

    # trace fiducials
    self.traceFiducials = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    self.traceFiducials.SetName("Trace Fiducial - ME: " + str(N))
    self.traceFiducials.GetDisplayNode().SetVisibility(0)
    shNode.SetItemParent(shNode.GetItemByDataNode(self.traceFiducials), self.folderID)

    # trace model
    self.traceModel = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode')
    self.traceModel.SetName("Trace Model - ME: " + str(N))
    self.traceModel.CreateDefaultDisplayNodes()
    self.traceModel.GetDisplayNode().SetAndObserveColorNodeID(slicer.util.getNode('Viridis').GetID())
    self.traceModel.GetDisplayNode().ScalarVisibilityOn()
    shNode.SetItemParent(shNode.GetItemByDataNode(self.traceModel), self.folderID)

    # observers
  
    # add fiducial every time the transform moves
    self.addObserver(self.translationTransform, slicer.vtkMRMLTransformNode.TransformModifiedEvent, self.onTransformModified)
    # update trace model every time the trace fiducials are modified 
    self.addObserver(self.traceFiducials, slicer.vtkMRMLMarkupsNode.PointAddedEvent,    self.updateModelFromFiducial)
    self.addObserver(self.traceFiducials, slicer.vtkMRMLMarkupsNode.PointModifiedEvent, self.updateModelFromFiducial)
    self.addObserver(self.traceFiducials, slicer.vtkMRMLMarkupsNode.PointRemovedEvent,  self.updateModelFromFiducial)


  def setAlphaOmegaChannelNode(self, channelNode):
    self.alphaOmegaChannelNode = channelNode
    self.cliNode = slicer.cli.run(slicer.modules.rootmeansquare, None, {'dataFileName': self.alphaOmegaChannelNode.GetChannelFullSavePath()})
    self.cliNode.SetAutoRun(1)
    self.addObserver(self.cliNode, slicer.vtkMRMLCommandLineModuleNode.StatusModifiedEvent, self.onCLIModified)

  def onCLIModified(self, caller, event):
    if self.cliNode.GetStatusString() == 'Completed':
      rmsValue = self.cliNode.GetParameterAsString('rootMeanSquare')
      fileName = os.path.basename(self.cliNode.GetParameterAsString('dataFileName'))
      fiducialLabels = vtk.vtkStringArray()
      self.traceFiducials.GetControlPointLabels(fiducialLabels)
      fiducialIndex = fiducialLabels.LookupValue("D = %.3f" % float(fileName[7:-3]))
      self.traceFiducials.SetNthControlPointDescription(fiducialIndex, rmsValue)
      self.cliNode.SetParameterAsString('dataFileName', self.alphaOmegaChannelNode.GetChannelFullSavePath())
      self.cliNode.InvokeEvent(slicer.vtkMRMLCommandLineModuleNode.AutoRunEvent)

  def onTransformModified(self, caller=None, event=None):
    # get current point
    currentPoint = [0.0] * 4
    matrix = vtk.vtkMatrix4x4()
    self.translationTransform.GetMatrixTransformToWorld(matrix)
    matrix.MultiplyPoint([0.0, 0.0, 0.0, 1.0], currentPoint)
    # add fiducial in current point with distance to target as name
    distanceToTargetNode = slicer.mrmlScene.GetNodeByID(self.translationTransform.GetTransformNodeID())
    distanceToTargetNode.GetMatrixTransformToParent(matrix)
    self.traceFiducials.AddFiducialFromArray(currentPoint[:3], "D = %.3f" % matrix.GetElement(2,3))

  def createTranslationTransform(self):
    index = np.array(np.unravel_index(self.trajectoryNumber, (3,3))) - 1  # from 0:8 to (-1,-1) : (1,1)
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

  def createMEModel(self):
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
    self.setDefaultOrientationToModel(model)
    return model

  def createTrajectoryLine(self):
    ls = vtk.vtkLineSource()
    ls.SetPoint1(0, -10, 0)
    ls.SetPoint2(0,  80, 0)
    ls.Update()
    model = slicer.modules.models.logic().AddModel(ls.GetOutput())
    model.CreateDefaultSequenceDisplayNodes()
    model.CreateDefaultDisplayNodes()
    model.SetDisplayVisibility(1)
    model.GetDisplayNode().SetColor(0.9,0.9,0.9)
    self.setDefaultOrientationToModel(model)
    return model

  def createTipFiducial(self):
    fiducialNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    fiducialNode.GetDisplayNode().SetTextScale(0)
    fiducialNode.GetDisplayNode().SetVisibility3D(1)
    fiducialNode.AddControlPointWorld(vtk.vtkVector3d([0, 0, 0]))
    fiducialNode.SetLocked(True)
    return fiducialNode

  def setDefaultOrientationToModel(self, model):
    # put in I-S axis
    vtkTransform = vtk.vtkTransform()
    vtkTransform.RotateWXYZ(90, 1, 0, 0)
    model.ApplyTransform(vtkTransform)


  def updateModelFromFiducial(self, caller=None, event=None):
    # get fiducial points to vtkPoints and description to vtkDoubleArray
    samplePoints = vtk.vtkPoints()
    valuesArray = []
    pos = [0.0] * 3
    for i in range(self.traceFiducials.GetNumberOfControlPoints()):
      if self.traceFiducials.GetNthControlPointDescription(i) in ['','nan']:
        pass
      else:
        self.traceFiducials.GetNthFiducialPosition(i,pos)
        samplePoints.InsertNextPoint(pos)
        valuesArray.append(float(self.traceFiducials.GetNthControlPointDescription(i)))
    if not valuesArray:
      return
    # array to vtk
    valuesMedian = np.median(valuesArray[:min(len(valuesArray),5)])
    vtkValuesArray = vtk.vtkDoubleArray()
    vtkValuesArray.SetName('values')
    for value in valuesArray:
      vtkValuesArray.InsertNextTuple((max((value-valuesMedian)/valuesMedian/2.0, 0.1),))
    # line source
    polyLineSource = vtk.vtkPolyLineSource()
    polyLineSource.SetPoints(samplePoints)
    polyLineSource.Update()
    # poly data
    polyData = polyLineSource.GetOutput()
    polyData.GetPointData().AddArray(vtkValuesArray)
    polyData.GetPointData().SetScalars(vtkValuesArray)
    # run tube filter
    tubeFilter = vtk.vtkTubeFilter()
    tubeFilter.SetInputData(polyData)
    tubeFilter.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    tubeFilter.SetNumberOfSides(16)
    tubeFilter.CappingOn()
    tubeFilter.Update()
    # smooth
    smoothFilter = vtk.vtkSmoothPolyDataFilter()
    smoothFilter.SetInputData(tubeFilter.GetOutput())
    smoothFilter.SetNumberOfIterations(2)
    smoothFilter.SetRelaxationFactor(0.5)
    smoothFilter.FeatureEdgeSmoothingOff()
    smoothFilter.BoundarySmoothingOn()
    smoothFilter.Update()
    # update
    self.traceModel.SetAndObservePolyData(smoothFilter.GetOutput())
    self.traceModel.GetDisplayNode().SetActiveScalarName('values')
    self.traceModel.GetDisplayNode().SetAutoScalarRange(False)
    self.traceModel.GetDisplayNode().SetScalarRange(0.0,1.0)
    self.traceModel.Modified()
    return 