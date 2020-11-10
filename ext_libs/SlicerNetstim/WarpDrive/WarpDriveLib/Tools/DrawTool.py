import vtk, slicer
import numpy as np


from ..Widgets.ToolWidget   import AbstractToolWidget
from ..Effects.DrawEffect import AbstractDrawEffect

from ..Helpers import GridNodeHelper, WarpDriveUtil

class DrawToolWidget(AbstractToolWidget):
    
  def __init__(self):
    toolTip = ''
    AbstractToolWidget.__init__(self, 'Draw', toolTip)


class DrawToolEffect(AbstractDrawEffect):


  def __init__(self, sliceWidget):
    AbstractDrawEffect.__init__(self, sliceWidget)


  def processEvent(self, caller=None, event=None):

    AbstractDrawEffect.processEvent(self, caller, event) 

    if event == 'LeftButtonReleaseEvent':

      sourceFiducial, targetFiducial = self.getSourceTargetFromDrawing()

      # reset
      self.resetPolyData()

      if not sourceFiducial: # return in case not created
        slicer.mrmlScene.RemoveNode(sourceFiducial)
        slicer.mrmlScene.RemoveNode(targetFiducial)
        return


      sourceFiducial.ApplyTransform(self.parameterNode.GetNodeReference("OutputGridTransform").GetTransformFromParent()) # undo current

      WarpDriveUtil.addCorrection(sourceFiducial, targetFiducial, 
                              spread=int(round(float(self.parameterNode.GetParameter("Spread")))),
                              referenceNode = self.parameterNode.GetNodeReference("InputNode"))   

      self.parameterNode.SetParameter("Update","true")
 

  def getSourceTargetFromDrawing(self):
    # get smplae distance
    sampleDistance = 1 #float(self.parameterNode.GetParameter("DrawSampleDistance"))

    # create curve from drawing and resample
    sourceCurve = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsCurveNode')
    sourceCurve.SetControlPointPositionsWorld(self.rasPoints)
    sourceCurve.ResampleCurveWorld(sampleDistance)  

    # get resampled points
    resampledPoints = vtk.vtkPoints()
    sourceCurve.GetControlPointPositionsWorld(resampledPoints)
    
    if resampledPoints.GetNumberOfPoints() <= 1:
      slicer.mrmlScene.RemoveNode(sourceCurve)
      return (None,)*2

    # get closest model sliced
    slicedModel, originalModel = self.sliceClosestModel(resampledPoints.GetPoint(0))

    if not slicedModel:
      slicer.mrmlScene.RemoveNode(sourceCurve)
      return (None,)*2
    
    # resample sourceCurve in sliced model with same amount of points
    targetCurve = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsCurveNode')
    targetCurve.GetDisplayNode().SetVisibility(0)
    targetCurve.SetControlPointPositionsWorld(resampledPoints)
    targetCurve.SetCurveTypeToShortestDistanceOnSurface(slicedModel)
    targetCurve.ResampleCurveSurface(sampleDistance, slicer.vtkMRMLModelNode().SafeDownCast(slicedModel), 0.0025)
    targetCurve.ResampleCurveWorld(targetCurve.GetCurveLengthWorld() / max((sourceCurve.GetNumberOfControlPoints() - 1), 1))
      
    # curve to fiducial
    sourceFiducial = self.curveToFiducial(sourceCurve)
    targetFiducial = self.curveToFiducial(targetCurve)

    # set name
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    modelParentName =  shNode.GetItemName(shNode.GetItemParent(shNode.GetItemByDataNode(originalModel)))
    targetFiducial.SetName(modelParentName + '_' + originalModel.GetName())

    # remove
    slicer.mrmlScene.RemoveNode(targetCurve)
    slicer.mrmlScene.RemoveNode(sourceCurve)
    slicer.mrmlScene.RemoveNode(slicedModel)

    return sourceFiducial, targetFiducial


  def copyControlPoints(self, sourceNode, targetNode):
    if targetNode.GetNumberOfControlPoints() == 0:
      label = '1'
    else:
      label = str( int(targetNode.GetNthFiducialLabel(targetNode.GetNumberOfControlPoints()-1)) + 1 )
    p = [0]*3
    for i in range(sourceNode.GetNumberOfControlPoints()):
      sourceNode.GetNthControlPointPosition(i,p)
      targetNode.AddFiducialFromArray(p, label)


  def curveToFiducial(self, curve):
    fiducial = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    fiducial.GetDisplayNode().SetVisibility(0)
    points = vtk.vtkPoints()
    curve.GetControlPointPositionsWorld(points)
    fiducial.SetControlPointPositionsWorld(points)
    return fiducial

  def sliceClosestModel(self, point):
    originalModel = None
    # set up plane
    normal = np.array([float(self.sliceLogic.GetSliceNode().GetName()==name) for name in ['Yellow','Green','Red']])
    plane = vtk.vtkPlane()
    plane.SetOrigin(point) # point in plane
    plane.SetNormal(normal)
    # set up cutter
    cutter = vtk.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetGenerateCutScalars(0)
    cutter.Update()
    # init point locator and output
    pointsLocator = vtk.vtkPointLocator() 
    globalMinDistance = 1000
    outPolyData = vtk.vtkPolyData()
    # iterate over models in scene
    nModels = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
    for i in range(nModels):
      model = slicer.mrmlScene.GetNthNodeByClass(i, 'vtkMRMLModelNode')
      polyData = model.GetPolyData()
      if model.GetDisplayNode() and model.GetDisplayNode().GetVisibility() and polyData.GetNumberOfCells() > 1 and model.GetName()!= 'auxSphereModel': # model visible and cells available
        cutter.SetInputData(polyData)
        cutter.Update()
        cutterOutput = cutter.GetOutput()
        if cutterOutput.GetNumberOfCells(): # model intersects with plane
          # get distance from input point to closest point in model
          pointsLocator.SetDataSet(cutterOutput)
          pointsLocator.BuildLocator()
          closestPoint = cutterOutput.GetPoint(pointsLocator.FindClosestPoint(point))
          localMinDistance = vtk.vtkMath().Distance2BetweenPoints(closestPoint, point)
          if localMinDistance < globalMinDistance: # new min
            outPolyData.DeepCopy(cutterOutput)
            globalMinDistance = localMinDistance
            originalModel = model
    # return in case no model found
    if not originalModel:
      return False, False
    # generate output
    triangulator = vtk.vtkContourTriangulator()
    triangulator.SetInputData(outPolyData)
    triangulator.Update()
    slicedModel = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode')
    slicedModel.SetAndObservePolyData(triangulator.GetOutput())
    slicedModel.CreateDefaultDisplayNodes()
    slicedModel.GetDisplayNode().SetVisibility(0)
    return slicedModel, originalModel

  def cleanup(self):
    AbstractDrawEffect.cleanup(self)

