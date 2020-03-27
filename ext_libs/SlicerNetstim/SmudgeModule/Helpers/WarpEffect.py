import slicer, vtk, qt
import numpy as np
import sys, os

from . import PointerEffect

import TransformsUtil





class WarpEffectTool():

  _instances = set()

  def __init__(self):
    self._instances.add(self)


  def removeRedoTransform(self, redoTransformID):
    if redoTransformID != "":
      slicer.mrmlScene.RemoveNode(slicer.util.getNode(self.parameterNode.GetParameter("redoTransformID")))
      self.parameterNode.SetParameter("redoTransformID","")

  def notityWarpModified(self, parameterNode):
    parameterNode.SetParameter("warpModified",str(1+int(parameterNode.GetParameter("warpModified"))))

  def cleanup(self):
    pass

  
  @classmethod
  def empty(cls):
    # clean instances and reset
    for inst in cls._instances:
      inst.cleanup()
    cls._instances = set()



#
# None Effect
#

class NoneEffect(PointerEffect.PointerEffectTool, WarpEffectTool):

  def __init__(self, parameterNode, sliceWidget):
    self.parameterNode = parameterNode
    PointerEffect.PointerEffectTool.__init__(self, sliceWidget)
    WarpEffectTool.__init__(self)



#
# SmudgeEffectTool
#

class SmudgeEffectTool(PointerEffect.CircleEffectTool, WarpEffectTool):

  def __init__(self, parameterNode, sliceWidget, auxTransformNode):

    self.parameterNode = parameterNode
    PointerEffect.CircleEffectTool.__init__(self, sliceWidget)
    WarpEffectTool.__init__(self)

    self.warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))
    
    # transform data
    self.auxTransformNode = auxTransformNode
    self.auxTransformSpacing = self.auxTransformNode.GetTransformFromParent().GetDisplacementGrid().GetSpacing()[0] # Asume isotropic!
    self.auxTransfromRASToIJK = TransformsUtil.TransformsUtilLogic().getTransformRASToIJK(self.auxTransformNode)  

    self.previousPoint = [0,0,0]   
    self.smudging = False


  def processEvent(self, caller=None, event=None):

    PointerEffect.CircleEffectTool.processEvent(self, caller, event)

    if event == 'LeftButtonPressEvent':
      self.smudging = True
      self.warpNode.SetAndObserveTransformNodeID(self.auxTransformNode.GetID())    
      self.auxTransformArray = slicer.util.array(self.auxTransformNode.GetID())
      xy = self.interactor.GetEventPosition()
      xyToRAS = self.sliceLogic.GetSliceNode().GetXYToRAS()
      self.previousPoint = xyToRAS.MultiplyDoublePoint( (xy[0], xy[1], 0, 1) )[0:3]
    elif event == 'LeftButtonReleaseEvent':
      self.smudging = False
      self.warpNode.HardenTransform()
      TransformsUtil.TransformsUtilLogic().emptyTransfrom(self.auxTransformNode)
      self.notityWarpModified(self.parameterNode)
      self.removeRedoTransform(self.parameterNode.GetParameter("redoTransformID"))


    elif event == 'MouseMoveEvent':
      if self.smudging:
        # create a sphere with redius
        r = int(round(float(self.parameterNode.GetParameter("radius")) / self.auxTransformSpacing))
        xx, yy, zz = np.mgrid[:2*r+1, :2*r+1, :2*r+1]
        sphereResult = (xx-r) ** 2 + (yy-r) ** 2 + (zz-r) ** 2
        sphereResult[r][r][r] = 1 # replace 0 by 1
        sphereLarge = sphereResult <= (r**2+1) # sphere that the mouse shows
        sphereSmall = sphereResult <= ((r * (1 - float(self.parameterNode.GetParameter("hardness")) / 100.0)) **2 + 1 ) # hardness amount
        sphereResult = 1.0 / sphereResult # invert
        # get value in the edge of the small sphere
        i1,i2,i3 = np.nonzero(sphereSmall)
        newMaxValue = sphereResult[i1[0]][i2[0]][i3[0]]
        # set same value inside the small sphere
        sphereResult[sphereSmall] = sphereSmall[sphereSmall] * newMaxValue
        # delete outside values 
        sphereResult = sphereResult * sphereLarge
        # set range to [0-1]
        newMinValue = sphereResult.min()
        sphereResult = (sphereResult - newMinValue) / (newMaxValue - newMinValue)
        # set force
        sphereResult = sphereResult * float(self.parameterNode.GetParameter("force")) / 100.0
        # get current IJK coord
        xy = self.interactor.GetEventPosition()
        xyToRAS = self.sliceLogic.GetSliceNode().GetXYToRAS()
        currentPoint = xyToRAS.MultiplyDoublePoint( (xy[0], xy[1], 0, 1) )[0:3]
        pos_i,pos_j,pos_k,aux = self.auxTransfromRASToIJK.MultiplyDoublePoint(currentPoint + (1,))
        k,j,i = int(round(pos_k)),int(round(pos_j)),int(round(pos_i))
        curr_index = slice(k-r,k+r+1), slice(j-r,j+r+1), slice(i-r,i+r+1)

        # apply to transform array
        self.auxTransformArray[curr_index] += np.stack([(sphereResult) * i for i in (np.array(self.previousPoint) - np.array(currentPoint))],3) # original

        # update view
        self.auxTransformNode.Modified()
        # update previous point
        self.previousPoint = currentPoint

  def cleanup(self):
    slicer.mrmlScene.RemoveNode(self.auxTransformNode)
    WarpEffectTool.cleanup(self)
    PointerEffect.CircleEffectTool.cleanup(self)


#
# Snap Effect Tool
#

class SnapEffectTool(PointerEffect.DrawEffectTool, WarpEffectTool):

    
  def __init__(self, parameterNode, sliceWidget):

    WarpEffectTool.__init__(self)
    PointerEffect.DrawEffectTool.__init__(self,sliceWidget)

    self.parameterNode = parameterNode
    self.warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))
    self.MNIBounds = [[-98.25, -134.25, -72.25], [-98.25, -134.25, 116.75], [-98.25, 98.75, -72.25], [-98.25, 98.75, 116.75], [98.75, -134.25, -72.25], [98.75, -134.25, 116.75], [98.75, 98.75, -72.25], [98.75, 98.75, 116.75]]
    
  def processEvent(self, caller=None, event=None):

    if event == 'LeftButtonReleaseEvent':
      self.actionState = None # stop drawing

      # create open curve from drawing
      curve1 = self.pointsToCurve(self.rasPoints)
      points1 = vtk.vtkPoints()
      curve1.GetControlPointPositionsWorld(points1)

      # get closest model to points
      modelNode = self.getClosestModel(points1)
      if not modelNode: # clean and continue
        slicer.mrmlScene.RemoveNode(curve1)
      else:
        # cut model with current plane
        slicedModel = self.cutModelWithSliceIntersection(modelNode)
        # get curve on sliced model
        points2 = self.samplePointsInModel(points1, slicedModel)
        curve2 = self.pointsToCurve(points2)
        # in case line is inside the model the curve will be the same. set curve type to shortest distance to surface and resample again
        curve2.SetCurveTypeToShortestDistanceOnSurface(slicedModel)
        curve2.ResampleCurveSurface(0.2, slicer.vtkMRMLModelNode().SafeDownCast(slicedModel), 0.0025)
        # get same number of points as other curve
        curve2.ResampleCurveWorld(curve2.GetCurveLengthWorld() / max((curve1.GetNumberOfControlPoints() - 1), 1) )
        curve2.GetControlPointPositionsWorld(points2)
        # calculate transform and apply
        transformNode = self.createTransform(points1, points2, modelNode)
        self.displayTransformFromPoints(transformNode, points1)
        self.warpNode.SetAndObserveTransformNodeID(transformNode.GetID())
        self.warpNode.HardenTransform()
        # cleanup
        slicer.mrmlScene.RemoveNode(curve1)
        slicer.mrmlScene.RemoveNode(curve2)
        slicer.mrmlScene.RemoveNode(transformNode)
        slicer.mrmlScene.RemoveNode(slicedModel)
        self.removeRedoTransform(self.parameterNode.GetParameter("redoTransformID"))
        self.notityWarpModified(self.parameterNode)

    PointerEffect.DrawEffectTool.processEvent(self, caller, event) 
    
  def samplePointsInModel(self, points, model, sampleDist = 0.2, maximumSearchRadius = 0.0025):
    auxCurve = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsCurveNode')
    auxCurve.SetControlPointPositionsWorld(points)
    constraintNode = slicer.vtkMRMLModelNode().SafeDownCast(model)
    auxCurve.ResampleCurveSurface(sampleDist, constraintNode, maximumSearchRadius)
    outPoints = vtk.vtkPoints()
    auxCurve.GetControlPointPositionsWorld(outPoints)
    slicer.mrmlScene.RemoveNode(auxCurve)
    return outPoints

  def getClosestModel(self, points):
    outModel = None
    ruler = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLAnnotationRulerNode')
    ruler.GetDisplayNode().SetVisibility(0)
    ruler.SetControlPoint(0,np.array(points.GetPoint(0)),0,0)
    models = slicer.mrmlScene.GetNodesByClass('vtkMRMLModelNode')
    minDistance = 10000.0
    for i in range(models.GetNumberOfItems()):
      model = models.GetItemAsObject(i)
      if bool(model.GetDisplayNode().GetVisibility()):
        pointsInModel = self.samplePointsInModel(points, model)
        ruler.SetControlPoint(1,np.array(pointsInModel.GetPoint(0)),0,0)
        if ruler.GetDistanceMeasurement() < minDistance:
          minDistance = ruler.GetDistanceMeasurement()
          outModel = model
    slicer.mrmlScene.RemoveNode(ruler)
    return outModel


  def createTransform(self, sourcePoints, targetPoints, modelNode):

    normal = np.array([float(self.sliceLogic.GetSliceNode().GetName()==name) for name in ['Yellow','Green','Red']])
    radius = float(self.parameterNode.GetParameter("radius"))
    bpoints = []
    for i in [0, int(sourcePoints.GetNumberOfPoints()/2), sourcePoints.GetNumberOfPoints()-1]:
      # get normal direction of curve in points
      sourcePoint = np.array(sourcePoints.GetPoint(i)) 
      direction = np.array(targetPoints.GetPoint(i)) - sourcePoint
      direction = direction / np.linalg.norm(direction)
      # add control point a radius away from the point in normal direction
      bpoints.append(sourcePoint + direction * radius)
      bpoints.append(sourcePoint - direction * radius)
      # add control point normal to the plane
      bpoints.append(sourcePoint + normal * radius)
      bpoints.append(sourcePoint - normal * radius)
      # add control points a radius away from first and last point of line
      if i == 0 or i == sourcePoints.GetNumberOfPoints() - 1:
        direction = sourcePoint - np.array(sourcePoints.GetPoint(max(1,i-1)))
        direction = direction / np.linalg.norm(direction)
        bpoints.append(sourcePoint + direction * radius)

    # get indexes of control points that are a radius away from all points in line
    for i in range(sourcePoints.GetNumberOfPoints()):
      sourcePoint = np.array(sourcePoints.GetPoint(i)) 
      keepIndex = [j for j in range(len(bpoints)) if np.linalg.norm(bpoints[j] - sourcePoint) > (radius - 0.5)]
    
    bpoints = [bpoints[i] for i in np.unique(keepIndex)] # only keep these points
    
    # addmni bounds to source and target points (so as to constrain the deformation)
    for bp in (bpoints + self.MNIBounds):
      sourcePoints.InsertNextPoint(*bp)
      targetPoints.InsertNextPoint(*bp)
    # create thin plate spline transfrom
    transform = vtk.vtkThinPlateSplineTransform()
    transform.SetSourceLandmarks(sourcePoints)
    transform.SetTargetLandmarks(targetPoints)
    transform.SetBasisToR() # 3D data
    transform.Inverse()
    transformNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLGridTransformNode')
    transformNode.SetAndObserveTransformFromParent(transform)
    return transformNode

  def delay(self):
    dieTime = qt.QTime().currentTime().addMSecs(500)
    while qt.QTime().currentTime() < dieTime:
      qt.QCoreApplication.processEvents(qt.QEventLoop.AllEvents, 100)

  def displayTransformFromPoints(self, transform, points):
    # aux fiducial node
    f = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    f.GetDisplayNode().SetVisibility(0)
    for i in range(points.GetNumberOfPoints()):
      f.AddFiducialFromArray(np.array(points.GetPoint(i)))
    if not transform.GetDisplayNode():
      transform.CreateDefaultDisplayNodes()
    transform.GetDisplayNode().SetVisibility(1)
    transform.GetDisplayNode().SetSliceIntersectionVisibility(1)
    transform.GetDisplayNode().SetVisibility3D(0)
    transform.GetDisplayNode().SetAndObserveGlyphPointsNode(f)
    self.delay()
    slicer.mrmlScene.RemoveNode(f)

  def pointsToCurve(self, points, sampleDist = 2):
    curve = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsCurveNode')
    curve.GetDisplayNode().SetVisibility(0)
    curve.SetControlPointPositionsWorld(points)
    curve.ResampleCurveWorld(sampleDist)
    return curve

  def cutModelWithSliceIntersection(self, model):
    plane = vtk.vtkPlane()
    plane.SetOrigin(self.rasPoints.GetPoint(0)) # point in plane
    normal = [float(self.sliceLogic.GetSliceNode().GetName()==name) for name in ['Yellow','Green','Red']]
    plane.SetNormal(normal)

    pd = model.GetPolyData()

    cutter = vtk.vtkCutter()
    cutter.SetInputData(pd)
    cutter.SetCutFunction(plane)
    cutter.SetGenerateCutScalars(0)
    cutter.Update()

    pd2 = cutter.GetOutput()

    triangulator = vtk.vtkContourTriangulator()
    triangulator.SetInputData(pd2)
    triangulator.Update()
    pd3 = triangulator.GetOutput()

    model3 = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode')
    model3.SetAndObservePolyData(pd3)
    model3.CreateDefaultDisplayNodes()
    model3.GetDisplayNode().SetVisibility(0)

    return model3

  def cleanup(self):
    WarpEffectTool.cleanup(self)
    PointerEffect.DrawEffectTool.cleanup(self)

    
  

  