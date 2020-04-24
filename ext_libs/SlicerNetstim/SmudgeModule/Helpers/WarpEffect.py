import slicer, vtk, qt
import numpy as np
import sys, os
from scipy import ndimage

from . import PointerEffect

import TransformsUtil
import SmudgeModule





class WarpEffectTool():

  _instances = set()

  def __init__(self):
    self._instances.add(self)
    self.parameterNode = SmudgeModule.SmudgeModuleLogic().getParameterNode()

  def createSphere(self, r):
    # create a sphere with redius
    xx, yy, zz = np.mgrid[:2*r+1, :2*r+1, :2*r+1]
    sphereResult = (xx-r) ** 2 + (yy-r) ** 2 + (zz-r) ** 2
    sphereResult[r][r][r] = 1 # replace 0 by 1
    sphereLarge = sphereResult <= (r**2+1) # sphere that the mouse shows
    sphereSmall = sphereResult <= ((r * float(self.parameterNode.GetParameter("hardness")) / 100.0) **2 + 1 ) # hardness amount
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
    return sphereResult

  def eventPositionToRAS(self):
    xy = self.interactor.GetEventPosition()
    xyToRAS = self.sliceLogic.GetSliceNode().GetXYToRAS()
    currentPoint = xyToRAS.MultiplyDoublePoint( (xy[0], xy[1], 0, 1) )[0:3]
    return currentPoint

  def getCurrentIndex(self, r, currentPoint, RASToIJK):
    # get current IJK
    pos_i,pos_j,pos_k,aux = RASToIJK.MultiplyDoublePoint(currentPoint + (1,))
    k,j,i = int(round(pos_k)),int(round(pos_j)),int(round(pos_i))
    # expand with radius
    currentIndex = slice(k-r,k+r+1), slice(j-r,j+r+1), slice(i-r,i+r+1)
    return currentIndex

  def applyChanges(self, transformNode):
    # apply to anchor points
    SmudgeModule.SmudgeModuleLogic().applyChangesToAnchorPoints(transformNode.GetTransformToParent())
    # get warp and harden transform
    warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))
    warpNode.HardenTransform()
    warpNode.InvokeEvent(slicer.vtkMRMLGridTransformNode.TransformModifiedEvent)
    # remove redo options
    SmudgeModule.SmudgeModuleLogic().removeRedoTransform()
    # save tool name
    self.parameterNode.SetParameter("lastOperation",self.toolName)


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

  def __init__(self, sliceWidget):
    WarpEffectTool.__init__(self)
    PointerEffect.PointerEffectTool.__init__(self, sliceWidget)
    



#
# SmudgeEffectTool
#

class SmudgeEffectTool(PointerEffect.CircleEffectTool, WarpEffectTool):

  def __init__(self, sliceWidget, auxTransformNode):

    WarpEffectTool.__init__(self)
    PointerEffect.CircleEffectTool.__init__(self, sliceWidget)
    
    self.toolName = 'smudge'

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
      # smooth
      sigma = float(self.parameterNode.GetParameter("sigma")) / self.auxTransformSpacing
      self.auxTransformArray[:] = np.stack([ndimage.gaussian_filter(self.auxTransformArray[:,:,:,i], sigma) for i in range(3)], 3).squeeze()
      # apply
      self.applyChanges(self.auxTransformNode)
      self.auxTransformArray[:] = np.zeros(self.auxTransformArray.shape)


    elif event == 'MouseMoveEvent':
      if self.smudging:

        r = int(round(float(self.parameterNode.GetParameter("radius")) / self.auxTransformSpacing))
        sphereResult = self.createSphere(r)
        currentPoint = self.eventPositionToRAS()
        currentIndex = self.getCurrentIndex(r, currentPoint, self.auxTransfromRASToIJK)

        # apply to transform array
        self.auxTransformArray[currentIndex] += np.stack([(sphereResult) * i for i in (np.array(self.previousPoint) - np.array(currentPoint))],3) # original


        # update view
        self.auxTransformNode.Modified()
        # update previous point
        self.previousPoint = currentPoint

  def cleanup(self):
    slicer.mrmlScene.RemoveNode(self.auxTransformNode)
    WarpEffectTool.cleanup(self)
    PointerEffect.CircleEffectTool.cleanup(self)



#
# BlurEffectTool
#

class BlurEffectTool(PointerEffect.CircleEffectTool, WarpEffectTool):

  def __init__(self, sliceWidget):

    WarpEffectTool.__init__(self)
    PointerEffect.CircleEffectTool.__init__(self, sliceWidget)
    
    self.toolName = 'blur'

    self.warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))
    self.warpRASToIJK = TransformsUtil.TransformsUtilLogic().getTransformRASToIJK(self.warpNode)
    self.warpSpacing = self.warpNode.GetTransformFromParent().GetDisplacementGrid().GetSpacing()[0] # Asume isotropic!

    self.blurContent = []
    self.currentIndex = []
    self.preview = False
    
  def processEvent(self, caller=None, event=None):

    PointerEffect.CircleEffectTool.processEvent(self, caller, event)

    if event =='LeftButtonDoubleClickEvent':
      self.preview = False
      self.transformArray[self.currentIndex] += self.blurContent
      # create aux transform to apply to anchor points
      size,origin,spacing = TransformsUtil.TransformsUtilLogic().getGridDefinition(self.warpNode)
      auxTransformNode = TransformsUtil.TransformsUtilLogic().emptyGridTransfrom(size,origin,spacing)
      auxTransformArray = slicer.util.array(auxTransformNode.GetID())
      auxTransformArray[self.currentIndex] = self.blurContent
      auxTransformNode.Modified()
      # apply and remove
      self.applyChanges(auxTransformNode)
      slicer.mrmlScene.RemoveNode(auxTransformNode)
      
    elif event == 'LeftButtonReleaseEvent':
      if self.preview:
        self.transformArray[self.currentIndex] -= self.blurContent
        self.warpNode.InvokeEvent(slicer.vtkMRMLGridTransformNode.TransformModifiedEvent)

    elif event == 'LeftButtonPressEvent':
      self.preview = True
      # get array 
      self.transformArray = slicer.util.array(self.warpNode.GetID())

      sigma = float(self.parameterNode.GetParameter("sigma")) / self.warpSpacing
      r = int(round(float(self.parameterNode.GetParameter("radius")) / self.warpSpacing))

      if r != int(round(float(self.parameterNode.GetParameter("maxRadius")) / self.warpSpacing)):
        # get shpere and index
        sphereResult = self.createSphere(r)
        currentPoint = self.eventPositionToRAS()
        self.currentIndex = self.getCurrentIndex(r, currentPoint, self.warpRASToIJK)   
        # gaussian filter for each component 
        self.blurContent =  np.stack([ndimage.gaussian_filter(self.transformArray[self.currentIndex + (slice(i,i+1),)], sigma) for i in range(3)], 3).squeeze()
        # substract original
        self.blurContent = self.blurContent - self.transformArray[self.currentIndex]
        # modulate result with the sphere
        self.blurContent = np.stack([self.blurContent[:,:,:,i] * sphereResult for i in range(3)], 3).squeeze()
      else: # maximum radius: take all warp field
        self.blurContent =  np.stack([ndimage.gaussian_filter(self.transformArray[:,:,:,i], sigma) for i in range(3)], 3).squeeze()
        self.blurContent = self.blurContent - self.transformArray
        self.currentIndex = tuple([slice(0,s) for s in self.blurContent.shape])
      # apply
      self.transformArray[self.currentIndex] += self.blurContent
      self.warpNode.InvokeEvent(slicer.vtkMRMLGridTransformNode.TransformModifiedEvent)
      

  def cleanup(self):
    WarpEffectTool.cleanup(self)
    PointerEffect.CircleEffectTool.cleanup(self)



#
# Snap Effect Tool
#

class SnapEffectTool(PointerEffect.DrawEffectTool, WarpEffectTool):

    
  def __init__(self, sliceWidget):

    WarpEffectTool.__init__(self)
    PointerEffect.DrawEffectTool.__init__(self,sliceWidget)

    self.toolName = 'snap'
    
    self.warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))
    size,origin,spacing = TransformsUtil.TransformsUtilLogic().getGridDefinition(self.warpNode)
    self.warpBounds = [[origin[0]+spacing[0]*size[0]*i,  origin[1]+spacing[1]*size[1]*j, origin[2]+spacing[2]*size[2]*k] for i in range(2) for j in range(2) for k in range(2)]
    
  def processEvent(self, caller=None, event=None):

    if event == 'LeftButtonReleaseEvent':
      self.actionState = None # stop drawing

      # create open curve from drawing
      curve1 = self.pointsToCurve(self.rasPoints)
      points1 = vtk.vtkPoints()
      curve1.GetControlPointPositionsWorld(points1)

      # get closest model to points
      modelNode = self.getClosestModel(points1)
      if (not modelNode) or (points1.GetNumberOfPoints() == 1): # clean and continue
        slicer.mrmlScene.RemoveNode(curve1)
      else:
        # cut model with current plane
        slicedModel = self.cutModelWithSliceIntersection(modelNode)
        # get curve on sliced model
        points2 = self.samplePointsInModel(points1, slicedModel)
        curve2 = self.pointsToCurve(points2)
        # in case line is inside the model the curve will be the same. set curve type to shortest distance to surface and resample again
        curve2.SetCurveTypeToShortestDistanceOnSurface(slicedModel)
        curve2.ResampleCurveSurface(2, slicer.vtkMRMLModelNode().SafeDownCast(slicedModel), 0.0025)
        # get same number of points as other curve
        curve2.ResampleCurveWorld(curve2.GetCurveLengthWorld() / max((curve1.GetNumberOfControlPoints() - 1), 1) )
        curve2.GetControlPointPositionsWorld(points2)
        # add drawings to hierarchy (before more points are added) TODO: better implement this
        self.addDrawingsToHierarchy(points1,points2,modelNode) 
        # calculate transform and apply
        transformNode = self.createTransform(points1, points2)
        self.displayTransformFromPoints(transformNode, points1)
        self.warpNode.SetAndObserveTransformNodeID(transformNode.GetID())
        self.applyChanges(transformNode)
        # cleanup
        slicer.mrmlScene.RemoveNode(curve1)
        slicer.mrmlScene.RemoveNode(curve2)
        slicer.mrmlScene.RemoveNode(slicedModel)
        slicer.mrmlScene.RemoveNode(transformNode)
        
    PointerEffect.DrawEffectTool.processEvent(self, caller, event) 
    

  def addDrawingsToHierarchy(self, points1, points2, modelNode):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    # add folder with nearest model name
    drawingsRootItem = int(self.parameterNode.GetParameter("drawingsRootItem"))
    modelParentName =  shNode.GetItemName(shNode.GetItemParent(shNode.GetItemByDataNode(modelNode)))
    subFolderID = shNode.CreateFolderItem(drawingsRootItem, modelParentName + '_' + modelNode.GetName())
    shNode.SetItemAttribute(subFolderID, 'drawing', '1')
    # create curve from points and add to folder
    for points,name in zip([points1, points2],['source','target']):
      auxCurve = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsCurveNode')
      auxCurve.SetControlPointPositionsWorld(points)
      auxCurve.SetName(name)
      auxCurve.SetLocked(1)
      auxCurve.GetDisplayNode().SetVisibility(name == 'source') # show source points
      auxCurve.GetDisplayNode().SetSelectedColor([1,0.5,0.5] if name=='source' else [0.8,0.6,0.5])
      shNode.SetItemParent(shNode.GetItemByDataNode(auxCurve), subFolderID)
      shNode.SetItemAttribute(shNode.GetItemByDataNode(auxCurve), 'drawing', '1')

  def samplePointsInModel(self, points, model, sampleDist = 2, maximumSearchRadius = 0.0025):
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


  def getControlPoints(self, sourcePoints, targetPoints):
    """
    Get points a radius away from the curve in normal direction
    """
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

    return bpoints

  def getAnchorPoints(self):
    # get drawings
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    drawingsRootItem = int(self.parameterNode.GetParameter("drawingsRootItem"))
    childrenIDs = vtk.vtkIdList()
    shNode.GetItemChildren(drawingsRootItem, childrenIDs, True)
    # init
    p=vtk.vtkPoints()
    sourcePoints=vtk.vtkPoints()
    targetPoints=vtk.vtkPoints()
    # iterete over source and target to get points
    for i in range(childrenIDs.GetNumberOfIds()):
      dataNode = shNode.GetItemDataNode(childrenIDs.GetId(i))
      if dataNode and isinstance(dataNode, slicer.vtkMRMLMarkupsCurveNode):
        dataNode.GetControlPointPositionsWorld(p)
        if dataNode.GetName() == 'source':
          sourcePoints.InsertPoints(sourcePoints.GetNumberOfPoints(),p.GetNumberOfPoints(),0,p)
        elif dataNode.GetName() == 'target':
          targetPoints.InsertPoints(targetPoints.GetNumberOfPoints(),p.GetNumberOfPoints(),0,p)
    
    return sourcePoints,targetPoints

  def createTransform(self, sourcePoints, targetPoints):
    # get control points
    bpoints = self.getControlPoints(sourcePoints, targetPoints)    
    # add warp bounds to source and target points (so as to constrain the deformation)
    for bp in (bpoints + self.warpBounds):
      sourcePoints.InsertNextPoint(*bp)
      targetPoints.InsertNextPoint(*bp)
    # add anchor points
    sourceAnchorPoints, targetAnchorPoints = self.getAnchorPoints()
    sourcePoints.InsertPoints(sourcePoints.GetNumberOfPoints(), sourceAnchorPoints.GetNumberOfPoints(), 0, sourceAnchorPoints)
    targetPoints.InsertPoints(targetPoints.GetNumberOfPoints(), targetAnchorPoints.GetNumberOfPoints(), 0, targetAnchorPoints)
    #create thin plate spline transfrom
    transform = vtk.vtkThinPlateSplineTransform()
    transform.SetSourceLandmarks(sourcePoints)
    transform.SetTargetLandmarks(targetPoints)
    transform.SetBasisToR() # 3D data
    transform.Inverse()
    transformNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLGridTransformNode')
    transformNode.SetAndObserveTransformFromParent(transform)
    # to grid transform
    outNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
    size,origin,spacing = TransformsUtil.TransformsUtilLogic().getGridDefinition(self.warpNode)
    referenceVolume = TransformsUtil.TransformsUtilLogic().createEmpyVolume(size,origin,spacing) # aux reference volume with specified resolution 
    # apply
    transformsLogic = slicer.modules.transforms.logic()
    transformsLogic.ConvertToGridTransform(transformNode, referenceVolume, outNode)
    transformNode.SetAndObserveTransformFromParent(outNode.GetTransformFromParent())
    # remove aux
    slicer.mrmlScene.RemoveNode(outNode)
    slicer.mrmlScene.RemoveNode(referenceVolume)
    # smooth
    sigma = float(self.parameterNode.GetParameter("sigma")) / spacing[0]
    a = slicer.util.array(transformNode.GetID())
    a[:] = np.stack([ndimage.gaussian_filter(a[:,:,:,i], sigma) for i in range(3)], 3).squeeze()
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
