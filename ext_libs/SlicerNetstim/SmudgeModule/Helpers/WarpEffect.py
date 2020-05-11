import slicer, vtk, qt
import numpy as np
import sys, os
from scipy import ndimage

from . import PointerEffect

import TransformsUtil
import SmudgeModule





class WarpEffectTool():

  _instances = set()

  def __init__(self, warpNode = True):
    self._instances.add(self)
    self.parameterNode = SmudgeModule.SmudgeModuleLogic().getParameterNode()
    if warpNode:
      self.warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))

  def createSphere(self, r):
    # create a sphere with redius
    xx, yy, zz = np.mgrid[:2*r+1, :2*r+1, :2*r+1]
    sphereResult = (xx-r) ** 2 + (yy-r) ** 2 + (zz-r) ** 2
    sphereResult[r][r][r] = 1 # replace 0 by 1
    sphereLarge = sphereResult <= (r**2+1) # sphere that the mouse shows
    currentEffect = self.parameterNode.GetParameter("currentEffect")
    sphereSmall = sphereResult <= ((r * float(self.parameterNode.GetParameter(currentEffect + "Hardness")) / 100.0) **2 + 1 ) # hardness amount
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
    if currentEffect == "Smudge":
      sphereResult = sphereResult * float(self.parameterNode.GetParameter("SmudgeForce")) / 100.0
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

  def applyChanges(self):
    # remove redo options
    SmudgeModule.SmudgeModuleLogic().removeRedoNodes()
    # flatten when there already 3 layers
    if TransformsUtil.TransformsUtilLogic().getNumberOfLayers(self.warpNode) == 3:
      TransformsUtil.TransformsUtilLogic().flattenTransform(self.warpNode, False)
    # harden transform
    self.warpNode.HardenTransform()
    self.warpNode.InvokeEvent(slicer.vtkMRMLGridTransformNode.TransformModifiedEvent)
    # save tool name
    self.parameterNode.SetParameter("lastOperation", self.parameterNode.GetParameter("currentEffect"))
    # update gui
    self.parameterNode.SetParameter("warpModified", str(int(self.parameterNode.GetParameter("warpModified"))+1))

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

class NoneEffectTool(PointerEffect.PointerEffectTool, WarpEffectTool):

  def __init__(self, sliceWidget):
    WarpEffectTool.__init__(self, warpNode = False)
    PointerEffect.PointerEffectTool.__init__(self, sliceWidget)
    



#
# SmudgeEffectTool
#

class SmudgeEffectTool(PointerEffect.CircleEffectTool, WarpEffectTool):

  def __init__(self, sliceWidget, auxTransformNode):

    WarpEffectTool.__init__(self)
    PointerEffect.CircleEffectTool.__init__(self, sliceWidget)
        
    # transform data
    self.auxTransformNode = auxTransformNode
    self.auxTransformSpacing = self.auxTransformNode.GetTransformFromParent().GetDisplacementGrid().GetSpacing()[0] # Asume isotropic!
    self.auxTransfromRASToIJK = TransformsUtil.TransformsUtilLogic().getTransformRASToIJK(self.auxTransformNode)  

    self.previousPoint = [0,0,0]   
    self.smudging = False
    self.outOfBounds = False


  def processEvent(self, caller=None, event=None):

    PointerEffect.CircleEffectTool.processEvent(self, caller, event)

    if event == 'LeftButtonPressEvent':
      self.smudging = True
      self.outOfBounds = False
      self.warpNode.SetAndObserveTransformNodeID(self.auxTransformNode.GetID())    
      self.auxTransformArray = slicer.util.array(self.auxTransformNode.GetID())
      xy = self.interactor.GetEventPosition()
      xyToRAS = self.sliceLogic.GetSliceNode().GetXYToRAS()
      self.previousPoint = xyToRAS.MultiplyDoublePoint( (xy[0], xy[1], 0, 1) )[0:3]
    elif event == 'LeftButtonReleaseEvent' and not self.outOfBounds:
      self.smudging = False
      # smooth
      if int(self.parameterNode.GetParameter("SmudgePostSmoothing")):
        sigma = float(self.parameterNode.GetParameter("SmudgeSigma")) / 100.0 * float(self.parameterNode.GetParameter("SmudgeRadius")) / self.auxTransformSpacing
        self.auxTransformArray[:] = np.stack([ndimage.gaussian_filter(self.auxTransformArray[:,:,:,i], sigma) for i in range(3)], 3).squeeze()
      # apply
      self.applyChanges()
      self.auxTransformArray[:] = np.zeros(self.auxTransformArray.shape)
      qt.QApplication.setOverrideCursor(qt.QCursor(qt.Qt.ArrowCursor))


    elif event == 'MouseMoveEvent':
      if self.smudging:

        r = int(round(float(self.parameterNode.GetParameter("SmudgeRadius")) / self.auxTransformSpacing))
        sphereResult = self.createSphere(r)
        currentPoint = self.eventPositionToRAS()
        currentIndex = self.getCurrentIndex(r, currentPoint, self.auxTransfromRASToIJK)

        # apply to transform array
        try:
          self.auxTransformArray[currentIndex] += np.stack([(sphereResult) * i for i in (np.array(self.previousPoint) - np.array(currentPoint))],3) # original
        except ValueError:
          qt.QMessageBox.warning(qt.QWidget(), '', 'Out of bounds. Try expanding the grid.')
          self.smudging = False
          self.outOfBounds = True
          self.cursorOn()

        # update view
        self.auxTransformNode.Modified()
        # update previous point
        self.previousPoint = currentPoint

  def cleanup(self):
    slicer.mrmlScene.RemoveNode(self.auxTransformNode)
    WarpEffectTool.cleanup(self)
    PointerEffect.CircleEffectTool.cleanup(self)



#
# SmoothEffectTool
#

class SmoothEffectTool(PointerEffect.CircleEffectTool, WarpEffectTool):

  def __init__(self, sliceWidget):

    WarpEffectTool.__init__(self)
    PointerEffect.CircleEffectTool.__init__(self, sliceWidget)
    
    self.transformArray = slicer.util.array(self.warpNode.GetID())
    self.warpRASToIJK = TransformsUtil.TransformsUtilLogic().getTransformRASToIJK(self.warpNode)
    self.warpSpacing = TransformsUtil.TransformsUtilLogic().getGridDefinition(self.warpNode)[2][0]

    self.smoothContent = []
    self.currentIndex = []
    self.preview = False
    
  def processEvent(self, caller=None, event=None):

    PointerEffect.CircleEffectTool.processEvent(self, caller, event)

    if self.preview and event in ['LeftButtonDoubleClickEvent','LeftButtonReleaseEvent']:
      # undo pressed opperation. sometimes double click is called before release and viceversa
      self.transformArray[self.currentIndex] -= self.smoothContent
      self.preview = False

    if event =='LeftButtonDoubleClickEvent':
      self.calculateSmoothContent()
      self.transformArray[self.currentIndex] += self.smoothContent
      self.applyChanges()
    elif event == 'LeftButtonReleaseEvent':
      self.warpNode.InvokeEvent(slicer.vtkMRMLGridTransformNode.TransformModifiedEvent)
      qt.QApplication.setOverrideCursor(qt.QCursor(qt.Qt.ArrowCursor))
    elif event == 'LeftButtonPressEvent':
      self.preview = True
      self.calculateSmoothContent()
      self.transformArray[self.currentIndex] += self.smoothContent
      self.warpNode.InvokeEvent(slicer.vtkMRMLGridTransformNode.TransformModifiedEvent)
      

  def calculateSmoothContent(self):
    sigma = float(self.parameterNode.GetParameter("SmoothSigma")) / self.warpSpacing
    r = int(round(float(self.parameterNode.GetParameter("SmoothRadius")) / self.warpSpacing))
    if int(self.parameterNode.GetParameter("SmoothUseRadius")):
      # get shpere and index
      sphereResult = self.createSphere(r)
      currentPoint = self.eventPositionToRAS()
      self.currentIndex = self.getCurrentIndex(r, currentPoint, self.warpRASToIJK)   
      # gaussian filter for each component 
      self.smoothContent =  np.stack([ndimage.gaussian_filter(self.transformArray[self.currentIndex + (slice(i,i+1),)], sigma) for i in range(3)], 3).squeeze()
      # substract original
      self.smoothContent = self.smoothContent - self.transformArray[self.currentIndex]
      # modulate result with the sphere
      self.smoothContent = np.stack([self.smoothContent[:,:,:,i] * sphereResult for i in range(3)], 3).squeeze()
    else: # maximum radius: take all warp field
      self.smoothContent =  np.stack([ndimage.gaussian_filter(self.transformArray[:,:,:,i], sigma) for i in range(3)], 3).squeeze()
      self.smoothContent = self.smoothContent - self.transformArray
      self.currentIndex = tuple([slice(0,s) for s in self.smoothContent.shape])

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
    
  def processEvent(self, caller=None, event=None):

    PointerEffect.DrawEffectTool.processEvent(self, caller, event) 

    if event == 'LeftButtonReleaseEvent':

      sampleDistance = float(self.parameterNode.GetParameter("DrawSampleDistance"))
      
      # create curve from drawing and resample
      sourceCurve = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsCurveNode')
      sourceCurve.SetControlPointPositionsWorld(self.rasPoints)
      sourceCurve.ResampleCurveWorld(sampleDistance)  

      # overwrite points with the resampled ones
      sourceCurve.GetControlPointPositionsWorld(self.rasPoints)

      # if only one point left exit
      if self.rasPoints.GetNumberOfPoints() <= 1:
        self.resetPolyData()
        slicer.mrmlScene.RemoveNode(sourceCurve)
        qt.QApplication.setOverrideCursor(qt.QCursor(qt.Qt.ArrowCursor))
        return

      # get closest model sliced
      slicedModel, originalModel = self.sliceClosestModel(self.rasPoints.GetPoint(0))

      # resample sourceCurve in sliced model with same amount of points
      targetCurve = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsCurveNode')
      targetCurve.GetDisplayNode().SetVisibility(0)
      targetCurve.SetControlPointPositionsWorld(self.rasPoints)
      targetCurve.SetCurveTypeToShortestDistanceOnSurface(slicedModel)
      targetCurve.ResampleCurveSurface(sampleDistance, slicer.vtkMRMLModelNode().SafeDownCast(slicedModel), 0.0025)
      targetCurve.ResampleCurveWorld(targetCurve.GetCurveLengthWorld() / max((sourceCurve.GetNumberOfControlPoints() - 1), 1))
      
      # curve to fiducial
      sourceFiducial = self.curveToFiducial(sourceCurve)
      targetFiducial = self.curveToFiducial(targetCurve)

      # compute warp
      landmarkWarp = self.computeWarp(sourceFiducial, targetFiducial)

      # visualize
      self.resetPolyData() # delete manual drawing
      sourceCurve.GetDisplayNode().SetSelectedColor(1,1,0)
      self.displayTransformFromFiducial(landmarkWarp, sourceFiducial)
      
      # apply
      self.warpNode.SetAndObserveTransformNodeID(landmarkWarp.GetID())
      self.applyChanges()

      # save target as fixed points
      self.addFiducialToHierarchy(targetFiducial, originalModel)

      # remove nodes
      slicer.mrmlScene.RemoveNode(sourceCurve)
      slicer.mrmlScene.RemoveNode(sourceFiducial)
      slicer.mrmlScene.RemoveNode(targetCurve)
      slicer.mrmlScene.RemoveNode(slicedModel)
      slicer.mrmlScene.RemoveNode(landmarkWarp)

      qt.QApplication.setOverrideCursor(qt.QCursor(qt.Qt.ArrowCursor))




  def addFiducialToHierarchy(self, fiducial, modelNode):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    # add folder with nearest model name
    modelParentName =  shNode.GetItemName(shNode.GetItemParent(shNode.GetItemByDataNode(modelNode)))
    fiducial.SetName(slicer.mrmlScene.GenerateUniqueName(modelParentName + '_' + modelNode.GetName()))
    fiducial.SetLocked(1)
    fiducial.GetDisplayNode().SetGlyphScale(2)
    fiducial.GetDisplayNode().SetVisibility(1)
    shNode.SetItemAttribute(shNode.GetItemByDataNode(fiducial), 'drawing', '1')
    for i in range(fiducial.GetNumberOfControlPoints()):
      fiducial.SetNthControlPointLabel(i,'')

  def delay(self):
    dieTime = qt.QTime().currentTime().addMSecs(500)
    while qt.QTime().currentTime() < dieTime:
      qt.QCoreApplication.processEvents(qt.QEventLoop.AllEvents, 100)

  def displayTransformFromFiducial(self, transform, fiducial):
    if not transform.GetDisplayNode():
      transform.CreateDefaultDisplayNodes()
    transform.GetDisplayNode().SetVisibility(1)
    transform.GetDisplayNode().SetVisibility2D(1)
    transform.GetDisplayNode().SetVisibility3D(0)
    transform.GetDisplayNode().SetAndObserveGlyphPointsNode(fiducial)
    self.delay()

  def getFixedPoints(self):
    points = vtk.vtkPoints()
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    nMarkups = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLMarkupsFiducialNode')
    for i in range(nMarkups):
      markupNode = slicer.mrmlScene.GetNthNodeByClass(i, 'vtkMRMLMarkupsFiducialNode')
      if 'drawing' in shNode.GetItemAttributeNames(shNode.GetItemByDataNode(markupNode)) and int(shNode.GetItemAttribute(shNode.GetItemByDataNode(markupNode),'drawing')):
        p = vtk.vtkPoints()
        markupNode.GetControlPointPositionsWorld(p)
        points.InsertPoints(points.GetNumberOfPoints(), p.GetNumberOfPoints(), 0, p)
    return points


  def computeWarp(self, sourceDrawing, targetDrawing):
    # points
    sourcePoints = vtk.vtkPoints()
    targetPoints = vtk.vtkPoints()
    # add drawing
    sourceDrawing.GetControlPointPositionsWorld(sourcePoints)
    targetDrawing.GetControlPointPositionsWorld(targetPoints)
    # add fixed points
    fixedPoints = self.getFixedPoints()
    sourcePoints.InsertPoints(sourcePoints.GetNumberOfPoints(), fixedPoints.GetNumberOfPoints(), 0, fixedPoints)
    targetPoints.InsertPoints(targetPoints.GetNumberOfPoints(), fixedPoints.GetNumberOfPoints(), 0, fixedPoints) 

    # source and target fiducials
    sourceFiducial = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    targetFiducial = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    sourceFiducial.SetControlPointPositionsWorld(sourcePoints)
    targetFiducial.SetControlPointPositionsWorld(targetPoints)

    # generate aux nodes
    size,origin,spacing = TransformsUtil.TransformsUtilLogic().getGridDefinition(self.warpNode)
    auxVolumeNode = TransformsUtil.TransformsUtilLogic().createEmpyVolume(size,origin,spacing)
    outWarp = TransformsUtil.TransformsUtilLogic().emptyGridTransform(size,origin,spacing)

    parameters = {}
    parameters["plmslc_landwarp_fixed_volume"] = auxVolumeNode.GetID()
    parameters["plmslc_landwarp_moving_volume"] = auxVolumeNode.GetID()
    parameters["plmslc_landwarp_fixed_fiducials"] = targetFiducial.GetID()
    parameters["plmslc_landwarp_moving_fiducials"] = sourceFiducial.GetID()
    parameters["plmslc_landwarp_output_vf"] = outWarp.GetID()
    parameters["plmslc_landwarp_rbf_type"] = "gauss"
    parameters["plmslc_landwarp_rbf_radius"] = float(self.parameterNode.GetParameter("DrawSpread"))
    parameters["plmslc_landwarp_stiffness"] = float(self.parameterNode.GetParameter("DrawStiffness"))

    cli = slicer.cli.run(slicer.modules.plastimatch_slicer_landwarp, None, parameters, wait_for_completion=True, update_display=False)

    slicer.mrmlScene.RemoveNode(auxVolumeNode)
    slicer.mrmlScene.RemoveNode(sourceFiducial)
    slicer.mrmlScene.RemoveNode(targetFiducial)

    return outWarp

   
  def curveToFiducial(self, curve):
    fiducial = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    fiducial.GetDisplayNode().SetVisibility(0)
    points = vtk.vtkPoints()
    curve.GetControlPointPositionsWorld(points)
    fiducial.SetControlPointPositionsWorld(points)
    return fiducial

  def sliceClosestModel(self, point):
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
      if model.GetDisplayNode() and model.GetDisplayNode().GetVisibility() and polyData.GetNumberOfCells() > 1: # model visible and cells available
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
    WarpEffectTool.cleanup(self)
    PointerEffect.DrawEffectTool.cleanup(self)
