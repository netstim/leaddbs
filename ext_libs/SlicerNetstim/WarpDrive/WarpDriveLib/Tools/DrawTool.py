import qt, vtk, slicer
import numpy as np


from ..Widgets.ToolWidget   import AbstractToolWidget
from ..Effects.DrawEffect import AbstractDrawEffect

from ..Helpers import GridNodeHelper

class DrawToolWidget(AbstractToolWidget):
  
  def __init__(self):
    toolTip = ''
    AbstractToolWidget.__init__(self, 'Draw', toolTip)

    # set up menu
    nearestModelAction = qt.QAction('To Nearest Model', self.effectButton)
    nearestModelAction.setCheckable(True)
    nearestModelAction.setChecked(True)
    fromNearestModelAction = qt.QAction('From Nearest Model', self.effectButton)
    fromNearestModelAction.setCheckable(True)
    twoLinesAction = qt.QAction('To Following Line', self.effectButton)
    twoLinesAction.setCheckable(True)
    actionsGroup = qt.QActionGroup(self.effectButton)
    actionsGroup.addAction(nearestModelAction)
    actionsGroup.addAction(fromNearestModelAction)
    actionsGroup.addAction(twoLinesAction)
    menu = qt.QMenu(self.effectButton)
    menu.addActions(actionsGroup.actions())
    self.effectButton.setMenu(menu)
    self.effectButton.setPopupMode(self.effectButton.DelayedPopup)


class DrawToolEffect(AbstractDrawEffect):


  def __init__(self, sliceWidget):
    AbstractDrawEffect.__init__(self, sliceWidget)
    self.sourceFiducial = None


  def processEvent(self, caller=None, event=None):

    AbstractDrawEffect.processEvent(self, caller, event) 

    if event == 'LeftButtonReleaseEvent':

      if self.sourceFiducial is None: # no previous drawing

        self.sourceFiducial = self.getFiducialFromDrawing()

        if self.sourceFiducial is None:
          self.resetDrawing()
          return

        if self.parameterNode.GetParameter("DrawMode") == 'To Nearest Model': # get target fiducial from nearest model 
          targetFiducial = self.getFiducialFromSlicedModel()
        elif self.parameterNode.GetParameter("DrawMode") == 'From Nearest Model': # get target fiducial from nearest model 
          targetFiducial = self.sourceFiducial
          self.sourceFiducial = self.getFiducialFromSlicedModel()
        elif self.parameterNode.GetParameter("DrawMode") == 'To Following Line': # return and wait for following drawing
          self.sourceFiducial.GetDisplayNode().SetVisibility(1)
          self.resetDrawing()
          return

      else: # use new drawing as target fiducial
        targetFiducial = self.getFiducialFromDrawing(nPoints = self.sourceFiducial.GetNumberOfControlPoints())
        targetFiducial.SetName(slicer.mrmlScene.GenerateUniqueName('drawing'))

      if targetFiducial is None:
        slicer.mrmlScene.RemoveNode(self.sourceFiducial)
        self.sourceFiducial = None
        self.resetDrawing()
        return

      self.sourceFiducial.ApplyTransform(self.parameterNode.GetNodeReference("OutputGridTransform").GetTransformFromParent()) # undo current

      self.setFiducialNodeAs("Source", self.sourceFiducial, targetFiducial.GetName(), self.parameterNode.GetParameter("Radius"))
      self.setFiducialNodeAs("Target", targetFiducial, targetFiducial.GetName(), self.parameterNode.GetParameter("Radius"))

      self.parameterNode.SetParameter("Update","true")
      self.sourceFiducial = None
      self.resetDrawing()

    elif event == 'RightButtonPressEvent' or (event == 'KeyPressEvent' and self.interactor.GetKeySym()=='Escape'):
      slicer.mrmlScene.RemoveNode(self.sourceFiducial)
      self.sourceFiducial = None

  def getFiducialFromDrawing(self, sampleDistance = 1, nPoints = None):

    # create curve from drawing and resample
    sourceCurve = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsCurveNode')
    sourceCurve.Copy(self.drawnCurveNode)
    sourceCurve.ResampleCurveWorld(sampleDistance)
    if nPoints is not None: # resample to get specified number of points
      sourceCurve.ResampleCurveWorld(sourceCurve.GetCurveLengthWorld() / max((nPoints - 1), 1))

    sourceFiducial = self.curveToFiducial(sourceCurve)

    slicer.mrmlScene.RemoveNode(sourceCurve)

    if sourceFiducial.GetNumberOfControlPoints() <= 1:
      slicer.mrmlScene.RemoveNode(sourceFiducial)
      return None
    else:
      return sourceFiducial

  def getFiducialFromSlicedModel(self, sampleDistance = 1):

    # get source fiducial points
    resampledPoints = vtk.vtkPoints()
    self.sourceFiducial.GetControlPointPositionsWorld(resampledPoints)
    # get closest model sliced
    slicedModel, originalModel = self.sliceClosestModel(resampledPoints.GetPoint(0))

    if not slicedModel:
      return None

    # resample sourceCurve in sliced model with same amount of points
    targetCurve = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsCurveNode')
    targetCurve.GetDisplayNode().SetVisibility(0)
    targetCurve.SetControlPointPositionsWorld(resampledPoints)
    targetCurve.SetCurveTypeToShortestDistanceOnSurface(slicedModel)
    targetCurve.ResampleCurveWorld(targetCurve.GetCurveLengthWorld() / max((resampledPoints.GetNumberOfPoints() - 1), 1))
      
    # curve to fiducial
    targetFiducial = self.curveToFiducial(targetCurve)

    # set name
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    modelParentName =  shNode.GetItemName(shNode.GetItemParent(shNode.GetItemByDataNode(originalModel)))
    targetFiducial.SetName(modelParentName + '_' + originalModel.GetName())

    # remove
    slicer.mrmlScene.RemoveNode(targetCurve)
    slicer.mrmlScene.RemoveNode(slicedModel)

    return targetFiducial


  def copyControlPoints(self, sourceNode, targetNode):
    if targetNode.GetNumberOfControlPoints() == 0:
      label = '1'
    else:
      label = str( int(targetNode.GetNthFiducialLabel(targetNode.GetNumberOfControlPoints()-1)) + 1 )
    p = np.zeros(3)
    for i in range(sourceNode.GetNumberOfControlPoints()):
      sourceNode.GetNthControlPointPosition(i,p)
      targetNode.AddFiducialFromArray(p, label)


  def curveToFiducial(self, curve):
    fiducial = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    fiducial.GetDisplayNode().SetGlyphTypeFromString('Sphere3D')
    fiducial.GetDisplayNode().SetGlyphScale(1)
    fiducial.GetDisplayNode().SetVisibility(0)
    fiducial.GetDisplayNode().SetPointLabelsVisibility(0)
    points = vtk.vtkPoints()
    curve.GetControlPointPositionsWorld(points)
    fiducial.SetControlPointPositionsWorld(points)
    return fiducial

  def sliceClosestModel(self, point):
    shNode = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)
    originalModel = None
    # set up plane
    sliceToRAS = self.sliceLogic.GetSliceNode().GetSliceToRAS()
    normal = np.array([sliceToRAS.GetElement(0,2), sliceToRAS.GetElement(1,2), sliceToRAS.GetElement(2,2)])
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
      if model.GetTransformNodeID():
        try:
          clonedItemID = slicer.modules.subjecthierarchy.logic().CloneSubjectHierarchyItem(shNode, shNode.GetItemByDataNode(model))
          tmpNode = shNode.GetItemDataNode(clonedItemID)
          tmpNode.SetAndObserveTransformNodeID(model.GetTransformNodeID())
        except:
          continue
        tmpNode.HardenTransform()
        polyData = tmpNode.GetPolyData()
      else:
        polyData = model.GetPolyData()
        tmpNode = None
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
      if tmpNode:
        slicer.mrmlScene.RemoveNode(tmpNode)
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
    slicer.mrmlScene.RemoveNode(self.sourceFiducial)
    self.sourceFiducial = None
    AbstractDrawEffect.cleanup(self)

