import qt, vtk, slicer
import numpy as np


from ..Widgets.ToolWidget   import AbstractToolWidget
from ..Effects.CircleEffect import AbstractCircleEffect

from ..Helpers import GridNodeHelper

class SmudgeToolWidget(AbstractToolWidget):
    
  def __init__(self):
    toolTip = 'Click and drag to deform the image'
    AbstractToolWidget.__init__(self, 'Smudge', toolTip)


class SmudgeToolEffect(AbstractCircleEffect):

  auxTransformNode = None
  auxTransfromRASToIJK = None

  def __init__(self, sliceWidget):
    AbstractCircleEffect.__init__(self, sliceWidget)
        
    # transform data
    if not type(self).auxTransformNode:
      size, origin, spacing, directionMatrix = GridNodeHelper.getGridDefinition(self.parameterNode.GetNodeReference("InputNode"))
      userSpacing = np.ones(3) * float(self.parameterNode.GetParameter("Spacing"))
      size = size * (spacing / userSpacing)
      type(self).auxTransformNode = GridNodeHelper.emptyGridTransform(size.astype(int), origin, userSpacing, directionMatrix) 
      type(self).auxTransfromRASToIJK = GridNodeHelper.getTransformRASToIJK(self.auxTransformNode)  

    # points
    self.interactionPoints = vtk.vtkPoints()

    self.previousPoint = np.zeros(3)
    self.smudging = False


  def processEvent(self, caller=None, event=None):

    AbstractCircleEffect.processEvent(self, caller, event)

    if event == 'LeftButtonPressEvent':
      self.parameterNode.GetNodeReference("OutputGridTransform").SetAndObserveTransformNodeID(self.auxTransformNode.GetID())    
      self.auxTransformArray = slicer.util.array(self.auxTransformNode.GetID())
      self.previousPoint = self.xyToRAS(self.interactor.GetEventPosition())
      self.interactionPoints.InsertNextPoint(self.previousPoint)
      self.smudging = True

    elif event == 'LeftButtonReleaseEvent' and self.smudging:
      self.smudging = False
      # resample
      self.resamplePoints()
      # get source and target
      sourceFiducial, targetFiducial = self.getSourceTargetFromPoints()
      # apply
      self.setFiducialNodeAs("Source", sourceFiducial, targetFiducial.GetName(), self.parameterNode.GetParameter("Radius"))
      self.setFiducialNodeAs("Target", targetFiducial, targetFiducial.GetName(), self.parameterNode.GetParameter("Radius"))
      self.parameterNode.SetParameter("Update","true")
      # reset
      self.parameterNode.GetNodeReference("OutputGridTransform").HardenTransform()
      self.interactionPoints = vtk.vtkPoints()
      self.auxTransformArray[:] = np.zeros(self.auxTransformArray.shape)
      self.auxTransformNode.Modified()
      # qt.QApplication.setOverrideCursor(qt.QCursor(qt.Qt.ArrowCursor))

    elif (event == 'RightButtonPressEvent' or (event == 'KeyPressEvent' and self.interactor.GetKeySym()=='Escape')) and self.smudging:
      self.cancelSmudging()

    elif event == 'MouseMoveEvent' and self.smudging:

      r = int(round(float(self.parameterNode.GetParameter("Radius")) / self.auxTransformNode.GetTransformFromParent().GetDisplacementGrid().GetSpacing()[0])) # Asume isotropic!
      sphereResult = self.createSphere(r)
      currentPoint = np.array(self.xyToRAS(self.interactor.GetEventPosition()))
      currentIndex = self.getCurrentIndex(r, currentPoint, self.auxTransfromRASToIJK)
      self.interactionPoints.InsertNextPoint(currentPoint)

      # apply to transform array
      try:
        self.auxTransformArray[currentIndex] += np.stack([(sphereResult) * i for i in (self.previousPoint - currentPoint)],3) # original
      except ValueError: # error when modifing outside the grid
        self.cancelSmudging()

      # update view
      self.auxTransformNode.Modified()
      # update previous point
      self.previousPoint = currentPoint


  def cancelSmudging(self):
    self.smudging = False
    self.interactionPoints = vtk.vtkPoints()
    self.auxTransformArray[:] = np.zeros(self.auxTransformArray.shape)
    self.auxTransformNode.Modified()

  def resamplePoints(self):
    # resample points 
    sourceCurve = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsCurveNode')
    sourceCurve.SetControlPointPositionsWorld(self.interactionPoints)
    sourceCurve.ResampleCurveWorld(1)  
    sourceCurve.GetControlPointPositionsWorld(self.interactionPoints)
    slicer.mrmlScene.RemoveNode(sourceCurve)

  def getSourceTargetFromPoints(self):
    # source points
    sourceFiducial = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    sourceFiducial.GetDisplayNode().SetGlyphTypeFromString('Sphere3D')
    sourceFiducial.GetDisplayNode().SetGlyphScale(1)
    sourceFiducial.GetDisplayNode().SetVisibility(0)
    sourceFiducial.SetControlPointPositionsWorld(self.interactionPoints)
    sourceFiducial.ApplyTransform(self.parameterNode.GetNodeReference("OutputGridTransform").GetTransformFromParent()) # undo current
    # target
    targetFiducial = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    targetFiducial.GetDisplayNode().SetGlyphTypeFromString('Sphere3D')
    targetFiducial.GetDisplayNode().SetGlyphScale(1)
    targetFiducial.GetDisplayNode().SetVisibility(0)
    targetFiducial.SetControlPointPositionsWorld(self.interactionPoints)
    targetFiducial.ApplyTransform(self.auxTransformNode.GetTransformToParent()) # apply smudge
    targetFiducial.SetName(slicer.mrmlScene.GenerateUniqueName('smudge'))

    return sourceFiducial, targetFiducial

  def createSphere(self, r):
    # create a sphere with redius
    xx, yy, zz = np.mgrid[:2*r+1, :2*r+1, :2*r+1]
    v = 0.5*r
    sphereResult = ((xx-r)/v) ** 2 + ((yy-r)/v) ** 2 + ((zz-r)/v) ** 2
    sphereResult = np.exp(-0.5 * sphereResult)
    sphereResult = sphereResult / sphereResult[r][r][r]
    return sphereResult

  def getCurrentIndex(self, r, currentPoint, RASToIJK):
    # get current IJK
    pos_i,pos_j,pos_k,aux = RASToIJK.MultiplyDoublePoint(np.append(currentPoint, 1))
    k,j,i = int(round(pos_k)),int(round(pos_j)),int(round(pos_i))
    # expand with radius
    currentIndex = slice(k-r,k+r+1), slice(j-r,j+r+1), slice(i-r,i+r+1)
    return currentIndex

  def cleanup(self):
    slicer.mrmlScene.RemoveNode(self.auxTransformNode)
    type(self).cleanAuxTransform()
    #WarpEffectTool.cleanup(self)
    AbstractCircleEffect.cleanup(self)

  @classmethod
  def cleanAuxTransform(cls):
    cls.auxTransformNode = None
    cls.auxTransfromRASToIJK = None
