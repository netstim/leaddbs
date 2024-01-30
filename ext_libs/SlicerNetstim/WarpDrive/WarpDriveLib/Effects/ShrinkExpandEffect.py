import vtk, qt, slicer
from math import sqrt, cos, sin

import numpy as np

from .CircleEffect import AbstractCircleEffect


class AbstractShrinkExpandEffect(AbstractCircleEffect):

  def __init__(self, sliceWidget):

    # keep a flag since events such as sliceNode modified
    # may come during superclass construction, which will
    # invoke our processEvents method
    self.initialized = False

    AbstractCircleEffect.__init__(self, sliceWidget)

    self.previewing = False

    self.transform = vtk.vtkThinPlateSplineTransform()
    self.transform.SetBasisToR()
    self.transform.Inverse()
    self.auxNodes = []

    self.initialized = True

  def cleanup(self):
    """
    call superclass to clean up actor
    """
    AbstractCircleEffect.cleanup(self)


  def processEvent(self, caller=None, event=None):
    """
    handle events from the render window interactor
    """

    if not self.initialized:
      return

    AbstractCircleEffect.processEvent(self, caller, event)

    if event == "LeftButtonPressEvent":
      self.previewing = True
      xy = self.interactor.GetEventPosition()
      point = self.xyToRAS(xy)
      self.initTransform(point)
    
    elif event == 'RightButtonPressEvent' or (event == 'KeyPressEvent' and self.interactor.GetKeySym()=='Escape'):
      self.previewing = False
      self.removeAuxNodes()

    self.sliceView.scheduleRender()

  def removeAuxNodes(self):
    while len(self.auxNodes):
      slicer.mrmlScene.RemoveNode(self.auxNodes.pop())

  def initTransform(self, centerPoint):
    sliceToRAS = self.sliceLogic.GetSliceNode().GetSliceToRAS()
    planeNormal   = np.array([sliceToRAS.GetElement(0,2), sliceToRAS.GetElement(1,2), sliceToRAS.GetElement(2,2)])
    inPlaneVector = np.array([sliceToRAS.GetElement(0,0), sliceToRAS.GetElement(1,0), sliceToRAS.GetElement(2,0)])

    ammount = float(self.parameterNode.GetParameter("ShrinkExpandAmmount")) / 100.0
    ammount = -ammount if self.parameterNode.GetParameter("ShrinkExpandMode") == "Shrink" else ammount
    radius = float(self.parameterNode.GetParameter("Radius"))

    sourcePoints, targetPoints = vtk.vtkPoints(), vtk.vtkPoints()
    point, transformedPoint = np.zeros(3), np.zeros(3)

    for deg in range(0,360,30):
      for rad,pts in zip([radius, radius*(1+ammount)],[sourcePoints, targetPoints]):
        rotationTransform = vtk.vtkTransform()  
        rotationTransform.Translate(centerPoint)      
        rotationTransform.RotateWXYZ(deg, planeNormal)
        rotationTransform.Translate(rad * inPlaneVector)
        rotationTransform.TransformPoint(point, transformedPoint)
        pts.InsertNextPoint(transformedPoint)

    sourceFiducialNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    sourceFiducialNode.GetDisplayNode().SetVisibility(0)
    sourceFiducialNode.SetControlPointPositionsWorld(sourcePoints)
    self.transform.SetSourceLandmarks(sourcePoints)
    self.transform.SetTargetLandmarks(targetPoints)
    transformNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
    transformNode.SetAndObserveTransformFromParent(self.transform)
    transformNode.CreateDefaultDisplayNodes()
    transformNode.GetDisplayNode().SetAndObserveGlyphPointsNode(sourceFiducialNode)
    transformNode.GetDisplayNode().SetVisibility(1)
    transformNode.GetDisplayNode().SetVisibility2D(1)
    transformNode.GetDisplayNode().SetVisibility3D(0)
    self.auxNodes.append(sourceFiducialNode)
    self.auxNodes.append(transformNode)

