import vtk, qt, slicer
from math import sqrt, cos, sin

from .CircleEffect import AbstractCircleEffect


class AbstractPointToPointEffect(AbstractCircleEffect):

  def __init__(self, sliceWidget):

    # keep a flag since events such as sliceNode modified
    # may come during superclass construction, which will
    # invoke our processEvents method
    self.initialized = False

    AbstractCircleEffect.__init__(self, sliceWidget)

    # interaction state variables
    self.actionState = None

    # initialization
    self.xyPoints = vtk.vtkPoints()
    self.rasPoints = vtk.vtkPoints()

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
      
    if event == "LeftButtonReleaseEvent":
      if self.actionState is None:
        xy = self.interactor.GetEventPosition()
        self.addPoint(self.xyToRAS(xy))
        self.initTransform()
        self.actionState = "placingPoint"
      elif self.actionState == "placingPoint":
        self.removeAuxNodes()
        self.actionState = None

    elif event == "MouseMoveEvent":
      if self.actionState == "placingPoint":
        p = vtk.vtkPoints()
        xy = self.interactor.GetEventPosition()
        p.InsertNextPoint(self.xyToRAS(xy))
        self.transform.SetTargetLandmarks(p)

    elif event == 'RightButtonPressEvent' or (event == 'KeyPressEvent' and self.interactor.GetKeySym()=='Escape'):
      self.removeAuxNodes()
      self.resetPoints()
      self.actionState = None

    # self.positionActors()
    self.sliceView.scheduleRender()

  def removeAuxNodes(self):
    while len(self.auxNodes):
      slicer.mrmlScene.RemoveNode(self.auxNodes.pop())

  def initTransform(self):
    sourceFiducialNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    sourceFiducialNode.GetDisplayNode().SetGlyphTypeFromString('Sphere3D')
    sourceFiducialNode.GetDisplayNode().SetGlyphScale(1)
    sourceFiducialNode.GetDisplayNode().SetVisibility(0)
    sourceFiducialNode.SetControlPointPositionsWorld(self.rasPoints)
    self.transform.SetSourceLandmarks(self.rasPoints)
    self.transform.SetTargetLandmarks(self.rasPoints)
    transformNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
    transformNode.SetAndObserveTransformFromParent(self.transform)
    transformNode.CreateDefaultDisplayNodes()
    transformNode.GetDisplayNode().SetVisibility(1)
    transformNode.GetDisplayNode().SetVisibility2D(1)
    transformNode.GetDisplayNode().SetVisibility3D(0)
    transformNode.GetDisplayNode().SetAndObserveGlyphPointsNode(sourceFiducialNode)
    self.auxNodes.append(sourceFiducialNode)
    self.auxNodes.append(transformNode)

  def resetPoints(self):
    """return the points to initial state with no points"""
    self.xyPoints.Reset()
    self.rasPoints.Reset()

  def addPoint(self,ras):
    """add a world space point"""
    p = self.rasPoints.InsertNextPoint(ras)


