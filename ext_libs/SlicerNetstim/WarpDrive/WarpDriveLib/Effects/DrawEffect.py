import vtk, qt, slicer
from math import sqrt, cos, sin

from .CircleEffect import AbstractCircleEffect


class AbstractDrawEffect(AbstractCircleEffect):

  def __init__(self, sliceWidget):

    # keep a flag since events such as sliceNode modified
    # may come during superclass construction, which will
    # invoke our processEvents method
    self.initialized = False

    AbstractCircleEffect.__init__(self, sliceWidget)

    self.drawnCurveNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsCurveNode")
    self.drawnCurveNode.GetDisplayNode().SetPropertiesLabelVisibility(0)
    self.drawnCurveNode.GetDisplayNode().SetLineThickness(1)
    self.drawnCurveNode.GetDisplayNode().SetGlyphScale(1)

    # interaction state variables
    self.activeSlice = None
    self.lastInsertSLiceNodeMTime = None
    self.actionState = None

    self.initialized = True

  def cleanup(self):
    """
    call superclass to clean up actor
    """
    AbstractCircleEffect.cleanup(self)
    slicer.mrmlScene.RemoveNode(self.drawnCurveNode)

  def processEvent(self, caller=None, event=None):
    """
    handle events from the render window interactor
    """

    if not self.initialized:
      return

    AbstractCircleEffect.processEvent(self, caller, event)

    # events from the interactor
    if event == "LeftButtonPressEvent":
      self.actionState = "drawing"
      xy = self.interactor.GetEventPosition()
      self.addPoint(self.xyToRAS(xy))
      self.abortEvent(event)
      
    elif event == "LeftButtonReleaseEvent":
      self.actionState = None

    elif event == "MouseMoveEvent":
      if self.actionState == "drawing":
        xy = self.interactor.GetEventPosition()
        self.addPoint(self.xyToRAS(xy))
        self.abortEvent(event)

    elif event == 'RightButtonPressEvent' or (event == 'KeyPressEvent' and self.interactor.GetKeySym()=='Escape'):
      self.resetDrawing()
      self.actionState = None

  def addPoint(self,ras):
    """add a world space point to the current outline"""
    resampleDistance = 1.0
    self.drawnCurveNode.AddControlPoint(vtk.vtkVector3d(ras))
    if self.drawnCurveNode.GetCurveLengthWorld() > resampleDistance:
      self.drawnCurveNode.ResampleCurveWorld(resampleDistance)

  def resetDrawing(self):
    self.drawnCurveNode.RemoveAllControlPoints()