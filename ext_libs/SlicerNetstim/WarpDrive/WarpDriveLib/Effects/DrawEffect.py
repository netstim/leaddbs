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


    # interaction state variables
    self.activeSlice = None
    self.lastInsertSLiceNodeMTime = None
    self.actionState = None

    # initialization
    self.xyPoints = vtk.vtkPoints()
    self.rasPoints = vtk.vtkPoints()
    self.polyData = self.createPolyData()

    self.mapper = vtk.vtkPolyDataMapper2D()
    self.actor = vtk.vtkActor2D()
    self.mapper.SetInputData(self.polyData)
    self.actor.SetMapper(self.mapper)
    property_ = self.actor.GetProperty()
    property_.SetColor(1,1,0)
    property_.SetLineWidth(1)
    self.renderer.AddActor2D( self.actor )
    self.actors.append( self.actor )

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

    # events from the interactor
    if event == "LeftButtonPressEvent":
      if not self.actionState:
        self.actionState = "drawing"
        xy = self.interactor.GetEventPosition()
        self.addPoint(self.xyToRAS(xy))
        self.abortEvent(event)
      
    elif event == "LeftButtonReleaseEvent":
      if self.actionState == "drawing":
        if self.rasPoints.GetNumberOfPoints() > 1:
          self.actionState = None
        else:
          self.actionState = "placingPoint"
          self.initTransform()
          # self.cursorOff()
      elif self.actionState == "placingPoint":
        self.removeAuxNodes()
        self.actionState = None

    elif event == "MouseMoveEvent":
      if self.actionState == "drawing":
        xy = self.interactor.GetEventPosition()
        self.addPoint(self.xyToRAS(xy))
        self.abortEvent(event)
      elif self.actionState == "placingPoint":
        p = vtk.vtkPoints()
        xy = self.interactor.GetEventPosition()
        p.InsertNextPoint(self.xyToRAS(xy))
        self.transform.SetTargetLandmarks(p)

    self.positionActors()

  def removeAuxNodes(self):
    while len(self.auxNodes):
      slicer.mrmlScene.RemoveNode(self.auxNodes.pop())

  def initTransform(self):
    sourceFiducialNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
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

  def positionActors(self):
    """
    update draw feedback to follow slice node
    """
    sliceLogic = self.sliceWidget.sliceLogic()
    sliceNode = sliceLogic.GetSliceNode()
    rasToXY = vtk.vtkTransform()
    rasToXY.SetMatrix( sliceNode.GetXYToRAS() )
    rasToXY.Inverse()
    self.xyPoints.Reset()
    rasToXY.TransformPoints( self.rasPoints, self.xyPoints )
    self.polyData.Modified()
    self.sliceView.scheduleRender()

  def createPolyData(self):
    """make an empty single-polyline polydata"""

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(self.xyPoints)

    lines = vtk.vtkCellArray()
    polyData.SetLines(lines)
    idArray = lines.GetData()
    idArray.Reset()
    idArray.InsertNextTuple1(0)

    polygons = vtk.vtkCellArray()
    polyData.SetPolys(polygons)
    idArray = polygons.GetData()
    idArray.Reset()
    idArray.InsertNextTuple1(0)

    return polyData


  def resetPolyData(self):
    """return the polyline to initial state with no points"""
    lines = self.polyData.GetLines()
    idArray = lines.GetData()
    idArray.Reset()
    idArray.InsertNextTuple1(0)
    self.xyPoints.Reset()
    self.rasPoints.Reset()
    lines.SetNumberOfCells(0)
    self.activeSlice = None

  def addPoint(self,ras):
    """add a world space point to the current outline"""
    # store active slice when first point is added
    sliceLogic = self.sliceWidget.sliceLogic()
    currentSlice = sliceLogic.GetSliceOffset()
    if not self.activeSlice:
      self.activeSlice = currentSlice

    # don't allow adding points on except on the active slice (where
    # first point was laid down)
    if self.activeSlice != currentSlice: return

    # keep track of node state (in case of pan/zoom)
    sliceNode = sliceLogic.GetSliceNode()
    self.lastInsertSliceNodeMTime = sliceNode.GetMTime()

    p = self.rasPoints.InsertNextPoint(ras)
    lines = self.polyData.GetLines()
    idArray = lines.GetData()
    idArray.InsertNextTuple1(p)
    idArray.SetTuple1(0, idArray.GetNumberOfTuples()-1)
    lines.SetNumberOfCells(1)

