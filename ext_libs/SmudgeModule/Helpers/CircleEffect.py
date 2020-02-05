import vtk, qt, slicer
from math import sqrt, cos, sin
from Helpers import Effect
from slicer.util import VTKObservationMixin


class CircleEffectTool(Effect.EffectTool, VTKObservationMixin):

  def __init__(self, sliceWidget):
    Effect.EffectTool.__init__(self, sliceWidget)
    VTKObservationMixin.__init__(self)

    self.addObserver(self.parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGlyph)
    
    self.rasToXY = vtk.vtkMatrix4x4()

    self.brush = vtk.vtkPolyData()
    self.createGlyph()

    self.mapper = vtk.vtkPolyDataMapper2D()
    self.mapper.SetInputData(self.brush)

    self.actor = vtk.vtkActor2D()
    self.actor.GetProperty().SetColor(.7, .7, 0)
    self.actor.SetMapper(self.mapper)
    self.actor.VisibilityOff()
    self.actors.append(self.actor) 

    self.renderer.AddActor2D(self.actor)   

  def createGlyph(self):
    """
    create a brush circle of the right radius in XY space
    - assume uniform scaling between XY and RAS which
      is enforced by the view interactors
    """
    polyData = self.brush
    radius = float(self.parameterNode.GetParameter("radius"))

    sliceNode = self.sliceWidget.sliceLogic().GetSliceNode()
    self.rasToXY.DeepCopy(sliceNode.GetXYToRAS())
    self.rasToXY.Invert()

    maximum, maxIndex = 0,0
    for index in range(3):
      if abs(self.rasToXY.GetElement(0, index)) > maximum:
        maximum = abs(self.rasToXY.GetElement(0, index))
        maxIndex = index
    point = [0, 0, 0, 0]
    point[maxIndex] = radius
    xyRadius = self.rasToXY.MultiplyPoint(point)
    
    xyRadius = sqrt( xyRadius[0]**2 + xyRadius[1]**2 + xyRadius[2]**2 )

    #if self.pixelMode:
    #  xyRadius = 0.01

    # make a circle paint brush
    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()
    polyData.SetPoints(points)
    polyData.SetLines(lines)
    PI = 3.1415926
    TWOPI = PI * 2
    PIoverSIXTEEN = PI / 16
    prevPoint = -1
    firstPoint = -1
    angle = 0
    while angle <= TWOPI:
      x = xyRadius * cos(angle)
      y = xyRadius * sin(angle)
      p = points.InsertNextPoint( x, y, 0 )
      if prevPoint != -1:
        idList = vtk.vtkIdList()
        idList.InsertNextId(prevPoint)
        idList.InsertNextId(p)
        polyData.InsertNextCell( vtk.VTK_LINE, idList )
      prevPoint = p
      if firstPoint == -1:
        firstPoint = p
      angle = angle + PIoverSIXTEEN

    # make the last line in the circle
    idList = vtk.vtkIdList()
    idList.InsertNextId(p)
    idList.InsertNextId(firstPoint)
    polyData.InsertNextCell( vtk.VTK_LINE, idList )

    # update
    self.sliceView.scheduleRender()

  def processEvent(self, caller=None, event=None):
    if event == 'MouseMoveEvent':
      xy = self.interactor.GetEventPosition()
      self.actor.SetPosition(xy)
      self.sliceView.scheduleRender()
    elif event == "EnterEvent":
      self.createGlyph() # update in case radius changed from slider
      self.actor.VisibilityOn()
    elif event == "LeaveEvent":
      self.actor.VisibilityOff()
    elif event == "KeyPressEvent":
      key = self.interactor.GetKeySym()
      if key == 'plus' or key == 'equal':
        self.scaleRadius(1.2)
      if key == 'minus' or key == 'underscore':
        self.scaleRadius(0.8) 

    if caller and caller.IsA('vtkMRMLSliceNode'):
      self.createGlyph()

  def scaleRadius(self,scaleFactor):
    radius = float(self.parameterNode.GetParameter("radius"))
    self.parameterNode.SetParameter( "radius", str(radius * scaleFactor) )
    
  def updateGlyph(self, caller, event):
    self.createGlyph()

  def cleanup(self):
    super(CircleEffectTool,self).cleanup()
