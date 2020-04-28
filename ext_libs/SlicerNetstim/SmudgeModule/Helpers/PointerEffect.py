import vtk, qt, slicer
from math import sqrt, cos, sin
from slicer.util import VTKObservationMixin

from . import Effect


class PointerEffectTool(Effect.EffectTool):

  def __init__(self, sliceWidget):
    Effect.EffectTool.__init__(self, sliceWidget)
    self.previousVisibleIDs = vtk.vtkIdList()
    self.previousBackgroundNodeID = None

  def processEvent(self, caller=None, event=None):

    if event == "KeyPressEvent":
      key = self.interactor.GetKeySym()
      if key.lower() == 's':
        if not self.previousVisibleIDs.GetNumberOfIds():
          # get all ids from subject hierarchy node
          shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
          shNode.GetItemChildren(shNode.GetSceneItemID(), self.previousVisibleIDs, True)
          ids = [self.previousVisibleIDs.GetId(i) for i in range(self.previousVisibleIDs.GetNumberOfIds())]
          for id in ids:
            data = shNode.GetItemDataNode(id)
            if isinstance(data, slicer.vtkMRMLModelNode) and data.GetDisplayNode() and data.GetDisplayNode().GetVisibility():
              data.GetDisplayNode().SetVisibility(0) # hide visible models
            else:
              self.previousVisibleIDs.DeleteId(id) # remove non visible models IDs
      elif key.lower() == 't':
        if self.previousBackgroundNodeID is None:
          compositeNode = slicer.app.layoutManager().sliceWidget('Red').sliceLogic().GetSliceCompositeNode()
          self.previousBackgroundNodeID = compositeNode.GetBackgroundVolumeID()
          self.previousForegroundOpacity = compositeNode.GetForegroundOpacity()
          slicer.util.setSliceViewerLayers(background = None)

    elif event == "KeyReleaseEvent":
      key = self.interactor.GetKeySym()
      if key.lower() == 's':
        shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
        for i in range(self.previousVisibleIDs.GetNumberOfIds()):
          shNode.GetItemDataNode(self.previousVisibleIDs.GetId(i)).GetDisplayNode().SetVisibility(1)
        self.previousVisibleIDs.Reset()
      if key.lower() == 't':
        slicer.util.setSliceViewerLayers(background = self.previousBackgroundNodeID, foregroundOpacity=self.previousForegroundOpacity)
        self.previousBackgroundNodeID = None


class CircleEffectTool(PointerEffectTool, VTKObservationMixin):

  def __init__(self, sliceWidget):
    PointerEffectTool.__init__(self, sliceWidget)
    VTKObservationMixin.__init__(self)

    self.addObserver(self.parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGlyph)
    
    self.rasToXY = vtk.vtkMatrix4x4()
    self.brush = vtk.vtkPolyData()
    self.updateGlyph([],[])
    self.mapper = vtk.vtkPolyDataMapper2D()
    self.mapper.SetInputData(self.brush)
    self.actor = vtk.vtkActor2D()
    self.actor.GetProperty().SetColor(.7, .7, 0)
    self.actor.SetMapper(self.mapper)
    self.actor.VisibilityOff()
    self.actors.append(self.actor) 
    self.renderer.AddActor2D(self.actor)


  def createGlyph(self, polyData, radius):
    """
    create a brush circle of the right radius in XY space
    - assume uniform scaling between XY and RAS which
      is enforced by the view interactors
    """
    #polyData = self.brush
    #radius = float(self.parameterNode.GetParameter("radius"))

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
      self.updateGlyph([],[]) # update in case radius changed from slider
      self.actor.VisibilityOn()
    elif event == "LeaveEvent":
      self.actor.VisibilityOff()
    elif event == 'LeftButtonPressEvent':
      self.cursorOff()
    elif event == 'LeftButtonReleaseEvent':
      self.cursorOn()
    elif event == "KeyPressEvent":
      key = self.interactor.GetKeySym()
      if key == 'plus' or key == 'equal':
        self.scaleRadius(1.2)
      if key == 'minus' or key == 'underscore':
        self.scaleRadius(0.8) 

    if caller and caller.IsA('vtkMRMLSliceNode'):
      self.updateGlyph([],[])

    PointerEffectTool.processEvent(self, caller, event)

  def scaleRadius(self,scaleFactor):
    effectRadius = self.parameterNode.GetParameter("currentEffect") + "Radius"
    radius = float(self.parameterNode.GetParameter(effectRadius))
    self.parameterNode.SetParameter(effectRadius, str(radius * scaleFactor) )
    
  def updateGlyph(self, caller, event):
    effectName = self.parameterNode.GetParameter("currentEffect")
    if effectName in ['Smudge','Blur']:
      r = float(self.parameterNode.GetParameter(effectName + "Radius"))
      self.createGlyph(self.brush, r)

  def cleanup(self):
    super(CircleEffectTool,self).cleanup()



class DrawEffectTool(PointerEffectTool):

  def __init__(self, sliceWidget):

    # keep a flag since events such as sliceNode modified
    # may come during superclass construction, which will
    # invoke our processEvents method
    self.initialized = False

    PointerEffectTool.__init__(self, sliceWidget)


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

    self.initialized = True

  def cleanup(self):
    """
    call superclass to clean up actor
    """
    super(DrawEffectTool,self).cleanup()


  def processEvent(self, caller=None, event=None):
    """
    handle events from the render window interactor
    """

    if not self.initialized:
      return

    # events from the interactor
    if event == "LeftButtonPressEvent":
      self.actionState = "drawing"
      xy = self.interactor.GetEventPosition()
      self.addPoint(self.xyToRAS(xy))
      self.abortEvent(event)
      self.cursorOff()
      
    elif event == "LeftButtonReleaseEvent":
      self.actionState = None
      self.cursorOn()
      self.resetPolyData()

    elif event == "MouseMoveEvent":
      if self.actionState == "drawing":
        xy = self.interactor.GetEventPosition()
        self.addPoint(self.xyToRAS(xy))
        self.abortEvent(event)

    PointerEffectTool.processEvent(self, caller, event)

    self.positionActors()

  def xyToRAS(self,xyPoint):
    """return r a s for a given x y"""
    sliceNode = self.sliceLogic.GetSliceNode()
    rast = sliceNode.GetXYToRAS().MultiplyPoint(xyPoint + (0,1,))
    return rast[:3]

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

