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

    if event == 'LeftButtonPressEvent' and self.parameterNode.GetParameter("currentEffect") != 'None':
      self.cursorOff()
    elif event == 'LeftButtonReleaseEvent' and self.parameterNode.GetParameter("currentEffect") != 'None':
      qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
      qt.QApplication.processEvents()

    elif event == "KeyPressEvent":
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

  sphere = None
  sphereModelNode = None

  def __init__(self, sliceWidget):
    PointerEffectTool.__init__(self, sliceWidget)
    VTKObservationMixin.__init__(self)

    # init class attributes
    if type(self).sphere is None:
      type(self).sphere = vtk.vtkSphereSource()
      self.sphere.SetPhiResolution(30)
      self.sphere.SetThetaResolution(30)
    if type(self).sphereModelNode is None:
      type(self).sphereModelNode = slicer.modules.models.logic().AddModel(self.sphere.GetOutput())
      self.sphereModelNode.GetDisplayNode().SetVisibility2D(True)
      self.sphereModelNode.GetDisplayNode().SetVisibility3D(False)
      self.sphereModelNode.GetDisplayNode().SetColor(0.7,0.7,0)
      
    self.effectName = self.parameterNode.GetParameter("currentEffect")

    self.addObserver(self.parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateSphere)
    self.updateSphere()
    self.sphereModelNode.GetDisplayNode().SetVisibility(False)

  def processEvent(self, caller=None, event=None):
    if event == 'MouseMoveEvent':
      xy = self.interactor.GetEventPosition()
      self.sphere.SetCenter(self.sliceWidget.sliceLogic().GetSliceNode().GetXYToRAS().MultiplyPoint(xy + (0, 1))[0:3])
      self.sphere.Update()
    elif event == "EnterEvent":
      self.updateSphere()
    elif event == "LeaveEvent":
      self.sphereModelNode.GetDisplayNode().SetVisibility(False)
    elif event == "KeyPressEvent":
      key = self.interactor.GetKeySym()
      if key == 'plus' or key == 'equal':
        self.scaleRadius(1.2)
      if key == 'minus' or key == 'underscore':
        self.scaleRadius(0.8) 

    PointerEffectTool.processEvent(self, caller, event)

  def scaleRadius(self,scaleFactor):
    radius = float(self.parameterNode.GetParameter(self.effectName + "Radius"))
    self.parameterNode.SetParameter(self.effectName + "Radius", str(radius * scaleFactor))

  def updateSphere(self,caller=None,event=None):
    if self.sphere:
      self.sphereModelNode.GetDisplayNode().SetVisibility((self.effectName=='Smooth' and int(self.parameterNode.GetParameter("SmoothUseRadius"))) or self.effectName=='Smudge')
      self.sphere.SetRadius(float(self.parameterNode.GetParameter(self.effectName + "Radius")))
      self.sphere.Update()

  def cleanup(self):
    slicer.mrmlScene.RemoveNode(self.sphereModelNode)
    type(self).cleanSphere()
    super(PointerEffectTool,self).cleanup()

  @classmethod
  def cleanSphere(cls):
    cls.sphere = None
    cls.sphereModelNode = None



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

    self.transform = vtk.vtkThinPlateSplineTransform()
    self.transform.SetBasisToR()
    self.transform.Inverse()
    self.auxNodes = []

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

    PointerEffectTool.processEvent(self, caller, event)

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
          self.cursorOff()
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

