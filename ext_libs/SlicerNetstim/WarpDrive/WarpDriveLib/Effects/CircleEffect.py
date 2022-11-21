import vtk, qt, slicer
from math import sqrt, cos, sin
from slicer.util import VTKObservationMixin

from .PointerEffect import AbstractPointerEffect


class AbstractCircleEffect(AbstractPointerEffect, VTKObservationMixin):

  sphere = None
  sphereModelNode = None

  def __init__(self, sliceWidget):
    AbstractPointerEffect.__init__(self, sliceWidget)
    VTKObservationMixin.__init__(self)

    # init class attributes
    if type(self).sphere is None:
      type(self).sphere = vtk.vtkSphereSource()
      self.sphere.SetPhiResolution(30)
      self.sphere.SetThetaResolution(30)
    if type(self).sphereModelNode is None:
      type(self).sphereModelNode = slicer.modules.models.logic().AddModel(self.sphere.GetOutput())
      self.sphereModelNode.SetName('auxSphereModel')
      self.sphereModelNode.GetDisplayNode().SetVisibility2D(True)
      self.sphereModelNode.GetDisplayNode().SetVisibility3D(False)
      self.sphereModelNode.GetDisplayNode().SetColor(0.7,0.7,0)
      
    self.addObserver(self.parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateSphere)
    self.updateSphere()
    self.sphereModelNode.GetDisplayNode().SetVisibility(False)

  def processEvent(self, caller=None, event=None):
    if event == 'MouseMoveEvent' or event == 'MouseWheelForwardEvent' or event == 'MouseWheelBackwardEvent':
      xy = self.interactor.GetEventPosition()
      self.sphere.SetCenter(self.xyToRAS(xy))
      self.sphere.Update()
    elif event == "EnterEvent":
      self.sphereModelNode.GetDisplayNode().SetVisibility(True)
      self.updateSphere()
    elif event == "LeaveEvent":
      self.sphereModelNode.GetDisplayNode().SetVisibility(False)
    elif event == "KeyPressEvent":
      key = self.interactor.GetKeySym()
      if key == 'plus' or key == 'equal':
        self.scaleRadius(1.2)
      if key == 'minus' or key == 'underscore':
        self.scaleRadius(0.8) 

    AbstractPointerEffect.processEvent(self, caller, event)

  def scaleRadius(self,scaleFactor):
    radius = float(self.parameterNode.GetParameter("Radius"))
    self.parameterNode.SetParameter("Radius", "%.2f" % (radius * scaleFactor))

  def updateSphere(self,caller=None,event=None):
    if self.sphere:
      self.sphere.SetRadius(float(self.parameterNode.GetParameter("Radius")))
      self.sphere.Update()

  def cleanup(self):
    slicer.mrmlScene.RemoveNode(self.sphereModelNode)
    self.removeObservers()
    type(self).cleanSphere()
    AbstractPointerEffect.cleanup(self)

  def xyToRAS(self,xyPoint):
    """return r a s for a given x y"""
    sliceNode = self.sliceLogic.GetSliceNode()
    rast = sliceNode.GetXYToRAS().MultiplyPoint(xyPoint + (0,1,))
    return rast[:3]

  @classmethod
  def cleanSphere(cls):
    cls.sphere = None
    cls.sphereModelNode = None

