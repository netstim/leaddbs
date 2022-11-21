import vtk, slicer
import numpy as np


from ..Widgets.ToolWidget   import AbstractToolWidget
from ..Effects.PointToPointEffect import AbstractPointToPointEffect

from ..Helpers import GridNodeHelper

class PointToPointToolWidget(AbstractToolWidget):
    
  def __init__(self):
    toolTip = ''
    AbstractToolWidget.__init__(self, 'PointToPoint', toolTip)


class PointToPointToolEffect(AbstractPointToPointEffect):


  def __init__(self, sliceWidget):
    AbstractPointToPointEffect.__init__(self, sliceWidget)


  def processEvent(self, caller=None, event=None):

    AbstractPointToPointEffect.processEvent(self, caller, event) 

    if event == 'LeftButtonReleaseEvent' and not self.actionState:

      # create source and target fiducials from points
      sourceFiducial, targetFiducial = self.getSourceTargetFromPoints()

      # reset
      self.resetPoints()

      sourceFiducial.ApplyTransform(self.parameterNode.GetNodeReference("OutputGridTransform").GetTransformFromParent()) # undo current

      self.setFiducialNodeAs("Source", sourceFiducial, targetFiducial.GetName(), self.parameterNode.GetParameter("Radius"))
      self.setFiducialNodeAs("Target", targetFiducial, targetFiducial.GetName(), self.parameterNode.GetParameter("Radius"))

      self.parameterNode.SetParameter("Update","true")
 

  def getSourceTargetFromPoints(self):
    # get clicked points from the vtk thinplate transform
    # source
    sourceFiducial = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    sourceFiducial.SetControlPointPositionsWorld(self.transform.GetSourceLandmarks())
    sourceFiducial.GetDisplayNode().SetGlyphTypeFromString('Sphere3D')
    sourceFiducial.GetDisplayNode().SetGlyphScale(1)
    sourceFiducial.GetDisplayNode().SetVisibility(0)
    # target
    targetFiducial = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    targetFiducial.SetControlPointPositionsWorld(self.transform.GetTargetLandmarks())
    targetFiducial.GetDisplayNode().SetGlyphTypeFromString('Sphere3D')
    sourceFiducial.GetDisplayNode().SetGlyphScale(1)
    targetFiducial.GetDisplayNode().SetVisibility(0)
    targetFiducial.SetName(slicer.mrmlScene.GenerateUniqueName('point'))
    return sourceFiducial, targetFiducial


  def cleanup(self):
    AbstractPointToPointEffect.cleanup(self)

