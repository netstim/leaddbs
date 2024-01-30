import vtk, slicer, ctk, qt
import numpy as np


from ..Widgets.ToolWidget   import AbstractToolWidget
from ..Effects.ShrinkExpandEffect import AbstractShrinkExpandEffect

from ..Helpers import GridNodeHelper

class ShrinkExpandToolWidget(AbstractToolWidget):
    
  def __init__(self):
    toolTip = ''
    AbstractToolWidget.__init__(self, 'ShrinkExpand', toolTip)

    # modes
    shrinkAction = qt.QAction('Shrink', self.effectButton)
    shrinkAction.setCheckable(True)
    shrinkAction.setChecked(True)
    expandAction = qt.QAction('Expand', self.effectButton)
    expandAction.setCheckable(True)
    actionsGroup = qt.QActionGroup(self.effectButton)
    actionsGroup.addAction(shrinkAction)
    actionsGroup.addAction(expandAction)

    # ammount
    ammountSlider = ctk.ctkSliderWidget()
    ammountSlider.singleStep = 5
    ammountSlider.minimum = 5
    ammountSlider.maximum = 95
    ammountSlider.value = 25
    ammountSlider.setToolTip("Ammount")

    ammountAction = qt.QWidgetAction(self.effectButton)
    ammountAction.setDefaultWidget(ammountSlider)
    # ammountAction.setToolTip("Ammount")

    settingsMenu = qt.QMenu(self.effectButton)
    settingsMenu.addAction(ammountAction)
    modeMenu = settingsMenu.addMenu("Mode")
    modeMenu.addActions(actionsGroup.actions())
    
    self.effectButton.setMenu(settingsMenu)
    self.effectButton.setPopupMode(self.effectButton.DelayedPopup)


class ShrinkExpandToolEffect(AbstractShrinkExpandEffect):


  def __init__(self, sliceWidget):
    AbstractShrinkExpandEffect.__init__(self, sliceWidget)


  def processEvent(self, caller=None, event=None):

    AbstractShrinkExpandEffect.processEvent(self, caller, event) 

    if event == 'LeftButtonReleaseEvent' and self.previewing:

      # create source and target fiducials from points
      sourceFiducial, targetFiducial = self.getSourceTargetFromPoints()

      sourceFiducial.ApplyTransform(self.parameterNode.GetNodeReference("OutputGridTransform").GetTransformFromParent()) # undo current

      if int(self.parameterNode.GetParameter("ModifiableCorrections")):
        self.modifyPreviousCorrections(sourceFiducial, targetFiducial)

      self.setFiducialNodeAs("Source", sourceFiducial, targetFiducial.GetName(), self.parameterNode.GetParameter("Radius"))
      self.setFiducialNodeAs("Target", targetFiducial, targetFiducial.GetName(), self.parameterNode.GetParameter("Radius"))

      self.parameterNode.SetParameter("Update","true")
 
      self.previewing = False
      self.removeAuxNodes()

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
    targetFiducial.SetName(slicer.mrmlScene.GenerateUniqueName(self.parameterNode.GetParameter("ShrinkExpandMode")))
    return sourceFiducial, targetFiducial


  def cleanup(self):
    AbstractShrinkExpandEffect.cleanup(self)

