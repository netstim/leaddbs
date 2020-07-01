import slicer
from slicer.util import VTKObservationMixin

from ..Widgets.ToolWidget import AbstractToolWidget
from ..Effects.PointerEffect import AbstractPointerEffect

class NoneToolWidget(AbstractToolWidget, VTKObservationMixin):
    
  def __init__(self):

    toolTip = 'None. Disable Effects.'
    AbstractToolWidget.__init__(self, 'None', toolTip)
    VTKObservationMixin.__init__(self)

    interactionNode = slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLInteractionNode')
    self.addObserver(interactionNode, interactionNode.InteractionModeChangedEvent, self.onInteractionModeChanged) 

  def onInteractionModeChanged(self, caller=None, event=None):
    interactionNode = slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLInteractionNode')
    if interactionNode.GetCurrentInteractionMode() != interactionNode.ViewTransform:
      self.effectButton.setChecked(True)

class NoneToolEffect(AbstractPointerEffect):

  def __init__(self, sliceWidget):
    AbstractPointerEffect.__init__(self, sliceWidget)
    
