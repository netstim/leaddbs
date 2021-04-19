import vtk, qt, slicer


class AbstractEffect():
  """
  One instance of this will be created per-view when the effect
  is selected.  It is responsible for implementing feedback and
  label map changes in response to user input.
  This class observes the editor parameter node to configure itself
  and queries the current view for background and label volume
  nodes to operate on.
  """

  def __init__(self,sliceWidget):

    # sliceWidget to operate on and convenience variables
    # to access the internals
    self.sliceWidget = sliceWidget
    self.sliceLogic = sliceWidget.sliceLogic()
    self.sliceView = self.sliceWidget.sliceView()
    self.interactor = self.sliceView.interactorStyle().GetInteractor()
    self.renderWindow = self.sliceWidget.sliceView().renderWindow()
    self.renderer = self.renderWindow.GetRenderers().GetItemAsObject(0)
    #self.editUtil = EditUtil.EditUtil()

    # optionally set by users of the class
    self.undoRedo = None

    # actors in the renderer that need to be cleaned up on destruction
    self.actors = []

    # the current operation
    self.actionState = None

    # set up observers on the interactor
    # - keep track of tags so these can be removed later
    # - currently all editor effects are restricted to these events
    # - make the observers high priority so they can override other
    #   event processors
    self.interactorObserverTags = []
    events = ( vtk.vtkCommand.LeftButtonPressEvent,
      vtk.vtkCommand.LeftButtonReleaseEvent,
      vtk.vtkCommand.MiddleButtonPressEvent,
      vtk.vtkCommand.MiddleButtonReleaseEvent,
      vtk.vtkCommand.RightButtonPressEvent,
      vtk.vtkCommand.RightButtonReleaseEvent,
      vtk.vtkCommand.LeftButtonDoubleClickEvent,
      vtk.vtkCommand.MouseMoveEvent,
      vtk.vtkCommand.KeyPressEvent,
      vtk.vtkCommand.KeyReleaseEvent,
      vtk.vtkCommand.EnterEvent,
      vtk.vtkCommand.LeaveEvent,
      vtk.vtkCommand.MouseWheelForwardEvent,
      vtk.vtkCommand.MouseWheelBackwardEvent)
    for e in events:
      tag = self.interactor.AddObserver(e, self.processEvent, 1.0)
      self.interactorObserverTags.append(tag)

    self.sliceNodeTags = []
    sliceNode = self.sliceLogic.GetSliceNode()
    tag = sliceNode.AddObserver(vtk.vtkCommand.ModifiedEvent, self.processEvent, 1.0)
    self.sliceNodeTags.append(tag)

    # spot for tracking the current cursor while it is turned off for paining
    self.savedCursor = None

  def processEvent(self, caller=None, event=None):
    """Event filter that lisens for certain key events that
    should be responded to by all events.
    Currently:
      '\\' - pick up paint color from current location (eyedropper)
    """
    if event == "KeyPressEvent":
      key = self.interactor.GetKeySym()
      if key.lower() == 's':
        return True
    return False

  def cursorOff(self):
    """Turn off and save the current cursor so
    the user can see the background image during editing"""
    qt.QApplication.setOverrideCursor(qt.QCursor(10))
    #self.savedCursor = self.sliceWidget.cursor
    #qt_BlankCursor = 10
    #self.sliceWidget.setCursor(qt.QCursor(qt_BlankCursor))

  def cursorOn(self):
    """Restore the saved cursor if it exists, otherwise
    just restore the default cursor"""
    qt.QApplication.restoreOverrideCursor()
    #if self.savedCursor:
    #  self.sliceWidget.setCursor(self.savedCursor)
    #else:
    #  self.sliceWidget.unsetCursor()

  def abortEvent(self,event):
    """Set the AbortFlag on the vtkCommand associated
    with the event - causes other things listening to the
    interactor not to receive the events"""
    # TODO: make interactorObserverTags a map to we can
    # explicitly abort just the event we handled - it will
    # be slightly more efficient
    for tag in self.interactorObserverTags:
      cmd = self.interactor.GetCommand(tag)
      cmd.SetAbortFlag(1)

  def cleanup(self):
    """clean up actors and observers"""
    for a in self.actors:
      self.renderer.RemoveActor2D(a)
    self.sliceView.scheduleRender()
    for tag in self.interactorObserverTags:
      self.interactor.RemoveObserver(tag)
    sliceNode = self.sliceLogic.GetSliceNode()
    for tag in self.sliceNodeTags:
      sliceNode.RemoveObserver(tag)

