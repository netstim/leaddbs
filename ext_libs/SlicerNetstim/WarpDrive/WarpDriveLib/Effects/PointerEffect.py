import vtk, qt, slicer

from .Effect import AbstractEffect

import WarpDrive

class AbstractPointerEffect(AbstractEffect):

  def __init__(self, sliceWidget):
    AbstractEffect.__init__(self, sliceWidget)
    self.parameterNode = WarpDrive.WarpDriveLogic().getParameterNode()
    self.previousVisibleIDs = vtk.vtkIdList()
    self.previousBackgroundNodeID = None
    self.previousTransformNodeID = None

  def processEvent(self, caller=None, event=None):

    # if event == 'LeftButtonPressEvent' and self.parameterNode.GetParameter("currentEffect") != 'None':
    #   self.cursorOff()
    # elif event == 'LeftButtonReleaseEvent' and self.parameterNode.GetParameter("currentEffect") != 'None':
    #   qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
    #   qt.QApplication.processEvents()

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
          compositeNode = slicer.app.layoutManager().sliceWidget(slicer.app.layoutManager().sliceViewNames()[0]).sliceLogic().GetSliceCompositeNode()
          self.previousBackgroundNodeID = compositeNode.GetBackgroundVolumeID()
          self.previousForegroundOpacity = compositeNode.GetForegroundOpacity()
          slicer.util.setSliceViewerLayers(background = None)
      elif key.lower() == 'o':
        if not self.previousTransformNodeID:
          self.previousTransformNodeID = self.parameterNode.GetNodeReference("InputNode").GetTransformNodeID()
          self.parameterNode.GetNodeReference("InputNode").SetAndObserveTransformNodeID(None)

    elif event == "KeyReleaseEvent":
      key = self.interactor.GetKeySym()
      if key.lower() == 's':
        shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
        for i in range(self.previousVisibleIDs.GetNumberOfIds()):
          shNode.GetItemDataNode(self.previousVisibleIDs.GetId(i)).GetDisplayNode().SetVisibility(1)
        self.previousVisibleIDs.Reset()
      elif key.lower() == 't':
        slicer.util.setSliceViewerLayers(background = self.previousBackgroundNodeID, foregroundOpacity=self.previousForegroundOpacity)
        self.previousBackgroundNodeID = None
      elif key.lower() == 'o':
        self.parameterNode.GetNodeReference("InputNode").SetAndObserveTransformNodeID(self.previousTransformNodeID)
        self.previousTransformNodeID = None
      

  def setFiducialNodeAs(self, type, fromNode, name, radius):
    toNode = self.parameterNode.GetNodeReference(type + "Fiducial")
    for i in range(fromNode.GetNumberOfControlPoints()):
      toNode.AddControlPoint(vtk.vtkVector3d(fromNode.GetNthControlPointPosition(i)), name)
      toNode.SetNthControlPointDescription(toNode.GetNumberOfControlPoints()-1, radius)
    slicer.mrmlScene.RemoveNode(fromNode)