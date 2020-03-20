import qt, vtk, slicer
import os
from slicer.util import VTKObservationMixin

class reducedToolbar(VTKObservationMixin):

  def __init__(self, parameterNode):


    VTKObservationMixin.__init__(self)

    self.parameterNode = parameterNode
    self.addObserver(self.parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateToolbarFromMRML)
    
    self.toolbar = qt.QToolBar()

    smw = slicer.util.mainWindow()
    interactionNode = slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLInteractionNode')
    layoutManager = slicer.app.layoutManager()
    affineNode = slicer.util.getNode(parameterNode.GetParameter("affineTransformID"))


    #
    # Layout
    #

    layoutFourUpAction = qt.QAction(smw)
    layoutFourUpAction.setIcon(qt.QIcon(qt.QPixmap(os.path.join(os.path.dirname(__file__),'Icons','LayoutFourUpView.png'))))
    layoutFourUpAction.connect('triggered(bool)', lambda t: layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpView))
    self.toolbar.addAction(layoutFourUpAction)

    layoutTabbedAction = qt.QAction(smw)
    layoutTabbedAction.setIcon(qt.QIcon(qt.QPixmap(os.path.join(os.path.dirname(__file__),'Icons','LayoutTabbedSliceView.png'))))
    layoutTabbedAction.connect('triggered(bool)', lambda t: layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutTabbedSliceView))

    self.toolbar.addAction(layoutTabbedAction)

    sliceIntersectionAction = qt.QAction(smw)
    sliceIntersectionAction.setIcon(qt.QIcon(qt.QPixmap(os.path.join(os.path.dirname(__file__),'Icons','SlicesCrosshair.png'))))
    sliceIntersectionAction.setCheckable(True)
    sliceIntersectionAction.connect('toggled(bool)', self.sliceIntersectionToggle)

    self.toolbar.addAction(sliceIntersectionAction)
    
    #
    # Window Level
    #

    windowLevelModeActions = qt.QActionGroup(smw)
    windowLevelModeActions.setExclusive(True)

    windowLevelAdjustModeAction = qt.QAction(smw)
    windowLevelAdjustModeAction.setText('Adjust')
    windowLevelAdjustModeAction.setCheckable(True)
    windowLevelAdjustModeAction.setChecked(True)
    windowLevelAdjustModeAction.connect('triggered(bool)', lambda t: interactionNode.SetAttribute('AdjustWindowLevelMode', 'Adjust'))

    windowLevelRegionModeAction = qt.QAction(smw)
    windowLevelRegionModeAction.setText('Select Region')
    windowLevelRegionModeAction.setCheckable(True)
    windowLevelRegionModeAction.connect('triggered(bool)', lambda t: interactionNode.SetAttribute('AdjustWindowLevelMode', 'Rectangle'))

    windowLevelCenteredRegionModeAction = qt.QAction(smw)
    windowLevelCenteredRegionModeAction.setText('Select Region - centered')
    windowLevelCenteredRegionModeAction.setCheckable(True)
    windowLevelCenteredRegionModeAction.connect('triggered(bool)', lambda t: interactionNode.SetAttribute('AdjustWindowLevelMode', 'RectangleCentered'))

    windowLevelModeActions.addAction(windowLevelAdjustModeAction)
    windowLevelModeActions.addAction(windowLevelRegionModeAction)
    windowLevelModeActions.addAction(windowLevelCenteredRegionModeAction)

    windowLevelMenu = qt.QMenu(smw)
    windowLevelMenu.addActions(windowLevelModeActions.actions())
    
    self.windowLevelAction = qt.QAction(smw)
    self.windowLevelAction.setIcon(qt.QIcon(qt.QPixmap('/Users/simon/Slicer/Base/QTGUI/Resources/Icons/MouseWindowLevelMode.png')))
    self.windowLevelAction.setCheckable(True)
    self.windowLevelAction.setMenu(windowLevelMenu)
    self.windowLevelAction.connect('toggled(bool)', lambda t: interactionNode.SetCurrentInteractionMode(5) if t else interactionNode.SetCurrentInteractionMode(2))

    self.addObserver(interactionNode, interactionNode.InteractionModeChangedEvent, self.onInteractionModeChanged) 

    self.toolbar.addAction(self.windowLevelAction)

    #
    # Warp visible in slice view
    #

    warpViewAction = qt.QAction(smw)
    warpViewAction.setIcon(qt.QIcon(qt.QPixmap(os.path.join(os.path.dirname(__file__),'Icons','GlyphIcon.png'))))
    warpViewAction.setCheckable(True)
    warpViewAction.connect('toggled(bool)', lambda t: affineNode.GetDisplayNode().SetVisibility(t))
    warpViewAction.connect('toggled(bool)', lambda t: affineNode.GetDisplayNode().SetVisibility2D(t))
    self.toolbar.addAction(warpViewAction)

    #
    # Space
    #

    empty = qt.QWidget()
    empty.setSizePolicy(qt.QSizePolicy.Expanding,qt.QSizePolicy.Preferred)
    self.toolbar.addWidget(empty)

    #
    # Subject
    #

    self.subjectNameLabel = qt.QLabel('Subject: ')    
    self.toolbar.addWidget(self.subjectNameLabel)

    #
    # Update
    #

    self.updateToolbarFromMRML()


  def onInteractionModeChanged(self, caller, event):
    interactionNode = slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLInteractionNode')
    self.windowLevelAction.setChecked(interactionNode.GetCurrentInteractionMode() == 5)

  def sliceIntersectionToggle(self, t):
    compositeCollection = slicer.mrmlScene.GetNodesByClass("vtkMRMLSliceCompositeNode")
    for i in range(compositeCollection.GetNumberOfItems()):
      compositeCollection.GetItemAsObject(i).SetSliceIntersectionVisibility(t)


  def updateToolbarFromMRML(self, caller=None, event=None):
    subjectN = int(self.parameterNode.GetParameter("subjectN"))
    subjectPaths = self.parameterNode.GetParameter("subjectPaths").split(' ')
    self.subjectNameLabel.text = 'Subject: ' + os.path.split(subjectPaths[subjectN])[-1]