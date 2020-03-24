import qt, vtk, slicer
from qt import QToolBar
import os
from slicer.util import VTKObservationMixin

import SmudgeModule
import ImportAtlas

class reducedToolbar(QToolBar, VTKObservationMixin):

  def __init__(self):

    QToolBar.__init__(self)
    VTKObservationMixin.__init__(self)

    self.parameterNode = SmudgeModule.SmudgeModuleLogic().getParameterNode()
    self.addObserver(self.parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateToolbarFromMRML)
    
    self.setWindowTitle(qt.QObject().tr("LeadDBS"))

    smw = slicer.util.mainWindow()
    interactionNode = slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLInteractionNode')
    layoutManager = slicer.app.layoutManager()
    affineNode = slicer.util.getNode(self.parameterNode.GetParameter("affineTransformID"))

    self.addObserver(interactionNode, interactionNode.InteractionModeChangedEvent, self.onInteractionModeChanged) 


    #
    # Layout
    #

    layoutFourUpAction = qt.QAction(smw)
    layoutFourUpAction.setIcon(qt.QIcon(qt.QPixmap(os.path.join(os.path.dirname(__file__),'Icons','LayoutFourUpView.png'))))
    layoutFourUpAction.connect('triggered(bool)', lambda t: layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpView))
    self.addAction(layoutFourUpAction)

    layoutTabbedAction = qt.QAction(smw)
    layoutTabbedAction.setIcon(qt.QIcon(qt.QPixmap(os.path.join(os.path.dirname(__file__),'Icons','LayoutTabbedSliceView.png'))))
    layoutTabbedAction.connect('triggered(bool)', lambda t: layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutTabbedSliceView))

    self.addAction(layoutTabbedAction)

    sliceIntersectionAction = qt.QAction(smw)
    sliceIntersectionAction.setIcon(qt.QIcon(qt.QPixmap(os.path.join(os.path.dirname(__file__),'Icons','SlicesCrosshair.png'))))
    sliceIntersectionAction.setCheckable(True)
    sliceIntersectionAction.connect('toggled(bool)', self.sliceIntersectionToggle)

    self.addAction(sliceIntersectionAction)
    
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
    self.windowLevelAction.setIcon(qt.QIcon(qt.QPixmap(os.path.join(os.path.dirname(__file__),'Icons','MouseWindowLevelMode.png'))))
    self.windowLevelAction.setCheckable(True)
    self.windowLevelAction.setMenu(windowLevelMenu)
    self.windowLevelAction.connect('toggled(bool)', lambda t: interactionNode.SetCurrentInteractionMode(5) if t else interactionNode.SetCurrentInteractionMode(2))

    self.addAction(self.windowLevelAction)

    #
    # Warp visible in slice view
    #

    warpViewAction = qt.QAction(smw)
    warpViewAction.setIcon(qt.QIcon(qt.QPixmap(os.path.join(os.path.dirname(__file__),'Icons','GlyphIcon.png'))))
    warpViewAction.setCheckable(True)
    warpViewAction.connect('toggled(bool)', lambda t: affineNode.GetDisplayNode().SetVisibility(t))
    warpViewAction.connect('toggled(bool)', lambda t: affineNode.GetDisplayNode().SetVisibility2D(t))
    self.addAction(warpViewAction)



    #
    # B <-> F slider
    #

    templateSlider = qt.QSlider(1)
    templateSlider.singleStep = 10
    templateSlider.minimum = 0
    templateSlider.maximum = 100
    templateSlider.value = 0
    templateSlider.setFixedWidth(120)
    templateSlider.connect('valueChanged(int)', lambda value: slicer.util.setSliceViewerLayers(foregroundOpacity = value / 100.0))
    self.addWidget(templateSlider)

    #
    # Atlas
    #

    self.atlasesComboBox = qt.QComboBox()
    self.atlasesComboBox.addItem('Import Atlas')
    self.atlasesComboBox.setMaximumWidth(350)
    self.atlasesComboBox.addItems(ImportAtlas.ImportAtlasLogic().getValidAtlases(self.parameterNode.GetParameter("MNIAtlasPath")))
    self.atlasesComboBox.connect('currentIndexChanged(int)', self.onAtlasChanged)
    self.atlasesComboBox.connect('activated(int)', lambda i: print(i))

    self.addWidget(self.atlasesComboBox)

    #
    # Space Separator
    #

    empty = qt.QWidget()
    empty.setSizePolicy(qt.QSizePolicy.Expanding,qt.QSizePolicy.Preferred)
    self.addWidget(empty)

    #
    # Subject
    #

    self.subjectNameLabel = qt.QLabel('Subject: ')    
    self.addWidget(self.subjectNameLabel)

    #
    # Save
    #
    self.saveButton = qt.QPushButton("Save and Next")
    self.saveButton.setFixedWidth(200)
    self.saveButton.setStyleSheet("background-color: green")
    #self.addWidget(self.saveButton)

    self.saveButton.connect("clicked(bool)", self.onSaveButton)

    #
    # Update
    #

    self.onAtlasChanged(1,'DISTAL Minimal (Ewert 2017)')
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

  def onAtlasChanged(self, index, atlasName = None):
    if not atlasName:
      atlasName = self.atlasesComboBox.itemText(index)
    if index != 0:
      atlasPath = os.path.join(self.parameterNode.GetParameter("MNIAtlasPath"), atlasName)
      folderID = ImportAtlas.ImportAtlasLogic().run(atlasPath)
    self.atlasesComboBox.setCurrentIndex(0)

  def onSaveButton(self):
    pass
    #self.SmudgeModuleWidget.exit()
    #SmudgeModuleLogic().applyChanges()

    # remove nodes
    #slicer.mrmlScene.RemoveNode(slicer.util.getNode(self.parameterNode.GetParameter("affineTransformID")))
    #slicer.mrmlScene.RemoveNode(slicer.util.getNode(self.parameterNode.GetParameter("warpID")))
    #layoutManager = slicer.app.layoutManager()
    #compositeNode = layoutManager.sliceWidget('Red').sliceLogic().GetSliceCompositeNode()
    #slicer.mrmlScene.RemoveNode(slicer.util.getNode(compositeNode.GetBackgroundVolumeID()))

    #nextSubjectN = int(self.parameterNode.GetParameter("subjectN"))+1
    #subjectPaths = self.parameterNode.GetParameter("subjectPaths").split(' ')
    
    #if nextSubjectN < len(subjectPaths):
    #  self.parameterNode.SetParameter("subjectN", str(nextSubjectN))
    #  self.parameterNode.SetParameter("subjectPath", subjectPaths[nextSubjectN])
    #  imageNode = SmudgeModuleLogic().loadSubjectData()
    #  SmudgeModuleLogic().initialize(imageNode)
    #  slicer.util.setSliceViewerLayers(background=imageNode.GetID())
    #  self.updateGuiFromMRML()
    #  self.toogleTools()
    #else:
    #  slicer.util.exit()



class reducedToolbarLogic(object):
  pass