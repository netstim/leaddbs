
import qt, ctk, vtk, slicer
import os
from slicer.util import VTKObservationMixin


import SmudgeModule
import TransformsUtil
from . import WarpEffect


class WarpAbstractEffect(VTKObservationMixin):
  

  def __init__(self, name, toolTip):
    VTKObservationMixin.__init__(self)
    self.parameterNode = SmudgeModule.SmudgeModuleLogic().getParameterNode()
    self.addObserver(self.parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGuiFromMRML)

    self.name = name

    # Tool Button
    effectPixmap = qt.QPixmap(os.path.join(os.path.split(__file__)[0] ,'Icons', self.name + '.png'))
    effectIcon = qt.QIcon(effectPixmap)
    self.effectButton = qt.QToolButton()
    self.effectButton.setToolButtonStyle(qt.Qt.ToolButtonTextUnderIcon)
    self.effectButton.setIcon(effectIcon)
    self.effectButton.setText(self.name)
    self.effectButton.setToolTip(toolTip)
    self.effectButton.setIconSize(effectPixmap.rect().size())
    self.effectButton.setAutoExclusive(True)
    self.effectButton.setCheckable(True)
    self.effectButton.setChecked(False)
    self.effectButton.setEnabled(False)
    self.effectButton.setSizePolicy(qt.QSizePolicy.MinimumExpanding,qt.QSizePolicy.Maximum)
    self.effectButton.connect('toggled(bool)', self.onEffectButtonToggle)
    self.effectButton.connect('clicked(bool)', self.onEffectButtonClicked)

    # Parameters Frame
    self.parametersFrame = qt.QFrame()
    self.parametersFrame.setLayout(qt.QFormLayout())
    self.parametersFrame.setVisible(False)

    # aux prev warp id
    self.prevWarpID = ""

  def onEffectButtonToggle(self):
    WarpEffect.WarpEffectTool.empty()
    self.parametersFrame.setVisible(self.effectButton.checked)

  def onEffectButtonClicked(self):
    self.parameterNode.SetParameter("currentEffect", self.name)
    interactionNode = slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLInteractionNode')
    interactionNode.SetCurrentInteractionMode(interactionNode.ViewTransform)

  def updateGuiFromMRML(self, caller=None, event=None):
    warpID = self.parameterNode.GetParameter("warpID")
    warpNode = slicer.util.getNode(warpID) if warpID != "" else None
    self.effectButton.enabled = bool(warpNode)
    if warpID != self.prevWarpID:
      self.resetEffect()
      self.prevWarpID = warpID
    return warpNode

  def sliceWidgets(self):
    for color in ['Red','Green','Yellow']:
      yield slicer.app.layoutManager().sliceWidget(color)


  def addEditButtonListeners(self, parent):
    for button in [parent.undoAllButton, parent.undoButton, parent.redoButton]:
      button.connect("pressed()", self.onEditButtonPressed)
      button.connect("released()", self.onEditButtonReleased)

  def onEditButtonPressed(self):
    WarpEffect.WarpEffectTool.empty()
    for sliceWidget in self.sliceWidgets():
      WarpEffect.NoneEffectTool(sliceWidget)

  def onEditButtonReleased(self):
    if self.effectButton.isChecked():
      self.effectButton.animateClick()

  def resetEffect(self):
    if self.effectButton.isChecked():
      WarpEffect.WarpEffectTool.empty()
      self.effectButton.animateClick()

#
# None
#

class NoneEffectParameters(WarpAbstractEffect):

  def __init__(self):

    toolTip = 'None. Disable Effects.'
    WarpAbstractEffect.__init__(self, 'None', toolTip)

    interactionNode = slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLInteractionNode')
    self.addObserver(interactionNode, interactionNode.InteractionModeChangedEvent, self.onInteractionModeChanged) 

  def onEffectButtonClicked(self):
    super().onEffectButtonClicked()
    for sliceWidget in self.sliceWidgets():
      WarpEffect.NoneEffectTool(sliceWidget)

  def updateGuiFromMRML(self, caller=None, event=None):
    self.effectButton.enabled = True

  def onInteractionModeChanged(self, caller=None, event=None):
    interactionNode = slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLInteractionNode')
    if interactionNode.GetCurrentInteractionMode() != interactionNode.ViewTransform:
      self.effectButton.setChecked(True)



#
# Smudge
#

class SmudgeEffectParameters(WarpAbstractEffect):

  def __init__(self):

    toolTip = 'Click and drag to deform the image'
    WarpAbstractEffect.__init__(self, 'Smudge', toolTip)

    # radius
    self.radiusSlider = ctk.ctkSliderWidget()
    self.radiusSlider.singleStep = 0.1
    self.radiusSlider.minimum = 5
    self.radiusSlider.maximum = float(self.parameterNode.GetParameter("maxRadius"))
    self.radiusSlider.decimals = 1
    self.radiusSlider.value = float(self.parameterNode.GetParameter("SmudgeRadius"))
    self.radiusSlider.setToolTip('Radius')
    self.parametersFrame.layout().addRow("Radius (mm):", self.radiusSlider)

    # hardness
    self.hardnessSlider = ctk.ctkSliderWidget()
    self.hardnessSlider.singleStep = 1
    self.hardnessSlider.minimum = 0
    self.hardnessSlider.maximum = 100
    self.hardnessSlider.decimals = 0
    self.hardnessSlider.value = float(self.parameterNode.GetParameter("SmudgeHardness"))
    self.hardnessSlider.setToolTip('Hardness')
    self.parametersFrame.layout().addRow("Hardness (%):", self.hardnessSlider)

    # force
    self.forceSlider = ctk.ctkSliderWidget()
    self.forceSlider.singleStep = 1
    self.forceSlider.minimum = 0
    self.forceSlider.maximum = 100
    self.forceSlider.decimals = 0
    self.forceSlider.value = float(self.parameterNode.GetParameter("SmudgeForce"))
    self.forceSlider.setToolTip('Force')
    self.parametersFrame.layout().addRow("Force (%):", self.forceSlider)

    # advanced
    advancedParametersGroupBox = ctk.ctkCollapsibleGroupBox()
    advancedParametersGroupBox.setTitle('Advanced')
    advancedParametersGroupBox.setLayout(qt.QFormLayout())
    advancedParametersGroupBox.collapsed = True
    self.parametersFrame.layout().addRow(advancedParametersGroupBox)
    
    # post smoothing
    self.postSmoothingCheckBox = qt.QCheckBox('')
    self.postSmoothingCheckBox.setChecked(int(self.parameterNode.GetParameter("SmudgePostSmoothing")))
    self.postSmoothingCheckBox.setToolTip('Enable to smooth the added warp once the operation finished.')
    advancedParametersGroupBox.layout().addRow("Post Smoothing:", self.postSmoothingCheckBox)

    # post smoothing value
    self.postSmoothingSlider = ctk.ctkSliderWidget()
    self.postSmoothingSlider.singleStep = 1
    self.postSmoothingSlider.minimum = 0
    self.postSmoothingSlider.maximum = 100
    self.postSmoothingSlider.decimals = 0
    self.postSmoothingSlider.value = float(self.parameterNode.GetParameter("SmudgeSigma"))
    self.postSmoothingSlider.setToolTip('Smoothing sigma as a percentage of the radius.')
    advancedParametersGroupBox.layout().addRow("Sigma (% Radius):", self.postSmoothingSlider)

    # expand edge
    self.expandGridSlider = ctk.ctkSliderWidget()
    self.expandGridSlider.singleStep = 1
    self.expandGridSlider.minimum = 0
    self.expandGridSlider.maximum = 100
    self.expandGridSlider.decimals = 0
    self.expandGridSlider.tracking = False
    self.expandGridSlider.value = float(self.parameterNode.GetParameter("expandGrid"))
    self.expandGridSlider.setToolTip('Expand the grid being modified. Useful when working with the edges of the image.')
    advancedParametersGroupBox.layout().addRow("Expand Grid (mm):", self.expandGridSlider) 

    # grid boounds
    self.gridBoundsCheckBox = qt.QCheckBox('')
    self.gridBoundsCheckBox.setChecked(False)
    self.gridBoundsCheckBox.setToolTip('Display grid bounds.')
    advancedParametersGroupBox.layout().addRow("Show grid bounds:", self.gridBoundsCheckBox)   


    self.radiusSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.hardnessSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.forceSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.postSmoothingCheckBox.connect('toggled(bool)', self.updateMRMLFromGUI)
    self.postSmoothingSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.expandGridSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.expandGridSlider.connect('valueChanged(double)', self.resetEffect)
    self.gridBoundsCheckBox.connect('toggled(bool)', self.onGridBoundsCheckBox)

  def onEffectButtonToggle(self):
    super().onEffectButtonToggle()
    if self.parameterNode.GetParameter("gridBoundsROIID") != "":
      gridBoundsROINode = slicer.util.getNode(self.parameterNode.GetParameter("gridBoundsROIID"))
      gridBoundsROINode.SetDisplayVisibility(self.gridBoundsCheckBox.checked and self.effectButton.checked)

  def onEffectButtonClicked(self):
    super().onEffectButtonClicked()
    warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))
    size,origin,spacing = self.getExpandedGrid()
    auxTransformNode = TransformsUtil.TransformsUtilLogic().emptyGridTransform(size, origin, spacing)
    for sliceWidget in self.sliceWidgets():
      WarpEffect.SmudgeEffectTool(sliceWidget, auxTransformNode)

  def updateGuiFromMRML(self, caller=None, event=None):
    warpNode = super().updateGuiFromMRML(caller,event)
    radius = float(self.parameterNode.GetParameter("SmudgeRadius"))
    self.radiusSlider.setValue( radius )
    if radius < self.radiusSlider.minimum or radius > self.radiusSlider.maximum:
      self.updateMRMLFromGUI()
    self.hardnessSlider.setValue(float(self.parameterNode.GetParameter("SmudgeHardness")))
    self.forceSlider.setValue(float(self.parameterNode.GetParameter("SmudgeForce")))
    self.postSmoothingSlider.setEnabled(self.postSmoothingCheckBox.checked)
    self.radiusSlider.maximum = float(self.parameterNode.GetParameter("maxRadius")) + float(self.parameterNode.GetParameter("expandGrid"))
    # roi bounds
    if warpNode and self.parameterNode.GetParameter("gridBoundsROIID") != "":
      gridBoundsROINode = slicer.util.getNode(self.parameterNode.GetParameter("gridBoundsROIID"))
      size,origin,spacing = TransformsUtil.TransformsUtilLogic().getGridDefinition(warpNode)
      sizeMM = [size[i]*spacing[i] for i in range(3)]
      gridBoundsROINode.SetRadiusXYZ([0.5*sizeMM[i]+float(self.parameterNode.GetParameter("expandGrid")) for i in range(3)])
    
  def updateMRMLFromGUI(self):
    self.parameterNode.SetParameter("SmudgeRadius", str(self.radiusSlider.value) )
    self.parameterNode.SetParameter("SmudgeHardness", str(self.hardnessSlider.value) )
    self.parameterNode.SetParameter("SmudgeForce", str(self.forceSlider.value) )
    self.parameterNode.SetParameter("SmudgePostSmoothing", str(int(self.postSmoothingCheckBox.isChecked())))
    self.parameterNode.SetParameter("SmudgeSigma", str(self.postSmoothingSlider.value) )
    self.parameterNode.SetParameter("expandGrid", str(self.expandGridSlider.value) )


  def getExpandedGrid(self):
    # create aux transform with same grid as warp
    warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))
    size,origin,spacing = TransformsUtil.TransformsUtilLogic().getGridDefinition(warpNode)
    # expand aux transform to deal with borders
    expandGrid = float(self.parameterNode.GetParameter("expandGrid"))
    origin = [o - expandGrid for o in origin]
    size = [int(round(s+expandGrid*2/spacing[0])) for s in size]
    return size, origin, spacing

  def onGridBoundsCheckBox(self):
    if self.parameterNode.GetParameter("gridBoundsROIID") == "" and self.parameterNode.GetParameter("warpID") != "":
      self.initROINode()
    gridBoundsROINode = slicer.util.getNode(self.parameterNode.GetParameter("gridBoundsROIID"))
    gridBoundsROINode.SetDisplayVisibility(self.gridBoundsCheckBox.checked)


  def initROINode(self):
    gridBoundsROINode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLAnnotationROINode')
    gridBoundsROINode.SetLocked(True)
    warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))
    size,origin,spacing = TransformsUtil.TransformsUtilLogic().getGridDefinition(warpNode)
    sizeMM = [size[i]*spacing[i] for i in range(3)]
    gridBoundsROINode.SetRadiusXYZ([0.5*sizeMM[i] for i in range(3)])
    gridBoundsROINode.SetXYZ([origin[i]+0.5*sizeMM[i] for i in range(3)])
    self.parameterNode.SetParameter("gridBoundsROIID", gridBoundsROINode.GetID())


#
# Draw
#

class DrawEffectParameters(WarpAbstractEffect):

  def __init__(self):

    toolTip = 'Draw a structure in the image to move it to the nearest visible model. When finished, corresponding fixed points will be added.'
    WarpAbstractEffect.__init__(self, 'Draw', toolTip)

    self.spreadSlider = ctk.ctkSliderWidget()
    self.spreadSlider.singleStep = 0.1
    self.spreadSlider.minimum = 5
    self.spreadSlider.maximum = 50
    self.spreadSlider.decimals = 1
    self.spreadSlider.value = float(self.parameterNode.GetParameter("DrawSpread"))
    self.spreadSlider.setToolTip('Specify up to how far away from the drawing the warp will be modified.')
    self.parametersFrame.layout().addRow("Spread (mm):", self.spreadSlider)

    self.spreadSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)

  def onEffectButtonClicked(self):
    super().onEffectButtonClicked()
    for sliceWidget in self.sliceWidgets():
      WarpEffect.SnapEffectTool(sliceWidget)

  def updateMRMLFromGUI(self):
    self.parameterNode.SetParameter("DrawSpread", str(self.spreadSlider.value))

#
# Smooth
#


class SmoothEffectParameters(WarpAbstractEffect):

  def __init__(self):

    toolTip = 'Smooth the warp field. Click and hold to preview, double-click to aply. Save current state to enable.'
    WarpAbstractEffect.__init__(self, 'Smooth', toolTip)

    # sigma
    self.sigmaSlider = ctk.ctkSliderWidget()
    self.sigmaSlider.singleStep = 1
    self.sigmaSlider.minimum = 0
    self.sigmaSlider.maximum = 30
    self.sigmaSlider.decimals = 0
    self.sigmaSlider.value = float(self.parameterNode.GetParameter("SmoothSigma"))
    self.sigmaSlider.setToolTip('Smoothing sigma')
    self.parametersFrame.layout().addRow("Sigma (mm):", self.sigmaSlider)

    # use radius / hole image
    self.useRadiusCheckbox = qt.QCheckBox()
    self.useRadiusCheckbox.setChecked(int(self.parameterNode.GetParameter("SmoothUseRadius")))
    self.useRadiusCheckbox.setToolTip('Use radius to focus editable area. Otherwise the hole warp will be modified.')
    self.parametersFrame.layout().addRow("Use Radius:", self.useRadiusCheckbox)

    # radius
    self.radiusSlider = ctk.ctkSliderWidget()
    self.radiusSlider.singleStep = 0.1
    self.radiusSlider.minimum = 5
    self.radiusSlider.maximum = float(self.parameterNode.GetParameter("maxRadius"))
    self.radiusSlider.decimals = 1
    self.radiusSlider.value = float(self.parameterNode.GetParameter("SmoothRadius"))
    self.radiusSlider.setToolTip('Radius')
    self.parametersFrame.layout().addRow("Radius (mm):", self.radiusSlider)

    # hardness
    self.hardnessSlider = ctk.ctkSliderWidget()
    self.hardnessSlider.singleStep = 1
    self.hardnessSlider.minimum = 0
    self.hardnessSlider.maximum = 100
    self.hardnessSlider.decimals = 0
    self.hardnessSlider.value = float(self.parameterNode.GetParameter("SmoothHardness"))
    self.hardnessSlider.setToolTip('Hardness')
    self.parametersFrame.layout().addRow("Hardness (%):", self.hardnessSlider)

    self.radiusSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.hardnessSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.sigmaSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.useRadiusCheckbox.connect('toggled(bool)', self.updateMRMLFromGUI)
    

  def onEffectButtonClicked(self):
    super().onEffectButtonClicked()
    for sliceWidget in self.sliceWidgets():
      WarpEffect.SmoothEffectTool(sliceWidget)

  def updateGuiFromMRML(self, caller=None, event=None):
    warpNode = super().updateGuiFromMRML(caller, event)
    warpNumberOfComponents = TransformsUtil.TransformsUtilLogic().getNumberOfLayers(warpNode)
    self.effectButton.enabled = warpNumberOfComponents == 1
    radius = float(self.parameterNode.GetParameter("SmoothRadius"))
    self.radiusSlider.setValue( radius )
    if radius < self.radiusSlider.minimum or radius > self.radiusSlider.maximum:
      self.updateMRMLFromGUI()
    self.hardnessSlider.setValue(float(self.parameterNode.GetParameter("SmoothHardness")))
    self.sigmaSlider.setValue(float(self.parameterNode.GetParameter("SmoothSigma")))
    # radius - hardness enables
    self.radiusSlider.setEnabled(self.useRadiusCheckbox.checked)
    self.hardnessSlider.setEnabled(self.useRadiusCheckbox.checked)


  def updateMRMLFromGUI(self):
    self.parameterNode.SetParameter("SmoothRadius", str(self.radiusSlider.value) )
    self.parameterNode.SetParameter("SmoothHardness", str(self.hardnessSlider.value) )
    self.parameterNode.SetParameter("SmoothSigma", str(self.sigmaSlider.value) )
    self.parameterNode.SetParameter("SmoothUseRadius", str(int(self.useRadiusCheckbox.checked)) )
