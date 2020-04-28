
import qt, ctk, vtk, slicer
import os
from slicer.util import VTKObservationMixin


import SmudgeModule
import TransformsUtil
from . import WarpEffect


class WarpAbstractEffect(VTKObservationMixin):
  

  def __init__(self, name):
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
    return warpNode

  def sliceWidgets(self):
    for color in ['Red','Green','Yellow']:
      yield slicer.app.layoutManager().sliceWidget(color)


  def addEditButtonListeners(self, parent):
    for button in [parent.undoAllButton, parent.undoButton, parent.redoButton, parent.overwriteButton]:
      button.connect("pressed()", self.onEditButtonPressed)
      button.connect("released()", self.onEditButtonReleased)

  def onEditButtonPressed(self):
    WarpEffect.WarpEffectTool.empty()
    for sliceWidget in self.sliceWidgets():
      WarpEffect.NoneEffectTool(sliceWidget)

  def onEditButtonReleased(self):
    if self.effectButton.isChecked():
      self.effectButton.animateClick()

#
# None
#

class NoneEffectParameters(WarpAbstractEffect):

  def __init__(self):

    WarpAbstractEffect.__init__(self, 'None')

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

    WarpAbstractEffect.__init__(self, 'Smudge')

    # radius
    self.radiusSlider = ctk.ctkSliderWidget()
    self.radiusSlider.singleStep = 0.1
    self.radiusSlider.minimum = 5
    self.radiusSlider.maximum = float(self.parameterNode.GetParameter("maxRadius"))
    self.radiusSlider.decimals = 1
    self.radiusSlider.value = float(self.parameterNode.GetParameter("SmudgeRadius"))
    self.parametersFrame.layout().addRow("Radius (mm):", self.radiusSlider)

    # hardness
    self.hardnessSlider = ctk.ctkSliderWidget()
    self.hardnessSlider.singleStep = 1
    self.hardnessSlider.minimum = 0
    self.hardnessSlider.maximum = 100
    self.hardnessSlider.decimals = 0
    self.hardnessSlider.value = float(self.parameterNode.GetParameter("SmudgeHardness"))
    self.parametersFrame.layout().addRow("Hardness (%):", self.hardnessSlider)

    # force
    self.forceSlider = ctk.ctkSliderWidget()
    self.forceSlider.singleStep = 1
    self.forceSlider.minimum = 0
    self.forceSlider.maximum = 100
    self.forceSlider.decimals = 0
    self.forceSlider.value = float(self.parameterNode.GetParameter("SmudgeForce"))
    self.parametersFrame.layout().addRow("Force (%):", self.forceSlider)

    # post smoothing
    self.postSmoothingCheckBox = qt.QCheckBox('')
    self.postSmoothingCheckBox.setChecked(int(self.parameterNode.GetParameter("SmudgePostSmoothing")))
    self.parametersFrame.layout().addRow("Post Smoothing:", self.postSmoothingCheckBox)

    # post smoothing value
    self.postSmoothingSlider = ctk.ctkSliderWidget()
    self.postSmoothingSlider.singleStep = 1
    self.postSmoothingSlider.minimum = 0
    self.postSmoothingSlider.maximum = 100
    self.postSmoothingSlider.decimals = 0
    self.postSmoothingSlider.value = float(self.parameterNode.GetParameter("SmudgeSigma"))
    self.parametersFrame.layout().addRow("Sigma (% Radius):", self.postSmoothingSlider)

    # force

    self.radiusSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.hardnessSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.forceSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.postSmoothingCheckBox.connect('toggled(bool)', self.updateMRMLFromGUI)
    self.postSmoothingSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)

  def onEffectButtonClicked(self):
    super().onEffectButtonClicked()
    warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))
    size,origin,spacing = TransformsUtil.TransformsUtilLogic().getGridDefinition(warpNode)
    auxTransformNode = TransformsUtil.TransformsUtilLogic().emptyGridTransfrom(size, origin, spacing)
    for sliceWidget in self.sliceWidgets():
      WarpEffect.SmudgeEffectTool(sliceWidget, auxTransformNode)

  def updateGuiFromMRML(self, caller=None, event=None):
    super().updateGuiFromMRML(caller,event)
    radius = float(self.parameterNode.GetParameter("SmudgeRadius"))
    self.radiusSlider.setValue( radius )
    if radius < self.radiusSlider.minimum or radius > self.radiusSlider.maximum:
      self.updateMRMLFromGUI()
    self.hardnessSlider.setValue(float(self.parameterNode.GetParameter("SmudgeHardness")))
    self.forceSlider.setValue(float(self.parameterNode.GetParameter("SmudgeForce")))
    
  def updateMRMLFromGUI(self):
    self.parameterNode.SetParameter("SmudgeRadius", str(self.radiusSlider.value) )
    self.parameterNode.SetParameter("SmudgeHardness", str(self.hardnessSlider.value) )
    self.parameterNode.SetParameter("SmudgeForce", str(self.forceSlider.value) )
    self.parameterNode.SetParameter("SmudgePostSmoothing", str(int(self.postSmoothingCheckBox.isChecked())))
    self.parameterNode.SetParameter("SmudgeSigma", str(self.postSmoothingSlider.value) )

#
# Draw
#

class DrawEffectParameters(WarpAbstractEffect):

  def __init__(self):

    WarpAbstractEffect.__init__(self, 'Draw')

    self.spreadSlider = ctk.ctkSliderWidget()
    self.spreadSlider.singleStep = 0.1
    self.spreadSlider.minimum = 5
    self.spreadSlider.maximum = 50
    self.spreadSlider.decimals = 1
    self.spreadSlider.value = float(self.parameterNode.GetParameter("DrawSpread"))
    self.parametersFrame.layout().addRow("Spread (mm):", self.spreadSlider)

    self.spreadSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)

  def onEffectButtonClicked(self):
    super().onEffectButtonClicked()
    for sliceWidget in self.sliceWidgets():
      WarpEffect.SnapEffectTool(sliceWidget)

  def updateMRMLFromGUI(self):
    self.parameterNode.SetParameter("DrawSpread", str(self.spreadSlider.value))

#
# Blur
#


class BlurEffectParameters(WarpAbstractEffect):

  def __init__(self):

    WarpAbstractEffect.__init__(self, 'Blur')

    # select transform
    transformSelectFrame = qt.QFrame()
    transformSelectFrame.setLayout(qt.QHBoxLayout())

    self.userModificationsRadioButton = qt.QRadioButton('User Modifications')
    self.CompleteRadioButton = qt.QRadioButton('Original + User Modifications')
    transformSelectFrame.layout().addWidget(self.userModificationsRadioButton)
    transformSelectFrame.layout().addWidget(self.CompleteRadioButton)
    self.parametersFrame.layout().addRow("Warp:", transformSelectFrame)

    # sigma
    self.sigmaSlider = ctk.ctkSliderWidget()
    self.sigmaSlider.singleStep = 1
    self.sigmaSlider.minimum = 0
    self.sigmaSlider.maximum = 30
    self.sigmaSlider.decimals = 0
    self.sigmaSlider.value = float(self.parameterNode.GetParameter("BlurSigma"))
    self.parametersFrame.layout().addRow("Sigma (mm):", self.sigmaSlider)

    # radius
    self.radiusSlider = ctk.ctkSliderWidget()
    self.radiusSlider.singleStep = 0.1
    self.radiusSlider.minimum = 5
    self.radiusSlider.maximum = float(self.parameterNode.GetParameter("maxRadius"))
    self.radiusSlider.decimals = 1
    self.radiusSlider.value = float(self.parameterNode.GetParameter("BlurRadius"))
    self.parametersFrame.layout().addRow("Radius (mm):", self.radiusSlider)

    # hardness
    self.hardnessSlider = ctk.ctkSliderWidget()
    self.hardnessSlider.singleStep = 1
    self.hardnessSlider.minimum = 0
    self.hardnessSlider.maximum = 100
    self.hardnessSlider.decimals = 0
    self.hardnessSlider.value = float(self.parameterNode.GetParameter("BlurHardness"))
    self.parametersFrame.layout().addRow("Hardness (%):", self.hardnessSlider)


    self.radiusSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.hardnessSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.sigmaSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.userModificationsRadioButton.connect('clicked(bool)', self.onUserModificationsRadioButton)
    self.CompleteRadioButton.connect('clicked(bool)', self.onCompleteRadioButton)

  def onCompleteRadioButton(self):
    # flatten if needed
    warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))
    if TransformsUtil.TransformsUtilLogic().getNumberOfLayers(warpNode) > 1:
      TransformsUtil.TransformsUtilLogic().flattenTransform(warpNode, True)
    # get array
    transformArray = slicer.util.array(warpNode.GetID())
    # init
    for sliceWidget in self.sliceWidgets():
      WarpEffect.BlurEffectTool(sliceWidget, transformArray)

  def onUserModificationsRadioButton(self):
    self.parameterNode.SetParameter("lastOperation","undoall") # disable undo last operation
    # flatten if needed
    warpNode = slicer.util.getNode(self.parameterNode.GetParameter("warpID"))
    hasCorrectNumberOFLayers = TransformsUtil.TransformsUtilLogic().getNumberOfLayers(warpNode) == 2
    isCorrectTransformType = isinstance(warpNode.GetTransformFromParent().GetConcatenatedTransform(0), slicer.vtkOrientedGridTransform)
    if TransformsUtil.TransformsUtilLogic().getNumberOfLayers(warpNode) > 2:
      TransformsUtil.TransformsUtilLogic().flattenTransform(warpNode, False)
    elif not isinstance(warpNode.GetTransformFromParent().GetConcatenatedTransform(0), slicer.vtkOrientedGridTransform):
      # get last transform (drawing)
      transformID = TransformsUtil.TransformsUtilLogic().removeLastLayer(slicer.util.getNode(self.parameterNode.GetParameter("warpID")))
      # transform to grid transform
      size,origin,spacing = TransformsUtil.TransformsUtilLogic().getGridDefinition(warpNode)
      outNode = TransformsUtil.TransformsUtilLogic().transformToGridTransform(slicer.util.getNode(transformID), size, origin, spacing)
      # re-apply
      warpNode.SetAndObserveTransformNodeID(outNode.GetID())
      warpNode.HardenTransform()
      # remove aux nodes
      slicer.mrmlScene.RemoveNode(outNode)
      slicer.mrmlScene.RemoveNode(slicer.util.getNode(transformID))
    # get array
    transformArray = TransformsUtil.TransformsUtilLogic().arrayFromGeneralTransform(warpNode, 0)
    # init
    for sliceWidget in self.sliceWidgets():
      WarpEffect.BlurEffectTool(sliceWidget, transformArray)


  def resetButtons(self):
    # uncheck radio buttons
    self.userModificationsRadioButton.setAutoExclusive(False)
    self.CompleteRadioButton.setAutoExclusive(False)
    self.userModificationsRadioButton.setChecked(False)
    self.CompleteRadioButton.setChecked(False)
    self.userModificationsRadioButton.setAutoExclusive(False)
    self.CompleteRadioButton.setAutoExclusive(False)

  def onEffectButtonClicked(self):
    super().onEffectButtonClicked()
    self.resetButtons()

  def updateGuiFromMRML(self, caller=None, event=None):
    warpNode = super().updateGuiFromMRML(caller, event)
    warpNumberOfComponents = TransformsUtil.TransformsUtilLogic().getNumberOfLayers(warpNode)
    self.userModificationsRadioButton.enabled = warpNumberOfComponents > 1
    self.CompleteRadioButton.enabled = warpNumberOfComponents == 1
    radius = float(self.parameterNode.GetParameter("BlurRadius"))
    self.radiusSlider.setValue( radius )
    if radius < self.radiusSlider.minimum or radius > self.radiusSlider.maximum:
      self.updateMRMLFromGUI()
    self.hardnessSlider.setValue(float(self.parameterNode.GetParameter("BlurHardness")))
    self.sigmaSlider.setValue(float(self.parameterNode.GetParameter("BlurSigma")))


  def updateMRMLFromGUI(self):
    self.parameterNode.SetParameter("BlurRadius", str(self.radiusSlider.value) )
    self.parameterNode.SetParameter("BlurHardness", str(self.hardnessSlider.value) )
    self.parameterNode.SetParameter("BlurSigma", str(self.sigmaSlider.value) )
