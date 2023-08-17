import qt, ctk, vtk, slicer
import os

import WarpDrive

import importlib

class AbstractToolWidget():
  
  toolEffect = []

  def __init__(self, name, toolTip):

    self.name = name

    # Tool Button
    effectPixmap = qt.QPixmap(os.path.join(os.path.split(WarpDrive.__file__)[0], 'Resources', 'Icons', self.name + '.png'))
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
    self.effectButton.setEnabled(True)
    self.effectButton.setSizePolicy(qt.QSizePolicy.MinimumExpanding,qt.QSizePolicy.Maximum)
    self.effectButton.connect('toggled(bool)', self.onEffectButtonToggle)
    self.effectButton.connect('clicked(bool)', self.onEffectButtonClicked)

  def onEffectButtonToggle(self):
    type(self).cleanEffects()

  def onEffectButtonClicked(self):
    # default interaction mode
    interactionNode = slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLInteractionNode')
    interactionNode.SetCurrentInteractionMode(interactionNode.ViewTransform)
    # get the respective effect class
    toolPkg = importlib.import_module('..Tools', 'WarpDriveLib.subpkg')
    toolMod = getattr(toolPkg, self.name + 'Tool')  
    toolCls = getattr(toolMod, self.name + 'ToolEffect')  
    # activate for each slice widget and save tool
    for sliceViewName in slicer.app.layoutManager().sliceViewNames():
      self.toolEffect.append(toolCls(slicer.app.layoutManager().sliceWidget(sliceViewName)))

  @classmethod
  def cleanEffects(cls):
    # cleanup created tools
    while cls.toolEffect:
      cls.toolEffect.pop().cleanup()
