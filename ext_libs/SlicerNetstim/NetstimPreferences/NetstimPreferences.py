import os
from abc import ABC
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin

#
# NetstimPreferences
#

class NetstimPreferences(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "NetstimPreferences Plugin"
    self.parent.categories = [""]
    self.parent.dependencies = [] 
    self.parent.contributors = ["Simon Oxenford (Charite Berlin)"]
    self.parent.helpText = ""
    self.parent.hidden = True

    # Additional initialization step after application startup is complete
    slicer.app.connect("startupCompleted()", setUpSettingsPanel)


def setUpSettingsPanel():
  if not slicer.app.commandOptions().noMainWindow:
    settingsPanel = NetstimPreferencesSettingsPanel()
    slicer.app.settingsDialog().addPanel("Lead-DBS", settingsPanel)


class NetstimPreferencesSettingsPanel(ctk.ctkSettingsPanel):
  def __init__(self):
    ctk.ctkSettingsPanel.__init__(self)
    self.ui = NetstimPreferencesSettingsUI(self)


class NetstimPreferencesSettingsUI:
  def __init__(self, parent):
      layout = qt.QFormLayout(parent)

      self.leadDBSPathButton = ctk.ctkDirectoryButton()
      self.leadDBSPathButton.setToolTip("Lead-DBS install directory")
      self.leadDBSPathButton.directoryChanged.connect(self.onLeadDBSPathChanged)
      layout.addRow("Lead-DBS Path: ", self.leadDBSPathButton)

      self.leadDBSSpaceComboBox = qt.QComboBox()
      self.leadDBSSpaceComboBox.setToolTip("Lead-DBS space")
      self.leadDBSSpaceComboBox.currentTextChanged.connect(lambda t: LeadDBSSpace().setValue(t))
      layout.addRow("Lead-DBS Space: ", self.leadDBSSpaceComboBox)

      self.useSmoothAtlasCheckBox = qt.QCheckBox()
      self.useSmoothAtlasCheckBox.checked = UseSmoothAtlas().getValue()
      self.useSmoothAtlasCheckBox.setToolTip("When checked, smoothed version will be used when loading atlases.")
      self.useSmoothAtlasCheckBox.connect("toggled(bool)", self.onUseSmoothAtlasCheckBoxToggled)
      layout.addRow("Use smooth atlases: ", self.useSmoothAtlasCheckBox)

      # initial set-up
      previousSpace = LeadDBSSpace().getValue()
      if previousSpace:
        self.leadDBSSpaceComboBox.addItems([previousSpace])
        self.leadDBSSpaceComboBox.setCurrentText(previousSpace)
      self.leadDBSPathButton.directory = LeadDBSPath().getValue()

  def onLeadDBSPathChanged(self):
    newDir = self.leadDBSPathButton.directory
    LeadDBSPath().setValue(newDir)
    if os.path.isfile(os.path.join(newDir,"lead.m")):
      import glob
      G = glob.glob(os.path.join(newDir, "templates", "space", "*", "atlases"))
      spaces = [os.path.basename(os.path.dirname(g)) for g in G]
      previousSpace = LeadDBSSpace().getValue()
      self.leadDBSSpaceComboBox.clear()
      self.leadDBSSpaceComboBox.addItems(spaces)
      if previousSpace in spaces:
        self.leadDBSSpaceComboBox.setCurrentText(previousSpace)
    else:
      qt.QMessageBox().warning(qt.QWidget(), "Error", "Invalid leaddbs path. Select leaddbs root install directory")

  def onUseSmoothAtlasCheckBoxToggled(self, checked):
    UseSmoothAtlas().setValue(checked)


class NetstimPreference(ABC):
  def __init__(self):
      self.rootKey = "NetstimPreferences/"
  
  def setValue(self, value):
    slicer.app.settings().setValue(self.rootKey + self.key, value)
  
  def getValue(self):
    return slicer.util.settingsValue(self.rootKey + self.key, self.default, converter=self.converter)

class LeadDBSPath(NetstimPreference):
  def __init__(self):
      super().__init__()
      self.key = "leadDBSPath"
      self.default = ""
      self.converter = str

class LeadDBSSpace(NetstimPreference):
  def __init__(self):
      super().__init__()
      self.key = "leadDBSSpace"
      self.default = "MNI152NLin2009bAsym"
      self.converter = str

class UseSmoothAtlas(NetstimPreference):
  def __init__(self):
      super().__init__()
      self.key = "useSmoothAtlas"
      self.default = True
      self.converter = slicer.util.toBool