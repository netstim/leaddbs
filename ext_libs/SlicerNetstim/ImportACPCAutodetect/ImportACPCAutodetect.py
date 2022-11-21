import logging
import os

import vtk

import slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin


#
# ImportACPCAutodetect
#

class ImportACPCAutodetect(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "ImportACPCAutodetect"  # TODO: make this more human readable by adding spaces
        self.parent.categories = [""]
        self.parent.dependencies = [] 
        self.parent.contributors = ["Simon Oxenford (Charite Berlin)"]
        self.parent.helpText = ""
        self.parent.hidden = True


#
# ImportACPCAutodetectLogic
#

class ImportACPCAutodetectLogic(ScriptedLoadableModuleLogic):
    """This class should implement all the actual
    computation done by your module.  The interface
    should be such that other python code can import
    this class and make use of the functionality without
    requiring an instance of the Widget.
    Uses ScriptedLoadableModuleLogic base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self):
        """
        Called when the logic class is instantiated. Can be used for initializing member variables.
        """
        ScriptedLoadableModuleLogic.__init__(self)

    def setDefaultParameters(self, parameterNode):
        """
        Initialize parameter node with default settings.
        """
        pass

    def loadPointsAndCreateTransform(self, filePath):
        try:
            AC,PC,MSP = self.importPointsWithScipy(filePath)
        except:
            AC,PC,MSP = self.importPointsWithH5PY(filePath)

        midlineNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode', 'ACPC_autodetect')
        midlineNode.AddControlPoint(AC, 'AC')
        midlineNode.AddControlPoint(PC, 'PC')
        midlineNode.AddControlPoint(MSP, 'MSP')

        acpcNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsLineNode', 'acpc')
        acpcNode.GetDisplayNode().SetVisibility(False)
        acpcNode.AddControlPoint(AC, 'AC')
        acpcNode.AddControlPoint(PC, 'PC')
        
        outputNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLinearTransformNode', 'ACPC_transform')
        parameters = {}
        parameters['centerVolume'] = True
        parameters['OutputTransform'] =  outputNode.GetID()   
        parameters['ACPC'] = acpcNode.GetID()
        parameters['Midline'] = midlineNode.GetID()
        slicer.cli.run(slicer.modules.acpctransform, None, parameters, wait_for_completion=True, update_display=False)

        slicer.mrmlScene.RemoveNode(acpcNode)

        return [midlineNode.GetID(), outputNode.GetID()]

    def importPointsWithScipy(self, filePath):
        from scipy import io
        load = io.loadmat(filePath)
        AC = load['AC'].flatten()
        PC = load['PC'].flatten()
        MSP = load['MSP'].flatten()
        return AC,PC,MSP

    def importPointsWithH5PY(self, filePath):
        try:
            import h5py
        except:
            slicer.util.pip_install('h5py')
        import h5py    
        with h5py.File(filePath, 'r') as ACPCFile:
            AC = ACPCFile['AC'][:].flatten()
            PC = ACPCFile['PC'][:].flatten()
            MSP = ACPCFile['MSP'][:].flatten()
        return AC,PC,MSP


#
# File Reader
#

class ImportACPCAutodetectFileReader:

  def __init__(self, parent):
    self.parent = parent

  def description(self):
    return 'Lead-DBS ACPC'

  def fileType(self):
    return 'LeadDBSACPC'

  def extensions(self):
    return ['Lead-DBS ACPC (*.mat)']

  def canLoadFile(self, filePath):
    # filename must be atlas_index.mat
    filename = os.path.split(filePath)[1]
    return filename == "ACPC_autodetect.mat"

  def load(self, properties):
    try:

      # Import data
      filePath = properties['fileName']
      logic = ImportACPCAutodetectLogic()
      loadedNodeIDs = logic.loadPointsAndCreateTransform(filePath)

    except Exception as e:
      logging.error('Failed to load file: '+str(e))
      import traceback
      traceback.print_exc()
      return False

    self.parent.loadedNodes = loadedNodeIDs
    return True


#
# ImportACPCAutodetectTest
#

class ImportACPCAutodetectTest(ScriptedLoadableModuleTest):
    """
    This is the test case for your scripted module.
    Uses ScriptedLoadableModuleTest base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def setUp(self):
        """ Do whatever is needed to reset the state - typically a scene clear will be enough.
        """
        slicer.mrmlScene.Clear()

    def runTest(self):
        """Run as few or as many tests as needed here.
        """
        self.setUp()
        self.test_ImportACPCAutodetect1()

    def test_ImportACPCAutodetect1(self):
        """ Ideally you should have several levels of tests.  At the lowest level
        tests should exercise the functionality of the logic with different inputs
        (both valid and invalid).  At higher levels your tests should emulate the
        way the user would interact with your code and confirm that it still works
        the way you intended.
        One of the most important features of the tests is that it should alert other
        developers when their changes will have an impact on the behavior of your
        module.  For example, if a developer removes a feature that you depend on,
        your test should break so they know that the feature is needed.
        """

        self.delayDisplay("Starting the test")

        self.delayDisplay('Test passed')
