import sys
from SlicerNetstimLib.util import NORASubject
import slicer

subject = NORASubject(sys.argv[1])
subject.createAnatToFrameTransform()
slicer.mrmlScene.Clear()
subject.createORScene()

slicer.util.exit()
