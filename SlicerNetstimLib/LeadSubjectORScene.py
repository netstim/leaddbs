import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__))))
from SlicerNetstimLib.util import LeadDBSSubject
import slicer

subject = LeadDBSSubject(sys.argv[2], sys.argv[1])
subject.createORScene()

slicer.util.exit()
