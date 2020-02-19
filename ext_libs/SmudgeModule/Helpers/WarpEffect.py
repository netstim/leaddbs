import PointerEffect
import FunctionsUtil
import slicer, vtk
import numpy as np


class WarpEffectTool():

  _instances = set()

  def __init__(self):
    self._instances.add(self)

    # deactivate window level adjust
    interactorStyle = self.sliceView.sliceViewInteractorStyle()
    interactorStyle.SetActionEnabled(interactorStyle.AdjustWindowLevelBackground, False)

  def cleanup(self):
    interactorStyle = self.sliceView.sliceViewInteractorStyle()
    interactorStyle.SetActionEnabled(interactorStyle.AdjustWindowLevelBackground, True)
  
  @classmethod
  def empty(cls):
    # clean instances and reset
    for inst in cls._instances:
      inst.cleanup()
    cls._instances = set()


#
# SmudgeEffectTool
#

class SmudgeEffectTool(PointerEffect.CircleEffectTool, WarpEffectTool):

  def __init__(self, parameterNode, sliceWidget):

    self.parameterNode = parameterNode
    PointerEffect.CircleEffectTool.__init__(self, sliceWidget)
    WarpEffectTool.__init__(self)

    self.transformNode = slicer.util.getNode(self.parameterNode.GetParameter("transformID"))
    
    # transform data
    self.auxTransformNode = slicer.util.getNode(self.parameterNode.GetParameter("auxTransformID"))
    self.auxTransformSpacing = self.auxTransformNode.GetTransformFromParent().GetDisplacementGrid().GetSpacing()[0] # Asume isotropic!
    self.auxTransfromRASToIJK = FunctionsUtil.getTransformRASToIJK(self.auxTransformNode)  
    self.auxTransformArray = slicer.util.array(self.auxTransformNode.GetID())

    self.previousPoint = [0,0,0]   
    self.smudging = False

  def processEvent(self, caller=None, event=None):

    PointerEffect.CircleEffectTool.processEvent(self, caller, event)

    if event == 'LeftButtonPressEvent':
      # get aux transform array in case it has been modified externalyy
      self.auxTransformArray = slicer.util.array(self.auxTransformNode.GetID())
      # clean redo transform
      redoTransformID = self.parameterNode.GetParameter("redoTransformID")
      if redoTransformID != "":
        slicer.mrmlScene.RemoveNode(slicer.util.getNode(self.parameterNode.GetParameter("redoTransformID")))
        self.parameterNode.SetParameter("redoTransformID","")
      self.smudging = True
      xy = self.interactor.GetEventPosition()
      xyToRAS = self.sliceLogic.GetSliceNode().GetXYToRAS()
      self.previousPoint = xyToRAS.MultiplyDoublePoint( (xy[0], xy[1], 0, 1) )[0:3]
    elif event == 'LeftButtonReleaseEvent':
      self.smudging = False
      self.transformNode.HardenTransform()
      FunctionsUtil.emptyTransform(self.auxTransformNode)
      self.transformNode.SetAndObserveTransformNodeID(self.parameterNode.GetParameter("auxTransformID"))
      self.parameterNode.SetParameter("currentLayer",str(int(self.parameterNode.GetParameter("currentLayer"))+1))

    elif event == 'MouseMoveEvent':
      if self.smudging:
        # create a sphere with redius
        r = int(round(float(self.parameterNode.GetParameter("radius")) / self.auxTransformSpacing))
        xx, yy, zz = np.mgrid[:2*r+1, :2*r+1, :2*r+1]
        sphereResult = (xx-r) ** 2 + (yy-r) ** 2 + (zz-r) ** 2
        sphereResult[r][r][r] = 1 # replace 0 by 1
        sphereLarge = sphereResult <= (r**2+1) # sphere that the mouse shows
        sphereSmall = sphereResult <= ((r * (1-float(self.parameterNode.GetParameter("blurr")) / 100.0)) **2 + 1 ) # Blurr amount
        sphereResult = 1.0 / sphereResult # invert
        # get value in the edge of the small sphere
        i1,i2,i3 = np.nonzero(sphereSmall)
        newMaxValue = sphereResult[i1[0]][i2[0]][i3[0]]
        # set same value inside the small sphere
        sphereResult[sphereSmall] = sphereSmall[sphereSmall] * newMaxValue
        # delete outside values 
        sphereResult = sphereResult * sphereLarge
        # set range to [0-1]
        newMinValue = sphereResult.min()
        sphereResult = (sphereResult - newMinValue) / (newMaxValue - newMinValue)
        # set hardness
        sphereResult = sphereResult * float(self.parameterNode.GetParameter("hardness")) / 100.0
        # get current IJK coord
        xy = self.interactor.GetEventPosition()
        xyToRAS = self.sliceLogic.GetSliceNode().GetXYToRAS()
        currentPoint = xyToRAS.MultiplyDoublePoint( (xy[0], xy[1], 0, 1) )[0:3]
        pos_i,pos_j,pos_k,aux = self.auxTransfromRASToIJK.MultiplyDoublePoint(currentPoint + (1,))
        k,j,i = int(round(pos_k)),int(round(pos_j)),int(round(pos_i))
        curr_index = slice(k-r,k+r+1), slice(j-r,j+r+1), slice(i-r,i+r+1)

        # apply to transform array
        self.auxTransformArray[curr_index] += np.stack([(sphereResult) * i for i in (np.array(self.previousPoint) - np.array(currentPoint))],3) # original

        # update view
        self.auxTransformNode.Modified()
        # update previous point
        self.previousPoint = currentPoint

  def cleanup(self):
    WarpEffectTool.cleanup(self)
    PointerEffect.CircleEffectTool.cleanup(self)


#
# Snap Effect Tool
#

class SnapEffectTool(PointerEffect.DrawEffectTool, WarpEffectTool):

    
  def __init__(self, parameterNode, sliceWidget):

    self.parameterNode = parameterNode
    PointerEffect.DrawEffectTool.__init__(self,sliceWidget)
    WarpEffectTool.__init__(self)

    self.markupNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
    
  def processEvent(self, caller=None, event=None):

  
    if event == 'LeftButtonReleaseEvent':

      modelNode = slicer.util.getNode(self.parameterNode.GetParameter("modelID"))
      childrenCollection = vtk.vtkCollection()
      modelNode.GetChildrenModelNodes(childrenCollection)
      pd = childrenCollection.GetItemAsObject(0).GetPolyData()

      sliceIndex = int(np.nonzero([self.sliceLogic.GetSliceNode().GetName()==name for name in ['Yellow','Green','Red']])[0])
      sliceValue = self.sliceLogic.GetSliceNode().GetSliceToRAS().MultiplyDoublePoint([0,0,0,1])[sliceIndex]

      pointIndex = [i for i in range(0,pd.GetNumberOfPoints(),5) if (pd.GetPoint(i)[sliceIndex] > sliceValue-0.2 and pd.GetPoint(i)[sliceIndex] < sliceValue+0.2)]

      points = vtk.vtkPoints()
      for i in pointIndex:
        points.InsertNextPoint(pd.GetPoint(i))

      pd2 = vtk.vtkPolyData()
      pd2.SetPoints(points)

      pointsLocator = vtk.vtkPointLocator()
      pointsLocator.SetDataSet(pd2)
      pointsLocator.BuildLocator()

      # get closest point of initial and last position
      #      closestPointId0 = pointsLocator.FindClosestPoint(self.rasPoints.GetPoint(0))
      #      closestPointIdend = pointsLocator.FindClosestPoint(self.rasPoints.GetPoint(self.rasPoints.GetNumberOfPoints()-1))
      #
      #      # go through points from initial to end in the two possible directions
      #      points_dir1 = vtk.vtkPoints()
      #      points_dir1.InsertNextPoint(points.GetPoint(closestPointId0))
      #      points_dir2 = vtk.vtkPoints()
      #      points_dir2.InsertNextPoint(points.GetPoint(closestPointId0))
      #
      #      idl = vtk.vtkIdList()
      #      pointsLocator.FindClosestNPoints(3,points.GetPoint(closestPointId0),idl)
      #
      #      # closest point is the same one. Append the other two to the other points respectivly
      #      points_dir1.InsertNextPoint(points.GetPoint(idl.GetId(1)))
      #      points_dir2.InsertNextPoint(points.GetPoint(idl.GetId(2)))
      #
      #      # add points
      #      lastID = idl.GetId(1)
      #      IDs = [closestPointId0, idl.GetId(1)]
      #
      #      while lastID != closestPointIdend:
      #        pointsLocator.FindClosestNPoints(3,points.GetPoint(lastID),idl)
      #        for i in range(1,3):
      #          if idl.GetId(i) not in IDs:
      #            IDs.append(idl.GetId(i))
      #            lastID = idl.GetId(i)
      #            
      #
      #      #print(closestPointId0)
      #      #print(idl.GetId(0))
      #      #rint(idl.GetId(1))
      #      print(points_dir1)
      #
      self.markupNode.RemoveAllMarkups()
      for i in range(self.rasPoints.GetNumberOfPoints()):
        closestPointId = pointsLocator.FindClosestPoint(self.rasPoints.GetPoint(i))
        self.markupNode.AddFiducialFromArray(points.GetPoint(closestPointId))

      #closestPointId0 = pointsLocator.FindClosestPoint(self.rasPoints.GetPoint(0))
      #closestPointIdend = pointsLocator.FindClosestPoint(self.rasPoints.GetPoint(self.rasPoints.GetNumberOfPoints()-1))

      #print(str(closestPointId0) + " " + str(closestPointIdend))
    
      #print(sliceIndex)

    PointerEffect.DrawEffectTool.processEvent(self, caller, event) 
    
    

  def cleanup(self):
    slicer.mrmlScene.RemoveNode(self.markupNode)
    WarpEffectTool.cleanup(self)
    PointerEffect.DrawEffectTool.cleanup(self)


    
  

  