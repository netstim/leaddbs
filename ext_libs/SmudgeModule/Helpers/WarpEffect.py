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

      pointIndex = [i for i in range(0,pd.GetNumberOfPoints()) if (pd.GetPoint(i)[sliceIndex] > sliceValue-0.2 and pd.GetPoint(i)[sliceIndex] < sliceValue+0.2)]

      points = vtk.vtkPoints()
      for i in pointIndex:
        points.InsertNextPoint(pd.GetPoint(i))

      pd2 = vtk.vtkPolyData()
      pd2.SetPoints(points)

      pointsLocator = vtk.vtkPointLocator()
      pointsLocator.SetDataSet(pd2)
      pointsLocator.BuildLocator()


      # get closest point of initial and last position
      closestPointId0 = pointsLocator.FindClosestPoint(self.rasPoints.GetPoint(0))
      closestPointIdend = pointsLocator.FindClosestPoint(self.rasPoints.GetPoint(self.rasPoints.GetNumberOfPoints()-1))
      
      # sort indexes
      idlist=vtk.vtkIdList()
      sortedPointIndex = [closestPointId0]
      pid = closestPointId0
      while len(sortedPointIndex) != points.GetNumberOfPoints():
        n=1
        while pid in sortedPointIndex:
          n+=1
          pointsLocator.FindClosestNPoints(n, points.GetPoint(sortedPointIndex[-1]), idlist )
          pid = idlist.GetId(n-1)
        sortedPointIndex.append(pid)

      # prepare filter for comparison
      import vtkSlicerSegmentComparisonModuleLogicPython
      pdf=vtkSlicerSegmentComparisonModuleLogicPython.vtkPolyDataDistanceHistogramFilter()
      pdf.SetInputReferencePolyData(self.polyData)

      ## one way
      
      sortedPoints1 = vtk.vtkPoints()
      for i in range(sortedPointIndex.index(closestPointIdend)+1):
        sortedPoints1.InsertNextPoint(points.GetPoint(sortedPointIndex[i]))
      
      polyLine = vtk.vtkPolyLine()
      polyLine.GetPointIds().SetNumberOfIds(sortedPoints1.GetNumberOfPoints())
      for i in range(sortedPoints1.GetNumberOfPoints()):
        polyLine.GetPointIds().SetId(i, i)
      
      cells = vtk.vtkCellArray()
      cells.InsertNextCell(polyLine)
      
      polyData1 = vtk.vtkPolyData()
      polyData1.SetPoints(sortedPoints1)
      polyData1.SetLines(cells)

      pointsLocator = vtk.vtkPointLocator()
      pointsLocator.SetDataSet(polyData1)
      pointsLocator.BuildLocator()
      midpt = self.rasPoints.GetPoint(int(self.rasPoints.GetNumberOfPoints()/2))
      closestPointIdmidle = pointsLocator.FindClosestPoint(midpt)
      pt = sortedPoints1.GetPoint(closestPointIdmidle)
      d1 = np.sqrt(np.sum((np.array(midpt)-np.array(pt))**2, axis=0))
      #print(d1)

      
      #m=slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode')
      #m.SetAndObservePolyData(polyData1)
      #m.CreateDefaultDisplayNodes()
      #m.GetDisplayNode().SetColor(0 ,0,1)
      #m.GetDisplayNode().SetSliceIntersectionVisibility(True)
      
      
      ## other way
      
      sortedPoints2 = vtk.vtkPoints()
      
      for i in range(sortedPointIndex.index(closestPointIdend),len(sortedPointIndex)):
        sortedPoints2.InsertNextPoint(points.GetPoint(sortedPointIndex[i]))
      
      #sortedPoints2.InsertNextPoint(points.GetPoint(sortedPointIndex[0]))
      
      polyLine = vtk.vtkPolyLine()
      polyLine.GetPointIds().SetNumberOfIds(sortedPoints2.GetNumberOfPoints())
      for i in range(sortedPoints2.GetNumberOfPoints()):
        polyLine.GetPointIds().SetId(i, i)
      
      cells = vtk.vtkCellArray()
      cells.InsertNextCell(polyLine)
      
      polyData2 = vtk.vtkPolyData()
      polyData2.SetPoints(sortedPoints2)
      polyData2.SetLines(cells)

      pointsLocator = vtk.vtkPointLocator()
      pointsLocator.SetDataSet(polyData2)
      pointsLocator.BuildLocator()
      midpt = self.rasPoints.GetPoint(int(self.rasPoints.GetNumberOfPoints()/2))
      closestPointIdmidle = pointsLocator.FindClosestPoint(midpt)
      pt = sortedPoints2.GetPoint(closestPointIdmidle)
      d2 = np.sqrt(np.sum((np.array(midpt)-np.array(pt))**2, axis=0))

      #print(pdf2>pdf1)
      
      #m=slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode')
      #if d1<d2:
      #  m.SetAndObservePolyData(polyData1)
      #else:
      #  m.SetAndObservePolyData(polyData2)
      #m.CreateDefaultDisplayNodes()
      #m.GetDisplayNode().SetColor(1,0,0)
      #m.GetDisplayNode().SetSliceIntersectionVisibility(True)

      if d1==d2: # can happen when midpoint is the same as start/end point
        sp = sortedPoints1 if sortedPoints1.GetNumberOfPoints() < sortedPoints2.GetNumberOfPoints() else sortedPoints2
      else:
        sp = sortedPoints1 if d1 < d2 else sortedPoints2
      

      self.markupNode.RemoveAllMarkups()
      #n=slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
      for i in range(sp.GetNumberOfPoints()):
      #  closestPointId = pointsLocator.FindClosestPoint(self.rasPoints.GetPoint(i))
        self.markupNode.AddFiducialFromArray(sp.GetPoint(i))

      #closestPointId0 = pointsLocator.FindClosestPoint(self.rasPoints.GetPoint(0))
      #closestPointIdend = pointsLocator.FindClosestPoint(self.rasPoints.GetPoint(self.rasPoints.GetNumberOfPoints()-1))

      #print(str(closestPointId0) + " " + str(closestPointIdend))
    
      #print(sliceIndex)

    PointerEffect.DrawEffectTool.processEvent(self, caller, event) 
    
    

  def cleanup(self):
    slicer.mrmlScene.RemoveNode(self.markupNode)
    WarpEffectTool.cleanup(self)
    PointerEffect.DrawEffectTool.cleanup(self)


    
  

  