import vtk
import os
import slicer
import SimpleITK as sitk
import sitkUtils
import numpy as np
import glob
import sys
from scipy import io

try:
  import h5py
except:
  slicer.util.pip_install('h5py')
  import h5py

try:
  import hdf5storage
except:
  slicer.util.pip_install('hdf5storage')
  import hdf5storage


def saveApprovedData(subjectPath):
  approvedFile = os.path.join(subjectPath,'ea_coreg_approved.mat')
  matfiledata = {}
  if os.path.isfile(approvedFile):
    try:
      # read file and copy data except for glanat
      with h5py.File(approvedFile,'r') as f:
        for k in f.keys():
          if k != 'glanat':
            keyValue = f[k][()]
            matfiledata[k] = keyValue
      # now add approved glanat
      matfiledata[u'glanat'] = np.array([2])

    except: # use other reader for .mat file
      f = io.loadmat(approvedFile)
      for k in f.keys():
        if k != 'glanat':
          keyValue = f[k]
          matfiledata[k] = keyValue
      matfiledata['glanat'] = np.array([[2]],dtype='uint8')
      io.savemat(approvedFile,matfiledata)
      return

  else:
    matfiledata[u'glanat'] = np.array([2])

  # save
  # for some reason putting subject path into hdf5storage.write doesnt work
  currentDir = os.getcwd()
  os.chdir(subjectPath)
  hdf5storage.write(matfiledata, '.', 'ea_coreg_approved.mat', matlab_compatible=True)
  os.chdir(currentDir)
