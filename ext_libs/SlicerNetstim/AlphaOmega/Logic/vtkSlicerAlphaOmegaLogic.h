/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// .NAME vtkSlicerAlphaOmegaLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the volumes


#ifndef __vtkSlicerAlphaOmegaLogic_h
#define __vtkSlicerAlphaOmegaLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"
#include <vtkVariantArray.h>


// MRML includes
class vtkMRMLTransformNode;
class vtkMRMLTableNode;
class vtkMRMLScriptedModuleNode;

// STD includes
#include <cstdlib>

#include "vtkSlicerAlphaOmegaModuleLogicExport.h"


/// \ingroup Slicer_QtModules_ExtensionTemplate
class VTK_SLICER_ALPHAOMEGA_MODULE_LOGIC_EXPORT vtkSlicerAlphaOmegaLogic :
  public vtkSlicerModuleLogic
{
public:

  static vtkSlicerAlphaOmegaLogic *New();
  vtkTypeMacro(vtkSlicerAlphaOmegaLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);

  void sayHelloWorld();
  void AOWaitForConnection();
  int AOIsConnected();
  int AODefaultStartConnection(const char* address);
  int AOCloseConnection();
  float AOGetDriveDepthMiliM();
  int AOTranslateNameToID(const char* channelName);
  int AOAddBufferChannel(int ChannelID, int BufferSizeMSec);
  int AOAddBufferChannel(const char* channelName, int BufferSizeMSec);
  int AOGetChannelData(int nChannelID, short* pData, int nData, int pDataCapture);
  int AOGetChannelData(const char* channelName, short* pData, int nData, int pDataCapture);
  int AOGetAlignedData(short* pData, int nData, int* pDataCapture, int* pChannels, int nChannels, unsigned long* pBeginTS);
  int AOClearBuffers();

  vtkMRMLScriptedModuleNode* getParameterNode();
  
protected:

  vtkSlicerAlphaOmegaLogic();
  virtual ~vtkSlicerAlphaOmegaLogic();
  virtual void SetMRMLSceneInternal(vtkMRMLScene* newScene);
  /// Register MRML Node classes to Scene. Gets called automatically when the MRMLScene is attached to this logic class.
  virtual void RegisterNodes() override;
  virtual void UpdateFromMRMLScene();
  virtual void OnMRMLSceneNodeAdded(vtkMRMLNode* node);
  virtual void OnMRMLSceneNodeRemoved(vtkMRMLNode* node);

  int TableIndex = 0;
  int StorageTableIndex = 0;

private:

  void createParameterNode();
  std::string GetModuleName();

  vtkSlicerAlphaOmegaLogic(const vtkSlicerAlphaOmegaLogic&); // Not implemented
  void operator=(const vtkSlicerAlphaOmegaLogic&); // Not implemented
};

#endif
