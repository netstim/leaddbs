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

// AlphaOmega Logic includes
#include "vtkSlicerAlphaOmegaLogic.h"

// MRML includes
#include <vtkMRMLScene.h>
#include "vtkMRMLTransformNode.h"
#include "vtkMRMLLinearTransformNode.h"
#include "vtkMRMLTableNode.h"
#include "vtkMRMLScriptedModuleNode.h"
#include "vtkMRMLAlphaOmegaChannelNode.h"

// VTK includes
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkMatrix4x4.h>
#include <vtkTable.h>
#include <vtkVariantArray.h>

// STD includes
#include <cassert>
#include <cmath> // NAN

// AlphaOmega SDK
#include "AOSystemAPI.h"
#include "AOTypes.h"
#include "StreamFormat.h"

// Windows
#include <Windows.h>

// sys
#include <sys/types.h>
#include <sys/stat.h>


//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerAlphaOmegaLogic);

//----------------------------------------------------------------------------
vtkSlicerAlphaOmegaLogic::vtkSlicerAlphaOmegaLogic()
{
}

//----------------------------------------------------------------------------
vtkSlicerAlphaOmegaLogic::~vtkSlicerAlphaOmegaLogic()
{
}

//----------------------------------------------------------------------------
void vtkSlicerAlphaOmegaLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
void vtkSlicerAlphaOmegaLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//-----------------------------------------------------------------------------
void vtkSlicerAlphaOmegaLogic::RegisterNodes()
{
  assert(this->GetMRMLScene() != nullptr);
  vtkMRMLScene *scene = this->GetMRMLScene();
  scene->RegisterNodeClass(vtkSmartPointer<vtkMRMLAlphaOmegaChannelNode>::New());
}

//---------------------------------------------------------------------------
void vtkSlicerAlphaOmegaLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerAlphaOmegaLogic
::OnMRMLSceneNodeAdded(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerAlphaOmegaLogic
::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerAlphaOmegaLogic::sayHelloWorld()
{
  std::cout << "Hello, World!\n";
}



//-----------------------------------------------------------------------------
void vtkSlicerAlphaOmegaLogic::AOWaitForConnection()
{
  for (int i=0; i<10; i++)
  {
    Sleep(1000);
    if (this->AOIsConnected()==eAO_CONNECTED)
      break;
  }
}


//-----------------------------------------------------------------------------
int vtkSlicerAlphaOmegaLogic::AOIsConnected()
{
  return (isConnected() == eAO_CONNECTED);
}

//-----------------------------------------------------------------------------
int vtkSlicerAlphaOmegaLogic::AODefaultStartConnection(const char* address)
{
	MAC_ADDR sysMAC = {0};
	sscanf(address, "%x:%x:%x:%x:%x:%x",
		&sysMAC.addr[0], &sysMAC.addr[1], &sysMAC.addr[2], &sysMAC.addr[3], &sysMAC.addr[4], &sysMAC.addr[5]);

  DefaultStartConnection(&sysMAC , 0);
  this->AOWaitForConnection();
  return this->AOIsConnected();
}

//-----------------------------------------------------------------------------
int vtkSlicerAlphaOmegaLogic::AOCloseConnection()
{
  return CloseConnection();
}


//-----------------------------------------------------------------------------
float vtkSlicerAlphaOmegaLogic::AOGetDriveDepthMiliM()
{
	int32 nDepthUm = 0;
	EAOResult eAORes = (EAOResult)GetDriveDepth(&nDepthUm);

	if (eAORes == eAO_OK)
	{
    return nDepthUm / 1000.0;
  }
  else
  {
    return NAN;
  }
}

//-----------------------------------------------------------------------------
int vtkSlicerAlphaOmegaLogic::AOTranslateNameToID(const char* channelName)
{
  int ChannelID = 0;
  TranslateNameToID((cChar*)channelName, strlen(channelName), &ChannelID);
  return ChannelID;
}

//-----------------------------------------------------------------------------
int vtkSlicerAlphaOmegaLogic::AOAddBufferChannel(int ChannelID, int BufferSizeMSec)
{
  return AddBufferChannel(ChannelID, BufferSizeMSec);
}

//-----------------------------------------------------------------------------
int vtkSlicerAlphaOmegaLogic::AOAddBufferChannel(const char* channelName, int BufferSizeMSec)
{
  return AddBufferChannel(this->AOTranslateNameToID(channelName), BufferSizeMSec);
}

//-----------------------------------------------------------------------------
int vtkSlicerAlphaOmegaLogic::AOGetChannelData(int nChannelID, short* pData, int nData, int pDataCapture)
{
  return GetChannelData(nChannelID,  pData, nData, &pDataCapture);
}

//-----------------------------------------------------------------------------
int vtkSlicerAlphaOmegaLogic::AOGetChannelData(const char* channelName, short* pData, int nData, int pDataCapture)
{
  return GetChannelData(this->AOTranslateNameToID(channelName),  pData, nData, &pDataCapture);
}

//-----------------------------------------------------------------------------
int vtkSlicerAlphaOmegaLogic::AOGetAlignedData(short* pData, int nData, int* pDataCapture, int* pChannels, int nChannels, unsigned long* pBeginTS)
{
  return GetAlignedData(pData, nData, pDataCapture, pChannels, nChannels, pBeginTS);
}

//-----------------------------------------------------------------------------
int vtkSlicerAlphaOmegaLogic::AOClearBuffers()
{
  return ClearBuffers();
}

//
// Parameter Node
//

//-----------------------------------------------------------------------------
vtkMRMLScriptedModuleNode* vtkSlicerAlphaOmegaLogic::getParameterNode()
{
  vtkMRMLScriptedModuleNode* node = vtkMRMLScriptedModuleNode::SafeDownCast(this->GetMRMLScene()->GetSingletonNode(this->GetModuleName().c_str(), "vtkMRMLScriptedModuleNode"));
  if (node != nullptr)
  {
    return node;
  }
  else
  {
    this->createParameterNode();
    return this->getParameterNode();
  }
}

//-----------------------------------------------------------------------------
void vtkSlicerAlphaOmegaLogic::createParameterNode()
{
  vtkNew<vtkMRMLScriptedModuleNode> node;
  this->GetMRMLScene()->AddNode(node.GetPointer());
  node->SetSingletonTag(this->GetModuleName().c_str());
  node->SetModuleName(this->GetModuleName().c_str());
  node->SetNodeReferenceID("DistanceToTargetTransform", "");
}

//-----------------------------------------------------------------------------
std::string vtkSlicerAlphaOmegaLogic::GetModuleName()
{
  std::string moduleNameFull = this->GetClassName(); // get module name from class
  return moduleNameFull.substr(9, moduleNameFull.size() - 14); // remove "vtkSlicer" and "Logic"
}

