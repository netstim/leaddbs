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
#include <vtkMRMLPlotSeriesNode.h>

// VTK includes
#include <vtkMRMLTableNode.h>
#include <vtkMultiThreader.h>
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
#include <mutex>

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
  this->MultiThreader = vtkMultiThreader::New();
}

//----------------------------------------------------------------------------
vtkSlicerAlphaOmegaLogic::~vtkSlicerAlphaOmegaLogic()
{
  this->TerminateGatherAlignedData();
  if (this->MultiThreader)
  {
    this->MultiThreader->Delete();
  }
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
int vtkSlicerAlphaOmegaLogic::AOGetAlignedData(short* pData, int nData, int* pChannels, int nChannels, int pDataCapture)
{
  int status = eAO_MEM_EMPTY;
  unsigned long* pBeginTS;
  while (status == eAO_MEM_EMPTY || pDataCapture == 0)
  {
    status =  GetAlignedData(pData, nData, &pDataCapture, pChannels, nChannels, pBeginTS);
  }
  return status;
}

//-----------------------------------------------------------------------------
int vtkSlicerAlphaOmegaLogic::AOClearBuffers()
{
  return ClearBuffers();
}

//
// Get Aligned Data
//


//----------------------------------------------------------------------------
static void *vtkSlicerAlphaOmegaLogic_ThreadFunction(vtkMultiThreader::ThreadInfo *genericData )
{
  ThreadInfoStruct *info = static_cast<ThreadInfoStruct*>(genericData);
  vtkSlicerAlphaOmegaLogic *self = static_cast<vtkSlicerAlphaOmegaLogic *>(info->UserData);
  self->ContinuousGatherAlignedData();
  return nullptr;
}

//----------------------------------------------------------------------------
void vtkSlicerAlphaOmegaLogic::ThreadedGatherAlignedData()
{
  this->ThreadActiveLock.lock();
  this->ThreadActive = true;
  this->ThreadActiveLock.unlock();

  this->ThreadID = this->MultiThreader->SpawnThread(
                            (vtkThreadFunctionType) &vtkSlicerAlphaOmegaLogic_ThreadFunction,
                            static_cast<void *>(this));
}

//----------------------------------------------------------------------------
void vtkSlicerAlphaOmegaLogic::TerminateGatherAlignedData()
{
  this->ThreadActiveLock.lock();
  this->ThreadActive = false;
  this->ThreadActiveLock.unlock();

  if (this->ThreadID >= 0)
  {
    this->MultiThreader->TerminateThread(this->ThreadID);
    this->ThreadID = -1;
  }
  
}

//----------------------------------------------------------------------------
void vtkSlicerAlphaOmegaLogic::ContinuousGatherAlignedData()
{
  vtkMRMLTableNode* tableNode = vtkMRMLTableNode::SafeDownCast(this->getParameterNode()->GetNodeReference("AlignedDataTable"));

  int numberOfChannels = this->GetMRMLScene()->GetNumberOfNodesByClass("vtkMRMLAlphaOmegaChannelNode");
  int* channelsID = new int[numberOfChannels];
  vtkMRMLAlphaOmegaChannelNode* AOChannelNode = nullptr;
  vtkFloatArray * channelPreviewSignalArray = nullptr;

  for (int i=0; i<numberOfChannels; i++)
  {
    AOChannelNode = vtkMRMLAlphaOmegaChannelNode::SafeDownCast(this->GetMRMLScene()->GetNthNodeByClass(i,"vtkMRMLAlphaOmegaChannelNode"));
    channelsID[i] = AOChannelNode->GetChannelID();

    channelPreviewSignalArray =  vtkFloatArray::New();
    channelPreviewSignalArray->SetName(AOChannelNode->GetChannelName().c_str());
    tableNode->GetTable()->AddColumn(channelPreviewSignalArray);

    vtkMRMLPlotSeriesNode* channelPreviewPlotSeriesNode = vtkMRMLPlotSeriesNode::SafeDownCast(this->GetMRMLScene()->AddNewNodeByClass("vtkMRMLPlotSeriesNode",AOChannelNode->GetChannelName()));
    channelPreviewPlotSeriesNode->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);
    channelPreviewPlotSeriesNode->SetMarkerStyle(vtkMRMLPlotSeriesNode::MarkerStyleNone);
    channelPreviewPlotSeriesNode->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeLine);
    channelPreviewPlotSeriesNode->SetAndObserveTableNodeID(tableNode->GetID());
    channelPreviewPlotSeriesNode->SetYColumnName(tableNode->GetColumnName(i));
    channelPreviewPlotSeriesNode->SetColor(0,0,0);

    AOChannelNode->InitializeChannelBuffer();
    AOChannelNode->InitializeSaveFile();
  }

  int bufferSize = AOChannelNode->GetChannelBufferSizeMiliSeconds() / 1000.0 * AOChannelNode->GetChannelSamplingRate(); // TODO: diferent buffersize
  short* dataBuffer = new short[bufferSize * numberOfChannels];
  int dataCapture = 0;
  int dataForEachChannel = 0;

  int active = true;
  float previousDistanceToTarget = AOChannelNode->GetDriveDistanceToTarget();
  float* newDataArray = new float[bufferSize];
  unsigned int tableNodeIndex = 0;
  unsigned int j,k;
  unsigned int previewSamples = AOChannelNode->GetChannelPreviewLengthMiliSeconds() / 1000.0 * AOChannelNode->GetChannelSamplingRate();
  tableNode->GetTable()->SetNumberOfRows(previewSamples);

  this->AOClearBuffers();

  while (active)
  {
    // Check distance to target
    if (!isnan(AOChannelNode->GetDriveDistanceToTarget()) && previousDistanceToTarget != AOChannelNode->GetDriveDistanceToTarget())
    {
      for (int i=0; i<numberOfChannels; i++)
      {
        AOChannelNode = vtkMRMLAlphaOmegaChannelNode::SafeDownCast(this->GetMRMLScene()->GetNthNodeByClass(i,"vtkMRMLAlphaOmegaChannelNode"));
        AOChannelNode->CloseSaveFile();
        AOChannelNode->InitializeSaveFile();
      }
      previousDistanceToTarget = AOChannelNode->GetDriveDistanceToTarget();
    }

    // Gather Data
    if (this->AOIsConnected())
    {
      this->AOGetAlignedData(dataBuffer, bufferSize * numberOfChannels, channelsID, numberOfChannels, dataCapture);
    }
    else
    {
      int sleepTimeMiliS = 100;
      Sleep(sleepTimeMiliS);
      dataCapture = sleepTimeMiliS / 1000.0 * AOChannelNode->GetChannelSamplingRate() * numberOfChannels;
      for (int i=0; i<numberOfChannels; i++)
      {
        for(j=0; j<(dataCapture/numberOfChannels); j++)
        {
          dataBuffer[(i*(dataCapture/numberOfChannels))+j] = j;
        }
      }
    }

    dataForEachChannel = dataCapture / numberOfChannels;

    // Add Data
    for (int i=0; i<numberOfChannels; i++)
    {
      for(j=0; j<dataForEachChannel; j++)
      {
        tableNode->GetTable()->SetValue((tableNodeIndex+j)%previewSamples, i, dataBuffer[(i*dataForEachChannel)+j]);
      }
      for(k=0; k<0.1*previewSamples; k++)
      {
        tableNode->GetTable()->SetValue((tableNodeIndex+j+k)%previewSamples, i, NAN);
      }
    }

    // Add Data
    for (int i=0; i<numberOfChannels; i++)
    {
      AOChannelNode = vtkMRMLAlphaOmegaChannelNode::SafeDownCast(this->GetMRMLScene()->GetNthNodeByClass(i,"vtkMRMLAlphaOmegaChannelNode"));
      for(j=0; j<dataForEachChannel; j++)
      {
        newDataArray[j] = dataBuffer[(i*dataForEachChannel)+j] * AOChannelNode->GetChannelBitResolution() / AOChannelNode->GetChannelGain();
      }
      AOChannelNode->SetNewDataArraySize(dataForEachChannel);
      AOChannelNode->AppendNewDataToSaveFile(newDataArray);
    }

    tableNodeIndex = (tableNodeIndex+j)%previewSamples;

    // Check to see if we should be shutting down
    this->ThreadActiveLock.lock();
    active = this->ThreadActive;
    this->ThreadActiveLock.unlock();
  }

  for (int i=0; i<numberOfChannels; i++)
  {
    AOChannelNode = vtkMRMLAlphaOmegaChannelNode::SafeDownCast(this->GetMRMLScene()->GetNthNodeByClass(i,"vtkMRMLAlphaOmegaChannelNode"));
    AOChannelNode->CloseSaveFile();
  }
  channelPreviewSignalArray->Delete();
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
  node->SetNodeReferenceID("AlignedDataTable", "");
  node->SetParameter("AlignedRunning", "false");
}

//-----------------------------------------------------------------------------
std::string vtkSlicerAlphaOmegaLogic::GetModuleName()
{
  std::string moduleNameFull = this->GetClassName(); // get module name from class
  return moduleNameFull.substr(9, moduleNameFull.size() - 14); // remove "vtkSlicer" and "Logic"
}

