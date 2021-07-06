
#include "vtkMRMLAlphaOmegaChannelNode.h"

// MRML includes
#include <vtkMRMLScene.h>

// VTK includes
#include <vtkNew.h>
#include <vtkCommand.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkMRMLTableNode.h>
#include <vtkMRMLPlotSeriesNode.h>
#include <vtkTable.h>
#include <vtksys/SystemTools.hxx>

// AlphaO2mega SDK
#include "AOSystemAPI.h"
#include "AOTypes.h"
#include "StreamFormat.h"

// STD includes
#include <cmath> // NAN
#include <algorithm> // replace

// Windows
#include <Windows.h>

// H5
#include "itk_H5Cpp.h"


static const int BUFFER_HEADER_SIZE = 7;
static const float MINIMUM_RECORDING_TIME_S = 4.0;

//------------------------------------------------------------------------------
vtkMRMLNodeNewMacro(vtkMRMLAlphaOmegaChannelNode);

//----------------------------------------------------------------------------
std::string vtkMRMLAlphaOmegaChannelNode::ChannelRootSavePath = "";

//----------------------------------------------------------------------------
float vtkMRMLAlphaOmegaChannelNode::DriveDistanceToTarget = NAN;

//----------------------------------------------------------------------------
std::mutex vtkMRMLAlphaOmegaChannelNode::H5Busy;

//----------------------------------------------------------------------------
vtkMRMLAlphaOmegaChannelNode::vtkMRMLAlphaOmegaChannelNode()
{
  this->SetHideFromEditors(false);
  this->MultiThreader = vtkMultiThreader::New();
  // Defaults for SPK Channel
  this->ChannelSamplingRate = 44000;
  this->ChannelGain = 20;
  this->ChannelBitResolution = 38.147;
  // Preview
  this->ChannelPreviewLengthMiliSeconds = 1000;
  this->ChannelPreviewTimeArray =  vtkFloatArray::New();
  this->ChannelPreviewTimeArray->SetName(const_cast<char *>("t"));
  this->ChannelPreviewSignalArray =  vtkFloatArray::New();
  this->ChannelPreviewSignalArray->SetName(const_cast<char *>("x"));
  this->ChannelPreviewTableNode = nullptr;
  this->ChannelPreviewPlotSeriesNode = nullptr;
}

//----------------------------------------------------------------------------
vtkMRMLAlphaOmegaChannelNode::~vtkMRMLAlphaOmegaChannelNode()
{
  if (this->MultiThreader)
  {
    this->TerminateGatherData();
    this->MultiThreader->Delete();
  }
  if (this->ChannelPreviewTimeArray)
  {
    this->ChannelPreviewTimeArray->Delete();
  }
  if (this->ChannelPreviewSignalArray)
  {
    this->ChannelPreviewSignalArray->Delete();
  }
}

//----------------------------------------------------------------------------
void vtkMRMLAlphaOmegaChannelNode::WriteXML(ostream& of, int nIndent)
{
  Superclass::WriteXML(of, nIndent);
  
  of <<" ChannelName=\"" << this->ChannelName << "\"";
  of <<" ChannelID=\"" << this->ChannelID << "\"";
  of <<" ChannelSamplingRate=\"" << this->ChannelSamplingRate << "\"";
  of <<" ChannelGain=\"" << this->ChannelGain << "\"";
  of <<" ChannelBitResolution=\"" << this->ChannelBitResolution << "\"";
  of <<" ChannelBufferSizeMiliSeconds=\"" << this->ChannelBufferSizeMiliSeconds << "\"";
  of <<" ChannelPreviewLengthMiliSeconds=\"" << this->ChannelPreviewLengthMiliSeconds << "\"";


}

//----------------------------------------------------------------------------
void vtkMRMLAlphaOmegaChannelNode::ReadXMLAttributes(const char** atts)
{
  vtkMRMLNode::ReadXMLAttributes(atts);

  const char* attName;
  const char* attValue;
  while (*atts != nullptr)
  {
    if (!strcmp(attName, "ChannelName"))
      {
      std::stringstream ss;
      ss << attValue;
      ss >> this->ChannelName;
      }
    else if (!strcmp(attName, "ChannelID"))
      {
      std::stringstream ss;
      ss << attValue;
      ss >> this->ChannelID;
      }
    else if (!strcmp(attName, "ChannelSamplingRate"))
      {
      std::stringstream ss;
      ss << attValue;
      ss >> this->ChannelSamplingRate;
      }      
    else if (!strcmp(attName, "ChannelSamplingRate"))
      {
      std::stringstream ss;
      ss << attValue;
      ss >> this->ChannelSamplingRate;
      }
    else if (!strcmp(attName, "ChannelBufferSizeMiliSeconds"))
      {
      std::stringstream ss;
      ss << attValue;
      ss >> this->ChannelBufferSizeMiliSeconds;
      }
    else if (!strcmp(attName, "ChannelPreviewLengthMiliSeconds"))
      {
      std::stringstream ss;
      ss << attValue;
      ss >> this->ChannelPreviewLengthMiliSeconds;
      }          
  }

}

//----------------------------------------------------------------------------
// Copy the node's attributes to this object.
// Does NOT copy: ID, FilePrefix, Name, VolumeID
void vtkMRMLAlphaOmegaChannelNode::Copy(vtkMRMLNode *anode)
{
  // TODO
  vtkMRMLAlphaOmegaChannelNode* node = vtkMRMLAlphaOmegaChannelNode::SafeDownCast(anode);
  if (!node)
    {
    vtkErrorMacro("Node copy failed");
    return;
    }

  int wasModified = this->StartModify();
  Superclass::Copy(anode);

  this->EndModify(wasModified);
}

//----------------------------------------------------------------------------
void vtkMRMLAlphaOmegaChannelNode::PrintSelf(ostream& os, vtkIndent indent)
{
  // TODO
  this->Superclass::PrintSelf(os, indent);

}


//----------------------------------------------------------------------------
// Clear all data from buffers
int vtkMRMLAlphaOmegaChannelNode::ClearAllChannelsBuffers()
{
  return ClearBuffers();
}

//----------------------------------------------------------------------------
// Get the names of all available channels
std::vector<std::string> vtkMRMLAlphaOmegaChannelNode::GetAllChannelsNames()
{
  // initialize output
  std::vector<std::string> channelsNames;

  // get channel count
	uint32 uChannelsCount = 0;
	GetChannelsCount(&uChannelsCount);

  // get information
  SInformation *pChannelsInfo = new SInformation[uChannelsCount];
	GetAllChannels(pChannelsInfo, uChannelsCount);

  channelsNames.push_back("");

  // fill vector. in case of error the uChannelCount = 0 and the vector is empty
  for (int i = 0; i<uChannelsCount; i++)
  {
      channelsNames.push_back(pChannelsInfo[i].channelName);
  }

  channelsNames.push_back("Test");
  channelsNames.push_back("Test2");
  channelsNames.push_back("Test3");

  return channelsNames;
}

//----------------------------------------------------------------------------
// Set Channel name and automatically set ID from SDK
int vtkMRMLAlphaOmegaChannelNode::SetChannelNameAndID(const char* name)
{
  this->ChannelName = name;
  return TranslateNameToID((cChar*)this->ChannelName.c_str(), this->ChannelName.length(), &this->ChannelID);
}

//----------------------------------------------------------------------------
void vtkMRMLAlphaOmegaChannelNode::SetChannelPreviewTableNode(vtkMRMLTableNode* tableNode)
{
  this->ChannelPreviewTableNode = tableNode;
  if (this->ChannelPreviewTableNode == nullptr)
  {
    return;
  }
  else
  {
    this->ChannelPreviewTableNode->Reset(nullptr);
    this->ChannelPreviewTableNode->GetTable()->AddColumn(this->ChannelPreviewTimeArray);
    this->ChannelPreviewTableNode->GetTable()->AddColumn(this->ChannelPreviewSignalArray);
    this->ChannelPreviewTableNode->Modified();
  }
  if (this->ChannelPreviewPlotSeriesNode == nullptr)
  {
    this->ChannelPreviewPlotSeriesNode = vtkMRMLPlotSeriesNode::SafeDownCast(this->GetScene()->AddNewNodeByClass("vtkMRMLPlotSeriesNode",this->GetChannelName()));
    this->ChannelPreviewPlotSeriesNode->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeScatter);
    this->ChannelPreviewPlotSeriesNode->SetMarkerStyle(vtkMRMLPlotSeriesNode::MarkerStyleNone);
    this->ChannelPreviewPlotSeriesNode->SetPlotType(vtkMRMLPlotSeriesNode::PlotTypeLine);
  }
  this->ChannelPreviewPlotSeriesNode->SetAndObserveTableNodeID(this->ChannelPreviewTableNode->GetID());
  this->ChannelPreviewPlotSeriesNode->SetXColumnName(this->ChannelPreviewTableNode->GetColumnName(0));
  this->ChannelPreviewPlotSeriesNode->SetYColumnName(this->ChannelPreviewTableNode->GetColumnName(1));
  this->ChannelPreviewPlotSeriesNode->SetColor(0,0,0);
}

//----------------------------------------------------------------------------
int vtkMRMLAlphaOmegaChannelNode::InitializeChannelBuffer()
{
  this->ChannelBufferSizeSamples = this->ChannelBufferSizeMiliSeconds / 1000.0 * this->ChannelSamplingRate;
  this->ChannelBufferSizeSamples += BUFFER_HEADER_SIZE;
  this->DataBuffer = new short[this->ChannelBufferSizeSamples];
  return AddBufferChannel(this->ChannelID, this->ChannelBufferSizeMiliSeconds);
}


//----------------------------------------------------------------------------
void vtkMRMLAlphaOmegaChannelNode::InitializeChannelPreview()
{
  this->ChannelPreviewLengthSamples = this->ChannelPreviewLengthMiliSeconds / 1000.0 * this->ChannelSamplingRate;
  // Initialize arrays
  this->ChannelPreviewTimeArray->SetNumberOfValues(this->ChannelPreviewLengthSamples);
  this->ChannelPreviewSignalArray->SetNumberOfValues(this->ChannelPreviewLengthSamples);
  for (int i=0; i<this->ChannelPreviewLengthSamples; i++)
  {
    this->ChannelPreviewTimeArray->SetValue(i, i/(float)this->ChannelSamplingRate);
    this->ChannelPreviewSignalArray->SetValue(i, NAN);
  }
  // Update table
  if (this->ChannelPreviewTableNode != nullptr)
  {
    this->SetChannelPreviewTableNode(this->ChannelPreviewTableNode);
  }
}


//----------------------------------------------------------------------------
void vtkMRMLAlphaOmegaChannelNode::InitializeSaveFile()
{
  vtkMRMLAlphaOmegaChannelNode::H5Busy.lock();

  std::string modifiedName = this->ChannelName;
  std::replace(modifiedName.begin(), modifiedName.end(), '/', '_');
  std::replace(modifiedName.begin(), modifiedName.end(), ' ', '_');

  std::vector<std::string> filesVector;
  filesVector.push_back(this->ChannelRootSavePath);
  filesVector.emplace_back("/");
  filesVector.push_back(modifiedName);
  this->ChannelFullSavePath = vtksys::SystemTools::JoinPath(filesVector);
  vtksys::SystemTools::MakeDirectory(this->ChannelFullSavePath.c_str());

  filesVector.emplace_back("/");
  filesVector.push_back(vtksys::SystemTools::GetCurrentDateTime("%H%M%S_") + std::to_string(this->DriveDistanceToTarget) + ".h5");
  this->ChannelFullSavePath = vtksys::SystemTools::JoinPath(filesVector);

  const hsize_t ndims = 1;

  // create file
  hid_t plist1 = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(plist1, H5F_CLOSE_STRONG);
  H5Pset_libver_bounds(plist1, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
  this->H5File = new H5::H5File(this->ChannelFullSavePath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist1);

  // save sampling rate
  H5::DataSpace samplingRateSpace(1, &ndims);
  H5::DataSet   samplingRateSet = this->H5File->createDataSet("sr", H5::PredType::NATIVE_INT, samplingRateSpace);
  samplingRateSet.write(&this->ChannelSamplingRate, H5::PredType::NATIVE_INT);
  samplingRateSet.close();

  // create dataspace
  hsize_t dims[ndims] = {0};
  hsize_t max_dims[ndims] = {H5S_UNLIMITED};
  hid_t file_space = H5Screate_simple(ndims, dims, max_dims);

  // create property list
  hid_t plist2 = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_layout(plist2, H5D_CHUNKED);
  hsize_t chunk_dims[ndims] = {1};
  H5Pset_chunk(plist2, ndims, chunk_dims);

  // create dataset
  H5Dcreate(this->H5File->getId(), "data", H5T_NATIVE_FLOAT, file_space, H5P_DEFAULT, plist2, H5P_DEFAULT);

  // initialize memory dataspace
  dims[0] = this->DataCapture;
  this->H5MemoryDataspace = H5Screate_simple(ndims, dims, NULL);

  // Single-Writer/Multiple-Reader
  H5Fstart_swmr_write(this->H5File->getId());

  // close
  H5Pclose(plist1);
  H5Pclose(plist2);
  H5Sclose(file_space);

  vtkMRMLAlphaOmegaChannelNode::H5Busy.unlock();
}

void vtkMRMLAlphaOmegaChannelNode::CloseSaveFile()
{
  vtkMRMLAlphaOmegaChannelNode::H5Busy.lock();
  
  H5Sclose(this->H5MemoryDataspace);
  this->H5File->close();
  
  vtkMRMLAlphaOmegaChannelNode::H5Busy.unlock();
}


//----------------------------------------------------------------------------
static void *vtkMRMLAlphaOmegaChannelNode_ThreadFunction(vtkMultiThreader::ThreadInfo *genericData )
{
  ThreadInfoStruct *info = static_cast<ThreadInfoStruct*>(genericData);
  vtkMRMLAlphaOmegaChannelNode *self = static_cast<vtkMRMLAlphaOmegaChannelNode *>(info->UserData);
  self->ContinuousGatherData();
  return nullptr;
}

//----------------------------------------------------------------------------
void vtkMRMLAlphaOmegaChannelNode::ThreadedGatherData()
{
  this->GatheringData = true;

  this->ThreadActiveLock.lock();
  this->ThreadActive = true;
  this->ThreadActiveLock.unlock();

  this->ThreadID = this->MultiThreader->SpawnThread(
                            (vtkThreadFunctionType) &vtkMRMLAlphaOmegaChannelNode_ThreadFunction,
                            static_cast<void *>(this));
}

//----------------------------------------------------------------------------
void vtkMRMLAlphaOmegaChannelNode::TerminateGatherData()
{
  this->GatheringData = false;

  // Signal the Thread that we are terminating.
  this->ThreadActiveLock.lock();
  this->ThreadActive = false;
  this->ThreadActiveLock.unlock();

  if (this->ThreadID >= 0)
  {
    this->MultiThreader->TerminateThread(this->ThreadID);
    this->ThreadID = -1;
    this->CloseSaveFile();
  }
  
}

//----------------------------------------------------------------------------
void vtkMRMLAlphaOmegaChannelNode::ContinuousGatherData()
{
  this->InitializeChannelBuffer();
  this->InitializeChannelPreview();
  this->InitializeSaveFile();

  this->ClearChannelBuffer();

  int active = true;
  float previousDistanceToTarget = this->DriveDistanceToTarget;

  while (active)
  {
    if (!isnan(this->DriveDistanceToTarget) && previousDistanceToTarget != this->DriveDistanceToTarget)
    {
      this->CloseSaveFile();
      this->InitializeSaveFile();
      previousDistanceToTarget = this->DriveDistanceToTarget;
    }

    this->GatherData();
    float* newDataArray = this->CreateNewDataArray();
    this->AppendNewDataToPreviewArray(newDataArray);
    this->AppendNewDataToSaveFile(newDataArray);
    
    // Check to see if we should be shutting down
    this->ThreadActiveLock.lock();
    active = this->ThreadActive;
    this->ThreadActiveLock.unlock();
  }
}

//----------------------------------------------------------------------------
void vtkMRMLAlphaOmegaChannelNode::ClearChannelBuffer()
{
  while (eAO_OK == GetChannelData(this->ChannelID, this->DataBuffer, this->ChannelBufferSizeSamples, &this->DataCapture) && this->DataCapture != 0){}
}

//----------------------------------------------------------------------------
void vtkMRMLAlphaOmegaChannelNode::GatherData()
{
  if (isConnected() == eAO_CONNECTED)
  {
    this->WaitAndGetChannelData();
  }
  else
  {
    this->SleepAndSimulateData();
  }
  this->NewDataArraySize = this->DataCapture-BUFFER_HEADER_SIZE;
}

//----------------------------------------------------------------------------
int vtkMRMLAlphaOmegaChannelNode::WaitAndGetChannelData()
{
  int status = eAO_MEM_EMPTY;
  while (status == eAO_MEM_EMPTY || this->DataCapture == 0)
  {
    status = GetChannelData(this->ChannelID, this->DataBuffer, this->ChannelBufferSizeSamples, &this->DataCapture);
  }
  return status;
}

//----------------------------------------------------------------------------
void vtkMRMLAlphaOmegaChannelNode::SleepAndSimulateData()
{
  int sleepTimeMiliS = 100;
  Sleep(sleepTimeMiliS);
  this->DataCapture = BUFFER_HEADER_SIZE + sleepTimeMiliS / 1000.0 * this->ChannelSamplingRate;
  for(int i=BUFFER_HEADER_SIZE; i<this->DataCapture; i++)
  {
    this->DataBuffer[i] = i;
  }
}

//----------------------------------------------------------------------------
float* vtkMRMLAlphaOmegaChannelNode::CreateNewDataArray()
{
  float* newDataArray = new float[this->NewDataArraySize];

  for(int i=0; i<this->NewDataArraySize; i++)
  {
    newDataArray[i] = this->DataBuffer[BUFFER_HEADER_SIZE+i] * this->ChannelBitResolution / this->ChannelGain;
  }
  return newDataArray;
}

//----------------------------------------------------------------------------
void vtkMRMLAlphaOmegaChannelNode::AppendNewDataToPreviewArray(float* newDataArray)
{
  if (this->ChannelPreviewSignalArray != nullptr)
  {
    // Add new values to array
    int i=0;
    for(; i<this->NewDataArraySize; i++)
    {
      this->ChannelPreviewSignalArray->SetValue((this->ChannelPreviewArrayIndex+i)%this->ChannelPreviewLengthSamples, newDataArray[i]);
    }
    // Update index for next iteration
    this->ChannelPreviewArrayIndex = (this->ChannelPreviewArrayIndex+i) % this->ChannelPreviewLengthSamples;
    // Add NANs to plot a gap
    for(i=0; i<(0.1*this->ChannelPreviewLengthSamples); i++)
    {
      this->ChannelPreviewSignalArray->SetValue((this->ChannelPreviewArrayIndex+i)%this->ChannelPreviewLengthSamples, NAN);
    }
  }
}

//----------------------------------------------------------------------------
void vtkMRMLAlphaOmegaChannelNode::AppendNewDataToSaveFile(float* newDataArray)
{
  vtkMRMLAlphaOmegaChannelNode::H5Busy.lock();

  const hsize_t ndims = 1;

  // Resize momory dataspace with current ammount of new data
  hsize_t dims[ndims] = {this->NewDataArraySize};  
  H5Sset_extent_simple(this->H5MemoryDataspace, ndims, dims, NULL);

  // Get current dimension of stored data
  H5::DataSet   H5DataSet = this->H5File->openDataSet("data");
  H5::DataSpace H5DataSpace = H5DataSet.getSpace();
  H5DataSpace.getSimpleExtentDims(dims, nullptr);

  hsize_t start[1] = {dims[0]}; // where to start appending
  hsize_t count[1] = {this->NewDataArraySize}; // ammount of new data
  dims[0] = dims[0] + this->NewDataArraySize; // new size

  H5Dset_extent(H5DataSet.getId(), dims); // Extend Data Set
  H5DataSpace = H5DataSet.getSpace(); // Get Space again (changed dimension)
  H5Sselect_hyperslab(H5DataSpace.getId(), H5S_SELECT_SET, start, NULL, count, NULL);
  H5Dwrite(H5DataSet.getId(), H5T_NATIVE_FLOAT, this->H5MemoryDataspace, H5DataSpace.getId(), H5P_DEFAULT, newDataArray);
  H5DataSpace.close();
  H5DataSet.close();
  this->H5File->flush(H5F_SCOPE_LOCAL);

  vtkMRMLAlphaOmegaChannelNode::H5Busy.unlock();
}