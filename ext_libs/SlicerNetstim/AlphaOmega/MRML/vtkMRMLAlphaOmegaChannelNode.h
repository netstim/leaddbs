
#ifndef __vtkMRMLAlphaOmegaChannelNode_h
#define __vtkMRMLAlphaOmegaChannelNode_h

#include "vtkSlicerAlphaOmegaModuleMRMLExport.h"

// VTK includes
#include <vtkMultiThreader.h>
#include <vtkFloatArray.h>
#include <vtkMRMLTableNode.h>
#include <vtkMRMLPlotSeriesNode.h>

// class vtkFloatArray;

// MRML includes
#include <vtkMRML.h>
#include <vtkMRMLNode.h>
#include <vtkNew.h>

// STD includes
#include <mutex>

// H5
#include "itk_H5Cpp.h"


class VTK_SLICER_ALPHAOMEGA_MODULE_MRML_EXPORT vtkMRMLAlphaOmegaChannelNode : public vtkMRMLNode
{
public:
  static vtkMRMLAlphaOmegaChannelNode *New();
  vtkTypeMacro(vtkMRMLAlphaOmegaChannelNode,vtkMRMLNode);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Create instance
  vtkMRMLNode* CreateNodeInstance() override;

  // Set node attributes from name/value pairs
  void ReadXMLAttributes( const char** atts) override;

  // Write this node's information to a MRML file in XML format.
  void WriteXML(ostream& of, int indent) override;

  // Copy the node's attributes to this object
  void Copy(vtkMRMLNode *node) override;

  // Get unique node XML tag name (like Volume, Model)
  const char* GetNodeTagName() override {return "AlphaOmegaChannel";};

  // AO
  int SetChannelNameAndID(const char* name);
  vtkGetMacro(ChannelName, std::string);
  vtkGetMacro(ChannelID, int);
  
  vtkGetMacro(DataCapture, int);

  static int ClearAllChannelsBuffers();
  static std::vector<std::string> GetAllChannelsNames();

  int WaitAndGetChannelData();
  void SleepAndSimulateData();
  void ClearChannelBuffer();
  void GatherData();

  // Thread
  vtkSetObjectMacro(MultiThreader, vtkMultiThreader);
  vtkGetObjectMacro(MultiThreader, vtkMultiThreader);

  vtkGetMacro(GatheringData, bool);

  void ContinuousGatherData();
  void ThreadedGatherData();
  void TerminateGatherData();

  // Channel Settings
  vtkSetMacro(ChannelSamplingRate, int);
  vtkGetMacro(ChannelSamplingRate, int);

  vtkSetMacro(ChannelGain, float);
  vtkGetMacro(ChannelGain, float);

  vtkSetMacro(ChannelBitResolution, float);
  vtkGetMacro(ChannelBitResolution, float);

  vtkSetMacro(ChannelBufferSizeMiliSeconds, int);
  vtkGetMacro(ChannelBufferSizeMiliSeconds, int);
  vtkGetMacro(ChannelBufferSizeSamples, int);
  int InitializeChannelBuffer();

  // Preview
  vtkSetMacro(ChannelPreviewLengthMiliSeconds, int);
  vtkGetMacro(ChannelPreviewLengthMiliSeconds, int);
  void InitializeChannelPreview();

  void SetChannelPreviewTableNode(vtkMRMLTableNode* tableNode);
  vtkGetObjectMacro(ChannelPreviewTableNode, vtkMRMLTableNode);

  void AppendNewDataToPreviewArray(float* newDataArray);

  // Save File
  static void SetChannelRootSavePath(std::string path){ChannelRootSavePath = path;};
  vtkGetMacro(ChannelRootSavePath, std::string);

  vtkGetMacro(ChannelFullSavePath, std::string);

  void InitializeSaveFile();
  void CloseSaveFile();
  void AppendNewDataToSaveFile(float* newDataArray);

  // Distance to target
  static void SetDriveDistanceToTarget(float dtt){DriveDistanceToTarget = dtt;};

protected:
  vtkMRMLAlphaOmegaChannelNode();
  ~vtkMRMLAlphaOmegaChannelNode() override;
  vtkMRMLAlphaOmegaChannelNode(const vtkMRMLAlphaOmegaChannelNode&);
  void operator=(const vtkMRMLAlphaOmegaChannelNode&);


protected:
  // AO
  int ChannelID{0};
  std::string ChannelName;
  int DataCapture{0};

  // Thread
  vtkMultiThreader *MultiThreader;
  bool GatheringData{false};

  // Channel Settings
  int ChannelSamplingRate;
  float ChannelGain;
  float ChannelBitResolution;
  int ChannelBufferSizeSamples;
  int ChannelBufferSizeMiliSeconds;

  // Signal Preview
  int ChannelPreviewLengthMiliSeconds;
  vtkMRMLTableNode *ChannelPreviewTableNode;
  vtkMRMLPlotSeriesNode *ChannelPreviewPlotSeriesNode;

  // Save File
  static std::string ChannelRootSavePath;
  std::string ChannelFullSavePath{""};
  
  // Distance to target
  static float DriveDistanceToTarget;

  
private:
  short* DataBuffer;
  float* CreateNewDataArray();
  int NewDataArraySize;

  int ChannelPreviewArrayIndex{0};
  vtkFloatArray *ChannelPreviewTimeArray;
  vtkFloatArray *ChannelPreviewSignalArray;

  // int ChannelBufferSizeSamples;
  int ChannelPreviewLengthSamples;

  std::mutex ThreadActiveLock;
  int ThreadActive{false};
  int ThreadID{-1};

  H5::H5File* H5File;
  hid_t H5MemoryDataspace;

  static std::mutex H5Busy;

};

#endif