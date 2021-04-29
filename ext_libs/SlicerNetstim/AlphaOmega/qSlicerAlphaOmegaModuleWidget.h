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

#ifndef __qSlicerAlphaOmegaModuleWidget_h
#define __qSlicerAlphaOmegaModuleWidget_h

#include <QListWidgetItem>
#include "vtkMRMLTransformNode.h"

// SlicerQt includes
#include "qSlicerAbstractModuleWidget.h"

#include "qSlicerAlphaOmegaModuleExport.h"

#include <vtkVariantArray.h>


class qSlicerAlphaOmegaModuleWidgetPrivate;
class vtkMRMLNode;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_ALPHAOMEGA_EXPORT qSlicerAlphaOmegaModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:

  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerAlphaOmegaModuleWidget(QWidget *parent=0);
  virtual ~qSlicerAlphaOmegaModuleWidget();

public slots:

  void onConnectPushButton();
  void onTestPushButton();
  void onDistanceToTargetTransformAdded(vtkMRMLNode * node);
  void onDistanceToTargetTransformModified(vtkMRMLNode * node);
  void onAlphaOmegaChannelNodeChanged(vtkMRMLNode * node);
  void updateGUIFromMRML();
  void updateChannelNodeFromGUI();
  void onChannelActiveStateChanged(int state);

protected:

  QScopedPointer<qSlicerAlphaOmegaModuleWidgetPrivate> d_ptr;
  virtual void setup();
  void populateChannelNamesComboBox();
  void setAndCreateRootSavePath();
  void setConnectingFeedback(bool connecting);
  void setChannelWidgetEnabled(bool enabled);

  
private slots:

  void onConnectionStatusModified(bool* deviceIsConnected);
  void onDistanceToTargetModified(float *distanceToTargetMiliM);

private:
  
  Q_DECLARE_PRIVATE(qSlicerAlphaOmegaModuleWidget);
  Q_DISABLE_COPY(qSlicerAlphaOmegaModuleWidget);

  bool modifyingGUI {false};

};

#endif
