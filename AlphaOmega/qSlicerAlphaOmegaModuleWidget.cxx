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

// Qt includes
#include <QApplication>
#include <QDebug>
#include <QListWidgetItem>
#include <QDir>
#include <QDirIterator>
#include <QDateTime>
#include <QTimer>

// SlicerQt includes
#include "qSlicerAlphaOmegaModuleWidget.h"
#include "ui_qSlicerAlphaOmegaModuleWidget.h"

// MRML
#include <vtkMRMLScene.h>
#include "vtkMRMLTransformNode.h"
#include "vtkMRMLTableNode.h"
#include "vtkMRMLScriptedModuleNode.h"
#include "vtkMRMLAlphaOmegaChannelNode.h"

// vtk
#include <vtkTable.h>
#include <vtkFloatArray.h>
#include <vtkMatrix4x4.h>

// STD
#include <cmath> // NAN
#include <algorithm> // find

// vtkSlicerLogic includes
#include "vtkSlicerAlphaOmegaLogic.h"

// thread
#include "qAlphaOmegaStatusThread.h"




//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerAlphaOmegaModuleWidgetPrivate: public Ui_qSlicerAlphaOmegaModuleWidget
{
  Q_DECLARE_PUBLIC(qSlicerAlphaOmegaModuleWidget);
protected:
  qSlicerAlphaOmegaModuleWidget* const q_ptr;
public:
  qSlicerAlphaOmegaModuleWidgetPrivate(qSlicerAlphaOmegaModuleWidget& object);
  vtkSlicerAlphaOmegaLogic*      logic()const;
  qAlphaOmegaStatusThread* alphaOmegaStatusThread = new qAlphaOmegaStatusThread;
  vtkMRMLScriptedModuleNode* parameterNode = nullptr;
  vtkWeakPointer<vtkMRMLAlphaOmegaChannelNode> CurrentChannelNode;
  QTimer* updateChannelsTablesTimer;
};

//-----------------------------------------------------------------------------
// qSlicerAlphaOmegaModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerAlphaOmegaModuleWidgetPrivate::qSlicerAlphaOmegaModuleWidgetPrivate(qSlicerAlphaOmegaModuleWidget& object)
  : q_ptr(&object)
{
}
//-----------------------------------------------------------------------------
vtkSlicerAlphaOmegaLogic* qSlicerAlphaOmegaModuleWidgetPrivate::logic()const
{
  Q_Q(const qSlicerAlphaOmegaModuleWidget);
  return vtkSlicerAlphaOmegaLogic::SafeDownCast(q->logic());
}

//-----------------------------------------------------------------------------
// qSlicerAlphaOmegaModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerAlphaOmegaModuleWidget::qSlicerAlphaOmegaModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerAlphaOmegaModuleWidgetPrivate(*this) )
{
}

//-----------------------------------------------------------------------------
qSlicerAlphaOmegaModuleWidget::~qSlicerAlphaOmegaModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerAlphaOmegaModuleWidget::setup()
{
  Q_D(qSlicerAlphaOmegaModuleWidget);
  d->setupUi(this);
  this->Superclass::setup();

  // UI
  connect(d->connectPushButton, SIGNAL(clicked()), this, SLOT(onConnectPushButton()));
  connect(d->testPushButton, SIGNAL(clicked()), this, SLOT(onTestPushButton()));
  connect(d->distanceToTargetTransformComboBox, SIGNAL(nodeAddedByUser(vtkMRMLNode *)), this, SLOT(onDistanceToTargetTransformAdded(vtkMRMLNode *)));
  connect(d->distanceToTargetTransformComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode *)), this, SLOT(onDistanceToTargetTransformModified(vtkMRMLNode *)));
  connect(d->alphaOmegaChannelComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode *)), this, SLOT(onAlphaOmegaChannelNodeChanged(vtkMRMLNode *)));

  // Channel Node
  connect(d->channelsNamesComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(updateChannelNodeFromGUI()));
  connect(d->samplingRateSpinBox, SIGNAL(valueChanged(int)), this, SLOT(updateChannelNodeFromGUI()));
  connect(d->gainSpinBox, SIGNAL(valueChanged(double)), this, SLOT(updateChannelNodeFromGUI()));
  connect(d->bitResolutionSpinBox, SIGNAL(valueChanged(double)), this, SLOT(updateChannelNodeFromGUI()));
  connect(d->bufferSizeSpinBox, SIGNAL(valueChanged(int)), this, SLOT(updateChannelNodeFromGUI()));
  connect(d->previewLengthSpinBox, SIGNAL(valueChanged(int)), this, SLOT(updateChannelNodeFromGUI()));
  connect(d->previewTableNodeComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode *)), this, SLOT(updateChannelNodeFromGUI()));
  connect(d->channelActiveCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onChannelActiveStateChanged(int)));

  // Thread
  connect(d->alphaOmegaStatusThread, SIGNAL(connectionStatusModified(bool*)),  this, SLOT(onConnectionStatusModified(bool*)));
  connect(d->alphaOmegaStatusThread, SIGNAL(distanceToTargetModified(float*)), this, SLOT(onDistanceToTargetModified(float*)));

  // Singleton Parameter Node
  d->parameterNode = d->logic()->getParameterNode();

  // Timer
  d->updateChannelsTablesTimer = new QTimer(this);
  QObject::connect(d->updateChannelsTablesTimer, SIGNAL(timeout()), this, SLOT(updateChannelsTables()));

}

//-----------------------------------------------------------------------------
void qSlicerAlphaOmegaModuleWidget::onConnectPushButton()
{
  Q_D(qSlicerAlphaOmegaModuleWidget);

  this->setConnectingFeedback(true);

  d->alphaOmegaStatusThread->terminate();

  if (d->connectPushButton->text() == "Connect")
  {
    d->logic()->AODefaultStartConnection(d->systemMACLineEdit->text().toLocal8Bit().constData());
  }
  else
  {
    d->logic()->AOCloseConnection();
  }

  // check connection
  if (d->logic()->AOIsConnected())
  {
    d->connectPushButton->setText("Disconnect");
    this->populateChannelNamesComboBox();
    this->setAndCreateRootSavePath();
    d->alphaOmegaStatusThread->start();
    d->updateChannelsTablesTimer->start(100);
  }
  else
  {
    d->connectPushButton->setText("Connect");
    d->distanceToTargetLabel->setText("-");
    d->channelsNamesComboBox->clear();
    d->updateChannelsTablesTimer->stop();
  }

  this->setConnectingFeedback(false);
}

//-----------------------------------------------------------------------------
void qSlicerAlphaOmegaModuleWidget::setConnectingFeedback(bool connecting)
{
  Q_D(qSlicerAlphaOmegaModuleWidget);
  d->connectPushButton->setEnabled(!connecting);
  if (connecting)
  {
    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
  }
  else
  {
    QApplication::setOverrideCursor(QCursor(Qt::ArrowCursor));
  }
  QApplication::processEvents();
}

//-----------------------------------------------------------------------------
void qSlicerAlphaOmegaModuleWidget::onTestPushButton()
{
  Q_D(qSlicerAlphaOmegaModuleWidget);
  d->logic()->sayHelloWorld();
  this->setAndCreateRootSavePath();
  this->populateChannelNamesComboBox();
  d->alphaOmegaStatusThread->start();
  d->updateChannelsTablesTimer->start(100);
}

//-----------------------------------------------------------------------------
void qSlicerAlphaOmegaModuleWidget::onConnectionStatusModified(bool* deviceIsConnected)
{
  Q_D(qSlicerAlphaOmegaModuleWidget);

  if (! *deviceIsConnected && d->connectPushButton->text() == "Disconnect")
  {
    d->connectPushButton->animateClick();
  }
}

//-----------------------------------------------------------------------------
void qSlicerAlphaOmegaModuleWidget::onDistanceToTargetModified(float* distanceToTargetMiliM)
{
  Q_D(qSlicerAlphaOmegaModuleWidget);

  d->distanceToTargetLabel->setText(QString::number(*distanceToTargetMiliM));

  if ((d->distanceToTargetTransformComboBox->currentNode() != nullptr) && (!isnan(*distanceToTargetMiliM)))
  {
    vtkMRMLTransformNode* distanceToTargetTransformNode = vtkMRMLTransformNode::SafeDownCast(d->distanceToTargetTransformComboBox->currentNode());
    vtkNew<vtkMatrix4x4> tempMatrix;
    tempMatrix->SetElement(2,3,*distanceToTargetMiliM);
    distanceToTargetTransformNode->SetMatrixTransformToParent(tempMatrix.GetPointer());
  }

  vtkMRMLAlphaOmegaChannelNode::SetDriveDistanceToTarget(*distanceToTargetMiliM);
}


//-----------------------------------------------------------------------------
void qSlicerAlphaOmegaModuleWidget::onDistanceToTargetTransformAdded(vtkMRMLNode * node)
{
  Q_D(qSlicerAlphaOmegaModuleWidget);
  node->SetName(d->logic()->GetMRMLScene()->GenerateUniqueName("Distance To Target").c_str());
}

//-----------------------------------------------------------------------------
void qSlicerAlphaOmegaModuleWidget::onDistanceToTargetTransformModified(vtkMRMLNode * node)
{
  Q_D(qSlicerAlphaOmegaModuleWidget);
  const char * nodeID = "";
  if (node != nullptr)  { nodeID = node->GetID(); }
  d->parameterNode->SetNodeReferenceID("DistanceToTargetTransform", nodeID);
}

//-----------------------------------------------------------------------------
void qSlicerAlphaOmegaModuleWidget::updateChannelsTables()
{
  Q_D(qSlicerAlphaOmegaModuleWidget);
  for (int i=0; i<d->logic()->GetMRMLScene()->GetNumberOfNodesByClass("vtkMRMLAlphaOmegaChannelNode"); i++)
  {
    vtkMRMLAlphaOmegaChannelNode* alphaOmegaChannelNode =  vtkMRMLAlphaOmegaChannelNode::SafeDownCast(d->logic()->GetMRMLScene()->GetNthNodeByClass(i,"vtkMRMLAlphaOmegaChannelNode"));
    vtkMRMLTableNode* channelTableNode = vtkMRMLTableNode::SafeDownCast(alphaOmegaChannelNode->GetChannelPreviewTableNode());
    if (channelTableNode != nullptr && alphaOmegaChannelNode->GetGatheringData())
    {
      channelTableNode->Modified();
    }
  }
}

//-----------------------------------------------------------------------------
void qSlicerAlphaOmegaModuleWidget::onAlphaOmegaChannelNodeChanged(vtkMRMLNode * node)
{
  Q_D(qSlicerAlphaOmegaModuleWidget);
  if (node == nullptr) {return;}

  qvtkReconnect(d->CurrentChannelNode, node, vtkCommand::ModifiedEvent, this, SLOT(updateGUIFromMRML()));
  d->CurrentChannelNode = node;

  this->updateGUIFromMRML();
}

//-----------------------------------------------------------------------------
void qSlicerAlphaOmegaModuleWidget::updateGUIFromMRML()
{
  Q_D(qSlicerAlphaOmegaModuleWidget);

  this->modifyingGUI = true;

  d->channelsNamesComboBox->setCurrentText(QString::fromUtf8(d->CurrentChannelNode->GetChannelName().c_str()));
  d->IDLabel->setText(QString::number(d->CurrentChannelNode->GetChannelID()));
  d->samplingRateSpinBox->setValue(d->CurrentChannelNode->GetChannelSamplingRate());
  d->gainSpinBox->setValue(d->CurrentChannelNode->GetChannelGain());
  d->bitResolutionSpinBox->setValue(d->CurrentChannelNode->GetChannelBitResolution());
  d->bufferSizeSpinBox->setValue(d->CurrentChannelNode->GetChannelBufferSizeMiliSeconds());
  d->previewLengthSpinBox->setValue(d->CurrentChannelNode->GetChannelPreviewLengthMiliSeconds());
  d->previewTableNodeComboBox->setCurrentNode(vtkMRMLTableNode::SafeDownCast(d->CurrentChannelNode->GetChannelPreviewTableNode()));
  d->channelActiveCheckBox->setChecked(d->CurrentChannelNode->GetGatheringData());
  d->channelActiveCheckBox->setEnabled(strcmp(d->channelsNamesComboBox->currentText().toLocal8Bit().constData(),"") != 0);

  this->modifyingGUI = false;
}


//-----------------------------------------------------------------------------
void qSlicerAlphaOmegaModuleWidget::updateChannelNodeFromGUI()
{
  Q_D(qSlicerAlphaOmegaModuleWidget);

  if (d->alphaOmegaChannelComboBox->currentNode() == nullptr || this->modifyingGUI) {return;}
  
  vtkMRMLAlphaOmegaChannelNode* alphaOmegaChannelNode =  vtkMRMLAlphaOmegaChannelNode::SafeDownCast(d->alphaOmegaChannelComboBox->currentNode());
  
  const char* newChannelName = d->channelsNamesComboBox->currentText().toLocal8Bit().constData();
  if (this->channelNameAlreadyInitialized(newChannelName) && strcmp(newChannelName,alphaOmegaChannelNode->GetChannelName().c_str()) != 0)
  {
    newChannelName = "";
    d->channelsNamesComboBox->setCurrentText(QString::fromUtf8(newChannelName));
  }

  alphaOmegaChannelNode->SetChannelNameAndID(newChannelName);
  alphaOmegaChannelNode->SetChannelSamplingRate(d->samplingRateSpinBox->value());
  alphaOmegaChannelNode->SetChannelGain(d->gainSpinBox->value());
  alphaOmegaChannelNode->SetChannelBitResolution(d->bitResolutionSpinBox->value());
  alphaOmegaChannelNode->SetChannelBufferSizeMiliSeconds(d->bufferSizeSpinBox->value());
  alphaOmegaChannelNode->SetChannelPreviewLengthMiliSeconds(d->previewLengthSpinBox->value());
  alphaOmegaChannelNode->SetChannelPreviewTableNode(vtkMRMLTableNode::SafeDownCast(d->previewTableNodeComboBox->currentNode()));

  this->updateGUIFromMRML();
}

//-----------------------------------------------------------------------------
bool qSlicerAlphaOmegaModuleWidget::channelNameAlreadyInitialized(const char* channelName)
{
  std::vector<std::string> initializedChannelsNames = this->getInitializedChannelsNames();
  if (std::find(initializedChannelsNames.begin(), initializedChannelsNames.end(), channelName) != initializedChannelsNames.end())
  {
    return true;
  }
  return false;
}

//-----------------------------------------------------------------------------
std::vector<std::string> qSlicerAlphaOmegaModuleWidget::getInitializedChannelsNames()
{
  Q_D(qSlicerAlphaOmegaModuleWidget);
  std::vector<std::string> initializedChannelsNames;
  for (int i=0; i<d->logic()->GetMRMLScene()->GetNumberOfNodesByClass("vtkMRMLAlphaOmegaChannelNode"); i++)
  {
    vtkMRMLAlphaOmegaChannelNode* alphaOmegaChannelNode =  vtkMRMLAlphaOmegaChannelNode::SafeDownCast(d->logic()->GetMRMLScene()->GetNthNodeByClass(i,"vtkMRMLAlphaOmegaChannelNode"));
    initializedChannelsNames.push_back(alphaOmegaChannelNode->GetChannelName());
  }
  return initializedChannelsNames;
}

//-----------------------------------------------------------------------------
void qSlicerAlphaOmegaModuleWidget::onChannelActiveStateChanged(int state)
{
  this->setChannelWidgetEnabled(!(bool)state);

  Q_D(qSlicerAlphaOmegaModuleWidget);
  if (d->alphaOmegaChannelComboBox->currentNode() == nullptr || this->modifyingGUI) {return;}
  vtkMRMLAlphaOmegaChannelNode* alphaOmegaChannelNode =  vtkMRMLAlphaOmegaChannelNode::SafeDownCast(d->alphaOmegaChannelComboBox->currentNode());

  if ((bool)state)
  {
    alphaOmegaChannelNode->ThreadedGatherData();
  }
  else
  {
    alphaOmegaChannelNode->TerminateGatherData();
  }
}

//-----------------------------------------------------------------------------
void qSlicerAlphaOmegaModuleWidget::setChannelWidgetEnabled(bool enabled)
{
  Q_D(qSlicerAlphaOmegaModuleWidget);

  d->channelsNamesComboBox->setEnabled(enabled);
  d->samplingRateSpinBox->setEnabled(enabled);
  d->gainSpinBox->setEnabled(enabled);
  d->bitResolutionSpinBox->setEnabled(enabled);
  d->bufferSizeSpinBox->setEnabled(enabled);
  d->previewLengthSpinBox->setEnabled(enabled);
  d->previewTableNodeComboBox->setEnabled(enabled);
}

//-----------------------------------------------------------------------------
void qSlicerAlphaOmegaModuleWidget::populateChannelNamesComboBox()
{
  Q_D(qSlicerAlphaOmegaModuleWidget);
  std::vector<std::string> channelsNames = vtkMRMLAlphaOmegaChannelNode::GetAllChannelsNames();
  for (int i = 0; i<channelsNames.size(); i++)
  {
    d->channelsNamesComboBox->addItem(QString::fromStdString(channelsNames[i]));
  }
}

//-----------------------------------------------------------------------------
int getNumberOfSubDirs(const char* path)
{
  int numberOfSubDirs = 0;
  QDirIterator it(path, QDir::Dirs | QDir::NoDotAndDotDot, QDirIterator::NoIteratorFlags);
  while (it.hasNext()) 
  {
    it.next();
    numberOfSubDirs += 1;
  }
  return numberOfSubDirs;
}

void qSlicerAlphaOmegaModuleWidget::setAndCreateRootSavePath()
{
  QString rootSavePath = QDir::cleanPath("C:/LeadORData/" + QDateTime::currentDateTime().toString("yyyyMMdd") + QDir::separator());
  QDir().mkpath(rootSavePath);

  int numberOfSubDirs = getNumberOfSubDirs(rootSavePath.toLocal8Bit().constData());
  rootSavePath = QDir::cleanPath(rootSavePath + QDir::separator() + QString::number(numberOfSubDirs) + QDir::separator());
  QDir().mkpath(rootSavePath);

  vtkMRMLAlphaOmegaChannelNode::SetChannelRootSavePath(rootSavePath.toLocal8Bit().constData());
}
