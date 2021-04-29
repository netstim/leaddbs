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
#include <vtkSlicerAlphaOmegaLogic.h>

// AlphaOmega includes
#include "qSlicerAlphaOmegaModule.h"
#include "qSlicerAlphaOmegaModuleWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerAlphaOmegaModulePrivate
{
public:
  qSlicerAlphaOmegaModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerAlphaOmegaModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerAlphaOmegaModulePrivate::qSlicerAlphaOmegaModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerAlphaOmegaModule methods

//-----------------------------------------------------------------------------
qSlicerAlphaOmegaModule::qSlicerAlphaOmegaModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerAlphaOmegaModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerAlphaOmegaModule::~qSlicerAlphaOmegaModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerAlphaOmegaModule::helpText() const
{
  return "This is a loadable module that can be bundled in an extension";
}

//-----------------------------------------------------------------------------
QString qSlicerAlphaOmegaModule::acknowledgementText() const
{
  return "This work was partially funded by NIH grant NXNNXXNNNNNN-NNXN";
}

//-----------------------------------------------------------------------------
QStringList qSlicerAlphaOmegaModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("John Doe (AnyWare Corp.)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerAlphaOmegaModule::icon() const
{
  return QIcon(":/Icons/AlphaOmega.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerAlphaOmegaModule::categories() const
{
  return QStringList() << "Netstim";
}

//-----------------------------------------------------------------------------
QStringList qSlicerAlphaOmegaModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerAlphaOmegaModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerAlphaOmegaModule
::createWidgetRepresentation()
{
  return new qSlicerAlphaOmegaModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerAlphaOmegaModule::createLogic()
{
  return vtkSlicerAlphaOmegaLogic::New();
}

//-----------------------------------------------------------------------------
QStringList qSlicerAlphaOmegaModule::associatedNodeTypes() const
{
  return QStringList() << "vtkMRMLAlphaOmegaChannelNode";
}