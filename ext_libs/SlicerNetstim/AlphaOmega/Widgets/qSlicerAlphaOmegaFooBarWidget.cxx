/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

// FooBar Widgets includes
#include "qSlicerAlphaOmegaFooBarWidget.h"
#include "ui_qSlicerAlphaOmegaFooBarWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_AlphaOmega
class qSlicerAlphaOmegaFooBarWidgetPrivate
  : public Ui_qSlicerAlphaOmegaFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerAlphaOmegaFooBarWidget);
protected:
  qSlicerAlphaOmegaFooBarWidget* const q_ptr;

public:
  qSlicerAlphaOmegaFooBarWidgetPrivate(
    qSlicerAlphaOmegaFooBarWidget& object);
  virtual void setupUi(qSlicerAlphaOmegaFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerAlphaOmegaFooBarWidgetPrivate
::qSlicerAlphaOmegaFooBarWidgetPrivate(
  qSlicerAlphaOmegaFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerAlphaOmegaFooBarWidgetPrivate
::setupUi(qSlicerAlphaOmegaFooBarWidget* widget)
{
  this->Ui_qSlicerAlphaOmegaFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerAlphaOmegaFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerAlphaOmegaFooBarWidget
::qSlicerAlphaOmegaFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerAlphaOmegaFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerAlphaOmegaFooBarWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerAlphaOmegaFooBarWidget
::~qSlicerAlphaOmegaFooBarWidget()
{
}
