//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      DefauleWidget.h.
//! @brief     DefauleWidget for load the view
//! @author    Hleo.
#ifndef DEFAULTWIDGET_H
#define DEFAULTWIDGET_H

#include <string>
#include <QWidget>
#include "QVTKWidget.h"
#include "vtkActor.h"
#include "vtkAxes.h"
#include "vtkContourFilter.h"
#include "vtkDataSet.h"
#include "vtkFloatArray.h"
#include "vtkGaussianSplatter.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkTubeFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"
#include "vtkVectorText.h"
#include "vtkElevationFilter.h"
#include "vtkTexture.h"
#include "vtkAxesActor.h"
#include "vtkTecplotReader.h"
#include "vtkExtractBlock.h"
#include "vtkCompositeDataGeometryFilter.h"
#include "vtkLookupTable.h"
#include "vtkScalarBarActor.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPointData.h"
using namespace std;
class vtkDataSet;
class vtkRenderer;
class vtkActor;
class vtkScalarBarActor;
class   DefaultWidget : public QVTKWidget
{
	Q_OBJECT
		typedef QVTKWidget Superclass;
public:
	DefaultWidget(QWidget* parent = NULL);
	virtual ~DefaultWidget();

	//! use vtk to load  flow field view
	void loadTecplotFile(QString fileName);

protected:
	vtkRenderer*    renderer;
	vtkMultiBlockDataSet* multiBlockDataSet;
	vtkSmartPointer<vtkActor> tecplotActor;
	vtkSmartPointer<vtkScalarBarActor> scalarBarActor;
	vtkSmartPointer<vtkAxesActor> axesActor;
	vtkSmartPointer<vtkActor> actor;
};
#endif
