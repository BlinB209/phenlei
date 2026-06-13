#include "DefaultWidget.h"
#include <QtGui>
#include <QtWidgets>
#include <QMessageBox>
#include <QString>
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
#include "vtkTextActor.h"
#include "vtkCaptionActor2D.h"
#include "vtkTextProperty.h"
#include "vtkProperty2D.h"

#ifdef WIN32
#if defined(VS) && (_MSC_VER >= 1600)
#pragma execution_character_set("utf-8")
#endif
#endif

#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

//! load init View Widget 
DefaultWidget::DefaultWidget(QWidget* parent)
	:QVTKWidget(parent)
{
	renderer = vtkRenderer::New();

	VTK_CREATE(vtkVectorText, text);
	text->SetText("风雷软件PHengLEI");
	
	VTK_CREATE(vtkPolyDataMapper, mapper);
	mapper->SetInputConnection(text->GetOutputPort());

	// Actor in scene
	actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->SetPosition(15, 15, 0);
	actor->GetProperty()->SetColor(255, 0.0, 0.0);
	actor->SetScale(4, 4, 4);
	renderer->AddActor(actor);

	VTK_CREATE(vtkAxes, axes);
	axes->SetOrigin(-2, -5, 0);
	axes->SetScaleFactor(3);
	axes->SymmetricOn();

	axesActor =	vtkSmartPointer<vtkAxesActor>::New();
	axesActor->SetOrigin(-2, -5, 0);
	axesActor->SetNormalizedShaftLength(1, 1, 1);
	axesActor->SetNormalizedTipLength(0.1, 0.1, 0.1);
	axesActor->SetTotalLength(15, 15, 15);
	axesActor->SetCylinderRadius(0.01);
	renderer->AddActor(axesActor);

	vtkCaptionActor2D* xAxisCaptionActor = axesActor->GetXAxisCaptionActor2D();
	vtkTextProperty* captionTextPropertyX = xAxisCaptionActor->GetCaptionTextProperty();
	captionTextPropertyX->SetColor(1, 0, 0);
	xAxisCaptionActor->SetWidth(0.5);
	xAxisCaptionActor->SetHeight(0.2);
	xAxisCaptionActor->GetTextActor()->SetTextScaleModeToNone();
	xAxisCaptionActor->GetCaptionTextProperty()->SetFontSize(20);

	vtkCaptionActor2D* yAxisCaptionActor = axesActor->GetYAxisCaptionActor2D();
	vtkTextProperty* captionTextPropertyY = yAxisCaptionActor->GetCaptionTextProperty();
	captionTextPropertyY->SetColor(0, 1, 0);
	yAxisCaptionActor->SetWidth(0.5);
	yAxisCaptionActor->SetHeight(0.2);
	yAxisCaptionActor->GetTextActor()->SetTextScaleModeToNone();
	yAxisCaptionActor->GetCaptionTextProperty()->SetFontSize(20);

	vtkCaptionActor2D* zAxisCaptionActor = axesActor->GetZAxisCaptionActor2D();
	vtkTextProperty* captionTextPropertyZ = zAxisCaptionActor->GetCaptionTextProperty();
	captionTextPropertyZ->SetColor(0, 0, 1);
	zAxisCaptionActor->SetWidth(0.5);
	zAxisCaptionActor->SetHeight(0.2);
	zAxisCaptionActor->GetTextActor()->SetTextScaleModeToNone();
	zAxisCaptionActor->GetCaptionTextProperty()->SetFontSize(20);

	tecplotActor = vtkSmartPointer<vtkActor>::New();
	scalarBarActor = vtkSmartPointer<vtkScalarBarActor>::New();

	renderer->SetBackground2(231 / 255.0, 234 / 255.0, 255 / 255.0);
	renderer->SetBackground(147 / 255.0, 149 / 255.0, 255 / 255.0);

	renderer->GradientBackgroundOn();
	GetRenderWindow()->AddRenderer(renderer);
}

DefaultWidget::~DefaultWidget()
{

}

void DefaultWidget::loadTecplotFile(QString fileName)
{

	vtkSmartPointer<vtkTecplotReader> tecplotReader = vtkSmartPointer<vtkTecplotReader>::New();
	tecplotReader->SetFileName(fileName.toLocal8Bit());
	tecplotReader->Update();

	multiBlockDataSet = tecplotReader->GetOutput();
	int num = multiBlockDataSet->GetNumberOfPoints();
	if (num==0)
	{
		QMessageBox::warning(this, QString("ERROR"), QString("File Error"));
	}
	
	vtkSmartPointer<vtkCompositeDataGeometryFilter> polyDataExtractFilter = vtkSmartPointer<vtkCompositeDataGeometryFilter>::New();
	polyDataExtractFilter->SetInputData(multiBlockDataSet);
	polyDataExtractFilter->Update();

	vtkPolyData* polyData = polyDataExtractFilter->GetOutput();
	QString value = "mach";
	polyData->GetPointData()->SetActiveScalars(value.toLocal8Bit());
	vtkDataArray* pointDataArray = polyData->GetPointData()->GetArray(value.toLocal8Bit().constData());
	double pointDataRange[2] = { 0 };
	double MaxValue = -1e30;
	double MinValue = 1e30;
	if (pointDataArray != 0)
	{
		pointDataArray->GetRange(pointDataRange);
		MinValue = min(pointDataRange[0], MinValue);
		MaxValue = max(pointDataRange[1], MaxValue);
	}
	if (MinValue > MaxValue)
		return;
	vtkSmartPointer<vtkLookupTable> lookupTable = vtkSmartPointer<vtkLookupTable>::New();
	lookupTable->SetTableRange(pointDataRange);
	lookupTable->SetHueRange(0.666,0);
	lookupTable->SetRange(MinValue, MaxValue);
	lookupTable->ForceBuild();

	vtkSmartPointer<vtkPolyDataMapper> tecplotMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	tecplotMapper->SetInputData(polyData);
	tecplotMapper->SetLookupTable(lookupTable);
	tecplotMapper->SetScalarRange(MinValue, MaxValue);
	
	tecplotActor->SetMapper(tecplotMapper);
	tecplotActor->GetProperty()->SetRepresentationToSurface();
	
	scalarBarActor->SetLookupTable(tecplotMapper->GetLookupTable());
	scalarBarActor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
	scalarBarActor->GetPositionCoordinate()->SetValue(0.25, 0.01);
	scalarBarActor->SetOrientationToHorizontal();
	scalarBarActor->SetWidth(0.5);
	scalarBarActor->SetHeight(0.1);
	scalarBarActor->SetTitle(value.toLocal8Bit());
	renderer->RemoveActor(actor);
	renderer->RemoveActor(axesActor);
	renderer->RemoveActor(tecplotActor);
	renderer->RemoveActor(scalarBarActor);
	renderer->AddActor(tecplotActor);
	renderer->AddActor(scalarBarActor);
	GetRenderWindow()->Render();
}
