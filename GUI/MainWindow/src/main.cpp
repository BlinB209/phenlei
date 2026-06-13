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
//! @file      main.cpp.
//! @brief     PHengLEI GUI main entrance.
//! @author    Hleo.
#include "GuiPHengLEI.h"
#include <QTextCodec>
#include <QtWidgets/QApplication>
#ifdef Q_OS_WIN
#include <Windows.h>
#include <MMSystem.h>
#endif
#include "vtkOutputWindow.h"
#include <QMessageBox>
#include <QTranslator>
#include <QFileInfo>
#include <QTextCodec>
#ifdef Q_OS_WIN
#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"" )
#endif

void initApplication()
{
	QTextCodec *codec = QTextCodec::codecForName("UTF-8");

	QTextCodec::setCodecForLocale(codec);

	qApp->setOrganizationName("CARDC");
	qApp->setApplicationName("PHengLEI");
}

int main(int argc, char *argv[])
{
	vtkOutputWindow::SetGlobalWarningDisplay(0);
	QApplication app(argc, argv);
	Q_INIT_RESOURCE(images);
	initApplication();

	GuiPHengLEI window;
	window.showMaximized();
	return app.exec();
}
