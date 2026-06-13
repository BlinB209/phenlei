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
//! @file      GUIPHengLEI.h.
//! @brief     GUIPHengLEI widget for load software interface
//! @author    Hleo.
#ifndef GUIPHENGLEI_H
#define GUIPHENGLEI_H

#include <QMainWindow>
#include <QPushButton>
#include "DefaultWidget.h"
#include "TextView.h"
QT_BEGIN_NAMESPACE
class QAction;
class QActionGroup;
class QLabel;
class QMenu;
QT_END_NAMESPACE

class GuiPHengLEI : public QMainWindow
{
	Q_OBJECT
private:
	GuiPHengLEI(GuiPHengLEI* q);
public:
	GuiPHengLEI();

protected:
#ifndef QT_NO_CONTEXTMENU
	void contextMenuEvent(QContextMenuEvent *event) override;
#endif //! QT_NO_CONTEXTMENU
	
	//! slots method
private slots:
	void newFile();
	//! import tecplotFile for load  flow field view
	void importTecplotFile();
	
private:
	//! create Action
	void createActions();
	//! create Menus
	void createMenus();

	bool containsChinese(const QString &target);
	//! MenuButton
	QMenu *fileMenu;
	QMenu *editMenu;
	QMenu *formatMenu;
	QMenu *viewMenu;
	QMenu *helpMenu;
	QMenu *importMenu;

	//! Action Button
	QAction *importTecplotAction;
	QActionGroup *alignmentGroup;
	QAction *newAction;
	QAction *openAction;
	QAction *saveAction;
	QAction *saveAsAction;
	QAction *exitAction;
	QAction *startAction;
	QAction *stopAction;
	QAction *cutAction;
	QAction *gridAction;
	QAction *airAction;
	QAction *flowAction;
	QAction *settingAction;
	QAction *leftAlignAction;
	QAction *rightAlignAction;
	QAction *justifyAction;
	QAction *centerAction;
	QAction *setLineSpacingAction;
	QAction *workAction;
	QAction *saveViewAction;
	QAction *dataAction;
	QAction *aboutAction;
	QAction *userAction;
	QAction *demoAction;
	QLabel *infoLabel;

	TextView*    textView;
	QTabWidget*  outputTab;
	DefaultWidget*     defaultWidget;
};
#endif
