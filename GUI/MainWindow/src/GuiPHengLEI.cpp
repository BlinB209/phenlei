#include "GuiPHengLei.h"
#include "DefaultWidget.h"
#include <QMap>
#include <QList>
#include <QAction>
#include <QIcon>
#include <QLabel>
#include <QMenu>
#include <QMenuBar>
#include <QContextMenuEvent>
#include <QMessageBox>
#include <QVBoxLayout>
#include <QStatusBar>
#include <QString>
#include <QSplitter>
#include <QFileDialog>
#if defined(VS) && (_MSC_VER >= 1600)
#pragma execution_character_set("utf-8")
#endif

GuiPHengLEI::GuiPHengLEI()
{
	QWidget *widget = new QWidget;
	setCentralWidget(widget);
	QWidget *topFiller = new QWidget;
	topFiller->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	defaultWidget = new DefaultWidget();

	textView = new TextView();
	textView->append(QString("\n"));
	textView->append(QString("欢迎使用风雷软件"));
	textView->setFont(QFont("宋体", 20));

	textView->setAlignment(Qt::AlignCenter);
	outputTab = new QTabWidget();
	outputTab->setObjectName("tabwgt_output");
	outputTab->addTab(textView, QString("解算器输出显示"));
	outputTab->resize(200, 200);

	QSplitter *splitter = new QSplitter(Qt::Vertical);
	splitter->addWidget(defaultWidget);
	splitter->addWidget(outputTab);

	QVBoxLayout *layout = new QVBoxLayout();
	layout->setMargin(5);
	layout->addWidget(splitter);
	widget->setLayout(layout);

	createActions();
	createMenus();

	QString message = tr("PHengLEI");
	statusBar()->showMessage(message);

	setWindowTitle(tr("PHengLEI"));
	resize(480, 320);
}


#ifndef QT_NO_CONTEXTMENU
void GuiPHengLEI::contextMenuEvent(QContextMenuEvent *event)
{
	QMenu menu(this);
	menu.exec(event->globalPos());
}
#endif 

void GuiPHengLEI::newFile()
{
	infoLabel->setText(tr("Invoked <b>File|New</b>"));
}

void GuiPHengLEI::createActions()
{
	importTecplotAction = new QAction(QString("导入"), this);
	connect(importTecplotAction, &QAction::triggered, this, &GuiPHengLEI::importTecplotFile);
	
	newAction = new QAction(QString("新建工程"), this);
	newAction->setStatusTip(QString("新建工程"));
	//connect(newAct, &QAction::triggered, this, &GUIPHengLEI::newFile);

	openAction = new QAction(QString("打开工程"), this);
	openAction->setStatusTip(QString("打开工程"));

	saveAction = new QAction(QString("保存工程"), this);
	saveAction->setStatusTip(QString("保存工程"));

	saveAsAction = new QAction(QString("另存为"), this);
	saveAsAction->setStatusTip(QString("另存为"));

	exitAction = new QAction(QString("退出"), this);
	exitAction->setStatusTip(QString("退出"));

	startAction = new QAction(QString("启动"), this);
	startAction->setStatusTip(QString("启动"));

	stopAction = new QAction(QString("停止"), this);
	stopAction->setStatusTip(QString("停止"));

	cutAction = new QAction(QString("清空"), this);
	cutAction->setStatusTip(QString("清空"));

	gridAction = new QAction(QString("网格视图"), this);
	gridAction->setStatusTip(QString("网格视图"));

	airAction = new QAction(QString("气动力视图"), this);
	airAction->setStatusTip(QString("气动力视图"));

	flowAction = new QAction(QString("流场视图"), this);
	flowAction->setStatusTip(QString("流场视图"));
	
	settingAction = new QAction(QString("设置"), this);
	settingAction->setStatusTip(QString("设置"));
	
	workAction = new QAction(QString("作业查询"), this);
	workAction->setStatusTip(QString("作业查询"));

	saveViewAction = new QAction(QString("保存视图"), this);
	saveViewAction->setStatusTip(QString("保存视图"));

	dataAction = new QAction(QString("数据处理"), this);
	dataAction->setStatusTip(QString("数据处理"));
	
	aboutAction = new QAction(QString("关于"), this);
	aboutAction->setStatusTip(QString("关于"));

	userAction = new QAction(QString("用户手册"), this);
	userAction->setStatusTip(QString("用户手册"));

	demoAction = new QAction(QString("算例库"), this);
	demoAction->setStatusTip(QString("算例库"));

}

void GuiPHengLEI::createMenus()
{
	importMenu = menuBar()->addMenu(QString("文件"));
	importMenu->addAction(importTecplotAction);

	fileMenu = menuBar()->addMenu(QString("工程"));
	fileMenu->addAction(newAction);
	fileMenu->addAction(openAction);
	fileMenu->addAction(saveAction);
	fileMenu->addAction(saveAsAction);
	fileMenu->addSeparator();
	fileMenu->addAction(exitAction);

	editMenu = menuBar()->addMenu(QString("计算"));
	editMenu->addAction(startAction);
	editMenu->addAction(stopAction);
	editMenu->addSeparator();
	editMenu->addAction(cutAction);
	editMenu->addSeparator();

	viewMenu = menuBar()->addMenu(QString("视图"));
	viewMenu->addAction(gridAction);
	viewMenu->addAction(airAction);
	viewMenu->addAction(flowAction);

	formatMenu = menuBar()->addMenu(QString("工具"));
	formatMenu->addAction(settingAction);
	formatMenu->addAction(workAction);
	formatMenu->addAction(saveViewAction);
	formatMenu->addAction(dataAction);

	helpMenu = menuBar()->addMenu(QString("帮助"));
	helpMenu->addAction(aboutAction);
	helpMenu->addAction(userAction);
	helpMenu->addAction(demoAction);
}

void GuiPHengLEI::importTecplotFile()
{
	QString fileName = QFileDialog::getOpenFileName(this, QString("打开流场结果文件"),
		NULL, QString("流场文件(*.plt)"));
	if (GuiPHengLEI::containsChinese(fileName))
	{
		QMessageBox::warning(this, QString("错误"), QString("工程目录不能含有中文，请重新选择工程目录"));
		return;
	}
	defaultWidget->loadTecplotFile(fileName);
}

bool GuiPHengLEI::containsChinese(const QString &target)
{
	return target.contains(QRegExp("[\\x4e00-\\x9fa5]"));
}
