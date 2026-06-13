
# 说明
&ensp;&ensp;&ensp;&ensp;风雷可视化(GUI)是一款基于QT和VTK开发的跨平台数值模拟可视化软件基础框架，本框架可集成求解器、前处理、后处理等模块，实现网格处理、参数设置、流场显示等功能。

## 源码文件夹结构
- 3rdparty:  编译所需要的三方库
- MainWindow： 风雷可视化框架源代码
- Resources： 资源文件

## 编译说明

&ensp;&ensp;&ensp;&ensp;项目采用cmake构建系统组织代码，可以直接使用cmake工具生成visual studio工程（建议使用VS2015），或者可以直接使用支持cmake系统的IDE（比如：vscode等）打开项目文件夹进行构建、编译。

## 运行说明
&ensp;&ensp;&ensp;&ensp;项目建议采用QT5.12和VTK7.1编译，运行时需要将VTK7.1版本的库文件拷贝至bin目录下（bin文件夹由程序编译后自动生成）。

## 特别说明

&ensp;&ensp;&ensp;&ensp;源码绝对路径中不要出现中文字符、空格以及特殊字符:(){}*/?|\等 要求cmake的最低版本为2.9。

## 开源说明
&ensp;&ensp;&ensp;&ensp;关于第三方开源库的声明如下：
&ensp;&ensp;&ensp;&ensp;风雷可视化框架使用了多种用于辅助开发的第三方开源库，比如利用cmake工具构建工程项目，QT用于搭建软件界面，VTK用于二维、三维图形绘制。
为方便用户使用，特在此准备了相关网址以及软件版本号，链接如下：

- CMAKE 2.9

​       网址：https://cmake.org/

- QT 5.12

​       网址：https://www.qt.io/

- VTK 7.1

​       网址：https://vtk.org/

## 致谢

&ensp;&ensp;&ensp;&ensp;感谢@肖翰为风雷团队作出的贡献。

## 联系我们

&ensp;&ensp;&ensp;&ensp;电子邮箱：phenglei@126.com。





