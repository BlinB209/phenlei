<div align=left><img src="PHengLEI\Documents/PHengLEI-LOGO.png" style="zoom: 15%;"/><br>

# 风雷软件（PHengLEI2506.v1818）

备注：主分支为master分支，更稳定，更多新功能请切换到ActiveBranch分支。

## 1.软件简介

&ensp;&ensp;&ensp;&ensp;风雷软件<sup>[1,2]</sup>（PHengLEI，Platform for Hybrid ENGineering simulation of flows）是中国空气动力研究与发展中心（CARDC）研发的面向流体工程的混合CFD平台。平台以面向对象的设计理念，采用C++语言编程。2020年12月，风雷软件正式面向全国开源，与其他开源CFD软件相比，风雷软件具有扩展能力强、开发难度低、计算效率高等特点。（更多介绍请阅读PHengLEI/Documents文件夹下的**《风雷软件应用与开发指南》**）

&ensp;&ensp;&ensp;&ensp;风雷软件更多动态和Demo请登录以下网址查看：

- 官网地址: https://www.cardc.cn/nnw/products.aspx?t=9
- 代码库地址：https://forge.osredm.com/PHengLEI/PHengLEI
- 算例库地址：https://forge.osredm.com/PHengLEI/PHengLEI-TestCases
- 论坛（常见问题、算例展示、技术分享）地址：https://osredm.com/forums/theme/38

&ensp;&ensp;&ensp;&ensp;风雷软件用户录制视频教程地址：

- 视频教程：https://www.bilibili.com/video/BV1eX4y1T7yW?from=search&seid=9482198996609923785

[1] 赵钟,等.风雷（PHengLEI）通用CFD软件设计[J]. 计算机工程与科学, 2020, 42(2): 210-219.( ZHAO Zhong, et al. Design of general CFD software PHengLEI [J]. Computer Engineering & Science, 2020, 42(2): 210-219. (in Chinese) )
[2] 赵钟,等.适用于任意网格的大规模并行CFD计算框架PHengLEI[J]. 计算机学报, 2019, 42(11):2368
-2383. ( ZHAO Zhong, et al. PHengLEI: A Large Scale Parallel CFD Framework for Arbitrary Grids [J]. Chinese Journal of Computers, 2019, 42(11): 2368-2383. (in Chinese) )

**声明：** 若用户将该软件用于学术研究或工程应用，须在相关的论文成果的显要位置处标注基于“风雷（PHengLEI）“软件，并引用“风雷（PHengLEI）”软件相关的参考文献（例如[1]和[2]）。

## 2.软件功能

&ensp;&ensp;&ensp;&ensp;风雷软件是一款结构/非结构通用CFD软件，计算范围覆盖低速、亚跨声速和高超声速。软件采用有限体积法求解定常/非定常的雷诺平均NS方程（RANS方程），集成了典型湍流模型，如SA、SST模型等；无粘项采用Roe、Vanleer、AUSM、Steger-Warming等格式；粘性项采用中心格式，时间推进采用LU-SGS或Block LU-SGS隐式方法求解；非定常计算时，采用双时间步方法。针对大规模问题，软件支持分区并行计算，并且使用多重网格技术加速收敛。同时，风雷软件也提供常用前/后置接口，如Gridgen、ICEM-CFD、FieldView、Tecplot、Ensight、ParaView等。


## 3.代码获取

1. 环境准备，安装git，官网地址：https://git-scm.com/;
2. 进入风雷代码库，点击右上角Fork按钮；
3. Fork完成后，会生成并跳转到新的仓库，复制新仓库地址，如https://www.osredm.com/p68217053/PHengLEI;
4. 在本地选择一个目录，右键打开git bash，输入命令进行代码克隆，如git clone https://www.osredm.com/p68217053/PHengLEI;
5. 输入用户名和密码，其中用户名是指上面命令中p开头的用户名，比如p68217053；
6. 项目克隆完成后，进入项目目录，默认分支为master分支，可通过命令切换到ActiveBranch分支，如git checkout ActiveBranch;
7. 切换到开发分支后，可输入git log查看日志，确认当前版本。

## 4.软件安装

&ensp;&ensp;&ensp;&ensp;风雷软件能够在Windows、Linux、Mac系统下运行，源代码采用C++语言编写，需要CMake软件构建项目，并行计算采用MPI库。因此，操作系统必须提供C++编译器、CMake2.8以上版本软件和MPI1.0或MPI2.0标准库。

### 4.1 Windows环境配置

&ensp;&ensp;&ensp;&ensp;Windows环境下所有必备软件按照默认步骤安装即可。

1. 安装Microsoft Visual Studio 2012以上版本；
2. 安装MPI库，推荐采用MSMPI；
3. 安装Cmake。

### 4.2 Linux环境配置

&ensp;&ensp;&ensp;&ensp;Linux环境配置的简要步骤如下：

1. 安装Cmake；
2. 安装MPICH3库；
3. 编译HDF5库；
4. 编译CGNS库；
5. 编译metis库和parmetis库。

**备注：**Linux环境配置的具体步骤请阅读《风雷软件应用与开发指南》。

## 5.如何贡献

&ensp;&ensp;&ensp;&ensp;针对基于风雷软件进行二次开发的用户，后续添加如何测试、提交Pull Request的步骤。

&ensp;&ensp;&ensp;&ensp;软件贡献者信息参见项目根目录下的Contributing文件（待增加）。

## 6.软件版权

&ensp;&ensp;&ensp;&ensp;风雷软件开源协议参见项目根目录下的LICENSE文件。

## 7.鸣谢

&ensp;&ensp;&ensp;&ensp;感谢所有参与风雷软件开发与推广的工作人员，也感谢所有支持风雷软件发展并提出宝贵意见和建议的广大用户。

## 8.联系我们

&ensp;&ensp;&ensp;&ensp;电子邮箱：phenglei@126.com。