关于第三方开源库的声明如下：

风雷开源软件使用了多种用于辅助的第三方开源库，比如利用Cmake工具构建工程项目，MPI支持基于消息传递的分布式并行编程模型，CGNS支持读入具有cgns格式的网格，利用风雷软件转换生成的网格文件，其数据格式为HDF5，支持并行读写。Metis和Parmetis实现对非结构网格的分区，利用Eigen实现一些简单的矢量或者矩阵运算，FFTW实现快速傅里叶变换，P3dfft实现"束分解"的并行傅里叶变换。Mgrid用于多重网格操作。Hypre和Yhamg是两种数学库。Sacado是一套C++应用的自动微分工具库；PETSc用于求解线性或者非线性方程组。

为方便用户使用，特在此准备了相关的源码压缩包，并提供相关网址链接如下：

- Cmake

​       网址：https://cmake.org/

- MPI

​       网址：https://www.mpich.org/

- CGNS

​       网址：http://cgns.github.io/

- HDF5

​       网址：https://www.hdfgroup.org/

- Metis

​       网址：http://glaros.dtc.umn.edu/gkhome/metis/metis/download

- Parmetis

​       网址：http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download

- Eigen

​       网址：https://eigen.tuxfamily.org/

- Tecio

​       网址：https://tecplot.com/products/tecio-library/

- FFTW

​       网址：http://fftw.org/

- P3dfft

​       网址：https://github.com/sdsc/p3dfft

- Mgrid

​       网址：https://github.com/mrklein/ParMGridGen

- Hypre

​       网址：https://github.com/hypre-space/hypre

- Yhamg

​       网址：https://gitee.com/e-level-parallel-algorithm/yhamg

- PETSc

​       网址：https://gitlab.com/petsc/petsc.git

- Sacado

​       网址：https://github.com/trilinos/Trilinos/tree/master/packages/sacado

下载的第三方库仅用于方便用户编译使用风雷软件，建议用户到官网自行下载并遵守其知识产权协议。 若有异议，可发送邮件至phenglei@126.com，通知开发团队删除。

