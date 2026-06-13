/**
 * @file: utilityInterfaces.h
 * @author: Ren Hu
 * @brief: fortran function interfaces
 * @date:   2019-11-11 10:56:28
 * @last Modified by:   lenovo
 * @last Modified time: 2019-11-26 14:18:23
 */
#ifndef UTILITY_UTILITYINTERFACES_H
#define UTILITY_UTILITYINTERFACES_H

#include "utilityType.h"

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************utility初始化*****************************************/
/**
 * @brief 初始化utility
 */

#define MY_FUNC0(...) MY_FUNC1((int *)0)
#define MY_FUNC1(...) MY_FUNC2(__VA_ARGS__, (char ***)0)
#define MY_FUNC2(...) __initUtility(__VA_ARGS__)
#define VAR_FUNC(_1, _2, _3, NAME, ...) NAME
#define GLUE(x, y) x y
#define initUtility(...) \
  GLUE(VAR_FUNC,( ,##__VA_ARGS__, MY_FUNC2, MY_FUNC1, MY_FUNC0)(__VA_ARGS__))

void __initUtility(int *argc, char ***argv);

/*****************************************进程间数据交换****************************************/

/**
 * @brief 获取进程号
 * @param[out]  pid 进程ID
 */
void getPid(int *pid);

/**
 * @brief 获取进程规模
 * @param[out]  commSize 进程总数
 */
void getCommsize(int *commSize);

/**
 * @brief 规约整形变量序列
 * @param[in][out]  data label类型变量
 * @param[in]  count 变量个数
 * @param[in]  op 操作类型
 */
void allReduceLabels(label *data, const label count, unsigned int op);

/**
 * @brief 规约浮点变量序列
 * @param[inout]  data scalar类型变量
 * @param[in]  count 变量个数
 * @param[in]  op 操作类型
 */
void allReduceScalars(scalar *data, const label count, unsigned int op);

/**
 * @brief 主进程广播整形序列
 * @param[in][out]  data 变量序列
 * @param[in]  count 变量个数
 */
void bcastLabels(label *data, const label count);

/**
 * @brief 主进程广播浮点型序列
 * @param[in][out]  data 变量序列
 * @param[in]  count 变量个数
 */
void bcastScalars(scalar *data, const label count);

/**
 * @brief 主进程收集整形序列
 * @param[in]  sdata 发送序列
 * @param[in][out]  rdata 接收序列
 * @param[in]  count 变量个数
 */
void gatherLabels(label *sdata, label *rdata, const label count);

/**
 * @brief 主进程收集浮点序列
 * @param[in]  sdata 发送序列
 * @param[in][out]  rdata 接收序列
 * @param[in]  count 变量个数
 */
void gatherScalars(scalar *sdata, scalar *rdata, const label count);

/**
 * @brief 主进程收集整数，计算最大值或最小值
 * @param[in] flag 可选参数 MAX，MIN
 * @param[in]  data 本进程值
 * @param[in][out]  result 所有进程最大值或最小值
 * @param[in]  count 变量个数
 */
void extremeLabelsInProcs(const char *flag, label *data, label *result,
                          const label count);

/**
 * @brief 主进程收集浮点数，计算最大值或最小值
 * @param[in] flag 可选参数 MAX，MIN
 * @param[in]  data 本进程值
 * @param[in][out]  result 所有进程最大值或最小值
 * @param[in]  count 变量个数
 */
void extremeScalarsInProcs(const char *flag, scalar *data, scalar *result,
                           const label count);

/*******************************************标准输出*******************************************/
/**
 * @brief 所有进程输出到特定文件
 * @param[in]  format 输出格式
 * @param[in]  ... 变量名，不定参数个数
 */
void par_std_out(const char *format, ...);
/**
 * @brief 所有进程输出到特定文件(fortran格式)
 * @param[in]  format 输出格式
 * @param[in]  ... 变量名，不定参数个数
 */
void par_std_out_(const char *format, ...);

/**
 * @brief 指定进程号输出到屏幕
 * @param[in]  pid 指定进程号
 * @param[in]  format 输出格式
 * @param[in]  ... 变量名，不定参数个数
 */
void proc_std_out(const label pid, const char *format, ...);
/**
 * @brief 指定进程号输出到屏幕(fortran格式)
 * @param[in]  pid 指定进程号
 * @param[in]  format 输出格式
 * @param[in]  ... 变量名，不定参数个数
 */
void proc_std_out_(const label pid, const char *format, ...);

/**
 * @brief 主核输出到屏幕
 * @param[in]  format 输出格式
 * @param[in]  ... 变量名，不定参数个数
 */
void master_std_out(const char *format, ...);
/**
 * @brief 主核输出到屏幕(fortran格式)
 * @param[in]  format 输出格式
 * @param[in]  ... 变量名，不定参数个数
 */
void master_std_out_(const char *format, ...);

// Fortran 字符串输出
// void fort_cout(const char* str, const int* len);
// void fort_pout(const char* str, const int* len);
// void fort_sout(const int* pid, const char* str, const int* len);

#ifdef __cplusplus
}
#endif

#endif  // UTILITY_UTILITYINTERFACES_H
