#ifndef RCMF_H_
#define RCMF_H_

#include "utilityType.h"

#ifdef __cplusplus
extern "C" {
#endif

/// \brief 输入为COO矩阵格式，输出排序数组，不改变原数据
/// \param nnz 非零元素个数
/// \param nRows 行数
/// \param row 行指针
/// \param col 列指针
/// \param postVetexOrder 行排序数组
/// \param postValOrder 非零元素排序数组
void rcmCOO_nowrite(const label nnz, const label nRows, const label *row,
                    const label *col, label *postVetexOrder,
                    label *postValOrder);

/// \brief 输入为COO矩阵格式，直接改变原数据
/// \param nnz 非零元素个数
/// \param nRows 行数
/// \param row 行指针
/// \param col 列指针
/// \param val 非零元素数组指针
void rcmCOO_rewrite(const label nnz, const label nRows, label *row, label *col,
                    scalar *val);

/// \brief 输入为CSR矩阵格式，输出排序数组，不改变原数据
/// \param nnz 非零元素个数
/// \param nRows 行数
/// \param rowsOffset 压缩行数组指针
/// \param col 非零元素数组列指针
/// \param postVetexOrder 行排序数组
/// \param postValOrder 非零元素排序数组
/// \param newRowsOffset 新的压缩行数组指针
void rcmCSR_nowrite(const label nnz, const label nRows, const label *rowsOffset,
                    const label *col, label *postVetexOrder,
                    label *postValOrder, label *newRowsOffset);

/// \brief 输入为CSR矩阵格式，直接改变原数据
/// \param nnz 非零元素个数
/// \param nRows 行数
/// \param rowsOffset 压缩行数组指针
/// \param col 非零元素数组列指针
/// \param val 非零元素数组指针
void rcmCSR_rewrite(const label nnz, const label nRows, label *rowsOffset,
                    label *col, scalar *val);

/// \brief 输入为LDU矩阵格式，输出排序数组，不改变原数据
/// \param nnz 上三角非零元素个数
/// \param nRows 行数
/// \param row 上三角非零元素行指针
/// \param col 上三角非零元素列指针
/// \param postVetexOrder 行排序数组
/// \param postEdgeOrder 上三角非零元素排序数组
void rcmLDU_nowrite(const label nnz, const label nRows, const label *row,
                    const label *col, label *postVetexOrder,
                    label *postEdgeOrder);

/// \brief 输入为LDU矩阵格式，直接改变原数据
/// \param nnz 上三角非零元素个数
/// \param nRows 行数
/// \param row 上三角非零元素行指针
/// \param col 上三角非零元素列指针
/// \param diagVal 对角线非零元素指针
/// \param upperVal 上三角非零元素指针
/// \param lowerVal 下三角非零元素指针
void rcmLDU_rewrite(const label nnz, const label nRows, label *row, label *col,
                    scalar *diagVal, scalar *upperVal, scalar *lowerVal);

/// \brief 输入为LDU矩阵格式，输出排序数组，不改变原数据
/// 快速模式，以上接口的速度在大规模矩阵时效率极低，建议使用此接口
/// \param nnz 上三角非零元素个数
/// \param nRows 行数
/// \param row 上三角非零元素行指针
/// \param col 上三角非零元素列指针
/// \param postVetexOrder 行排序数组
void rcmLDU_nowrite_fast(const label32 nnz, const label32 nRows,
                         const label32 *row, const label32 *col,
                         label32 *postVetexOrder);

#ifdef __cplusplus
}
#endif

#endif
