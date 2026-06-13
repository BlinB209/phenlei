/*! \file sparseMatrix.hpp
 *  \author Zhao Chengpeng (chengpeng_zhao@foxmail.com)
 *  \date 2022-01-26
 *  \modified 2022-08-01
 *  \brief SparseMatrix class，用来储存稀疏矩阵的基类
 */

#ifndef MATRIX_DIRECTSOLVERMATRIX_SPARSEMATRIX_HPP
#define MATRIX_DIRECTSOLVERMATRIX_SPARSEMATRIX_HPP

#include <complex>

#include "unap.hpp"
#include "unapMPI.hpp"

namespace UNAP {

/*!
 * \class SparseMatrix
 * \brief Base class for compressed sparse matrix storage
 * 可以是CSR、CSC、COO格式，不负责分配空间，但负责指针的销毁
 */
template <typename Cmpt>
class SparseMatrix {
 protected:
  label dim_, nnz_;
  bool symm_;
  label *ptr_, *ind_;
  Cmpt* val_;

  /*! 默认构造 */
  SparseMatrix();
  /*! 此构造函数将指针初始化为NULL */
  SparseMatrix(label n, label nnz, bool symm = true);
  /*! 浅拷贝 */
  SparseMatrix(label n, label nnz, label* ptr, label* ind, Cmpt* val,
               bool symm = true);
  /*! 基类析构函数有必要为虚函数 */
  virtual ~SparseMatrix();

 public:
  /*! 矩阵大小（默认行列相等） */
  label size() const { return dim_; }
  /*! 非零元个数 */
  label nnz() const { return nnz_; }
  /*! 矩阵大小，左值引用 */
  label& size() { return dim_; }
  /*! 非零元个数，左值引用 */
  label& nnz() { return nnz_; }

  /*! const pointer, ptr与ind存放索引，val存放值 */
  const label* ptr() const { return ptr_; }
  const label* ind() const { return ind_; }
  const Cmpt* val() const { return val_; }
  /*! const reference, 返回矩阵某个位置的索引或值 */
  const label& ptr(label i) const { return ptr_[i]; }
  const label& ind(label i) const { return ind_[i]; }
  const Cmpt& val(label i) const { return val_[i]; }
  /*! pointer reference, ptr与ind存放索引，val存放值 */
  label*& ptr() { return ptr_; }
  label*& ind() { return ind_; }
  Cmpt*& val() { return val_; }
  /*! reference, 返回矩阵某个位置的索引或值 */
  label& ptr(label i) { return ptr_[i]; }
  label& ind(label i) { return ind_[i]; }
  Cmpt& val(label i) { return val_[i]; }

  /*! 对称矩阵时为true */
  bool symm() const { return symm_; }
  /*! 设置矩阵的对称性（模式，非数值） */
  void setSymm(bool symm = true) { symm_ = symm; }
  /*!
   * 对矩阵进行重新排序
   *\param iperm 原矩阵的编号i是新矩阵的编号iperm[i]
   *\param perm  新矩阵的编号j是原矩阵的编号perm[j]
   */
  void reorder(const std::vector<label>& iperm, const std::vector<label>& perm);
  /*! 输出大小和非零元个数，虚函数 */
  virtual void print() const;
  /*! 释放矩阵空间 */
  void free();

  /*!
   * 从稀疏矩阵A提取波前阵F11部分的一个块（参考DistributeMatrix格式），虚函数
   * \param F output，存放矩阵值的数组，头指针指向矩阵F的top-left位置
   * \param ldF leading dimension
   * \param row top-left位置在稀疏矩阵A中的对应行标
   * \param nrRows 预提取的该block的行数
   * \param col top-left位置在稀疏矩阵A中的对应列标
   * \param nrCols 预提取的该block的列数
   */
  virtual void extractF11Block(Cmpt* F, label ldF, label row, label nrRows,
                               label col, label nrCols) const {}
  /*!
   * 从稀疏矩阵A提取波前阵F12部分的一个块（参考DistributeMatrix格式），虚函数
   * \param F output，存放矩阵值的数组，头指针指向矩阵F的top-left位置
   * \param ldF leading dimension
   * \param row top-left位置在稀疏矩阵A中的对应行标
   * \param nrRows 预提取的该block的行数
   * \param col top-left位置在稀疏矩阵A中的对应列标
   * \param nrCols 预提取的该block的列数
   * \param upd 存有该波前阵的贡献块在A中对应标号的数组指针
   */
  virtual void extractF12Block(Cmpt* F, label ldF, label row, label nrRows,
                               label col, label nrCols,
                               const label* upd) const {}
  /*!
   * 从稀疏矩阵A提取波前阵F21部分的一个块（参考DistributeMatrix格式），虚函数
   * \param F output，存放矩阵值的数组，头指针指向矩阵F的top-left位置
   * \param ldF leading dimension
   * \param row top-left位置在稀疏矩阵A中的对应行标
   * \param nrRows 预提取的该block的行数
   * \param col top-left位置在稀疏矩阵A中的对应列标
   * \param nrCols 预提取的该block的列数
   * \param upd 存有该波前阵的贡献块在A中对应标号的数组指针
   */
  virtual void extractF21Block(Cmpt* F, label ldF, label row, label nrRows,
                               label col, label nrCols,
                               const label* upd) const {}
};

/*! sorts both indices and values at the same time */
template <typename Cmpt>
void sortIndicesValues(label* ind, Cmpt* val, label begin, label end) {
  if (end > begin) {
    label left = begin + 1, right = end;
    label pivot = (begin + (end - begin) / 2);
    std::swap(ind[begin], ind[pivot]);
    std::swap(val[begin], val[pivot]);
    pivot = ind[begin];
    while (left < right) {
      if (ind[left] <= pivot)
        left++;
      else {
        while (left < --right && ind[right] >= pivot) {
        }
        std::swap(ind[left], ind[right]);
        std::swap(val[left], val[right]);
      }
    }
    left--;
    std::swap(ind[begin], ind[left]);
    std::swap(val[begin], val[left]);
    sortIndicesValues<Cmpt>(ind, val, begin, left);
    sortIndicesValues<Cmpt>(ind, val, right, end);
  }
}

} /* End namespace UNAP */

#endif /* end of include guard: MATRIX_DIRECTSOLVERMATRIX_SPARSEMATRIX_HPP */
