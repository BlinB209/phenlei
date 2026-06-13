/*! \file csrMatrix.hpp
 *  \author Zhao Chengpeng (chengpeng_zhao@foxmail.com)
 *  \date 2022-01-27
 *  \modified 2022-11-08
 *  \brief Block-row distributed compressed sparse row matrix
 */

#ifndef MATRIX_DIRECTSOLVERMATRIX_CSRMATRIX_HPP
#define MATRIX_DIRECTSOLVERMATRIX_CSRMATRIX_HPP

#include "sparseMatrix.hpp"

namespace UNAP {
template <typename Cmpt>
class SPMVBuffers {
 public:
  bool initialized_;
  SPMVBuffers() : initialized_(false) {}
  std::vector<label> sranks_, rranks_, soff_, roffs_, sind_;
  std::vector<Cmpt> sbuf_, rbuf_;
  /*! for each off-diagonal entry spmv_prbuf stores the
   *  corresponding index in the receive buffer
   */
  std::vector<label> prbuf_;
};

/*!
 * \class CSRMatrix
 * \brief CSR格式的分布式稀疏矩阵，按行分块
 */
template <typename Cmpt>
class CSRMatrix : public SparseMatrix<Cmpt> {
 private:
  Communicator* comm_;
  std::vector<label> vtxdist_;
  label beginRow_;
  label localRows_;
  label localnnz_;

  std::vector<label> offDiagStart_;
  mutable SPMVBuffers<Cmpt> spmvBufs_;
  void setupSpMVBuffers() const;

 public:
  /*! 构造空矩阵 */
  CSRMatrix(Communicator* otherComm);
  /*! 调用基类构造函数 */
  CSRMatrix(label n, label nnz, bool symm = true);
  /*!
   * 通过三个数组直接构造CSR矩阵
   * \param rowPtr 反映当前进程矩阵各行元素数量的数组，长localRows_+1
   * \param colInd 当前进程矩阵的列标数组，长localnnz_
   * \param values 当前进程矩阵的数值数组，长localnnz_
   * \param dist 各进程矩阵首行行标数组，长nProcs+1
   * \param symm  矩阵结构是否对称（默认结构对称）
   * \param comm utilities中的Communicator
   * \param deepCopy 深拷贝时为true，浅拷贝为false
   */
  CSRMatrix(label* rowPtr, label* colInd, Cmpt* values, label* dist, bool symm,
            Communicator* comm, bool deepCopy = false);

  /*!
   * 通过COO（分布式）构建CSR矩阵
   * \param beginRow  即起始行行标（全局）
   * \param n  即localRows
   * \param nnz 即localnnz
   * \param cooI  全局行标，需按行顺序依次排列
   * \param cooJ  全局列标，每一行的元素可不按顺序排
   * \param cooV  对应矩阵位置数值
   * \param comm  utilities中的Communicator
   * \param symm  矩阵结构是否对称（默认结构对称）
   */
  CSRMatrix(label beginRow, label n, label nnz, label* cooI, label* cooJ,
            Cmpt* cooV, Communicator* comm, bool symm = true);

  /*!
   * 构建一个串行CSR矩阵
   * \param rowNum 矩阵总行数
   * \param rowPtr 反映矩阵各行元素数量的数组，长Rows_+1
   * \param colInd 矩阵的列标数组，长nnz_
   * \param values 矩阵的数值数组，长nnz_
   * \param symm  矩阵结构是否对称（默认结构对称）
   * \param deepCopy 深拷贝时为true，浅拷贝为false
   */
  CSRMatrix(label rowNum, label* rowPtr, label* colInd, Cmpt* values, bool symm,
            bool deepCopy = false);

  //由串行生成分布式CSR矩阵
  CSRMatrix(const CSRMatrix<Cmpt>& Aseq, Communicator* comm);

  /*! const reference，返回各进程矩阵首行行标vector，长nProcs+1 */
  const std::vector<label>& vtxdist() const { return vtxdist_; }
  /*! reference，返回各进程矩阵首行行标vector，长nProcs+1 */
  std::vector<label>& vtxdist() { return vtxdist_; }
  /*! const reference，返回进程p矩阵首行行标 */
  const label& vtxdist(label p) const { return vtxdist_[p]; }
  /*! reference，返回进程p矩阵首行行标 */
  label& vtxdist(label p) { return vtxdist_[p]; }
  /*! value，当前进程分块矩阵的行数 */
  label lRows() const { return localRows_; }
  /*! value，当前进程分块矩阵的起始行 */
  label beginRow() const { return beginRow_; }
  /*! value，当前进程分块矩阵的末行 */
  label endRow() const { return beginRow_ + localRows_; }
  /*! value，当前进程分块矩阵的非零元个数 */
  label lnnz() const { return localnnz_; }
  /*! reference，当前进程分块矩阵的行数 */
  label& lRows() { return localRows_; }
  /*! reference，当前进程分块矩阵的起始行 */
  label& beginRow() { return beginRow_; }
  /*!reference，当前进程分块矩阵的非零元个数 */
  label& lnnz() { return localnnz_; }
  /*! pointer，utilities中的Communicator */
  Communicator* comm() const { return comm_; }
  /*! 打印矩阵信息，调试用 */
  virtual void print() const /*override*/;

  /*! 析构函数，基类的析构函数负责释放空间 */
  ~CSRMatrix() {}
  /*! 区分对角块与非对角部分 */
  void splitOffDiag();
  /*! 稀疏矩阵向量乘 */
  void spmv(const Cmpt* x, Cmpt* y) const;
  /*! iterativeRefinement中的backward error */
  scalar maxScaledRes(const Cmpt* x, const Cmpt* b) const;
  /*! 1范数 */
  scalar norm1() const;
  /*! gathers the matrix to 1 process */
  CSRMatrix<Cmpt>* gather() const;
  /*! 选主元排序后矩阵结构变为非对称，需要增加额外非零元 */
  void getSymmPattern();
};
/*!
 * 从MatrixMarket格式文件中读取数据，构建CSRMatrix，不同进程保留各自数据
 * \param A output，CSR格式的分布式稀疏矩阵
 * \param filename，文件所在完整路径
 * \param zeroBased，矩阵位置编号是否从0开始，默认为false
 */
template <typename Cmpt>
void readMatrixMarket(CSRMatrix<Cmpt>& A, const char* filename,
                      bool zeroBased = false);

/*!
 * \class CSRGraph
 * \brief 图的CSR格式表示。可以是分布式，按行分块
 * 图包含的边中可能有不属于此图的节点（非对角块）
 */
class CSRGraph {
 private:
  std::vector<label> ptr_, ind_;

 public:
  /*! 默认构造函数 */
  CSRGraph() {}
  /*! 构造函数，根据点和边的数量初始化vector */
  CSRGraph(label nrVert, label nrEdge) : ptr_(nrVert + 1), ind_(nrEdge) {}
  /*! 构造函数，由已有vector构造 */
  CSRGraph(const std::vector<label>& ptr, const std::vector<label>& ind)
      : ptr_(ptr), ind_(ind) {}
  /*! 点的数量 */
  label size() const { return ptr_.size() - 1; }
  /*! 点的数量 */
  label vertices() const { return ptr_.size() - 1; }
  /*! 边的数量 */
  label edges() const { return ind_.size(); }
  /*! const pointer，当前进程图各节点边数量的数组 */
  const label* ptr() const { return ptr_.data(); }
  /*! const pointer，当前进程图的各节点所有边对应另一端的节点编号，按节点顺序 */
  const label* ind() const { return ind_.data(); }
  /*! pointer，当前进程图各节点边数量的数组 */
  label* ptr() { return ptr_.data(); }
  /*! pointer，当前进程图的各节点所有边对应另一端的节点编号，按节点顺序 */
  label* ind() { return ind_.data(); }
  /*! const reference，当前进程节点i对应的边数量 */
  const label& ptr(label i) const { return ptr_[i]; }
  /*! const reference，第i条边对应另一端的节点编号 */
  const label& ind(label i) const { return ind_[i]; }
  /*! reference，当前进程节点i对应的边数量 */
  label& ptr(label i) { return ptr_[i]; }
  /*! reference，第i条边对应另一端的节点编号 */
  label& ind(label i) { return ind_[i]; }
  /*!
   * Apply permutation in local subgraphs（仅对角块）
   * 非对角元素与distributed separators有关，不会被重新编号
   * size of order and iorder: this->size() == high-low
   * order has elements in the range [low, high)
   * iorder had elements in the range [0, high-low)
   */
  void permuteLocal(const std::vector<label>& order,
                    const std::vector<label>& iorder, label low, label high);
  /*! 各节点对应边的另一端节点编号从小到大顺序排列 */
  void sortRows();
};

} /* End namespace UNAP */

#endif /* end of include guard: MATRIX_DIRECTSOLVERMATRIX_CSRMATRIX_HPP */
