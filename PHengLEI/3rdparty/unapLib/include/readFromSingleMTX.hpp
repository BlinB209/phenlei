#ifndef READFROMSINGLEMTX_HPP
#define READFROMSINGLEMTX_HPP

#include "unapLduMatrix.hpp"

namespace UNAP {

namespace TOOLS {
template <typename T>
inline static bool compareTri(const Triplet<T> &a, const Triplet<T> &b) {
  if (a.r_ != b.r_) return (a.r_ < b.r_);
  return (a.c_ < b.c_);
}
template <typename T>
inline static bool compareTriCol(const Triplet<T> &a, const Triplet<T> &b) {
  if (a.c_ != b.c_) return (a.c_ < b.c_);
  return (a.r_ < b.r_);
}
class CompareOrder {
 public:
  std::vector<label> postVetexOrder_;
  CompareOrder(const std::vector<label> &order) : postVetexOrder_(order) {}
  inline bool operator()(const label &i, const label &j);
};
}  // namespace TOOLS

/*! 对Martrix
 * Market格式矩阵进行重新排序（串行），使得上三角按照行优先、下三角按照列优先的顺序排列，以便后续构建拓扑，该函数会生成fileName.new新矩阵文件，mtx文件格式细节详见<https://sparse.tamu.edu/>
 *  \param[in] fileName mtx文件所在完整路径
 */
void reorderingMTXFile(const char *fileName);

/*! 对Martrix
 * Market格式矩阵进行RCM减带宽重排（串行），同时重新排序，使得上三角按照行优先、下三角按照列优先的顺序排列，以便后续构建拓扑，该函数会生成fileName.new新矩阵文件，mtx文件格式细节详见<https://sparse.tamu.edu/>
 *  \param[in] fileName A.mtx文件所在完整路径
 *  \return postVetexOrder
 * 返回RCM排序后旧编号对应的新编号数组，用来后续对向量b进行重排
 */
std::vector<label> rcmReorderingMTXFile(const char *filename);

/*! 读取排序后的mtx矩阵文件（单文件），并构建unap的矩阵拓扑（并行）
 *  \param[out] LduA unap分布式矩阵，包括对角部分和interfaces
 *  \param[in] fileName A.mtx文件所在完整路径（不需要添加.new后缀）
 *  \param[in] zero_based mtx文件编号是否是从0开始
 *  \param[in] symm
 * 矩阵是否对称(默认结构对称），对称的话文件会少存一个三角的信息
 */
void constructMatrixFromMTX(LduMatrix &LduA, const char *fileName,
                            bool zero_based, bool symm);

/*! 读取线性方程组的右端项b
 *  \param[out] b unap中的向量，右端项分布式存储
 *  \param[in] fileName b.mtx文件所在完整路径
 */
void constructVectorFromMTX(scalarVector &b, const char *fileName);

/*! 读取线性方程组的右端项b之前利用RCM排序结果对b进行重排
 *  \param[out] b unap中的向量，右端项分布式存储，与排序后的Ldu矩阵对应
 *  \param[in] fileName b.mtx文件所在完整路径
 *  \param[in] postVetexOrder  rcmReorderingMTXFile函数的返回值
 *  \return distGlobalOrder 该进程的解向量所对应原矩阵的编号数组（全局）
 */
std::vector<label> constructVectorFromMTX(scalarVector &b, const char *fileName,
                                          std::vector<label> postVetexOrder);

}  // namespace UNAP
#endif
