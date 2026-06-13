#ifndef FORTRANINTERFACE_HPP
#define FORTRANINTERFACE_HPP

#include "unap.hpp"

namespace UNAP {
#ifdef __cplusplus
extern "C" {
#endif

/// \brief 初始化MPI环境 返回communicator
/// \param CommPtr communicator pointer
void comminit_(label64 *CommPtr);

/// \brief get current processor ID and size of communication domain
/// \param CommPtr communicator pointer
/// \param rank current processor ID
/// \param size size of communication domain
void commgetmyidsize_(label64 *CommPtr, label *rank, label *size);

/// \brief create a lduMatrix from giving parameters
/// \param APtrPtr value to store the address of the lduMatrix
/// \param rowSize number of rows in matrix
/// \param upperSize size of upper part in matrix
/// \param lowerAddr row index array for upper part in matrix
/// \param upperAddr column index array for upper part in matrix
/// \param lower data of lower part
/// \param diag data of diagonal part
/// \param upper data of upper part
/// \param commPtr communicator pointer
void ldumatrixcreat_(label64 *APtrPtr, label *rowSize, label *upperSize,
                     label *lowerAddr, label *upperAddr, scalar *lower,
                     scalar *diag, scalar *upper, label64 *commPtr);

/// \brief convert a coo matrix to ldu matrix
/// \param APtrPtr value to store the address of the lduMatrix
/// \param dataPtr coefficients array in coo
/// \param rowsPtr row index array
/// \param columnPtr column index array
/// \param rowSize number of rows
/// \param size size of non-zeros
/// \param symm if matrix is symmetric or not
/// \param commPtr communicator pointer
void coo2ldumatrixcreat_(label64 *APtrPtr, const scalar *dataPtr,
                         const label *rowsPtr, const label *columnPtr,
                         const label *rowSize, const label *size,
                         const label *symm, label64 *commPtr);

/// \brief convert a csr matrix to ldu matrix
/// \param APtrPtr value to store the address of the lduMatrix
/// \param dataPtr coefficients array in coo
/// \param rowsPtr compressed row index array
/// \param columnPtr column index array
/// \param rowSize number of rows
/// \param size size of non-zeros
/// \param symm if matrix is symmetric or not
/// \param commPtr communicator pointer
void csr2ldumatrixcreat_(label64 *APtrPtr, const scalar *dataPtr,
                         const label *compRowsPtr, const label *columnPtr,
                         const label *rowSize, const label *size,
                         const label *symm, label64 *commPtr);

/// \brief matrix interfaces create
/// \param APtrPtr value to store the address of the lduMatrix
/// \param nNeiProcs number of neighbor processors
/// \param destRank array of neighbor processor IDs
/// \param locPosition accumulated number of coefficients per neighbor processor
/// \param faceCells array of local row index for each coefficient
/// \param data coefficients array
void matrixinterfacescreat_(label64 *APtrPtr, const label *nNeiProcs,
                            const label *destRank, const label *locPosition,
                            const label *faceCells, const scalar *data);

/// \brief solve Ax=b by PCG
/// \param x unknowns
/// \param APtrPtr value to store the address of the lduMatrix
/// \param b rhs
/// \param rowSize number of rows
/// \param precond Preconditioner type
/// \param tol absolute tolerance
/// \param relTol relative tolerance
/// \param maxIter maximum number of iterations allowed
/// \param minIter minimum number of iterations
/// \param num_iterations number of iterations generated
/// \param final_res_norm final residual generated
void pcgsolversolve_(scalar *x, label64 *APtrPtr, scalar *b, label *rowSize,
                     label *precond, scalar *tol, scalar *relTol,
                     label *maxIter, label *minIter, label *num_iterations,
                     scalar *final_res_norm);

/// \brief solve Ax=b by PBiCGStab
/// \param x unknowns
/// \param APtrPtr value to store the address of the lduMatrix
/// \param b rhs
/// \param rowSize number of rows
/// \param precond Preconditioner type
/// \param tol absolute tolerance
/// \param relTol relative tolerance
/// \param maxIter maximum number of iterations allowed
/// \param minIter minimum number of iterations
/// \param num_iterations number of iterations generated
/// \param final_res_norm final residual generated
void pbicgstabsolversolve_(scalar *x, label64 *APtrPtr, scalar *b,
                           label *rowSize, label *precond, scalar *tol,
                           scalar *relTol, label *maxIter, label *minIter,
                           label *num_iterations, scalar *final_res_norm);

/// \brief solve Ax=b by Multigrid
/// \param x unknowns
/// \param APtrPtr value to store the address of the lduMatrix
/// \param b rhs
/// \param rowSize number of rows
/// \param aggl Agglomeration type
/// \param Smoother Smoother type
/// \param tol absolute tolerance
/// \param relTol relative tolerance
/// \param maxIter maximum number of iterations allowed
/// \param minIter minimum number of iterations
/// \param num_iterations number of iterations generated
/// \param final_res_norm final residual generated
/// \param faceAreaPtr array of face area
void mgsolversolve_(scalar *x, label64 *APtrPtr, scalar *b, label *rowSize,
                    label *aggl, label *Smoother, scalar *tol, scalar *relTol,
                    label *maxIter, label *minIter, label *num_iterations,
                    scalar *final_res_norm, scalar *faceAreaPtr);

/// \brief reorder COO matrix
/// \param val pointer of coefficients array
/// \param row row index
/// \param col column index
/// \param rowSizePtr size of rows
/// \param sizePtr size of non-zeros
/// \param commPtr communication domain
void reordercoo_(scalar *val, label *row, label *col, label *rowSizePtr,
                 label *sizePtr, label64 *commPtr);

/// \brief reorder upper part as row-wise
/// \param row row index
/// \param col column index
/// \param rowSizePtr size of rows
/// \param sizePtr size of non-zeros
/// \param newOrder position index in the new upper array
/// \param commPtr communication domain
void reorderuface__(label *row, label *col, label *rowSizePtr, label *sizePtr,
                    label *newOrder, label64 *commPtr);

/// \brief reorder lower part as row-wise
/// \param row row index
/// \param col column index
/// \param rowSizePtr size of rows
/// \param sizePtr size of non-zeros
/// \param newOrder position index in the new lower array
/// \param commPtr communication domain
void reorderlface__(label *row, label *col, label *rowSizePtr, label *sizePtr,
                    label *newOrder, label64 *commPtr);

/// \brief reorder values by newOrder
/// \param val pointer of values array
/// \param newOrder new ordering of values
/// \param sizePtr size of non-zeros
/// \param commPtr communication domain
void reordervalue__(scalar *val, label *newOrder, label *sizePtr,
                    label64 *commPtr);

/// more general version as adopted in compass
/// \brief construct the Topology(inner) of matrix
/// \param APtrPtr value to store the address of the lduMatrix
/// \param rowSizePtr number of rows
/// \param rowsPtr row index
/// \param colsPtr column index
/// \param sizePtr number of non-zeros
/// \param commPtr communication domain
void contruct_sw_matrix__(label64 *APtrPtr, const label *rowSizePtr,
                          const label *rowsPtr, const label *colsPtr,
                          const label *sizePtr, label64 *commPtr);

/// \brief construct the Topology(inter-processor) of matrix
/// \param APtrPtr value to store the address of the lduMatrix
/// \param nNeiProcsPtr number of neighbor processors
/// \param destRankPtr array of neighbor processor IDs
/// \param offDiagRowsPtr array of local row index for each coefficient
/// \param offDiagStartsPtr accumulated number of coefficients per neighbor
/// processor
void contruct_sw_matrix_interfaces__(label64 *APtrPtr,
                                     const label *nNeiProcsPtr,
                                     const label *destRankPtr,
                                     const label *offDiagRowsPtr,
                                     const label *offDiagStartsPtr);

#ifdef SW_SLAVE
/// \brief construct a mlb iterator from giving matrix
/// \param APtrPtr value to store the address of the lduMatrix
void construct_mlb_iterator__(label64 *APtrPtr);
#endif

/// \brief fill the data(inner) for giving matrix
/// \param APtrPtr value to store the address of the lduMatrix
/// \param diagPtr diagonal coefficients array
/// \param upperPtr upper coefficients array
/// \param lowerPtr lower coefficients array
void fill_sw_matrix_coefficients__(label64 *APtrPtr, const scalar *diagPtr,
                                   const scalar *upperPtr,
                                   const scalar *lowerPtr);

/// \brief fill the data(inter-processor) for giving matrix
/// \param APtrPtr value to store the address of the lduMatrix
/// \param offDiagCoeffs inter-processor coefficients array
void fill_sw_matrix_interfaces_coefficients__(label64 *APtrPtr,
                                              const scalar *offDiagCoeffs);

/// \brief construct a multigrid solver from giving parameters
/// \param mgPtrPtr value to store the address of the multigrid solver
/// \param APtrPtr value to store the address of the lduMatrix
/// \param AgglPtr value to store the address of the Agglomeration
/// \param faceAreaPtr coefficients array to generate coarse matrix, could be
/// face areas or upper coefficients \param smTypePtr Smoother type \param
/// maxLevelsPtr maximum number of coarse levels \param rowSizeCoarsestPtr
/// number of rows in the coarsest level
void contruct_solver_mg__(label64 *mgPtrPtr, label64 *APtrPtr, label64 *AgglPtr,
                          const scalar *faceAreaPtr, const label *smTypePtr,
                          const label *maxLevelsPtr,
                          const label *rowSizeCoarsestPtr);

#ifdef SW_SLAVE
/// \brief generate mlb iterator in coarse level matrix
/// \param AgglPtr value to store the address of the Agglomeration
void mg_coarse_mlb__(label64 *agglPtr);
#endif

/// \brief set maximum number of iterations in multigrid solver
/// \param solverPtrPtr value to store the address of the multigrid solver
/// \param maxIterPtr maximum number of iterations
void sw_solver_mg_set_maxiter__(label64 *solverPtrPtr, const label *maxIterPtr);

/// \brief set minimum number of iterations in multigrid solver
/// \param solverPtrPtr value to store the address of the multigrid solver
/// \param maxIterPtr minimum number of iterations
void sw_solver_mg_set_miniter__(label64 *solverPtrPtr, const label *minIterPtr);

/// \brief set absolute tolerance in multigrid solver
/// \param solverPtrPtr value to store the address of the multigrid solver
/// \param tolPtr absolute tolerance
void sw_solver_mg_set_tol__(label64 *solverPtrPtr, const scalar *tolPtr);

/// \brief set relative tolerance in multigrid solver
/// \param solverPtrPtr value to store the address of the multigrid solver
/// \param reltolPtr relative tolerance
void sw_solver_mg_set_reltol__(label64 *solverPtrPtr, const scalar *reltolPtr);

/// \brief set number of pre smooth in multigrid solver
/// \param solverPtrPtr value to store the address of the multigrid solver
/// \param numPtr number of pre smooth
void sw_solver_mg_set_npresweeps__(label64 *solverPtrPtr, const label *numPtr);

/// \brief set number of post smooth in multigrid solver
/// \param solverPtrPtr value to store the address of the multigrid solver
/// \param numPtr number of post smooth
void sw_solver_mg_set_npostsweeps__(label64 *solverPtrPtr, const label *numPtr);

/// \brief set number of smooth in finest level in multigrid solver
/// \param solverPtrPtr value to store the address of the multigrid solver
/// \param numPtr smooth in finest level
void sw_solver_mg_set_nfinestsweeps__(label64 *solverPtrPtr,
                                      const label *numPtr);

/// \brief solve Ax=b by Multigrid
/// \param mgPtrPtr value to store the address of the multigrid solver
/// \param APtrPtr value to store the address of the lduMatrix
/// \param xPtr unknowns
/// \param bPtr rhs
/// \param res_normPtr final residual generated
void sw_solve_mg__(label64 *mgPtrPtr, label64 *APtrPtr, scalar *xPtr,
                   scalar *bPtr, scalar *res_normPtr);

/// \brief construct a pbicgstab solver from giving parameters
/// \param solverPtrPtr value to store the address of the solver
/// \param commPtr value to store the address of the communicator
void contruct_solver_pbicgstab__(label64 *solverPtrPtr, label64 *commPtr);

/// \brief set maximum number of iterations in pbicgstab solver
/// \param solverPtrPtr value to store the address of the solver
/// \param maxIterPtr maximum number of iterations
void sw_solver_pbicgstab_set_maxiter__(label64 *solverPtrPtr,
                                       const label *maxIterPtr);

/// \brief set minimum number of iterations in pbicgstab solver
/// \param solverPtrPtr value to store the address of the solver
/// \param minIterPtr minimum number of iterations
void sw_solver_pbicgstab_set_miniter__(label64 *solverPtrPtr,
                                       const label *minIterPtr);

/// \brief set absolute tolerance in pbicgstab solver
/// \param solverPtrPtr value to store the address of the solver
/// \param tolPtr absolute tolerance
void sw_solver_pbicgstab_set_tol__(label64 *solverPtrPtr, const scalar *tolPtr);

/// \brief set relative tolerance in pbicgstab solver
/// \param solverPtrPtr value to store the address of the solver
/// \param reltolPtr relative tolerance
void sw_solver_pbicgstab_set_reltol__(label64 *solverPtrPtr,
                                      const scalar *reltolPtr);

/// \brief set Preconditioner in pbicgstab solver
/// \param solverPtrPtr value to store the address of the solver
/// \param APtrPtr value to store the address of the matrix
/// \param precondTypePtr Preconditioner type
void sw_solver_pbicgstab_set_precond__(label64 *solverPtrPtr, label64 *APtrPtr,
                                       const label *precondTypePtr);

/// \brief solve Ax=b by pbicgstab solver
/// \param solverPtrPtr value to store the address of the solver
/// \param APtrPtr value to store the address of the matrix
/// \param xPtr unknowns
/// \param bPtr rhs
/// \param res_normPtr final residual generated
void sw_solve_pbicgstab__(label64 *solverPtrPtr, label64 *APtrPtr, scalar *xPtr,
                          scalar *bPtr, scalar *res_normPtr);

/// \brief free the matrix
/// \param APtrPtr value to store the address of the matrix
void sw_matrix_destroy__(label64 *APtrPtr);

/// \brief free solver multigrid
/// \param solverPtrPtr value to store the address of the solver
void sw_solver_destroy_mg__(label64 *solverPtrPtr);

/// \brief free solver pbicgstab
/// \param solverPtrPtr value to store the address of the solver
void sw_solver_destroy_pbicgstab__(label64 *solverPtrPtr);

/// \brief create the Topology(inter-processor) of matrix
/// \param APtrPtr value to store the address of the lduMatrix
/// \param offDiagRows array of local row index for each coefficient
/// \param offDiagCols array of local column index for each coefficient
/// \param offDiagPids array of neighbor processor index for each coefficient
/// \param cellNumPtr number of rows in current processor
/// \param faceNumPtr number of inter-processor coefficients in current
/// processor \param postOrders new position array for current position after
/// sorting by neighbor processors
void createinterfaces__(label64 *APtrPtr, label64 *offDiagRows,
                        label64 *offDiagCols, label *offDiagPids,
                        label *cellNumPtr, label *faceNumPtr,
                        label *postOrders);

/// \brief print vector into files
/// \param data values
/// \param size number of values
/// \param name output file name
/// \param commPtr value to store the address of the communicator
void printvector__(scalar *data, label *size, char *name, label64 *commPtr);

#ifdef __cplusplus
}
#endif

}  // namespace UNAP

#endif  // FORTRANINTERFACE_HPP
