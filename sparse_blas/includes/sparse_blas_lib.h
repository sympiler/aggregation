//
// Created by kazem on 10/9/19.
//

#ifndef FUSION_SPARSEBLASLIB_H
#define FUSION_SPARSEBLASLIB_H

#include "def.h"
namespace sym_lib {

 ///////////////////////// SPTRSV

 ///
 /// Forward-substitution on lower triangular matrix L
 /// \param n rank of matrix L
 /// \param Lp row pointers of matrix L
 /// \param Li column indices of matrix L
 /// \param Lx nonzero values of matrix L
 /// \param x initial RHS and outputs as result
 void
 sptrsv_csr(int n, int *Lp, int *Li, double *Lx, double *x);


 void
 sptrsv_csr_levelset(int n, const int *Lp, const int *Li, const double *Lx,
                     double *x,
                     int levels, const int *levelPtr, const int *levelSet);

 void sptrsv_csr_levelset_seq(int n, const int *Lp, const int *Li, const double *Lx,
                              double *x,
                              int levels, const int *levelPtr, const int *levelSet);


 void
 sptrsv_csr_lbc(int n, int *Lp, int *Li, double *Lx, double *x,
                int level_no, int *level_ptr,
                int *par_ptr, int *partition);

 /**
  * @brief perform SpTrsv_CSR based on grouping and LBC
  * @param n number of rows
  * @param Lp Row pointer in the CSR format
  * @param Li Index array in the CSR format
  * @param Lx Val array in the CSR format
  * @param x  right hand side
  * @param level_no number of levels
  * @param level_ptr  pointer to one coarsen level
  * @param par_ptr  pointer to one partition
  * @param partition pointer to array storing node Index
  * @param groupPtr Pointer to the starting location of one group
  * @param groupSet Pointer to the column indices in one group
  */
 void sptrsv_csr_group_lbc(int n, int *Lp, int *Li, double *Lx, double *x,
                              int level_no, int *level_ptr,
                              int *par_ptr, int *partition,  int *groupPtr, int *groupSet);

 void sptrsv_csr_lbc_seq(int n, int *Lp, int *Li, double *Lx, double *x,
                            int level_no, int *level_ptr,
                            int *par_ptr, int *partition);


 /////////////////////////// ICholesky0

 /// Takes the lower part of a SPD matrix and factorize it with incomplete Cholesky0
 /// \param n
 /// \param val
 /// \param rowPtr
 /// \param rowIdx
 void spic0_csr(int n, double *val, int *rowPtr, int *rowIdx);
 void spic0_csr_lbc(int n, double *val, int *rowPtr, int *colIdx,
                    int level_no, int *level_ptr,
                    int *par_ptr, int *partition);


}


#endif //FUSION_SPARSEBLASLIB_H
