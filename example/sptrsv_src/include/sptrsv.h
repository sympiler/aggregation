//
// Created by kazem on 10/9/19.
//

#ifndef FUSION_SPARSEBLASLIB_H
#define FUSION_SPARSEBLASLIB_H

#include "aggregation/def.h"
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


;



    ///=============================================================================
    ///============================= HDAGG SPARSE BLASS ============================
    ///=============================================================================
    ///\Description It is a parallel left looking Sparse Triangular Solve which uses wavefronts schedule
    ///\input n Number of iterations or node
    ///\input Lp the pointer array in CSC version
    ///\input Li the index array in CSC version
    ///\input Lx the value array in CSC version
    ///\input levels number of levels in the DAg
    ///\input levelPtr the pointer array in CSC format
    /// that point to starting and ending point of nodes in a level
    ///\input LevelSet the array that store nodes sorted based on their level
    ///\inout x the output
    void sptrsv_csr_levelset(int n, const int *Lp, const int *Li, const double *Lx,
                             int levels, const int *levelPtr, const int *levelSet,
                             double *x);

    /*
     * It is left looking Sparse Triangular Solve
     * @param n Number of iterations or node
     * @param Lp the pointer array in CSC version
     * @param Li the index array in CSC version
     * @param Lx the value array in CSC version
     * @return x the output
     * @param levels number of levels in the DAg
     * @param levelPtr the pointer array in CSC format
     * that point to starting and ending point of nodes in a level
     * @param LevelSet the array that store nodes sorted based on their level
     * @param groupPtr the array pointer for groups
     * @param groupSet the array set for groups. Nodes are sorted based on their group
     */
    void sptrsv_csr_group_levelset(int *Lp, int *Li, double *Lx, double *x,
                                   int level_no, int *level_ptr,
                                   int *levelSet, int *groupPtr, int *groupSet);

}


#endif //FUSION_SPARSEBLASLIB_H
