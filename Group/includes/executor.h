//
// Created by labuser (Bangtian Liu) on 9/15/20.
//

#ifndef LBC_LIB_EXECUTOR_H
#define LBC_LIB_EXECUTOR_H
#include <algorithm>
#include <omp.h>
#include <cstring>
#include <functional>

namespace group_cols
{
    int fs_csr_executor_sgroup(int n, int *Lp, int *Li, double *Lx, double *b, double *x, int *groupPtr, int *groupSet,  int ngroup,
                               int levels, int *levelPtr, int *levelSet);
    /**
     * @brief profile serial implementation of triangualr solver based on CSR format
     * @param n number of rows
     * @param Lp pointers to starting address of one row
     * @param Li index array of non-zeros
     * @param flops number of floating operations for serial code
     * @param access_nnz  the total memory access operations on sparse matrix and x array
     * @param reuse_nnz   the reused memory access across two consecutive iterations
     */
    void fs_csr_stat(int n, int *Lp, int *Li,  int &flops, int &access_nnz, int &reuse_nnz);

    /**
     * brief  profile parallel implementation of triangular solver based on CSR
     * @param Lp row pointer in the CSR format
     * @param Li index array in the CSR format
     * @param groupPtr Pointer to the starting location of one group
     * @param groupSet Pointer to the column indices in one group
     * @param levels number of levels
     * @param levelPtr Pointer to the starting location of one level
     * @param levelSet Pointer to index array of one level
     * @param lcost store the maximum difference between nodes for each level
     */
    void fs_csr_levelset_stat(int *Lp, int *Li, int *groupPtr, int *groupSet,
                              int levels, int *levelPtr, int *levelSet, int *lcost);
    /**
     * @brief similar to previous function, but works for coarsening level method
     * @param n
     * @param Lp
     * @param Li
     * @param level_no
     * @param level_ptr
     * @param par_ptr
     * @param partition
     * @param lcost
     */
    void sptrsv_csr_lbc_stat(int n, int *Lp, int *Li,
                             int level_no, int *level_ptr,
                             int *par_ptr, int *partition, int *lcost);

}

#endif //LBC_LIB_EXECUTOR_H
