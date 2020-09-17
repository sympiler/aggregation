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

    void fs_csr_stat(int n, int *Lp, int *Li,  int &flops, int &access_nnz, int &reuse_nnz);


    void fs_csr_levelset_stat(int *Lp, int *Li, int *groupPtr, int *groupSet,
                              int levels, int *levelPtr, int *levelSet, int *lcost);



}

#endif //LBC_LIB_EXECUTOR_H
