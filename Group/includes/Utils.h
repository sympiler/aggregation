//
// Created by labuser (Bangtian Liu) on 9/15/20.
//

#ifndef LBC_LIB_UTILS_H
#define LBC_LIB_UTILS_H

#include <cstdlib>
#include <vector>
#include <deque>
#include <cstring>

namespace group_cols
{
    // Makes an edge inside dependence graph
    inline void connect(int v, int w, std::vector<std::vector<int>> &DAG){
        DAG[v].push_back( w );
    }

    int buildLevelSet_CSC_Queue(int n, int nnz, int *Lp, int *Li, int *&levelPtr,
                                int *&levelSet);


    void fs_csr_inspector_dep(int ngroup, int *groupPtr, int *groupSet, int *gInv, int *Lp, int *Li, std::vector<std::vector<int>> &DAG);


    void rhsInit_csr(int n, int *Ap, int *Ai, double *Ax, double *b);


}
#endif //LBC_LIB_UTILS_H
