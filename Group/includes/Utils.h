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

    /**
     * @brief dependence inspector
     * @param n
     * @param Lp
     * @param Li
     * @param DAG
     */
    void fs_csr_inspector_dep(int n, int *Lp, int *Li, std::vector<std::vector<int>> &DAG);

    void fs_csr_inspector_dep(int ngroup, int *groupPtr, int *groupSet, int *gInv, int *Lp, int *Li, std::vector<std::vector<int>> &DAG);


    void rhsInit_csr(int n, int *Ap, int *Ai, double *Ax, double *b);

    /**
     * @brief Convert the input DAG to a group-version DAG by applying the group information(ngroup, groupPtr, groupSet, groupInv)
     * @param DAG  the input DAG, which indicates the dependence between each row or column from SpMat.
     * @param groupPtr the pointer to the starting addree of one group
     * @param groupSet the pointer to the index array
     * @param groupInv mapping the column Id to group Id
     * @param ngroup number of groups
     * @return the grouped version DAG, in which each node presents a couple of columns in one group rather than one column
     */
    std::vector<std::vector<int>> Group_DAG(std::vector<std::vector<int>> DAG, int *groupPtr, int *groupSet, int *groupInv, int ngroup);

    /**
     * @brief used for verifying the correctness of the generated DAG. if circle exists (A->B, B->A), the DAG is wrong
     * @param DAG
     * @return
     */
    bool detectDAGCircle(std::vector<std::vector<int>> DAG);


}
#endif //LBC_LIB_UTILS_H
