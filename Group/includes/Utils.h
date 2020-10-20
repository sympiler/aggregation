//
// Created by labuser (Bangtian Liu) on 9/15/20.
//

#ifndef LBC_LIB_UTILS_H
#define LBC_LIB_UTILS_H

#include <cstdlib>
#include <vector>
#include <deque>
#include <cstring>

namespace sym_lib
{
    // Makes an edge inside dependence graph
    inline void connect(int v, int w, std::vector<std::vector<int>> &DAG){
     DAG[v].push_back( w );
    }



    /**
     * @brief dependence inspector for no-grouping version code by taking SpMat in CSR format as Input
     * @param n Number of Rows
     * @param Lp Row pointer in the CSR format
     * @param Li Index array in the CSR format
     * @param DAG In which, each node represents one rows, the edge represents the dependence between two rows
     */
    void fs_csr_inspector_dep(int n, int *Lp, int *Li, std::vector<std::vector<int>> &DAG);
    /**
     * @brief dependence inspector for grouping version code by taking SpMat in CSR format and Grouping Information as Input
     * @param ngroup Number of groups
     * @param groupPtr Pointer to the starting location of one group
     * @param groupSet Pointer to the column indices in one group
     * @param gInv mapping column idx to group idx
     * @param Lp Row pointer in the CSR format
     * @param Li Index array in the CSR format
     * @param DAG In which, each node represents one group of rows, the edge represents the dependence between two groups
     */
    void fs_csr_inspector_dep(int ngroup, int *groupPtr, int *groupSet, int *gInv, int *Lp, int *Li, std::vector<std::vector<int>> &DAG);

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
 * @brief Inspector the dependence between Groups
 * @param ngroup
 * @param groupPtr
 * @param groupSet
 * @param gInv
 * @param Lp
 * @param Li
 * @param DAG
 */
    void fs_ic0csc_inspector_group(int ngroup, int *groupPtr, int *groupSet, int *gInv, int *Lp, int *Li, std::vector<std::vector<int>> &DAG);


    /**
  * @brief used for verifying the correctness of the generated DAG. if circle exists (A->B, B->A), the DAG is wrong
  * @param DAG  the input DAG
  * @return
  */
    void detectDAGCircle(std::vector<std::vector<int>> DAG);


}
#endif //LBC_LIB_UTILS_H
