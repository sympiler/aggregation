//
// Created by behrooz on 2021-07-11.
//

#ifndef LBC_LIB_GLC_H
#define LBC_LIB_GLC_H


#include <iomanip>
#include "sparse_inspector.h"
#include "lbc_utils.h"
#include "utils.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <cmath>
#include <iostream>
#include "sparse_io.h"
#include "test_utils.h"
#include <omp.h>
#include "metis_interface.h"

//#include "sparse_blas_lib.h"
#include "SuperNodalTools.h"
#include "BCSCMatrix.h"
//#include "BLAS.h"

#ifdef SPMP

#include <LevelSchedule.hpp>
#include <synk/barrier.hpp>

#endif


namespace HDAGG{
    enum Kernel { SpTrSv_LL, SpTrSv_RL, SpICh0_LL, SpICh0_RL, SpICh0_UL, SpILU0_UL, General };
    struct twoDIteration{
        int i;
        int k;
        int j;
    };

#ifdef SPMP

    void Convert_LBC_CSR_to_SpMP(sym_lib::CSR* LBC_A, SpMP::CSR* SpMP_A);

#endif
    bool sortbyfirst(const twoDIteration& a, const twoDIteration& b);

    bool sortbysec(const twoDIteration& a, const twoDIteration& b);

    /*
     * @brief Computing Data Dependency Graph of the kernel and store it in CSC format
     * @param n Number of Nodes in the DAG
     * @param DAG_ptr The pointer array in CSC format
     * @param DAG_set The index array in CSC format. It stores the child of each node
     * @return LevelPtr the pointer array in CSC format that store pointers to the nodes inside a level
     * @return LevelSet the nodes are sorted based on level in this array
    */
    sym_lib::timing_measurement computingLevelSet_CSC(int n, const std::vector<int>& DAG_ptr, const std::vector<int>& DAG_set,
                                             std::vector<int>& LevelPtr, std::vector<int>& LevelSet, int& nlevels);

    ///\Description: This function compute the post order, then it will group chains at the end it going to generate
    // a DAG which the DAG-node id's are based on post order, so if we just sort them in the coarse level ptr, we going
    // to execute the DAG based on its post-order traversal
    ///\param n: number of nodes in the dependency DAG
    ///\param DAG_ptr: The pointer array of the DAG in CSC format
    ///\param DAG_set: The index array of the DAG in CSC format
    ///\param groupPtr: number of groups or chains in the DAG
    ///\param post_order: The post order of the DAG traversal using dfs. Note that the groupPtr can group chains based
    // on their post order
    ///\param group_DAG_ptr: The DAG pointer array after changing the id's
    ///\param group_DAG_set: The DAG index array after changing the id's
    void computePostOrderAndGroup(int n, const int* DAG_ptr, const int* DAG_set,
                                  int& ngroups, std::vector<int>& groupPtr, std::vector<int>& post_order,
                                  std::vector<int>& group_DAG_ptr, std::vector<int>& group_DAG_set);

    ///\Description: This function compute the post order, then it will group chains at the end it going to generate
    // a DAG which the DAG-node id's are based on post order, so if we just sort them in the coarse level ptr, we going
    // to execute the DAG based on its post-order traversal
    ///\input n: number of nodes in the dependency DAG
    ///\input DAG_ptr: The pointer array of the DAG in CSC format
    ///\input DAG_set: The index array of the DAG in CSC format
    ///\output post_order: The post order of the DAG traversal using dfs. Note that the groupPtr can group chains based
    void computePostOrder(int n, const int* DAG_ptr, const int* DAG_set,
                          std::vector<int>& post_order);

    int build_levelSet_CSC(size_t n, const int *Lp, const int *Li,
                           int *levelPtr, int *levelSet, int* node_to_level);


    /*
     * @brief assign the nodes' level to each node and store it in Node2Level
     * @param LevelPtr The pointer array in CSC format
     * @param LevelSet The index array in CSC format
     * @param nlevels Number of levels
     * @param num_nodes
     * @return Node2Level array that stores the data
     */
    sym_lib::timing_measurement computingNode2Level(const std::vector<int>& LevelPtr, const std::vector<int>& LevelSet, int nlevels,
                                           int num_nodes, std::vector<int>& Node2Level);


    ///\Description: This function finds
    ///\param n: number of nodes in the dependency DAG
    ///\param lC: The pointer array of the DAG in CSC format
    ///\param lR: The index array of the DAG in CSC format
    ///\param dfsLevel: The starting level
    ///\param ubLevel: number of nodes in the dependency DAG
    ///\param n: number of nodes in the dependency DAG
    ///\param n: number of nodes in the dependency DAG
    ///\param n: number of nodes in the dependency DAG
    ///\param n: number of nodes in the dependency DAG
    void parallel_cc(int n, const int *lC, const int *lR, int lbLevel, int ubLevel,
                     const int *node2Level, int numNodes, const int *nodes,
                     int *node2partition);


    ///\description: Given the WM vectr and levelset information and chunk sizes drove by the HDAGG load balance
    ///unit it will provide the merged level schedule
    ///\input num_threads: The number of cores that required for the computation
    ///\input n: number of nodes in the DAG
    ///\input DAG_ptr: The DAG pointer array stored in the CSC format
    ///\input DAG_set: The DAG index array stored in CSC format
    ///\input level_ptr: The pointer array that shows the nodes in each level
    ///\input level_set: The index array sorted by the nodes in each level
    ///\input node_to_level: The array that maps node id to its respective level
    ///\input WM: The vector that shows how levels are going to be merged
    ///\input Chunk_Sizes: The schedule inside each level
    ///\output coarse_level_ptr: The output pointer array to shows CCs in each level
    ///\output coarse_part_ptr: The output pointer that shows nodes inside each CCs
    ///\output coarse_node_ptr: Nodes sorted based on their CCs
    ///\input apply_serial_schedule: apply the Chunk_Sizes schedule to each level if set to true
    ///\input sort: sort the nodes based on their ids
    void computeSchedule(int num_threads,
                         int n,
                         const int* DAG_ptr,
                         const int* DAG_set,
                         const int* level_ptr,
                         const int* level_set,
                         const int* node_to_level,
                         std::vector<int>& WM,
                         std::vector<std::vector<int>> Chunk_Sizes,
                         int& coarse_level_no,
                         std::vector<int>& coarse_level_ptr,
                         std::vector<int>& coarse_part_ptr,
                         std::vector<int>& coarse_node_ptr,
                         bool apply_serial_schedule,
                         bool postOrder);


    ///\description: This unit find the CCs in each iteration of the HDAGG algorithm
    /// and it will generate the cost of each CC and their positions
    ///\input lb_level: The starting level that merging happens from it (inclusive)
    ///\input computed_up_to: The level that we already want to merge levels with
    ///\input ub_level: The upper level that we want to check whether merging worth it or not
    ///\input n: number of nodes in the DAG
    ///\input DAG_ptr: The DAG pointer array stored in the CSC format
    ///\input DAG_set: The DAG index array stored in CSC format
    ///\input level_ptr: The pointer array that shows the nodes in each level
    ///\input level_set: The index array sorted by the nodes in each level
    ///\input node_to_level: The array that maps node id to its respective level
    ///\input_Output node2partition: It shows that each node belong to which partition up to computed_up_to level
    ///\output part_cnt: Number of CCs so far
    ///\output partition_cost: Cost of computing each CCs
    ///\output node_cost: Cost of each individual nodes
    ///\output total_cost: The total cost of computing these number of nodes in all the cores
    void computeCC(int lb_level, int computed_up_to,
                   int ub_level,
                   int n,
                   const int* DAG_ptr,
                   const int* DAG_set,
                   const int* level_ptr,
                   const int* level_set,
                   const int* node_to_level,
                   std::vector<int>& node2partition,
                   int& part_cnt,
                   std::vector<double>& partition_cost,
                   const std::vector<double>& node_cost,
                   double& total_cost);


    ///\description: This is the main function which greedily search for the best merging window
    ///\input n: Number of nodes inside the current computational DAG
    ///\input nnz: Number of edges inside the DAG
    ///\input DAG_ptr_not_prune: The input computational DAG pointer array which is prune or unprune
    ///\input DAG_set_not_prune: The input computational DAG index array which is prune or unprune
    ///\input node_cost: The cost of each node (floating point operations) or nnz touched
    ///\input cores: Number of cores
    ///\output coarse_level_no: Number of merged levels
    ///\output coarse_level_ptr: The output pointer array to shows CCs in each level
    ///\output coarse_part_ptr: The output pointer that shows nodes inside each CCs
    ///\output coarse_node_ptr: Nodes sorted based on their CCs
    ///\input sparsification: Should we sparsify the DAG or not?
    ///\input sort: Sort the serial's chunk based on the node id?
    ///\return WM: it will return the WM vector for analysis
    std::vector<int> HDAGG(int n, int nnz,
                           std::vector<int> & DAG_ptr_not_prune, std::vector<int> & DAG_set_not_prune,
                           const std::vector<double> & node_cost,
                           int cores,
                           int & coarse_level_no,
                           std::vector<int> & coarse_level_ptr,
                           std::vector<int> & coarse_part_ptr,
                           std::vector<int> & coarse_node_ptr,
                           bool parallelLevelset = false,
                           bool PostOrder = false,
                           bool bin_pack = true);


    std::vector<int> GLC_V2(int n, int nnz,
                         std::vector<int> & DAG_ptr_not_prune, std::vector<int> & DAG_set_not_prune,
                         const std::vector<double> & node_cost,
                         int cores,
                         int & coarse_level_no,
                         std::vector<int> & coarse_level_ptr,
                         std::vector<int> & coarse_part_ptr,
                         std::vector<int> & coarse_node_ptr,
                         bool parallelLevelset = false,
                         bool PostOrder = false,
                         bool bin_pack = true);

    ///\description: This function convert CN groups into Supernode - if they are convertible, it simply change them
    //to nodes
    ///\input n: Number of nodes inside the current computational DAG
    ///\input num_supernode: Number of supernodes
    ///\input group_ptr: the CN group pointer array
    ///\input group_set: The nodes inside the group set array
    ///\input coarse_level_no: Number of merged levels
    ///\input supernode_ptr shows the supernodes
    ///\input_output coarse_level_ptr: The output pointer array to shows CCs in each level
    ///\input_output coarse_part_ptr: The output pointer that shows nodes inside each CCs
    ///\input_output coarse_node_ptr: Nodes sorted based on their CCs
    void convertCNScheduleToSuperNode(int n,
                                      int num_supernode,
                                      int* group_ptr,
                                      int* group_set,
                                      int* supernode_ptr,
                                      int & coarse_level_no,
                                      std::vector<int> & coarse_level_ptr,
                                      std::vector<int> & coarse_part_ptr,
                                      std::vector<int> & coarse_node_ptr);


    double aggressiveSupernodeGrouping(int num_nodes,
                                       int* group_ptr,
                                       int* group_set,
                                       int* DAG_ptr,
                                       int* DAG_set,
                                       int & coarse_level_no,
                                       std::vector<int> & coarse_level_ptr,
                                       std::vector<int> & coarse_part_ptr,
                                       std::vector<int> & coarse_node_ptr);

    ///\Description: This function get the grouped node and generate a Coarsen DAG based on group_ptr and group_set
    ///\input n: Number of nodes inside the DAG
    ///\input ngroups: Number of groups - normally the size of group_ptr is ngroups + 1
    ///\input group_ptr: Pointer that separate groups sorted in group_set array (like CSC/CSR format)
    ///\input group_set: Nodes sorted based on their group in this array
    ///\input DAG_ptr: The original DAG pointer array (in CSC format)
    ///\input DAG_set: The original DAG index array (in CSC format)
    ///\output group_DAG_ptr: The pointer array in CSC format for the grouped DAG
    ///\output group_DAG_set: The index array in CSC format for grouped DAG (shows the child of each node)
    void buildGroupDAG(const int& n,
                       const int& ngroups, const int* group_ptr, const int* group_set,
                       const int* DAG_ptr, const int* DAG_set,
                       std::vector<int>& group_DAG_ptr, std::vector<int>& group_DAG_set);


    ///\Description: Compute elimination tree as an array of parents
    ///\input n: Number of nodes in the matrix A
    ///\input Ap: the pointer array in CSC format
    ///\input Ai: the index array in CSC format
    ///\output parents: The array that shows the parent of each node in elimination tree
    void computeParents(int& n, const int* Ap, const int* Ai, std::vector<int>& parents);

    ///\Description: This Function get the input DAG in CSC format and generates the etree
    ///\input A: The input matrix in CSC format
    ///\output nedges: Number of edges in etree
    ///\output tree_ptr: The pointer array in CSC format
    ///\output tree_set: The index array in CSC format
    void buildETree(sym_lib::CSR* A, int& nedges, std::vector<int>& tree_ptr, std::vector<int>& tree_set);


    ///\Description: Given the array in CSR format, it will compute the none-zero pattern in a row
    ///\input n: The number of rows or columns
    ///\input Ap: The input pointer array in CSR format
    ///\input Ai: The input index array in CSR format
    ///\input k: The row that none-zero pattern should be computed
    ///\input parent: The parent of each node in elimination tree is stored in here
    ///\output s: The stack required for storing then nonzero pattern from s[top] to s[n-1] (should have the size of n)
    ///\internal w: is the array to mark the visited and unvisited nodes (should have the size of n)
    int computeNoneZeroPatternOfRowK(const int &n, const int* Ap, const int* Ai,
                                     const int k, const int* parent, int* s, std::vector<bool>& w);
    ///\Description This Function get the input DAG in CSC format and generate a tree based grouping
    ///\Description Make sure to prune the DAG before passing it to this function
    ///\input n: Number of nodes
    ///\input DAG_ptr: array pointer in CSC format
    ///\input DAG_set: The index array pointer
    ///\output group_ptr: The pointer array that shows the group of nodes
    ///\output group_set: The index array that which is sorted based on groups
    ///\input apply_size_restriction: Apply size restriction
    bool treeBasedGrouping(int n, std::vector<int>& DAG_ptr, std::vector<int>& DAG_set,
                           int& ngroups, std::vector<int>& group_ptr,
                           std::vector<int>& group_set, bool Lfactor = false);


    ///\Description This Function get the input DAG in CSC format and generate a chain based grouping
    ///\Description Make sure to prune the DAG before passing it to this function
    ///\input n: Number of nodes
    ///\input DAG_ptr: array pointer in CSC format
    ///\input DAG_set: The index array pointer
    ///\output group_ptr: The pointer array that shows the group of nodes
    ///\output group_set: The index array that which is sorted based on groups
    ///\input apply_size_restriction: Apply size restriction
    bool chainGrouping(int n, std::vector<int>& DAG_ptr, std::vector<int>& DAG_set,
                           int& ngroups, std::vector<int>& group_ptr,
                           std::vector<int>& group_set);

    ///\Description This Function get the input DAG in CSC format and generate a tree based grouping
    ///\Description Make sure to prune the DAG before passing it to this function
    ///\input n: Number of nodes
    ///\input DAG_ptr: array pointer in CSC format
    ///\input DAG_set: The index array pointer
    ///\output group_ptr: The pointer array that shows the group of nodes
    ///\output group_set: The index array that which is sorted based on groups
    ///\input apply_size_restriction: Apply size restriction
    bool treeBasedGroupingBFS(int n, std::vector<int>& DAG_ptr, std::vector<int>& DAG_set,
                              int& ngroups, std::vector<int>& group_ptr,
                              std::vector<int>& group_set, bool Lfactor = false);


    ///\Description This Function unpacks the nodes in the group in each big w-partitions and then apply
    /// an ordering to that pack. For now, only sort is supported.
    ///\input n: for checking purposes.
    ///\input final_level_no: final number of coarsen levels
    ///\input final_level_ptr: Pointer to the part array that shows parts inside a level
    ///\input_output final_part_ptr: Pointer to the part node pointer where nodes are sorted based on w-partitions
    ///\input_output final_node_ptr: Array that stores nodes based on w-partitions
    ///\input group_ptr: The pointer array that shows the group of nodes
    ///\input group_set: The index array that which is sorted based on groups
    ///\input apply_sort_order: Apply a sort order on each part (nodes which are going to be computed by a single core)
    ///\input apply_post_order: Apply a post order of the DAG on each part
    void ungroupingScheduleAndApplyOrdering(int n, int final_level_no, std::vector<int>& final_level_ptr,
                                            std::vector<int>& final_part_ptr, std::vector<int>& final_node_ptr,
                                            std::vector<int>& group_ptr, std::vector<int>& group_set,
                                            const int* DAG_ptr = NULLPNTR, const int* DAG_set = NULLPNTR,
                                            bool apply_sort_order = true, bool apply_post_order = false);


    ///\Description This function applies a partial sparsification to the input DAG. The sparsification is drived
    /// from sparsifying intel paper 2014. Child in a DAG must be stored in order
    ///\input n: Number of nodes
    ///\input nnz: Number of edges
    ///\input DAG_ptr_not_prune: The input DAG pointer in CSC format
    ///\input DAG_set_not_prune: The input DAG index in CSC format
    ///\output DAG_ptr: The pruned-DAG pointer in CSC format
    ///\output DAG_set: The pruned-DAG index in CSC format
    void partialSparsification(int n, int nnz, const int* DAG_ptr_not_prune, const int* DAG_set_not_prune,
                               std::vector<int>& DAG_ptr, std::vector<int>& DAG_set, bool cliqueSimplification=true);


    ///\Description This function compute cost per node of the DAG based on the kernel. It's simply the number of nonezeros
    ///Touched. For SpTrSv LL version, pass the CSR format of the matrix and for RL version pass the CSC format. For SpICH0
    /// LL version please pass the CSC and CSR format and for RL please pass CSC version. For UL version please pass CSR format of the matrix.
    ///\input nodes: Number of nodes in the DAG (normally equal to the group size or number of columns/rows
    ///\input CSC_Lp: The pointer array in CSC format
    ///\input CSC_Li: The index array in CSC format
    ///\input CSR_Lp: The pointer array in CSR format
    ///\input CSR_Li: The index array in CSR format
    ///\input kernel: The index array in CSC/CSR format
    ///\input group_ptr: The pointer array that shows group in group_set
    ///\input group_set: The nodes sorted based on the group
    ///\input grouped: whether to do a grouping or not
    ///\output cost: The return cost per node
    void costComputation(int nodes, const int* CSC_Lp, const int* CSC_Li, const int* CSR_Lp, const int* CSR_Li,
                            Kernel kernel, const int* group_ptr, const int* group_set, bool grouped,
                            std::vector<double>& cost);

    ///\Description This function creates a 2D tiling of the two outer most loops
    ///\input n: number of rows/columns
    ///\input nnz: Number of none zero elements inside the matrix
    ///\input Lp: The pointer array
    ///\input Li: The index array
    ///\input final_level_no: final number of coarsen levels
    ///\input final_level_ptr: Pointer to the part array that shows parts inside a level
    ///\input final_part_ptr: Pointer to the part node pointer where nodes are sorted based on w-partitions
    ///\input final_node_ptr: Array that stores nodes based on w-partitions
    ///\input T1_size: The index array in CSC/CSR format
    ///\output final_tuple_ptr: The pointer to the pair array
    ///\output final_iter_ptr: The i and j of the outer and inner for loops
    void twoDReordering(int n,
                        int nnz,
                        int* Lp,
                        int* Li,
                        int final_level_no,
                        std::vector<int>& final_level_ptr,
                        std::vector<int>& final_part_ptr,
                        std::vector<int>& final_node_ptr,
                        int T1_size,
                        std::vector<int>& final_tuple_ptr,
                        std::vector<twoDIteration>& final_iter_ptr);


    ///\Description This function create groups and the DAG of this group from final partition
    /// Note that it will clear the elements in output vectors at the beginning.
    ///\input n: Number of nodes that final partitions are representing (size of final_node_ptr)
    ///\input final_level_no: final number of levels
    ///\input final_level_ptr: final level ptr
    ///\input final_part_ptr: final part ptr
    ///\input final_node_ptr: final nodes sorted based on their parts
    ///\input orig_DAG_ptr: Original DAG pointer array in CSC format
    ///\input orig_DAG_set: Original DAG index array in CSC format
    ///\input ngroups: number of parts
    ///\output group_ptr: pointer to group_set array that shows the boundary of each group
    ///\output group_set: Nodes sorted based on their groups
    ///\output DAG_ptr: CSC format pointer array
    ///\output DAG_set: CSC format index array
    void getFinalScheduleDAG(int n, int final_level_no, const int* final_level_ptr,
                                const int* final_part_ptr, const int* final_node_ptr,
                                const int* orig_DAG_ptr, const int* orig_DAG_set,
                                int& ngroups, std::vector<int>& group_ptr, std::vector<int>& group_set,
                                std::vector<int>& DAG_ptr, std::vector<int>& DAG_set);


    /// building coarsened levels based on DAG sparsity
    /// \param n IN dimension
    /// \param Lp IN pointers to columns
    /// \param Li IN: row indices
    /// \param cores IN: number of cores
    /// \param coarse_level_no OUT
    /// \param coarse_level_ptr OUT
    /// \param coarse_part_ptr OUT
    /// \param coarse_node_ptr OUT
    /// \return
    int build_coarsened_level_parallel(const int n, const int* Lp, const int* Li,
                                    const int cores,
                                    int& coarse_level_no,
                                    std::vector<int>& coarse_level_ptr,
                                    std::vector<int>& coarse_part_ptr,
                                    std::vector<int>& coarse_node_ptr);

#ifdef SPMP

    ///\Description This is a parallel Levelset function which works on
    /// a DAG in CSC version and uses SPMP levelset implementation
    ///\input A: The full matrix A in SpMP CSR format
    ///\input nthreads: Number of parallel threads
    ///\output level_ptr: Pointer to the Levelset array
    ///\output level_set: Nodes in the levelset are sort based on their level
    ///\return nlevels: Number of levels => if == -1 means the input graph has cycle
    int levelsetCSRParallel_SpMP(SpMP::CSR* A, int nthreads, std::vector<int>& level_ptr, std::vector<int>& level_set);

#endif



}

#endif //LBC_LIB_GLC_H
