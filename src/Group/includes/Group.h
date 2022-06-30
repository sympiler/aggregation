//
// Created by labuser (Bangtian Liu) on 9/15/20.
//

#ifndef LBC_LIB_GROUP_H
#define LBC_LIB_GROUP_H
#include <cstring>
#include <vector>
#include <cstdlib>
#include <unordered_set>
#include <iostream>
#include <set>
#include <algorithm>
#include <math.h>
//#define NDEBUG
#include <assert.h>
#include <Utils.h>
#include <def.h>

namespace sym_lib
{

    class group {
    private:
        int *status; // status for each column to help to do grouping
        int *next; // pointer to its next column of one column/row
        int *pre; // pointer to its previous column of one column/row
        int *mP; // Pointer array in CSR/CSC
        int *mI; // Index array in CSR/CSC
        int ncol; // number of columns
        bool *visited; // indicate whether one column/row is visited or not, used for grouping

        int *child; // indicates its first dependent column/row based on first off-diagonal element

    public:
        /**
         * @brief construction function
         * @param n number of rows/columns
         * @param p  row/column pointer in the CSR/CSC format
         * @param i index array in the CSR/CSC format
         */
         group(int n, int *p, int *i) : ncol(n), mP(p), mI(i) {
         status = (int *)malloc(sizeof(int)*n);
         memset(status, 0, sizeof(int)*n);

         next = (int *)malloc(sizeof(int)*n);
         memset(next, 0, sizeof(int)*n);

         pre = (int *)malloc(sizeof(int)*n);
         memset(pre, 0, sizeof(int)*n);

         visited = (bool *)malloc(sizeof(bool)*n);
         memset(visited, 0, sizeof(bool)*n);

         child =(int *)malloc(sizeof(int)*ncol);
         memset(child, 0, sizeof(int)*ncol);


        }


        /**
         * @brief Grouping consecutive columns with dependence of sptrsv_csc
         * @param groupPtr
         * @param groupSet
         * @param ngroup
         * @param groupInv
         */
        void inspection_spicocsc_v1(int *groupPtr, int *groupSet, int &ngroup, int *groupInv);

        void inspection_spicocsc_v2(int *groupPtr, int *groupSet, int &ngroup, int *groupInv);

        /**
      * @brief Grouping consecutive columns with dependence for sptrsv_csr
      * @param groupPtr Pointer to the starting location of one group
      * @param groupSet Pointer to the column indices in one group
      * @param ngroup Number of groups
      * @param groupInv  mapping the column Id to group Id
      */
        void inspection_sptrsvcsr_v1(int *groupPtr, int *groupSet, int &ngroup, int *groupInv);

        /**
         * @brief Grouping consecutive columns with dependence and also includes empty columns/rows for sptrsv_csr
         * @param groupPtr Pointer to the starting location of one group
         * @param groupSet Pointer to the column indices in one group
         * @param ngroup Number of groups
         * @param groupInv mapping the column Id to group Id
         */
        void inspection_sptrsvcsr_v2(int *groupPtr, int *groupSet, int &ngroup, int *groupInv);


        /**
         * @brief group every blksize rows/columns
         * @param n  number of rows/cols
         * @param groupPtr the pointer to the starting address of one group
         * @param groupSet the pointer to index array
         * @param ngroup number of columns
         * @param ginv mapping column idx to group idx
         * @param blksize  parameter for grouping
           **/
        void NaiveGrouping(int n, int *groupPtr, int *groupSet,  int &ngroup, int *ginv, int blksize=1);
    };

    //============================================ GreedyGrouping ================================================
    class GreedyGrouping{
    private:
        //It will maintain a work-set for grouping maximum of group_max_size elements;
        int max_search_space_size;
        // Minimum cost that a group should have, any group with less than this cost will be merged into nearest column with lowest cost
        double minimum_group_cost_threshold;
        bool merge_less_th;
        // a 2D array where the first dimension is the node number and the second dimension is the set of columns
        // it presents.
        std::vector<int> group_set;
        std::vector<std::vector<int>> group_children;
        //store the group index of each node
        std::vector<int> node_group_idx;
        //Both of these variables should be equal at the end of the grouping
        int total_num_child;
        int group_nnz;
        int group_size;
        // pointers to the CSR and CSC format of the input matrix
        int ptr_size;
        int* csc_ptr;
        int* csc_idx;
        // Candidates for grouping
        std::unordered_set<int> candidate_pool;
        std::vector<int> current_candidates;
    public:
        /*
         * @brief: constructor of the greedyGrouping class
         * @param: cols_rows: the number of column/row of the matrix
         * @param: csc_ptr: sparse matrix csc format column ptr
         * @param: csc_idx: sparse matrix csc format row indices
         * @param: max_search_space_size: Whenever the algorithm see the possible choices,
         * it only stores the choices that are inside the range of [group_head, group_head + max_search_space_size)
         * minimum_group_cost_threshold
         * @param: minimum_group_cost_threshold:  The algorithm will merge groups with costs less than this threshold
         */
        GreedyGrouping(int cols_rows, int* csc_ptr, int* csc_idx,
                       int max_search_space_size, double minimum_group_cost_threshold, bool merge_less_th);
        ~GreedyGrouping();
        /*
         * @brief: Start grouping (find consecutive columns in a connected component
         */
        void startGrouping();
        /*
         * @brief: give the closest candidate for grouping to node with index group_head
         */
        int getBestCandidate(int group_head);

        /*
         * @brief: give the closest candidate for grouping to node with index group_head
         * It doesn't have any sorting function in it because I couldn't make it work
         * for big matrices
         */
        int naive_getBestCandidate(int group_head);

        /*
         * @brief: Add candidates (the children of node) to the candidates pool
         */
        void addChildesToCandidate(int node);
        /*
         * @brief Shrinking the group_set, group_children, group_parent space
         * Also, it will do a couple of simple check to somehow make sure the grouping was correct
         */
        void finalizeGroupingProcedure();

        /*
         * @brief Get grouped CSC format: if we assume all the grouped nodes as a single node
         * we will have a new sparse matrix. This function will return the CSC format of
         * That sparse matrix DAG. Note that you should allocate the
         * memory needed for these pointers before calling it.
         * The size of group_csc_ptr should be number_of_group + 1
         * and the size of group_csc_idx should be equal to grouped matrix nnz
         * @param group_csc_ptr: same as a normal column ptr in CSC format
         * @param group_csc_idx: same as a normal row idx in CSC format
         * @param group_csc_ptr_size: This is for checking purposes
         * @param group_csc_idx_size: Same as group_size, it is for checking purposes
         */
        void getGroupedCSCFormat(int* group_csc_ptr, int* group_csc_idx,
                                 int group_csc_ptr_size, int group_csc_idx_size);
        /*
         * @brief Get grouped conversion to original matrix. since we only group consecutive nodes
         * the size of group_range is 2 * group_size, and each two consecutive elements in this group_range array
         * shows the range of node indices inside a group
         * @param group_range: The array in which the group range are stored
         * for example the range of group[2] = [ group_range[2 * 2], group_range[2 * 2 + 1])
         * @param group_range_size = It is for checking purposes to make sure that you are
         * allocating the right amount of data for group_range.
         */
        void getGroupedRangeArray(int* group_range, int group_range_size);

        /*
         * @brief get the number of nnz in the group sparse matrix
         */
        int getGroupNNZ(){return this->group_nnz;}

        /*
         * @brief get the number of group in the group sparse matrix
         * (number of columns or rows in the aggregated sparse matrix
         */
        int getGroupSize(){return this->group_size;}

        /*
         * @brief It will return the number of nodes that are going to be touched by processing this group
         * @param: The cost of group_idx is going the be returned
         */
        double getGroupCost(int group_idx);

        /*
         * @brief Merged group into nearest group with lowest cost
         */
        void mergeSmallGroups();

    };
    //============================================ Collection Grouping ================================================
    class CollectionGrouping{
    private:
        //Output DAG variables
        int grouped_DAG_num_nodes;
        int grouped_DAG_num_edges;
        std::vector<int> group_set;
        std::vector<int> group_set_ptr;
        std::vector<int> node_to_group;
        std::vector<int> grouped_DAG_nodes;
        std::vector<int> grouped_DAG_edges;
        timing_measurement grouping_time;
        timing_measurement cmp_DAG_creation_time;
        //Input DAG variables
        int input_DAG_num_nodes;
        int input_DAG_num_edges;

        int* input_DAG_nodes;
        int* input_DAG_edges;
        //Safety Check
        bool grouping_happened;
    public:
        CollectionGrouping(int num_nodes, int* DAG_nodes, int* DAG_edges);
        void startCNGrouping();
        void startSingleChildGrouping();
        int getGroupedDAGNumEdges() const{return grouped_DAG_num_edges;};
        int getGroupedDAGNumNodes() const{return grouped_DAG_num_nodes;};
        timing_measurement getGroupingTime(){return grouping_time;}
        timing_measurement getCMPDAGCreationTime(){return cmp_DAG_creation_time;}
        void getGroupedDAG(std::vector<int>& DAG_nodes, std::vector<int>& DAG_edges);
        void getGroups(std::vector<int>& groupSet, std::vector<int>& groupPtr);
        ~CollectionGrouping() = default;
    };
}





#endif //LBC_LIB_GROUP_H
