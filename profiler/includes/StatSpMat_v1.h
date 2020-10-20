//
// Created by labuser (Bangtian Liu) on 10/19/20.
//

#ifndef LBC_LIB_STATSPMAT_V1_H
#define LBC_LIB_STATSPMAT_V1_H

#include <def.h>
#include <Group.h>
#include <Utils.h>
#include <omp.h>
#include <executor.h>
#include <algorithm>
#include <numeric>
#include <cmath>


namespace sym_lib
{
    typedef enum
    {
        SpTrsv_CSR,
        SpTrsv_CSC,
        SpInChol_CSC
    } SpKerType; // the type for sparse kernel

    class StatSpMat_v1{
        StatSpMat_v1(CSR *L, SpKerType kerType, int num_threads, int *levelPtr, int *partPtr, int *nodePtr, int levelNo,
                     int partNo);

        int nnz; // number of non-zeros in SpMat
        int n;  // number of row or columns
        int npart; // number of partitions
        int ngroup; // number of groups when grouping is enabled
        int nnz_access; // total number of non-zeros which is accessed
        int nFlops; // Number of flops for one Sparse Kernel
        int nnz_reuse; // reused nnz across ierations
        double AverageMaxDiff; // Maximal difference per (coarsened) level. (nnz cost)
        double VarianceMaxDiff; // Variance difference per (coarsened) level. (nnz cost)
        long long int SumMaxDiff; // Sumimum of Maximal difference per (coarsened) level. (nnz cost)

        SpKerType spkernel;

        double NnzPerRows; // number of non-zeros per row/cols
        int numofcores; // number of parallelism
        double density;

    public:

        StatSpMat_v1(){};

        /**
         * @brief do profiling for coarsening method by taking the matrix L in CSR format and output of coarsening method directly
         * @param L   Triangular Sparse Matrix in CSR format
         * @param kerType kerType specifies the sparse kernel
         * @param num_threads
         * @param levelPtr The pointer to levelset
         * @param partPtr The pointer to one partitions
         * @param nodePtr The points to one node, which is one row for CSR
         * @param levelNo The number of coarsen levels
         * @param partNo The number of partitions
         * @param levelSetNo The number of levels from levelset method
         */
        StatSpMat_v1(CSC *L, SpKerType kerType, int num_threads, int *levelPtr, int *partPtr, int *nodePtr, int levelNo, int partNo);


        /**
         * @brief do profiling for coarsening method combined with grouping by taking the matrix L in CSR format and output of coarsening and grouping method directly
         * @param L  Triangular Sparse Matrix in CSR format
         * @param kerType kerType specifies the sparse kernel
         * @param num_threads
         * @param levelPtr The pointer to levelset
         * @param partPtr  The pointer to one partitions
         * @param nodePtr The points to one node, which is one group
         * @param groupPtr  Pointer to the starting location of one group
         * @param groupSet  Pointer to the column indices in one group
         * @param levelNo  Number of Coarsen Levels
         * @param partNo Number of Partitioning
         * @param levelSetNo  Number of Levels
         * @param groupNo Number of groups
         */
        StatSpMat_v1(CSC *L, SpKerType kerType, int num_threads, int *levelPtr, int *partPtr, int *nodePtr, int *groupPtr, int *groupSet,
                     int levelNo, int partNo, int groupNo);



        void Setup(CSC *L, SpKerType kerType);

        /**
         * @brief print the collected metrics
         */
        void PrintData();

    };

}


#endif //LBC_LIB_STATSPMAT_V1_H
