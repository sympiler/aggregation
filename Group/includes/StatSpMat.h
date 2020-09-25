//
// Created by labuser (Bangtian Liu) on 9/15/20.
//

#ifndef LBC_LIB_STATSPMAT_H
#define LBC_LIB_STATSPMAT_H

#include <def.h>
#include <Group.h>
#include <Utils.h>
#include <omp.h>
#include <executor.h>
#include <algorithm>
#include <numeric>
#include <cmath>


namespace group_cols
{
    typedef enum
    {
        SpTrsv_CSR,
        SpTrsv_CSC,
        SpInChol_CSC
    } SpKerType; // the type for sparse kernel

    class StatSpMat{
        int nnz; // number of non-zeros in SpMat
        int n;  // number of row or columns
        int ngroup; // number of groups when grouping is enabled
        int nnz_access; // total number of non-zeros which is accessed
        int nFlops; // Number of flops for one Sparse Kernel
        double AverageMaxDiff; // Maximal difference per (coarsened) level. (nnz cost)
        double VarianceMaxDiff; // Variance difference per (coarsened) level. (nnz cost)
        long long int SumMaxDiff; // Sumimum of Maximal difference per (coarsened) level. (nnz cost)
        int numofcores; // number of parallelism


        double t_serial; // running time for serial code
        double t_level; // running time for levelset parallel method
        double t_lbc; // running time for lbc parallel method
        double t_group_level; // runing time for level parallel method combined with grouping method

        int nlevels; // number of levels
        int num_sys; // number of synchronization
        double NnzPerRows; // number of non-zeros per row/cols
        int nnz_reuse; // reused nnz across ierations
        double nnzPerLevels; // number of non-zeros per (coarsen) levels
        double averParallelism; // average parallelism
        SpKerType spkernel;

    public:
        /**
         * @brief do profiling for levelset method by taking the matrix L in CSR format directly
         * @param L Triangular Sparse Matrix in CSR format
         * @param kerType specifies the sparse kernel
         * @param num_threads number of threads
         * @param blksize the grouping parameter for group method, no grouping can be enabled by setting blksize to 1
         */
        StatSpMat(CSR *L, SpKerType kerType, int num_threads, int blksize);

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
        StatSpMat(CSR *L, SpKerType kerType, int num_threads, int *levelPtr, int *partPtr, int *nodePtr, int levelNo, int partNo, int levelSetNo);

        /**
         *
         * @param L SpMat in CSR format
         * @param DAG depdence information from L
         * @param kerType kerType specifies the sparse kernel
         * @param num_threads number of
         * @param blksize parameters for grouping method, enble by setting blksize>1
         */
        StatSpMat(CSR *L, std::vector<std::vector<int>> DAG, SpKerType kerType, int num_threads, int blksize=1);


        void Setup(CSR *L, SpKerType kerType);

        StatSpMat(CSC *L, SpKerType kerType, int num_threads);

        /**
         * @brief print the collected metrics
         */
        void PrintData();

        void set_seq_time(double serial_time){
            t_serial=serial_time;
        }

        void set_level_time(double  level_time){
            t_level=level_time;
        }

        void set_lbc_time(double lbc_time){
            t_lbc=lbc_time;
        }

        void set_glevel_time(double glevel_time){
            t_group_level=glevel_time;
        }

    };

}




#endif //LBC_LIB_STATSPMAT_H
