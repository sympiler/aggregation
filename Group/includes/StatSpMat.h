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
//#define CSV_LOG
//
//#ifdef CSV_LOG
//#define PRINT_CSV(x) std::cout <<(x)<<","
//#else
//#define PRINT_CSV(x)
//#endif

    typedef enum
    {
        SpTrsv_CSR,
        SpTrsv_CSC,
        SpInChol_CSC
    } SpKerType;

    class StatSpMat{
        int nnz; // number of non-zeros in SpMat
        int n;  // number of row or columns
        int ngroup; // number of groups when grouping is enabled
        int nnz_access; // total number of non-zeros which is accessed
        int nFlops; // Number of flops for one Sparse Kernel
        double AverageMaxDiff; // Maximal difference per (coarsened) level. (nnz cost)
        double VarianceMaxDiff; // Variance difference per (coarsened) level. (nnz cost)
        long long int SumMaxDiff; // Sumimum of Maximal difference per (coarsened) level. (nnz cost)
        int numofcores;


        double t_serial;
        double t_level;
        double t_lbc;
        double t_group_level;

        int nlevels;
        int num_sys;
        double NnzPerRows;
        int nnz_reuse;
        double nnzPerLevels;
        double averParallelism;
        SpKerType spkernel;

    public:
        StatSpMat(CSR *L, SpKerType kerType, int num_threads, int blksize);

        StatSpMat(CSR *L, SpKerType kerType, int num_threads, int lparm, int divrate);


        StatSpMat(CSC *L, SpKerType kerType, int num_threads);

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
