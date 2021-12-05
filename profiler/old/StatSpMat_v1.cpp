//
// Created by labuser (Bangtian Liu) on 10/19/20.
//

#include <StatSpMat_v1.h>
#include <lbc.h>
#include <sparse_inspector.h>

#include <sparse_blas_lib.h>

namespace sym_lib{
    StatSpMat_v1::StatSpMat_v1(CSR *L, SpKerType kerType, int num_threads, int *levelPtr, int *partPtr, int *nodePtr, int levelNo, int partNo)
    {

    }

    StatSpMat_v1::StatSpMat_v1(CSC *L, SpKerType kerType, int num_threads, int *levelPtr, int *partPtr, int *nodePtr, int *groupPtr, int *groupSet, int levelNo, int partNo, int groupNo) {

    }


    void StatSpMat_v1::Setup(CSC *L, SpKerType kerType) {
     this->n = L->n;
     this->nnz = L->nnz;
     this->spkernel = kerType;
     this->NnzPerRows = L->nnz*1.0/L->n;
     this->density = (this->nnz)*1.0/(this->n * this->n);
     this->numofcores = omp_get_max_threads();

     ic0_csc_stat(L->n, L->p, L->i, this->nFlops, this->nnz_access,this->nnz_reuse);
    }


    void StatSpMat_v1::PrintData() {
     PRINT_CSV(this->n);
     PRINT_CSV(this->nnz);
     PRINT_CSV(this->density);
     PRINT_CSV(this->nFlops);
     PRINT_CSV(this->nnz_access);
     PRINT_CSV(this->nFlops*1.0/this->nnz_access);
    }



}
