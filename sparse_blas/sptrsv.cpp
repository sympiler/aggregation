//
// Created by george on 2019-10-09.
//
#include <omp.h>

#include "sparse_blas_lib.h"

namespace sym_lib {

    void sptrsv_csr(int n, int *Lp, int *Li, double *Lx, double *x) {
     int i, j;
     for (i = 0; i < n; i++) {
      for (j = Lp[i]; j < Lp[i + 1] - 1; j++) {
       x[i] -= Lx[j] * x[Li[j]];
      }
      x[i] /= Lx[Lp[i + 1] - 1];
     }
    }



    void sptrsv_csr_levelset(int n, const int *Lp, const int *Li, const double *Lx,
                             double *x,
                             int levels, const int *levelPtr,
                             const int *levelSet) {
#pragma omp parallel
     {
      for (int l = 0; l < levels; l++) {
//#pragma omp parallel for default(shared) schedule(auto)
#pragma omp for schedule(auto)
       for (int k = levelPtr[l]; k < levelPtr[l + 1]; ++k) {
        int i = levelSet[k];
        for (int j = Lp[i]; j < Lp[i + 1] - 1; j++) {
         x[i] -= Lx[j] * x[Li[j]];
        }
        x[i] /= Lx[Lp[i + 1] - 1];
       }
      }


     };

    }

    void sptrsv_csr_levelset_seq(int n, const int *Lp, const int *Li, const double *Lx,
                                 double *x,
                                 int levels, const int *levelPtr,
                                 const int *levelSet) {
     for (int l = 0; l < levels; l++) {
      for (int k = levelPtr[l]; k < levelPtr[l + 1]; ++k) {
       int i = levelSet[k];
       for (int j = Lp[i]; j < Lp[i + 1] - 1; j++) {
        x[i] -= Lx[j] * x[Li[j]];
       }
       x[i] /= Lx[Lp[i + 1] - 1];
      }
     }
    }


    void sptrsv_csr_lbc_seq(int n, int *Lp, int *Li, double *Lx, double *x,
                            int level_no, int *level_ptr,
                            int *par_ptr, int *partition) {
//#pragma omp parallel
     {
      for (int i1 = 0; i1 < level_no; ++i1) {
//#pragma omp  for schedule(auto)
//                printf("len=%d\n", level_ptr[i1+1]-level_ptr[i1]);
       for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1) {
        for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1) {
         int i = partition[k1];
         for (int j = Lp[i]; j < Lp[i + 1] - 1; j++) {
          x[i] -= Lx[j] * x[Li[j]];
         }
         x[i] /= Lx[Lp[i + 1] - 1];
        }
       }
      }
     };
    }

    void sptrsv_csr_lbc(int n, int *Lp, int *Li, double *Lx, double *x,
                        int level_no, int *level_ptr,
                        int *par_ptr, int *partition) {
#pragma omp parallel
     {
      for (int i1 = 0; i1 < level_no; ++i1) {
#pragma omp  for schedule(auto)
       for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1) {
        for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1) {
         int i = partition[k1];
         for (int j = Lp[i]; j < Lp[i + 1] - 1; j++) {
          x[i] -= Lx[j] * x[Li[j]];
         }
         x[i] /= Lx[Lp[i + 1] - 1];
        }
       }
      }
     };
    }

    void sptrsv_csr_group_lbc(int n, int *Lp, int *Li, double *Lx, double *x,
                              int level_no, int *level_ptr,
                              int *par_ptr, int *partition,  int *groupPtr, int *groupSet) {
#pragma omp parallel
     {
      for (int i1 = 0; i1 < level_no; ++i1) {
#pragma omp  for schedule(auto)
       for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1) {
        for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1) {
         int p = partition[k1];

         for (int k = groupPtr[p]; k < groupPtr[p+1]; ++k) {
          int i = groupSet[k];
          for (int j = Lp[i]; j < Lp[i + 1] - 1; j++) {
           x[i] -= Lx[j] * x[Li[j]];
          }
          x[i] /= Lx[Lp[i + 1] - 1];
         }
        }
       }
      }
     };
    }
    ///=============================================================================
    ///============================= HDAGG SPARSE BLASS ============================
    ///=============================================================================
    //=========================== Left Looking SpTrSv ==========================
    void sptrsv_csr_levelset(int n, const int *Lp, const int *Li, const double *Lx,
                             int levels, const int *levelPtr, const int *levelSet,
                             double *x)
    {
#pragma omp parallel
      {
        for (int l = 0; l < levels; l++)
        {
#pragma omp for schedule(auto)
          for (int k = levelPtr[l]; k < levelPtr[l + 1]; ++k)
          {
            int i = levelSet[k];
            for (int j = Lp[i]; j < Lp[i + 1] - 1; j++)
            {
              x[i] -= Lx[j] * x[Li[j]];//S1
            }
            x[i] /= Lx[Lp[i + 1] - 1]; //S2
          }
        }
      }
    }

    void sptrsv_csr_group_levelset(int *Lp, int *Li, double *Lx, double *x,
                                   int level_no, int *level_ptr, int *level_set,
                                   int *groupPtr, int *groupSet)
    {
#pragma omp parallel
      {
        for (int i1 = 0; i1 < level_no; ++i1)
        {
#pragma omp for schedule(auto)
          for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1)
          {
            int group_idx = level_set[j1];
            for (int k = groupPtr[group_idx]; k < groupPtr[group_idx + 1]; ++k)
            {
              int i = groupSet[k];
              for (int j = Lp[i]; j < Lp[i + 1] - 1; j++)
              {
                x[i] -= Lx[j] * x[Li[j]];
              }
              x[i] /= Lx[Lp[i + 1] - 1];
            }
          }
        }
      }
    }
}