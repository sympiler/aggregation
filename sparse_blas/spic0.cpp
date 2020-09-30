//
// Created by kazem on 2020-04-23.
//
#include <cmath>
#include <cassert>

namespace sym_lib {

 void spic0_csr(int n, double *Lx, int *Lp, int *Li)  {
  for (int i = 0; i < n; i++) {
   assert(Lx[Lp[i+1]-1] >0);
   Lx[Lp[i+1]-1] = sqrt(Lx[Lp[i+1]-1]);//S1

   for (int m = Lp[i]; m < Lp[i + 1]-1; m++) {
    Lx[m] = Lx[m] / Lx[Lp[i+1]-1];//S2
   }

   for (int m = Lp[i] ; m < Lp[i + 1]-1; m++) {
    for (int k = Lp[Li[m]]; k < Lp[Li[m] + 1]; k++) {
     for (int l = m; l < Lp[i + 1]-1; l++) {
       if (Li[l] == Li[k] && Li[l + 1] <= Li[k]) {
        Lx[k] -= Lx[m] * Lx[l]; //S3
       }
     }
    }
   }
  }
 }

 void spic0_csr_lbc(int n, double *Lx, int *Lp, int *Li,
                  int level_no, int *level_ptr,
                  int *par_ptr, int *partition) {
  for (int i1 = 0; i1 < level_no; ++i1) {
#pragma omp parallel
   {
#pragma omp  for schedule(auto)
    for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1) {
     for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1) {
      int i = partition[k1];
      Lx[Lp[i+1]-1] = sqrt(Lx[Lp[i+1]-1]);//S1

      for (int m = Lp[i]; m < Lp[i + 1]-1; m++) {
       Lx[m] = Lx[m] / Lx[Lp[i+1]-1];//S2
      }

      for (int m = Lp[i] ; m < Lp[i + 1]-1; m++) {
       for (int k = Lp[Li[m]]; k < Lp[Li[m] + 1]; k++) {
        for (int l = m; l < Lp[i + 1]-1; l++) {
         if (Li[l] == Li[k] && Li[l + 1] <= Li[k]) {
          Lx[k] -= Lx[m] * Lx[l]; //S3
         }
        }
       }
      }
     }}}}//LBC outermost
 }


}