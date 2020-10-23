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

    void ic0_csc(int n, double *val, int * colPtr, int *rowIdx)
    {
     int i, k,l,m;
     double temp;
     for (i = 0; i < n ; i++){
      temp = val[colPtr[i]];
      val[colPtr[i]] = val[colPtr[i]]/sqrt(temp);//S1

      for (m = colPtr[i] + 1; m < colPtr[i+1]; m++){
       val[m] = val[m] / val[colPtr[i]];//S2
      }

      for (m = colPtr[i] + 1; m < colPtr[i+1]; m++) { // j
       for (k = colPtr[rowIdx[m]] ; k < colPtr[rowIdx[m]+1]; k++){ // i
        for ( l = m; l < colPtr[i+1] ; l++){
         if (rowIdx[l] == rowIdx[k] ){
//                            if(rowIdx[l+1] <= rowIdx[k]){
          val[k] -= val[m]* val[l]; //S3   a(i,j) = a(i,j) - a(i,k) * a(j,k)
//                            }
         }
         else if(rowIdx[l] > rowIdx[k])
          break;
        }
       }
      }
     }
    }


    void spico_csc_levelset(int n, const int *Lp, const int *Li,  double *Lx,
                            int levels, const int *levelPtr,
                            const int *levelSet)
    {
#pragma omp parallel
     {
      for (int s = 0; s < levels; s++) {
#pragma omp for schedule(auto)
       for (int k1 = levelPtr[s]; k1 < levelPtr[s + 1]; ++k1) {
        int i = levelSet[k1];
        double temp = Lx[Lp[i]];
        Lx[Lp[i]] = Lx[Lp[i]]/sqrt(temp);//S1

        for (int m = Lp[i] + 1; m < Lp[i+1]; m++){
         Lx[m] = Lx[m] / Lx[Lp[i]];//S2
        }

        for (int m = Lp[i] + 1; m < Lp[i+1]; m++) {
         for (int k = Lp[Li[m]] ; k < Lp[Li[m]+1]; k++){
          for ( int l = m; l < Lp[i+1] ; l++){
           if (Li[l] == Li[k] ){
//                            if(rowIdx[l+1] <= rowIdx[k]){
#pragma omp atomic
            Lx[k] -= Lx[m]* Lx[l]; //S3
//                            }
           }
           else if(Li[l] > Li[k])
            break;
          }
         }
        }
       }
      }
     };
    }

    void spico_csc_group_levelset(int n, const int *Lp, const int *Li,  double *Lx,
                                  int levels, const int *levelPtr,
                                  const int *levelSet, int *groupPtr, int *groupSet)
    {
#pragma omp parallel
     {
      for (int l = 0; l < levels; ++l){
#pragma omp for schedule(auto)
       for (int li = levelPtr[l]; li < levelPtr[l + 1]; ++li) {
        int lidx = levelSet[li];
        for (int k1 = groupPtr[lidx]; k1 < groupPtr[lidx + 1]; ++k1)
        {
         int i = groupSet[k1];
         double temp = Lx[Lp[i]];
         Lx[Lp[i]] = Lx[Lp[i]]/sqrt(temp);//S1

         for (int m = Lp[i] + 1; m < Lp[i+1]; m++){
          Lx[m] = Lx[m] / Lx[Lp[i]];//S2
         }

         for (int m = Lp[i] + 1; m < Lp[i+1]; m++) {
          for (int k = Lp[Li[m]] ; k < Lp[Li[m]+1]; k++){
           for ( int l = m; l < Lp[i+1] ; l++){
            if (Li[l] == Li[k] ){
//                            if(rowIdx[l+1] <= rowIdx[k]){
#pragma omp atomic
             Lx[k] -= Lx[m]* Lx[l]; //S3
//                            }
            }
            else if(Li[l] > Li[k])
             break;
           }
          }
         }
        }
       }
      }
     };
    }

    void spic0_csc_group_lbc(int n, int *Lp, int *Li, double *Lx,
                             int level_no, int *level_ptr,
                             int *par_ptr, int *partition,  int *groupPtr, int *groupSet)
    {
#pragma omp parallel
     {
      for (int i1 = 0; i1 < level_no; ++i1) {
#pragma omp  for schedule(auto)
       for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1) {
        for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1) {
         int p = partition[k1];

         for (int k = groupPtr[p]; k < groupPtr[p+1]; ++k) {
          int i = groupSet[k];
          double temp = Lx[Lp[i]];
          Lx[Lp[i]] = Lx[Lp[i]]/sqrt(temp);//S1

          for (int m = Lp[i] + 1; m < Lp[i+1]; m++){
           Lx[m] = Lx[m] / Lx[Lp[i]];//S2
          }

          for (int m = Lp[i] + 1; m < Lp[i+1]; m++) {
           for (int k = Lp[Li[m]] ; k < Lp[Li[m]+1]; k++){
            for ( int l = m; l < Lp[i+1] ; l++){
             if (Li[l] == Li[k] ){
//                            if(rowIdx[l+1] <= rowIdx[k]){
#pragma omp atomic
              Lx[k] -= Lx[m]* Lx[l]; //S3
//                            }
             }
             else if(Li[l] > Li[k])
              break;

            }
           }
          }
         }
        }
       }
      }
     };
    }

    void spic0_csr_lbc(int n, double *Lx, int *Lp, int *Li,
                       int level_no, int *level_ptr,
                       int *par_ptr, int *partition) {

#pragma omp parallel
     {
      for (int i1 = 0; i1 < level_no; ++i1) {
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
     };
    }

    void spico_csc_lbc(int n, double *Lx, int *Lp, int *Li,
                       int level_no, int *level_ptr,
                       int *par_ptr, int *partition)
    {
#pragma omp parallel
     {
      for (int i1 = 0; i1 < level_no; ++i1) {
#pragma omp  for schedule(auto)
       for (int j1 = level_ptr[i1]; j1 < level_ptr[i1 + 1]; ++j1) {
        for (int k1 = par_ptr[j1]; k1 < par_ptr[j1 + 1]; ++k1) {
         int i = partition[k1];
         double temp = Lx[Lp[i]];
         Lx[Lp[i]] = Lx[Lp[i]]/sqrt(temp);//S1

         for (int m = Lp[i] + 1; m < Lp[i+1]; m++){
          Lx[m] = Lx[m] / Lx[Lp[i]];//S2
         }

         for (int m = Lp[i] + 1; m < Lp[i+1]; m++) {
          for (int k = Lp[Li[m]] ; k < Lp[Li[m]+1]; k++){
           for ( int l = m; l < Lp[i+1] ; l++){
            if (Li[l] == Li[k] ){
//                            if(rowIdx[l+1] <= rowIdx[k]){
#pragma omp atomic
             Lx[k] -= Lx[m]* Lx[l]; //S3
//                            }
            }
           }
          }
         }
        }
       }
      }
     };
    }



}