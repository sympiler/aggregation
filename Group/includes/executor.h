//
// Created by labuser (Bangtian Liu) on 9/15/20.
//

#ifndef LBC_LIB_EXECUTOR_H
#define LBC_LIB_EXECUTOR_H

namespace group_cols
{
    int fs_csr_executor_sgroup(int n, int *Lp, int *Li, double *Lx, double *b, double *x, int *groupPtr, int *groupSet,  int ngroup,
                               int levels, int *levelPtr, int *levelSet) {
        if (!Lp || !Li || !x) return (0);
//#pragma omp parallel
        {
            for (int l = 0; l < levels; ++l) {
//#pragma omp for schedule(auto)
                for (int li = levelPtr[l]; li < levelPtr[l + 1]; ++li) {
                    int lidx = levelSet[li];
                    for (int k1 = groupPtr[lidx]; k1 < groupPtr[lidx + 1]; ++k1) {
                        int j1 = groupSet[k1];
                        double tmp = b[j1];
                        for (int j = Lp[j1]; j < Lp[j1 + 1] - 1; ++j) {
                            tmp -= Lx[j] * x[Li[j]];
                        }
                        x[j1] = tmp / Lx[Lp[j1 + 1] - 1];
                    }
                }
            }
        };
    }

}

#endif //LBC_LIB_EXECUTOR_H
