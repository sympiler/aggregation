//
// Created by kazem on 9/18/19.
//

#ifndef CROSSKERNEL_TRIANGULAR_CSC_H
#define CROSSKERNEL_TRIANGULAR_CSC_H


/*
 ****** Serial implementation
 */
int lsolve(int n, int *Lp, int *Li, double *Lx, double *x) {
    int p, j;
    if (!Lp || !Li || !x) return (0);                     /* check inputs */
    for (j = 0; j < n; j++) {
        x[j] /= Lx[Lp[j]];
        for (p = Lp[j] + 1; p < Lp[j + 1]; p++) {
            x[Li[p]] -= Lx[p] * x[j];
        }
    }
    return (1);
}


/*
 ****** Parallel
 */
int lsolvePar(int n, int *Lp, int *Li, double *Lx, double *x,
              int levels, int *levelPtr, int *levelSet, int chunk) {
    if (!Lp || !Li || !x) return (0);                     /* check inputs */
    for (int l = 0; l < levels; ++l) {
        int li = 0;
#pragma omp parallel for \
   default(shared) private(li)  \
   schedule(auto)
        for (li = levelPtr[l]; li < levelPtr[l + 1]; ++li) {
            int j = levelSet[li];
            x[j] /= Lx[Lp[j]];
            for (int p = Lp[j] + 1; p < Lp[j + 1]; p++) {
                double tmp = Lx[p] * x[j];
                int idx = Li[p];
#pragma omp atomic
                x[idx] -= tmp;
            }
        }
    }
    return (1);
}


void LLcross(int n, const int *Lp, const int *Li, const double *Lx, double *x, const int *levelPtr, const int *levels, int n_lev) {
    int i, j, p, index;
    int temp, base;

    for (i = 0; i < n_lev; i++) {
#pragma omp parallel default(shared) private(index, j, temp, base, p)
        {
#pragma omp for schedule(auto)
            for (j = levelPtr[i]; j < levelPtr[i + 1]; j++) {
                index = levels[j];

                int kern_num = index / n;
                base = kern_num * n;
                temp = index - base;
                if (base != 0)
                    x[index] += x[index - n];

                x[index] /= Lx[Lp[temp]];
                for (p = Lp[temp] + 1; p < Lp[temp + 1]; p++) {
                    int idx = Li[p] + base;
                    double tmp = Lx[p] * x[index];

#pragma omp atomic
                    x[idx] -= tmp;
                }
            }
        }
    }
}

/*
 ****** Parallel H2
 */
int lsolveParH2(int n, int *Lp, int *Li, double *Lx, double *x,
                int levels, int *levelPtr, int *levelSet,
                int parts, int *parPtr, int *partition,
                int chunk) {
    if (!Lp || !Li || !x) return (0);                     /* check inputs */
    for (int i1 = 0; i1 < levels; ++i1) {
#pragma omp parallel //shared(lValues)//private(map, contribs)
        {
#pragma omp  for schedule(auto)
            for (int j1 = levelPtr[i1]; j1 < levelPtr[i1 + 1]; ++j1) {
                for (int k1 = parPtr[j1]; k1 < parPtr[j1 + 1]; ++k1) {
                    int j = partition[k1];
                    x[j] /= Lx[Lp[j]];
                    for (int p = Lp[j] + 1; p < Lp[j + 1]; p++) {
                        double tmp = Lx[p] * x[j];
                        int idx = Li[p];
#pragma omp atomic
                        x[idx] -= tmp;
                    }
                }
            }
        }
    }
    return (1);
}

int llcrossParH2(int n, int *Lp, int *Li, double *Lx, double *x,
                int levels, int *levelPtr, int *levelSet,
                int parts, int *parPtr, int *partition,
                int chunk) {
//    if (!Lp || !Li || !x) return (0);                     /* check inputs */
    int i, j, p, index;
    int temp, base;
    for (int i1 = 0; i1 < levels; ++i1) {
#pragma omp parallel default(shared) private(index, j, temp, base, p)//shared(lValues)//private(map, contribs)
        {
#pragma omp  for schedule(auto)
            for (int j1 = levelPtr[i1]; j1 < levelPtr[i1 + 1]; ++j1) {
                for (int k1 = parPtr[j1]; k1 < parPtr[j1 + 1]; ++k1) {
                    index = partition[k1];

                    int kern_num = index / n;
                    base = kern_num * n;
                    temp = index - base;
                    if (base != 0)
                        x[index] += x[index - n];

                    x[index] /= Lx[Lp[temp]];
                    for (p = Lp[temp] + 1; p < Lp[temp + 1]; p++) {
                        int idx = Li[p] + base;
                        double tmp = Lx[p] * x[index];
#pragma omp atomic
                        x[idx] -= tmp;
                    }
                }
            }
        }
    }
    return (1);
}

#endif //CROSSKERNEL_TRIANGULAR_CSC_H
