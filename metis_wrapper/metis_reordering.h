//
// Created by george on 2019-09-12.
//

#ifndef METIS_REORDERING_H
#define METIS_REORDERING_H


/**
 * A does not get deallocated!!
 *
 * @param A
 * @param new_A
 * @param n_par
 */
void metis_reordering(CSC *A, CSC *&new_A, int n_par) {
    int n = A->ncol;
    int nnz = A->nzmax;
    int *Ap = A->p;
    int *Ai = A->i;
    double *Ax = A->x;

    int *perm = new int[n]();
    Metis_partitioning(A, perm, n_par);

    int *count = new int[n_par]();
    for(int i = 0; i < n; i++)
        count[perm[i]]++;
    int **temp = new int*[n_par](); // column #
    for(int i = 0; i < n_par; i++) {
        temp[i] = new int[count[i]]();
    }

    int *counter = new int[n_par]();
    for(int i = 0; i < n; i++) {
        int par = perm[i];
        temp[par][counter[par]] = i;
        counter[par]++;
    }
    delete[]perm;
    delete[]counter;

    new_A = new CSC;
    new_A->ncol = new_A->nrow = A->ncol;
    new_A->nzmax = A->nzmax;
    new_A->stype = A->stype;
    new_A->sorted = A->sorted;
    new_A->packed = A->packed;
    new_A->nz = A->nz;
    new_A->p = new int[n+1]();
    new_A->i = new int[nnz]();
    new_A->x = new double[nnz]();

    new_A->p[0] = 0;
    int ptr_count = 0;
    int row_count = 0;
    for(int i = 0; i < n_par; i++) {
        for(int j = 0; j < count[i]; j++) {
            int col = temp[i][j];
            new_A->p[ptr_count+1] = new_A->p[ptr_count] + Ap[col+1] - Ap[col];
            ptr_count++;

            for(int k = Ap[col]; k < Ap[col+1]; k++) {
                new_A->i[row_count] = Ai[k];
                new_A->x[row_count] = Ax[k];
                row_count++;
            }
        }
    }

    for(int i = 0; i < n_par; i++)
        delete[](temp[i]);
    delete[]temp;
    delete[]count;
}

#endif //METIS_REORDERING_H
