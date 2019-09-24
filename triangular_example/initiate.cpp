//
// Created by George Huang on 2019-07-08.
//

#include "initiate.h"

void process_matrix(FILE *mf, CSC *&A) {
    int M;
    if (!mf) exit(1);

    MM_typecode mcode;
    if (mm_read_banner(mf, &mcode) != 0) {
        std::cerr << "Error processing matrix banner\n";
        fclose(mf);
        exit(1);
    }

    A = new CSC;

    int n, nz;
    if (mm_read_mtx_crd_size(mf, &n, &M, &nz) != 0) exit(1);
    load_matrix(mf, n, nz, A->i, A->p, A->x);

    A->nrow = A->ncol = n;
    A->nzmax = nz;
    A->packed = TRUE;
    A->sorted = TRUE;
    A->stype = -1;
    A->nz = nullptr;
}


void load_matrix(FILE *f, int n, int nz, int *&Ai, int *&Ap, double *&Ax) {
    int *J = new int[nz]();

    Ai = new int[nz]();
    Ap = new int[nz + 1]();
    Ax = new double[nz]();

    // Copy matrix data into COO format
    for (int i = 0; i < nz; i++) {
        if(fscanf(f, "%d %d %lg\n", (Ai) + i, &J[i], (Ax) + i) == EOF)
            exit(1);
        (Ai)[i]--;
        J[i]--;
    }
    to_csc(n, nz, J, Ap);

    delete[]J;
}


void to_csc(int n, int nz, int *J, int *Lp) {
    if (!J || !Lp) {
        fprintf(stderr, "Error converting to CSC format\n");
        exit(1);
    }

    int i = 0, j;
    int index = 0, cur;
    for (; i < nz; i++) {
        Lp[index] = i;
        cur = J[i];
        for (j = i + 1; j < nz; j++) {
            if (J[j] != cur)
                break;
            else
                i++;
        }
        index += 1;
    }
    Lp[n] = nz;
}


void rhsInit(int n, int *Lp, int *Li, double *Lx, double *b) {
    for (int c = 0; c < n; ++c) {
        for (int cc = Lp[c]; cc < Lp[c + 1]; ++cc) {
            b[Li[cc]] += Lx[cc];
        }
    }
}


void load_vector_sparse(FILE *f, int n, int vz, double *b, int **index) {
    double value;
    int temp;

    *index = new int[vz]();
    for (int i = 0; i < vz; i++) {
        if(fscanf(f, "%d %d %lg\n", (*index) + i, &temp, &value) == EOF)
            exit(1);
        (*index)[i]--;
        b[(*index)[i]] = value;
    }
}


void load_vector_dense(FILE *f, int n, double *b) {
    for (int i = 0; i < n; i++)
        if(fscanf(f, " %lg\n", b + i) == EOF)
            exit(1);
}

int allocateLC(BCSC *L, int sw) {
    int sNo = L->nsuper;
    if (sw) {
        L->super = new int[sNo + 1]();
//  L->col2Sup = new int[L->n](); //TODO HACK
        L->p = new size_t[L->ssize + 1]();
        L->pi = new size_t[sNo + 1]();
        L->i_ptr = new size_t[L->xsize + 1](); // index pointers
        L->i = new int[sNo + 1]();//Nothing for now
        // L->px = new int[L->xsize]();
        L->s = new int[L->xsize]();//Index values
        //L->sParent = new int[sNo](); //TODO HACK
        // L->x = new double[L->xsize]();
        L->is_ll = TRUE;
        L->xtype = CHOLMOD_REAL;
        L->is_super = TRUE;

    } else {
        delete[]L->super;
        delete[]L->col2Sup;
        delete[]L->p;
        delete[]L->pi;
        delete[]L->i;
        delete[]L->i_ptr;
        // delete []L->px;
        delete[]L->s;
        delete[]L->sParent;
        delete[]L->Parent;
        delete[]L->IPerm;
        delete[]L->Perm;
//  delete []L->x;
    }
    return 0;
}
//
//int allocateAC(CSC *A, int nrow, int nnz, int sytpe, int sw) {
//    if (sw) {
//        A->nrow = A->ncol = nrow;
//        A->nzmax = nnz;
//        A->stype = sytpe;
//        A->xtype = CHOLMOD_REAL;//TODO removed later
//        A->packed = TRUE; // Always
//        A->p = new int[nrow + 1]();
//        A->i = new int[nnz]();
//        A->x = new double[nnz]();
//        A->nz = nullptr;
//    } else {
//        delete[]A->p;
//        delete[]A->i;
//        delete[]A->x;
//    }
//
//    return 0;
//}


int transpose(CSC *L, CSC *U) {
//    int n = L->nrow;
//    int nnz = L->nzmax;
//    int *Lp = L->p;
//    int *Li = L->i;
//    double *Lx = L->x;
//
//    auto *Up = new int[n+1]();
//    auto *Ui = new int[nnz]();
//    auto *Ux = new double[nnz]();
//
//    auto *rowCnt = new int[n]();
//    for(int i = 0; i < n; i++)
//        rowCnt[Li[i]]++;
//
//
//    delete[]rowCnt;
//
//    U->nrow = U->ncol = n;
//    U->nzmax = nnz;
//    U->p = Up;
//    U->i = Ui;
//    U->x = Ux;

    return 0;
}
