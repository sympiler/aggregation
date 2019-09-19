//
// Created by Kazem on 7/3/19.
//

#include <fstream>
#include <sstream>
#include "Util.h"
#include "metis_wrapper.h"


int main(int argc, char *argv[]) {
    std::string f1 = argv[1];
    int *colA, *rowA;
    double *valL;
    double *valA;
    int maxSupWid, maxCol;
    size_t n, nnzA;
    //std::vector<profilingInfo> piArray;
    if (!readMatrix(f1, n, nnzA, colA, rowA, valA))
        return -1;

    CSC *A = new CSC;
    A->ncol = A->nrow = n;
    A->p = colA;
    A->i = rowA;
    A->x = valA;
    A->nzmax = nnzA;
    A->packed = TRUE;
    A->sorted = TRUE;
    A->xtype = CHOLMOD_REAL;
    A->stype = -1;
    A->nz = nullptr;
    int *perm = new int[n]();
    int k = strtol(argv[2], nullptr, 10);

    // Uncomment the following for ordering
    // Metis_oredering(A, perm);

    // Uncomment the following for partitining
    Metis_partitioning(A,perm,k);

    for (int i = 0; i < n; ++i) {
        std::cout << perm[i] << "\n";
    }
    delete[]perm;
    allocateAC(A, 0, 0, 0, FALSE);
    return 0;
}