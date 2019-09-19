//
// Created by Kazem on 7/3/19.
//

#ifndef PROJECT_METIS_WRAPPER_H
#define PROJECT_METIS_WRAPPER_H

#include <iostream>
#include "metis.h"
#include "utils/initiate.h"
#include "Transpose.h"


int Metis_oredering
        (
                /* ---- input ---- */
                CSC *A,/* matrix to order and analyze */
                /* ---- output ---- */
                int *Lperm
        ) {
    CSC *ATrans;
    int status = 0;
    unsigned long nnzFull = A->nzmax * 2;//Symmetric case
    ATrans = ptranspose(A, 0, NULL, NULL, 0, status);

    //Making the graph for passing it to metis, it should have
    //both upper and lower parts
    //allocateAC(AFull,ncol,nnzFull,0,TRUE);
    idx_t options1[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options1);
    int ncol = A->ncol;
    idx_t *AFullp = new idx_t[ncol + 1]();
    idx_t *AFulli = new idx_t[nnzFull]();
    idx_t ncolIDXT = ncol;
    idx_t *weigt = new idx_t[ncol];
    idx_t *LpermIDX = new idx_t[ncol];
    idx_t *ILpermIDX = new idx_t[ncol];
    for (int i = 0; i < ncol; ++i) {
        LpermIDX[i] = 0;
        ILpermIDX[i] = 0;
        weigt[i] = 1;
    }
    AFullp[0] = 0;
    int ind = 1;
    for (int i = 0; i < ncol; ++i) {
        int nnzOfCurCol = ATrans->p[i + 1] - ATrans->p[i] - 1;
        nnzOfCurCol += A->p[i + 1] - A->p[i] - 1;
        assert(nnzOfCurCol >= 0);
        AFullp[i + 1] = (long int) AFullp[i] + nnzOfCurCol;
        //copying Upper part, ignoring diagonal
        int base = AFullp[i];
        for (int j = ATrans->p[i], k = 0; j < ATrans->p[i + 1] - 1; ++j, ++k) {
            AFulli[base + k] = (long int) ATrans->i[j];
        }
        //copying L part
        base += ATrans->p[i + 1] - ATrans->p[i] - 1;
        for (int j = A->p[i] + 1, k = 0; j < A->p[i + 1]; ++j, ++k) {
            AFulli[base + k] = (long int) A->i[j];
        }
    }

    int retMet = METIS_NodeND(&ncolIDXT, AFullp, AFulli, NULL, options1,
                              LpermIDX, ILpermIDX);
    assert(retMet == METIS_OK);
    if (retMet != METIS_OK) {
        std::cout << " " << retMet << "\n";
        return 0;
    }
    for (int i = 0; i < ncol; ++i) {
        Lperm[i] = LpermIDX[i];
        //std::cout<<Lperm[i];
    }
    allocateAC(ATrans, ATrans->nrow, ATrans->nzmax, ATrans->stype, false);
    METIS_Free(AFullp);
    METIS_Free(AFulli);
    METIS_Free(weigt);
    METIS_Free(LpermIDX);
    METIS_Free(ILpermIDX);

    return 1;
}

int Metis_partitioning
        (
                /* ---- input ---- */
                CSC *A,/* matrix to order and analyze */
                /* ---- output ---- */
                int *Lpart,
                int k
        ) {
    CSC *ATrans;
    int status = 0;
    unsigned long nnzFull = A->nzmax * 2;//Symmetric case
    ATrans = ptranspose(A, 0, nullptr, nullptr, 0, status);

    //Making the graph for passing it to metis, it should have
    //both upper and lower parts
    //allocateAC(AFull,ncol,nnzFull,0,TRUE);
    idx_t options1[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options1);
    int ncol = A->ncol;
    auto *AFullp = new idx_t[ncol + 1]();
    auto *AFulli = new idx_t[nnzFull]();
    idx_t ncolIDXT = ncol;
    auto *weigt = new idx_t[ncol];
    auto *partIDX = new idx_t[ncol];
    idx_t nconIDX = 1;
    idx_t nparsIDX = k;
    idx_t objvalIDX = 0;
    for (int i = 0; i < ncol; ++i) {
        partIDX[i] = 0;
        weigt[i] = 1;
    }
    AFullp[0] = 0;
    int ind = 1;
    for (int i = 0; i < ncol; ++i) {
        int nnzOfCurCol = ATrans->p[i + 1] - ATrans->p[i] - 1;
        nnzOfCurCol += A->p[i + 1] - A->p[i] - 1;
        AFullp[i + 1] = (long int) AFullp[i] + nnzOfCurCol;
        //copying Upper part, ignoring diagonal
        int base = AFullp[i];
        for (int j = ATrans->p[i], k = 0; j < ATrans->p[i + 1] - 1; ++j, ++k) {
            AFulli[base + k] = (long int) ATrans->i[j];
        }
        //copying L part
        base += ATrans->p[i + 1] - ATrans->p[i] - 1;
        for (int j = A->p[i] + 1, k = 0; j < A->p[i + 1]; ++j, ++k) {
            AFulli[base + k] = (long int) A->i[j];
        }
    }


    idx_t retMet = METIS_PartGraphKway(&ncolIDXT, &nconIDX, AFullp, AFulli, nullptr, nullptr, nullptr, &nparsIDX, nullptr,
                                     nullptr, options1, &objvalIDX, partIDX);
    assert(retMet == METIS_OK);
    if (retMet != METIS_OK) {
        std::cout << " " << retMet << "\n";
        return 0;
    }
    for (int i = 0; i < ncol; ++i) {
        Lpart[i] = partIDX[i];
    }
    allocateAC(ATrans, ATrans->nrow, ATrans->nzmax, ATrans->stype, false);
    METIS_Free(AFullp);
    METIS_Free(AFulli);
    METIS_Free(weigt);
    METIS_Free(partIDX);

    return 1;
}

#endif //PROJECT_METIS_WRAPPER_H
