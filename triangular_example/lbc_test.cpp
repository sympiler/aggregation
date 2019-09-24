//
// Created by kazem on 9/18/19.
//


#include <iostream>
#include <chrono>
#include <algorithm>
#include <omp.h>
#include "Transpose.h"

#include "mergeGraph.h"

#include "metis_wrapper.h"
#include "Triangular_CSC.h"
#include "Inspection_Level.h"
#include "lbc_csc.h"
#include "Util.h"
#include "initiate.h"

void spmv_csc(int n, int *Ap, int *Ai, double *Ax, double *x, double *y) {
    int p, j;
    double *temp = x;
    for (j = 0; j < n; j++) {
        for (p = Ap[j]; p < Ap[j + 1]; p++) {
            y[Ai[p]] += Ax[p] * *temp;
        }
        temp++;
    }
}

#define CPUTIME (SuiteSparse_time ( ))
#define CSC_SER
#define CSC_LVL
#define CSC_LBC
#define NUM_TEST 9
#undef DEBUG
//#define FLOPCNT
#define METIS 1

int testTriangular(size_t n, const double *x, double epsilon = 1e-9);

int compareDouble(const void *a, const void *b) {
    if(*(double*)a < *(double*)b)
        return -1;
    else if(*(double*)a == *(double*)b)
        return 0;
    else
        return 1;
}

double get_median(double *list, int num) {
    qsort(list, num, sizeof(double), compareDouble);
    return list[num / 2];
}

bool check(int n, double *x) {
   double max = 0;
   for(int i = 0; i < n; i++) {
       if(std::fabs(x[i] - 1.0) > max)
           max = std::fabs(x[i] - 1.0);
   }
   return max <= 1e-9;
}

std::chrono::time_point<std::chrono::system_clock> start, end;
std::chrono::duration<double> elapsed_seconds;


void test_LL(const CSC *A, const double *b1, const double *b2, int inner_part, int level_param, int div_rate);


int main(int argc, char *argv[]) {
    std::string file = argv[1];

    FILE *fp = fopen(file.c_str(), "r");

    CSC *Amat = nullptr;
    process_matrix(fp, Amat);

    int numThread = (int) strtol(argv[2], nullptr, 10);
    int chunk =     (int) strtol(argv[3], nullptr, 10); // not used
    int innerPart = (int) strtol(argv[4], nullptr, 10);//Inner parts

    /** minimum # levels to consider **/
    int levelParam =  (int) strtol(argv[5], nullptr, 10);// level distance
    int blasThreads = (int) strtol(argv[6], nullptr, 10); // not used
    int divRate =     (int) strtol(argv[7], nullptr, 10);

    /**
     * METIS reordering of matrix
     */
    int *Perm = new int[Amat->ncol]();
    start = std::chrono::system_clock::now();
    Metis_oredering(Amat, Perm);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    std::cout << "METIS: " << elapsed_seconds.count() << "\n";
    start = std::chrono::system_clock::now();
    int status = 0;
    CSC *A1 = ptranspose(Amat, 2, Perm, nullptr, 0, status);
    CSC *A2 = ptranspose(A1, 2, nullptr, nullptr, 0, status);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    std::cout << "Transpose: " << elapsed_seconds.count() << "\n";
    delete[]Perm;
    allocateAC(Amat, 0, 0, 0, false);
    allocateAC(A1, 0, 0, 0, false);

    /**
     * Initialize rhs for testing
     */
    auto *b1 = new double[A2->ncol]();
    auto *b2 = new double[A2->ncol]();
    rhsInit(A2->ncol, A2->p, A2->i, A2->x, b1);
    spmv_csc(A2->ncol, A2->p, A2->i, A2->x, b1, b2);

    /** Test LLx = b **/
    omp_set_num_threads(numThread);
    test_LL(A2, b1, b2, innerPart, levelParam, divRate);

    allocateAC(A2, 0, 0, 0, false);
    delete[]b1;
    delete[]b2;
    return 0;
}


void test_LL(const CSC *A, const double *b1, const double *b2, int inner_part, int level_param, int div_rate) {
    int n = A->ncol;
    int nnz = A->nzmax;
    int *Ap = A->p;
    int *Ai = A->i;
    double *Ax = A->x;

    /** one kernel **/
    double duration;
    auto *x1 = new double[n]();
#ifdef CSC_SER
    std::cout << "SER1: ";
    for(int i = 0; i < NUM_TEST; i++) {
        std::memcpy(x1, b1, sizeof(double) * n);
        start = std::chrono::system_clock::now();
        lsolve(n, Ap, Ai, Ax, x1);
        end = std::chrono::system_clock::now();
        elapsed_seconds = (end - start);
        duration = elapsed_seconds.count();
        if(check(n, x1))
            std::cout << duration << ", ";
        else
            std::cerr << "##Serial, ";
    }
    std::cout << "\n";
#endif
#ifdef CSC_LVL
    std::cout << "LVL1: ";
    int *levelPtr, *levelSet;
    int levels = buildLevelSet_CSC(n, Ap, Ai, levelPtr, levelSet);
//    std::cout << "\n" << levels << "\n";
    for(int i = 0; i < NUM_TEST; i++) {
        std::memcpy(x1, b1, sizeof(double) * n);
        start = std::chrono::system_clock::now();
        lsolvePar(n, Ap, Ai, Ax, x1, levels, levelPtr, levelSet, 0);
        end = std::chrono::system_clock::now();
        elapsed_seconds = (end - start);
        duration = elapsed_seconds.count();
        if(check(n, x1))
            std::cout << duration << ", ";
        else
            std::cerr << "##Levelset, ";
    }
    std::cout << "\n";
    delete[]levelPtr;
    delete[]levelSet;
#endif
#ifdef CSC_LBC
    std::cout << "LBC1: ";
    int *HlevelPtr, *HlevelSet, *parPtr, *partition;
    int nLevels = 0, nPar = 0;
    auto *nodeCost = new double[n]();
    for(int i = 0; i < n; i++)
        nodeCost[i] = 1;
    getCoarseLevelSet_DAG_CSC03(n, Ap, Ai, nLevels, HlevelPtr, HlevelSet, nPar, parPtr, partition,
            inner_part, level_param, div_rate, nodeCost);
//    std::cout << "\n" << nLevels << "\n";
    delete[]nodeCost;
    for(int i = 0; i < NUM_TEST; i++) {
        std::memcpy(x1, b1, sizeof(double) * n);
        start = std::chrono::system_clock::now();
        lsolveParH2(n, Ap, Ai, Ax, x1, nLevels, HlevelPtr, HlevelSet, nPar, parPtr, partition, 0);
        end = std::chrono::system_clock::now();
        elapsed_seconds = (end - start);
        duration = elapsed_seconds.count();
        if(check(n, x1))
            std::cout << duration << ", ";
        else
            std::cerr << "##LBC_ver, ";
    }
    std::cout << "\n";
    delete[]HlevelPtr;
    delete[]partition;
    delete[]parPtr;
#endif
    delete[]x1;

    /** two kernels **/
    std::cout << "============== LLx = b =============\n";

     /**
     * Merge two L graphs here
     */
     int *nLp, *nLi;
     int **Lps = new int*[2];
     int **Lis = new int*[2];
     Lps[0] = Lps[1] = Ap;
     Lis[0] = Lis[1] = Ai;

     merge_graph(2, n, Lps, Lis, nLp, nLi);

    int i;
    auto *x2 = new double[2 * n]();
#ifdef CSC_SER
    std::cout << "SER2: ";
    for(i = 0; i < NUM_TEST; i++) {
        std::memcpy(x2, b2, sizeof(double) * n);
        start = std::chrono::system_clock::now();
        lsolve(n, Ap, Ai, Ax, x2);
        lsolve(n, Ap, Ai, Ax, x2);
        end = std::chrono::system_clock::now();
        elapsed_seconds = (end - start);
        duration = elapsed_seconds.count();
        if(check(n, x2))
            std::cout << duration << ", ";
        else
            std::cerr << "##Serial, ";
    }
    std::cout << "\n";
#endif
#ifdef CSC_LVL
    std::cout << "LVL2: ";
    int *nlevelPtr, *nlevelSet;
    int nlevels = buildLevelSet_CSC(2*n, nLp, nLi, nlevelPtr, nlevelSet);
//    std::cout << "\n" << nlevels << "\n";
    for(i = 0; i < NUM_TEST; i++) {
        std::memcpy(x2, b2, sizeof(double) * n);
        std::memset(x2+n, 0, sizeof(double) * n);
        start = std::chrono::system_clock::now();
        LLcross(n, Ap, Ai, Ax, x2, nlevelPtr, nlevelSet, nlevels);
        end = std::chrono::system_clock::now();
        elapsed_seconds = (end - start);
        duration = elapsed_seconds.count();
        if(check(n, x2+n))
            std::cout << duration << ", ";
        else
            std::cerr << "##Levelset, ";
    }
    std::cout << "\n";
    delete[]nlevelPtr;
    delete[]nlevelSet;
#endif
#ifdef CSC_LBC
    int *nHlevelPtr, *nHlevelSet, *nparPtr, *npartition;
    int nnLevels = 0, nnPar = 0;
    nodeCost = new double[2*n]();
    for(i = 0; i < 2*n; i++)
        nodeCost[i] = 1;
    getCoarseLevelSet_DAG_CSC03(2*n, nLp, nLi, nnLevels, nHlevelPtr, nHlevelSet, nnPar, nparPtr, npartition,
                                inner_part, level_param, div_rate, nodeCost);
    delete[]nodeCost;
    std::cout << "LBC2: ";
    // std::cout << "\n" << nnLevels << "\n";
    for(i = 0; i < NUM_TEST; i++) {
        std::memcpy(x2, b2, sizeof(double) * n);
        std::memset(x2+n, 0, sizeof(double) * n);
        start = std::chrono::system_clock::now();
        llcrossParH2(n, Ap, Ai, Ax, x2, nnLevels, nHlevelPtr, nHlevelSet, nnPar, nparPtr, npartition, 0);
        end = std::chrono::system_clock::now();
        elapsed_seconds = (end - start);
        duration = elapsed_seconds.count();
        if(check(n, x2+n))
            std::cout << duration << ", ";
        else
            std::cerr << "##LBC_ver, ";
    }
    std::cout << "\n";
    delete[]nHlevelPtr;
    delete[]npartition;
    delete[]nparPtr;
#endif
    delete[]x2;
    delete[]nLp;
    delete[]nLi;
}


//int main(int argc, char *argv[]) {
//
//    std::string f1 = argv[1];
//    int *colA;
//    int *rowA;
//    double *valL;
//    double *valA;
//    int maxSupWid, maxCol;
//    size_t n, nnzA, ncol;
//
//    std::chrono::time_point<std::chrono::system_clock> start, end;
//    std::chrono::duration<double> elapsed_seconds;
//    double durationSym = 0, duration3 = 0, duration2 = 0, duration1 = 0;
//    long totalIntraContribNNZs = 0, totalInterCotribNNZs = 0, numOfEdgeCuts = 0;
//    int numberOfIntraCore = 0, numberOfInterCore = 0;
//
//    if (!readMatrix(f1, n, nnzA, colA, rowA, valA))
//        return -1;
//
//    int numThread = (int) strtol(argv[2], nullptr, 10);
//    int chunk =     (int) strtol(argv[3], nullptr, 10);
//    int innerPart = (int) strtol(argv[4], nullptr, 10);//Inner parts
//
//    /** minimum # levels to consider **/
//    int levelParam =  (int) strtol(argv[5], nullptr, 10);// level distance
//    int blasThreads = (int) strtol(argv[6], nullptr, 10);
//    int divRate =     (int) strtol(argv[7], nullptr, 10);
//
////    cout << "Inputs:\n" << f1 << "," << numThread << "," << chunk << "," << innerPart << "," <<
////         levelParam << "," << blasThreads << "," << divRate << ";\n";
//    /*
//     * Calling Cholesky to generate blocked triangular matrix
//     */
//
//    omp_set_num_threads(numThread);
////    mkl_set_num_threads(numThread);
//
//    // MKL_Set_Num_Threads(1);
////    MKL_Domain_Set_Num_Threads(blasThreads, MKL_DOMAIN_BLAS);
//
//    double *timingChol = new double[4 + numThread]();//for time measurement
//    double orderingTime = 0;
//    int status = 0;
//    CSC *Amat = new CSC;
//    Amat->nzmax = nnzA;
//    Amat->ncol = Amat->nrow = n;
//    Amat->stype = -1;
//    Amat->xtype = CHOLMOD_REAL;
//    Amat->packed = TRUE;
//    Amat->p = colA;
//    Amat->i = rowA;
//    Amat->x = valA;
//    Amat->nz = NULL;
//    Amat->sorted = TRUE;
//    ncol = Amat->ncol;
//
//    //Ordering
//
//    start = std::chrono::system_clock::now();
//#ifdef GIVEN
//    //pastix_data_t **pastix_data;
//    L->ordering = CHOLMOD_METIS;
//    for (int l = 0; l < A->nrow; ++l) {
//     Lperm[l] = inPerm[l];
//    }
//
//#elif METIS
//    int *Perm = new int[n]();
//    start = std::chrono::system_clock::now();
//    Metis_oredering(Amat, Perm);
//    end = std::chrono::system_clock::now();
//    elapsed_seconds = end - start;
//    std::cout << "METIS: " << elapsed_seconds.count() << "\n";
//#else
//    double info[20]={0};
//    double Control[2];
//    Control [0] = 10; //TODO check later //AMD_Dense
//    Control [1] = TRUE; //AMD_AGGRESSIVE
//    L->ordering = CHOLMOD_AMD;
//    amd_order(ncol,A->p,A->i,Lperm,NULL,info);
//#endif
//
//#ifdef VERIFY
//    auto checkOrder = new bool[ncol]();
//    for (int i = 0; i < ncol; ++i) checkOrder[i] = false;
//    for (int i = 0; i < ncol; ++i) {
//     checkOrder[Perm[i]] = true;
//    }
//    for (int i = 0; i < ncol; ++i) {
//     assert(checkOrder[i] == true);
//    }
//    delete checkOrder;
//#endif
//
//    start = std::chrono::system_clock::now();
//    CSC *A1 = ptranspose(Amat, 2, Perm, NULL, 0, status);
//    CSC *A2 = ptranspose(A1, 2, NULL, NULL, 0, status);
//    end = std::chrono::system_clock::now();
//    elapsed_seconds = end - start;
//    std::cout << "Transpose: " << elapsed_seconds.count() << "\n";
//    delete[]Perm;
//#if 0
//    for (int i = 0; i < n; ++i) {
//     for (int j = A2->p[i]; j < A2->p[i+1]; ++j) {
//      std::cout<<A2->i[j]<<";";
//     }
//     std::cout<<"\n";
//    }
//#endif
//
//    allocateAC(Amat, 0, 0, 0, FALSE);
//    allocateAC(A1, 0, 0, 0, FALSE);
///*
// * ********************* Triangular Solve
// */
//
//
//    /**
//     * Merge two L graphs here
//     */
//     int *nLp = new int[2 * n + 1]();
//     int *nLi = new int[2 * nnzA + n]();
//     int new_i_counter = 0;
//     size_t i = 0;
////     nLp[0] = 0;
//     for(i = 0; i < n; i++) {
//         nLp[i] = int(A2->p[i] + (i));
//         for(int j = A2->p[i]; j < A2->p[i+1]; j++) {
//             nLi[new_i_counter] = int(A2->i[j]);
//             new_i_counter ++;
//         }
//         nLi[new_i_counter] = int(i + n);
//         new_i_counter++;
//     }
//     for(; i < 2 * n; i++) {
//         nLp[i] = int(A2->p[i-n] + nnzA + n);
//         for(int j = A2->p[i-n]; j < A2->p[i-n+1]; j++) {
//             nLi[new_i_counter] = int(A2->i[j] + n);
//             new_i_counter++;
//         }
//     }
//     nLp[2*n] = int(2 * nnzA + n);
//     assert(new_i_counter == 2 * nnzA + n);
//
//
//#ifdef FLOPCNT
//    //***************Serial
//    int *ia = new int[n + 1];
//    int *ja = new int[L->xsize];
//    double *a = new double[L->xsize];
//    bcsc2csc(n, L->nsuper, L->p, L->s, L->i_ptr, L->super, valL, ia, ja, a);
//    unsigned long counts=0;
//    rhsInit(n,ia,ja,a,x);
//    counts = flopCoutLSolve(n,ia,ja,a,x);
//    std::cout<<L->xsize<<";"<<counts<<";";
//    delete []ia;
//    delete []ja;
//    delete []a;
//#endif
////Running LBC here
//
//    int *HLevelPtr = NULL, *HLevelSet = NULL, *parPtr = NULL,
//            *partition = NULL;
//    int *levelPtr = NULL, *levelSet = NULL;
//    int nLevels = 0, nPar = 0, levels = 0;
//    double *nodeCost;
//    int iterNo = 5;
//
//    //Computing node cost
//    nodeCost = new double[2*n];
//    for (int s = 0; s < n; ++s)
//        nodeCost[s] = 1;
//    start = std::chrono::system_clock::now();
//    int avg = getCoarseLevelSet_DAG_CSC03(2*n, nLp, nLi,
//                                          nLevels, HLevelPtr,
//                                          HLevelSet, nPar,
//                                          parPtr, partition,
//                                          innerPart, levelParam, divRate, nodeCost);
//    end = std::chrono::system_clock::now();
//    elapsed_seconds = end - start;
//    duration1 = elapsed_seconds.count();
//    std::cout << "LBC: " << avg << "," << duration1 << ",\n";
//    delete[]nodeCost;
//    delete[]levelPtr;
//    delete[]levelSet;
//
//    auto *ones = new double[n]();
//    for(int i = 0; i < n; i++)
//        ones[i] = 1.0;
//    auto *x = new double[n]();
//#ifdef CSC_SER
//    //***************Serial
//    std::cout << "Serial: ";
//    for (int l = 0; l < iterNo; ++l) {
//        std::memset(x, 0, sizeof(double) * n);
//        auto *temp = new double[n]();
//        spmv_csc(n, A2->p, A2->i, A2->x, ones, temp);
//        spmv_csc(n, A2->p, A2->i, A2->x, temp, x);
//        delete[]temp;
//        start = std::chrono::system_clock::now();
//        lsolve(n, A2->p, A2->i, A2->x, x);
//        lsolve(n, A2->p, A2->i, A2->x, x);
//        end = std::chrono::system_clock::now();
//        elapsed_seconds = (end - start);
//        duration1 = elapsed_seconds.count();
//
//        double max = 0;
//        for(int i = 0; i < n; i++) {
//            if(std::fabs(x[i] - 1.0) > max)
//                max = std::fabs(x[i] - 1.0);
//        }
////        if (!testTriangular(n, x))
//        if(max > 1e-8)
//            std::cout << "##serial,";
//        else
//            std::cout << duration1 << ",";
//    }
//    delete[]x;
//    cout << "\n";
//#endif
//#ifdef CSC_LVL
//    x = new double[2 * n]();
//    //****************Parallel CSC
//    levels = buildLevelSet_CSC(n, A2->p, A2->i,
//                               levelPtr, levelSet);
//    int *nlevelPtr, *nlevelSet;
//    int nlevels = buildLevelSet_CSC(2*n, nLp, nLi, nlevelPtr, nlevelSet);
//
//    std::cout << "Levelset: ";
//    for (int l = 0; l < iterNo; ++l) {
//        std::memset(x, 0, sizeof(double) * 2 * n);
//        auto *temp = new double[n]();
//        spmv_csc(n, A2->p, A2->i, A2->x, ones, temp);
//        spmv_csc(n, A2->p, A2->i, A2->x, temp, x);
//        delete[]temp;
//        start = std::chrono::system_clock::now();
////        lsolvePar(n, A2->p, A2->i, A2->x, x, levels, levelPtr, levelSet, chunk);
//        LLcross(n, A2->p, A2->i, A2->x, x, nlevelPtr, nlevelSet, nlevels);
//        end = std::chrono::system_clock::now();
//        elapsed_seconds = end - start;
//        duration2 = elapsed_seconds.count();
//
//        double max = 0;
//        for(int i = 0; i < n; i++) {
//            if(std::fabs(x[i+n] - 1.0) > max)
//                max = std::fabs(x[i+n] - 1.0);
//        }
////        if (!testTriangular(n, x))
//        if(max > 1e-8)
//            std::cout << "##LevelSet,";
//        else
//            std::cout << duration2 << ",";
//    }
//    cout << "\n";
//    delete[]levelPtr;
//    delete[]levelSet;
//    delete[]x;
//#endif
//#ifdef CSC_LBC
//    x = new double[2 * n]();
//    //****************Parallel H2 CSC
//    std::cout << "LBC time: ";
//    for (int l = 0; l < iterNo; ++l) {
//        std::memset(x, 0, sizeof(double) * 2 * n);
//        auto *temp = new double[n]();
//        spmv_csc(n, A2->p, A2->i, A2->x, ones, temp);
//        spmv_csc(n, A2->p, A2->i, A2->x, temp, x);
//        delete[]temp;
//        start = std::chrono::system_clock::now();
//        llcrossParH2(n, A2->p, A2->i, A2->x, x, nLevels, HLevelPtr, HLevelSet,
//                    nPar, parPtr, partition, chunk);
//        end = std::chrono::system_clock::now();
//        elapsed_seconds = (end - start);
//        duration2 = elapsed_seconds.count();
//
//        double max = 0;
//        for(int i = 0; i < n; i++) {
//            if(std::fabs(x[i+n] - 1.0) > max)
//                max = std::fabs(x[i+n] - 1.0);
//        }
////        if (!testTriangular(n, x))
//        if(max > 1e-8)
//            std::cout << "##HlevelSet,";
//        else
//            std::cout << duration2 << ",";
//    }
//    cout << "\n";
//    delete[]x;
//#endif
//
//#ifdef CSC_LBC_SER
//    //****************Parallel H2 CSC
//    std::cout << "Serial LBC: ";
//    omp_set_num_threads(1);
//    for (int l = 0; l < iterNo; ++l) {
//        rhsInit(n, A2->p, A2->i, A2->x, x);
//        start = std::chrono::system_clock::now();
//        lsolveParH2(n, A2->p, A2->i, A2->x, x, nLevels, HLevelPtr, HLevelSet,
//                    nPar, parPtr, partition, chunk);
//        //lsolvePar2(n,col,row,val,x);
//        end = std::chrono::system_clock::now();
//        elapsed_seconds = end - start;
//        duration2 = elapsed_seconds.count();
//        if (!testTriangular(n, x))
//            std::cout << "##HlevelSetSerial,";
//        else
//            std::cout << duration2 << ",";
//    }
//    cout << nLevels << "\n";
//#endif
//
//
//#ifdef BCSC_TRNG
//    //*************** BCSC
//    int *levelbPtr, *levelbSet, blevels=0;
//    int *col2sup = L->col2Sup;
//    int supNo=0, newNNZ=0;
//    int *sup2col = L->super;
//    size_t *newCol = L->p;
//    int *newRow = L->s;
//    double *newVal = valL;
//    size_t *rowP = L->i_ptr;
//    size_t nBlocks = L->nsuper;
//    double max_tmp=0.0;
//    int nChild=NULL;
//    int ATreeHeight = getTreeHeightBruteForce(nBlocks,L->sParent);
//
//
//   /* for (int k = 0; k < nBlocks; ++k) {
//     cout<<L->sParent[k]<<",";
//    }
//    cout<<"\n";*/
//   /* int *cutPerLevel = getInitialCut_DAG(newCol,rowP,newRow,nBlocks,sup2col,
//                                         col2sup);*/
//
//   /* std::vector<int> empty;
//    int *xi = new int[2*nBlocks]();
//    int *marked = new int[nBlocks]();
//    marked[25]=false;
//    int tmp = dfs_BCSC_CC(nBlocks, 1, newCol, rowP, newRow, sup2col,
//                       col2sup,marked,
//                       nBlocks,xi,xi+nBlocks,empty,NULL);
//    for(int i=tmp; i<nBlocks; ++i){
//     std::cout<<xi[i]<<";";
//    }
//    std::cout<<"\n";*/
//
//    std::cout<<f1<<","<<levelParam<<","<<divRate<<","<<n<<","<< ATreeHeight <<","
//             <<nBlocks<<","<<nnz <<","<<durationSym<<","<<orderingTime<<",,";
//    int iterno=5;
//    //*************** Serial Blocked
//    for (int j = 0; j < iterno; ++j) {
//     rhsInitBlocked(n, nBlocks, newCol, newRow, rowP, newVal, x);
//     start = std::chrono::system_clock::now();
//     blockedLsolve(n, newCol, newRow, newVal, nnz, rowP, col2sup, sup2col, nBlocks, x);
//     end = std::chrono::system_clock::now();
//     elapsed_seconds = end - start;
//     duration2 = elapsed_seconds.count();
//     std::cout << duration2 << ",";
//    }
//    std::cout<<"*:";
//
//#if 0
//    for (int j = 0; j < blevels; ++j) {
//           for (int i = levelbPtr[j]; i < levelbPtr[j+1]; ++i) {
//               std::cout<<levelbSet[i]<<",";
//           }
//           std::cout<<"\n";
//       }
//#endif
//    //*************** Parallel Blocked
//    //omp_set_num_threads(1);
//    /*blevels = buildLevelSet_BCSC(n,nnz,col,rowP,newRow,nBlocks,
//
//                                sup2col,col2sup,levelbPtr,levelbSet);*/
//    levelbPtr = new int[nBlocks]();
//    levelbSet = new int[nBlocks]();
//    blevels = getLevelSet(nBlocks,L->sParent,levelbPtr,levelbSet);
//
//#if 0
//    for (int i = 0; i < blevels; ++i) {
//     int cnode = levelbPtr[i+1] - levelbPtr[i];
//     //assert(cnode == cutPerLevel[i]);
//     std::cout<<cnode<<";";
//    }
//    std::cout<<"\n";
//#endif
//
//    for (int j = 0; j < iterno; ++j) {
//     rhsInitBlocked(n,nBlocks,newCol,newRow,rowP,newVal,x);
//     start = std::chrono::system_clock::now();
//     //blockedLsolve(n,newCol,newRow,newVal,nnz,rowP,col2sup,sup2col,supNo,x);
//     leveledBlockedLsolve(n,newCol,newRow,newVal,nnz,rowP,col2sup,sup2col,nBlocks,
//                          x,blevels,levelbPtr,levelbSet,chunk);
//     end = std::chrono::system_clock::now();
//     elapsed_seconds = end-start;
//     duration2=elapsed_seconds.count();
//     if(testTriangular(n, x))
//      std::cout<<duration2<<",";
//     else
//      std::cout<<"H1 failed,";
//    }
//
//    std::cout<<"*:";
//    //*************** Parallel Blocked H2
//    //omp_set_num_threads(1);
//    for (int j = 0; j < iterno; ++j) {
//     rhsInitBlocked(n, nBlocks, newCol, newRow, rowP, newVal, x);
//     start = std::chrono::system_clock::now();
//     //blockedLsolve(n,newCol,newRow,newVal,nnz,rowP,col2sup,sup2col,supNo,x);
//     H2LeveledBlockedLsolve(n, newCol, newRow, newVal, nnz, rowP, col2sup, sup2col, nBlocks,
//                            x, nLevels, HLevelPtr, HLevelSet,
//                            nPar, parPtr, partition, chunk);
//     end = std::chrono::system_clock::now();
//     elapsed_seconds = end - start;
//     duration2 = elapsed_seconds.count();
//     if (testTriangular(n, x))
//      std::cout << duration2 << ",";
//     else
//      std::cout << "H2 failed,";
//    }
//    std::cout<<"*:";
//
//    //*************** Parallel Blocked H2 Peeled
//    //omp_set_num_threads(1);
//   /* for (int j = 0; j < iterno; ++j) {
//     rhsInitBlocked(n, nBlocks, newCol, newRow, rowP, newVal, x);
//     start = std::chrono::system_clock::now();
//     //blockedLsolve(n,newCol,newRow,newVal,nnz,rowP,col2sup,sup2col,supNo,x);
//     H2LeveledBlockedLsolve_Peeled(n, newCol, newRow, newVal, nnz, rowP, col2sup, sup2col, nBlocks,
//                                   x, nLevels, HLevelPtr, HLevelSet,
//                                   nPar, parPtr, partition, chunk,numThread);
//     end = std::chrono::system_clock::now();
//     elapsed_seconds = end - start;
//     duration2 = elapsed_seconds.count();
//     if (testTriangular(n, x))
//      std::cout << duration2 << ",";
//     else
//      std::cout << "H2 failed,";
//    }
//    std::cout<<"*:";*/
//   //TODO HACK
//    delete []L->super;
//    delete []L->sParent;
//    delete []L->s;
//    delete []L->col2Sup;
//    delete []L->p;
//    delete []L->i_ptr;
//#endif
//
//
//#if DEBUG > 0
//    for (int i = n-10; i < n; ++i) {
//               std::cout<<i<<":\n";
//               for (int m = colL[i],cnt=0; m < colL[i+1]; ++m, ++cnt) {
//                   if(!std::isfinite(valL[m])) {
//                       std::cout << "Error in colA "<< i;
//                       return -1;
//                   }
//                   if(rowL[li_ptr[i]+cnt] >= i )
//                       std::cout<<valL[m]<<",";
//               }
//               std::cout<<"\n";
//           }
//           std::cout<<"\n";
//#endif
//// delete []col2sup;
//#ifdef PRUNE
//    delete []prunePtr; delete []pruneSet;
//#endif
//    if (HLevelPtr != NULL)
//        delete[]HLevelPtr;
//    if (HLevelPtr != NULL)
//        delete[]HLevelSet;
//    if (parPtr != NULL)
//        delete[]parPtr;
//    if (partition != NULL)
//        delete[]partition;
//    //delete []contribs;
//    //delete []map;
//// delete[]valL;
//    //delete []colL;
//    //delete []li_ptr;
//    delete[]timingChol;
//    allocateAC(A2, 0, 0, 0, FALSE);
//}
//
//
///*
// * Testing lower triangular solve to make sure all x elements are ONE.
// * Warning: This works only for when RHS does not have zero
// */
//
//int testTriangular(size_t n, const double *x, double epsilon) {//Testing
//    int test = 0;
//    for (int i = 0; i < n; ++i) {
//        if (std::abs(1 - x[i]) < epsilon) {
//            test++;
//        } /*else{
//   std::cout<<i<<" : "<<1-x[i]<<";";
//  }*/
//        //else
//        // cout<<i<<";";
//    }
//    if (n - test > 0) {
//        return false;
//    }
//    return true;
//}
