//
// Created by kazem on 9/18/19.
//


#include <iostream>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <cholmod.h>
#include <cholmod_function.h>
#include <omp.h>
#include "Transpose.h"

#include "mkl.h"
#include "metis_wrapper.h"

#include "Triangular_CSC.h"
#include "Inspection_Level.h"
#include "lbc_csc.h"
#include "Util.h"



void rhsInit(int n, int *Lp, int *Li, double *Lx, double *b) {
 for(int i=0; i<n; i++) b[i]=0;
 for (int c = 0; c < n; ++c) {
  for (int cc = Lp[c]; cc < Lp[c + 1]; ++cc) {
   b[Li[cc]] += Lx[cc];
  }
 }
}
#define CPUTIME (SuiteSparse_time ( ))
#define CSC_TRNG
#undef DEBUG
//#define FLOPCNT
#define METIS 1

using namespace std;

int testTriangular(size_t n, const double *x, double epsilon=1e-9);

int main(int argc, char *argv[]) {

 std::string f1 = argv[1];
 int *colA;
 int *rowA;
 double *valL;
 double *valA;
 int maxSupWid, maxCol;
 size_t n, nnzA, ncol;

 std::chrono::time_point<std::chrono::system_clock> start, end;
 std::chrono::duration<double> elapsed_seconds;
 double durationSym = 0, duration3 = 0, duration2 = 0, duration1 = 0;
 long totalIntraContribNNZs = 0, totalInterCotribNNZs = 0, numOfEdgeCuts = 0;
 int numberOfIntraCore = 0, numberOfInterCore = 0;

 if (!readMatrix(f1, n, nnzA, colA, rowA, valA))
  return -1;

 int numThread = atoi(argv[2]);
 int chunk = atoi(argv[3]);
 int innerPart = atoi(argv[4]);//Inner parts
 int levelParam = atoi(argv[5]);// level distance
 int blasThreads = atoi(argv[6]);
 int divRate = atoi(argv[7]);

 cout<<f1<<","<<numThread<<","<<chunk<<","<<innerPart<<","<<
     levelParam<<","<<blasThreads<<","<<divRate<<";";
 /*
  * Calling Cholesky to generate blocked triangular matrix
  */

 omp_set_num_threads(numThread);
 mkl_set_num_threads(numThread);

 // MKL_Set_Num_Threads(1);
 MKL_Domain_Set_Num_Threads(blasThreads, MKL_DOMAIN_BLAS);


 int *Perm = new int[n]();

 double *timingChol = new double[4 + numThread]();//for time measurement
 double orderingTime = 0;
 int status = 0;
 CSC *Amat = new CSC;
 Amat->nzmax = nnzA;
 Amat->ncol = Amat->nrow = n;
 Amat->stype = -1;
 Amat->xtype = CHOLMOD_REAL;
 Amat->packed = TRUE;
 Amat->p = colA;
 Amat->i = rowA;
 Amat->x = valA;
 Amat->nz = NULL;
 Amat->sorted = TRUE;
 ncol = Amat->ncol;

 //Ordering

 start = std::chrono::system_clock::now();
#ifdef GIVEN
 //pastix_data_t **pastix_data;
 L->ordering = CHOLMOD_METIS;
 for (int l = 0; l < A->nrow; ++l) {
  Lperm[l] = inPerm[l];
 }

#elif METIS
 Metis_oredering(Amat,Perm);
#else
 double info[20]={0};
 double Control[2];
 Control [0] = 10; //TODO check later //AMD_Dense
 Control [1] = TRUE; //AMD_AGGRESSIVE
 L->ordering = CHOLMOD_AMD;
 amd_order(ncol,A->p,A->i,Lperm,NULL,info);
#endif
 end = std::chrono::system_clock::now();
 elapsed_seconds = end - start;
 orderingTime = elapsed_seconds.count();
//printf("ddddd %f ddd \n",orderingTime);
 elapsed_seconds = end - start;
 durationSym = elapsed_seconds.count();
#ifdef VERIFY
 auto checkOrder = new bool[ncol]();
 for (int i = 0; i < ncol; ++i) checkOrder[i] = false;
 for (int i = 0; i < ncol; ++i) {
  checkOrder[Perm[i]] = true;
 }
 for (int i = 0; i < ncol; ++i) {
  assert(checkOrder[i] == true);
 }
 delete checkOrder;
#endif


 CSC *A1 = ptranspose(Amat, 2, Perm, NULL, 0, status);
 CSC *A2 = ptranspose(A1, 2, NULL, NULL, 0, status);
#if 0
 for (int i = 0; i < n; ++i) {
  for (int j = A2->p[i]; j < A2->p[i+1]; ++j) {
   std::cout<<A2->i[j]<<";";
  }
  std::cout<<"\n";
 }
#endif

 allocateAC(Amat, 0, 0, 0, FALSE);
 allocateAC(A1, 0, 0, 0, FALSE);
/*
 * ********************* Triangular Solve
 */


/* std::chrono::time_point<std::chrono::system_clock> start, end;
 std::chrono::duration<double> elapsed_seconds;
 double duration4 = 0 ,duration3 = 0, duration2=0, duration1=0;*/

 double *x = new double[n]();

#ifdef FLOPCNT
 //***************Serial
 int *ia = new int[n + 1];
 int *ja = new int[L->xsize];
 double *a = new double[L->xsize];
 bcsc2csc(n, L->nsuper, L->p, L->s, L->i_ptr, L->super, valL, ia, ja, a);
 unsigned long counts=0;
 rhsInit(n,ia,ja,a,x);
 counts = flopCoutLSolve(n,ia,ja,a,x);
 std::cout<<L->xsize<<";"<<counts<<";";
 delete []ia;
 delete []ja;
 delete []a;
#endif
//Running LBC here

 int *HLevelPtr = NULL, *HLevelSet = NULL, *parPtr = NULL,
   *partition =NULL;
 int *levelPtr = NULL, *levelSet = NULL;
 int nLevels=0, nPar=0, levels=0;
 double *nodeCost;
 int iterNo=5;

 //Computing node cost
 nodeCost = new double[n];
 int *xi = new int[2*n]();
 for (int s = 0; s < n; ++s) {
  //nodeCost[s]=A2->p[s+1]-A2->p[s];
  nodeCost[s]=1;
 }
 delete []xi;
 start = std::chrono::system_clock::now();
 int avg = getCoarseLevelSet_DAG_CSC03(n, A2->p, A2->i,
                                       nLevels, HLevelPtr,
                                       HLevelSet, nPar,
                                       parPtr, partition,
                                       innerPart, levelParam, divRate, nodeCost);
 end = std::chrono::system_clock::now();
 elapsed_seconds = end - start;
 duration1 = elapsed_seconds.count();
 std::cout <<avg<<"," <<duration1 << ",";
 delete[]nodeCost;
#ifdef CSC_TRNG
 //***************Serial
 for (int l = 0; l < iterNo; ++l) {
  rhsInit(n, A2->p, A2->i, A2->x, x);
  start = std::chrono::system_clock::now();
  lsolve(n, A2->p, A2->i, A2->x, x);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  duration1 = elapsed_seconds.count();
  if (!testTriangular(n, x))
   std::cout << "##serial,";
  else
   std::cout << duration1 << ",";
 }
 cout<< ";;";

 //****************Parallel CSC
 for (int l = 0; l < iterNo; ++l) {
  levels = buildLevelSet_CSC(n, A2->nzmax, A2->p, A2->i,
                             levelPtr, levelSet);
  rhsInit(n, A2->p, A2->i, A2->x, x);
  start = std::chrono::system_clock::now();
  lsolvePar(n, A2->p, A2->i, A2->x, x, levels, levelPtr, levelSet, chunk);
  //lsolvePar2(n,col,row,val,x);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  duration2 = elapsed_seconds.count();
  if (!testTriangular(n, x) )
   std::cout << "##LevelSet,";
  else
   std::cout << duration2 << ",";
 }
 cout<< levels << ";;";

 //****************Parallel H2 CSC
 for (int l = 0; l < iterNo; ++l) {
  rhsInit(n,A2->p,A2->i,A2->x,x);
  start = std::chrono::system_clock::now();
  lsolveParH2(n,A2->p,A2->i,A2->x,x,nLevels,HLevelPtr,HLevelSet,
              nPar,parPtr,partition, chunk);
  //lsolvePar2(n,col,row,val,x);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  duration2=elapsed_seconds.count();
  if(!testTriangular(n,x))
   std::cout<<"##HlevelSet,";
  else
   std::cout <<duration2<<",";
 }
 cout<< nLevels <<";;";

 //****************Parallel H2 CSC
 //omp_set_num_threads(1);
 for (int l = 0; l < iterNo; ++l) {
  rhsInit(n,A2->p,A2->i,A2->x,x);
  start = std::chrono::system_clock::now();
  lsolveParH2(n,A2->p,A2->i,A2->x,x,nLevels,HLevelPtr,HLevelSet,
              nPar,parPtr,partition, chunk);
  //lsolvePar2(n,col,row,val,x);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  duration2=elapsed_seconds.count();
  if(!testTriangular(n,x))
   std::cout<<"##HlevelSetSerial,";
  else
   std::cout<<duration2<<",";
 }
 cout<< nLevels <<";;";
#endif


#ifdef BCSC_TRNG
 //*************** BCSC
 int *levelbPtr, *levelbSet, blevels=0;
 int *col2sup = L->col2Sup;
 int supNo=0, newNNZ=0;
 int *sup2col = L->super;
 size_t *newCol = L->p;
 int *newRow = L->s;
 double *newVal = valL;
 size_t *rowP = L->i_ptr;
 size_t nBlocks = L->nsuper;
 double max_tmp=0.0;
 int nChild=NULL;
 int ATreeHeight = getTreeHeightBruteForce(nBlocks,L->sParent);


/* for (int k = 0; k < nBlocks; ++k) {
  cout<<L->sParent[k]<<",";
 }
 cout<<"\n";*/
/* int *cutPerLevel = getInitialCut_DAG(newCol,rowP,newRow,nBlocks,sup2col,
                                      col2sup);*/

/* std::vector<int> empty;
 int *xi = new int[2*nBlocks]();
 int *marked = new int[nBlocks]();
 marked[25]=false;
 int tmp = dfs_BCSC_CC(nBlocks, 1, newCol, rowP, newRow, sup2col,
                    col2sup,marked,
                    nBlocks,xi,xi+nBlocks,empty,NULL);
 for(int i=tmp; i<nBlocks; ++i){
  std::cout<<xi[i]<<";";
 }
 std::cout<<"\n";*/

 std::cout<<f1<<","<<levelParam<<","<<divRate<<","<<n<<","<< ATreeHeight <<","
          <<nBlocks<<","<<nnz <<","<<durationSym<<","<<orderingTime<<",,";
 int iterno=5;
 //*************** Serial Blocked
 for (int j = 0; j < iterno; ++j) {
  rhsInitBlocked(n, nBlocks, newCol, newRow, rowP, newVal, x);
  start = std::chrono::system_clock::now();
  blockedLsolve(n, newCol, newRow, newVal, nnz, rowP, col2sup, sup2col, nBlocks, x);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  duration2 = elapsed_seconds.count();
  std::cout << duration2 << ",";
 }
 std::cout<<"*:";

#if 0
 for (int j = 0; j < blevels; ++j) {
        for (int i = levelbPtr[j]; i < levelbPtr[j+1]; ++i) {
            std::cout<<levelbSet[i]<<",";
        }
        std::cout<<"\n";
    }
#endif
 //*************** Parallel Blocked
 //omp_set_num_threads(1);
 /*blevels = buildLevelSet_BCSC(n,nnz,col,rowP,newRow,nBlocks,

                             sup2col,col2sup,levelbPtr,levelbSet);*/
 levelbPtr = new int[nBlocks]();
 levelbSet = new int[nBlocks]();
 blevels = getLevelSet(nBlocks,L->sParent,levelbPtr,levelbSet);

#if 0
 for (int i = 0; i < blevels; ++i) {
  int cnode = levelbPtr[i+1] - levelbPtr[i];
  //assert(cnode == cutPerLevel[i]);
  std::cout<<cnode<<";";
 }
 std::cout<<"\n";
#endif

 for (int j = 0; j < iterno; ++j) {
  rhsInitBlocked(n,nBlocks,newCol,newRow,rowP,newVal,x);
  start = std::chrono::system_clock::now();
  //blockedLsolve(n,newCol,newRow,newVal,nnz,rowP,col2sup,sup2col,supNo,x);
  leveledBlockedLsolve(n,newCol,newRow,newVal,nnz,rowP,col2sup,sup2col,nBlocks,
                       x,blevels,levelbPtr,levelbSet,chunk);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  duration2=elapsed_seconds.count();
  if(testTriangular(n, x))
   std::cout<<duration2<<",";
  else
   std::cout<<"H1 failed,";
 }

 std::cout<<"*:";
 //*************** Parallel Blocked H2
 //omp_set_num_threads(1);
 for (int j = 0; j < iterno; ++j) {
  rhsInitBlocked(n, nBlocks, newCol, newRow, rowP, newVal, x);
  start = std::chrono::system_clock::now();
  //blockedLsolve(n,newCol,newRow,newVal,nnz,rowP,col2sup,sup2col,supNo,x);
  H2LeveledBlockedLsolve(n, newCol, newRow, newVal, nnz, rowP, col2sup, sup2col, nBlocks,
                         x, nLevels, HLevelPtr, HLevelSet,
                         nPar, parPtr, partition, chunk);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  duration2 = elapsed_seconds.count();
  if (testTriangular(n, x))
   std::cout << duration2 << ",";
  else
   std::cout << "H2 failed,";
 }
 std::cout<<"*:";

 //*************** Parallel Blocked H2 Peeled
 //omp_set_num_threads(1);
/* for (int j = 0; j < iterno; ++j) {
  rhsInitBlocked(n, nBlocks, newCol, newRow, rowP, newVal, x);
  start = std::chrono::system_clock::now();
  //blockedLsolve(n,newCol,newRow,newVal,nnz,rowP,col2sup,sup2col,supNo,x);
  H2LeveledBlockedLsolve_Peeled(n, newCol, newRow, newVal, nnz, rowP, col2sup, sup2col, nBlocks,
                                x, nLevels, HLevelPtr, HLevelSet,
                                nPar, parPtr, partition, chunk,numThread);
  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  duration2 = elapsed_seconds.count();
  if (testTriangular(n, x))
   std::cout << duration2 << ",";
  else
   std::cout << "H2 failed,";
 }
 std::cout<<"*:";*/
//TODO HACK
 delete []L->super;
 delete []L->sParent;
 delete []L->s;
 delete []L->col2Sup;
 delete []L->p;
 delete []L->i_ptr;
#endif


#if DEBUG > 0
 for (int i = n-10; i < n; ++i) {
            std::cout<<i<<":\n";
            for (int m = colL[i],cnt=0; m < colL[i+1]; ++m, ++cnt) {
                if(!std::isfinite(valL[m])) {
                    std::cout << "Error in colA "<< i;
                    return -1;
                }
                if(rowL[li_ptr[i]+cnt] >= i )
                    std::cout<<valL[m]<<",";
            }
            std::cout<<"\n";
        }
        std::cout<<"\n";
#endif
// delete []col2sup;
#ifdef PRUNE
 delete []prunePtr; delete []pruneSet;
#endif
 if (HLevelPtr != NULL)
  delete[]HLevelPtr;
 if (HLevelPtr != NULL)
  delete[]HLevelSet;
 if (parPtr != NULL)
  delete[]parPtr;
 if (partition != NULL)
  delete[]partition;
 //delete []contribs;
 //delete []map;
// delete[]valL;
 //delete []colL;
 //delete []li_ptr;
 delete[]timingChol;
 allocateAC(A2, 0, 0, 0, FALSE);
}


/*
 * Testing lower triangular solve to make sure all x elements are ONE.
 * Warning: This works only for when RHS does not have zero
 */

int testTriangular(size_t n, const double *x, double epsilon) {//Testing
 int test=0;
 for (int i = 0; i < n; ++i) {
  if(std::abs(1-x[i])<epsilon){
   test++;
  } /*else{
   std::cout<<i<<" : "<<1-x[i]<<";";
  }*/
  //else
  // cout<<i<<";";
 }
 if(n-test>0){
  return false;
 }
 return true;
}
