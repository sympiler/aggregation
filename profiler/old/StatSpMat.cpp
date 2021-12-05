//
// Created by labuser (Bangtian Liu) on 9/15/20.
//

#include <StatSpMat.h>
#include <lbc.h>
#include <sparse_inspector.h>

#include <sparse_blas_lib.h>

namespace sym_lib{
    StatSpMat::

    StatSpMat(CSR *L, SpKerType kerType, int num_threads, int blksize)
    {
     omp_set_num_threads(num_threads);

     Setup(L, kerType);

     //profile the serial code
     fs_csr_stat(L->n, L->p, L->i, this->nFlops, this->nnz_access, this->nnz_reuse);

     std::vector<std::vector<int>> DAG;
     DAG.resize(L->n);

     fs_csr_inspector_dep(L->n, L->p, L->i, DAG);


     int *groupPtr = (int *)malloc(sizeof(int)*(L->n+1));
     memset(groupPtr, 0, sizeof(int)*(1+L->n));
     int *groupSet = (int *)malloc(sizeof(int)*L->n);
     memset(groupSet, 0, sizeof(int)*L->n);
     int *groupInv = (int *)malloc(sizeof(int)*L->n);
     memset(groupInv, 0, sizeof(int)*L->n);

     int ngroup;
     group g(L->n, L->p, L->i);
     g.inspection_sptrsvcsr_v1(groupPtr, groupSet, ngroup, groupInv);
//        NaiveGrouping(L->n,  groupPtr, groupSet, ngroup, groupInv, blksize);
     this->ngroup = ngroup;
     this->npart = ngroup;

     /**
      * apply grouping information to the DAG and generated a smaller DAG.
      */
     auto gDAG=Group_DAG(DAG, groupPtr, groupSet, groupInv, ngroup);


     size_t count=0;
     for (int j = 0; j < gDAG.size(); ++j) {
//            DAG[j].erase(std::unique(DAG[j].begin(), DAG[j].end()), DAG[j].end());
      count+=gDAG[j].size();
     }

//        detectDAGCircle(DAG);

     int *gv, *gedg;
     gv = new int[L->n+1]();
     gedg = new int[count+L->n]();
     int *levelPtr = new int[L->n+1]();
     int *levelSet = new int[L->n]();

     long int cti,edges=0;
     for(cti = 0, edges = 0; cti < ngroup; cti++){
      gv[cti] = edges;
      gedg[edges++] = cti;
      for (int ctj = 0; ctj < gDAG[cti].size(); ctj++) {
       gedg[edges++] = gDAG[cti][ctj];
      }
     }
     gv[cti] = edges;


     this->nlevels = buildLevelSet_CSC_Queue(ngroup, 0, gv, gedg, levelPtr, levelSet);

     this->num_sys = this->nlevels;


     this->nnzPerLevels = this->nnz * 1.0 /this->nlevels;
     this->averParallelism = this->ngroup * 1.0 / this->nlevels;

     std::vector<int> lcost;
     lcost.resize(this->nlevels);

     fs_csr_levelset_stat(L->p, L->i, groupPtr, groupSet, this->nlevels, levelPtr, levelSet, lcost.data());

     this->SumMaxDiff=0;
     for (auto &cost: lcost) {
      this->SumMaxDiff +=cost;
     }

     this->AverageMaxDiff = this->SumMaxDiff * 1.0 / lcost.size();

     double accum=0.0;
     std::for_each (std::begin(lcost), std::end(lcost), [&](const double d) {
         accum += (d - this->AverageMaxDiff) * (d - this->AverageMaxDiff);
     });
     this->VarianceMaxDiff = std::sqrt(accum/(lcost.size()-1));

     free(groupPtr);
     free(groupSet);
     free(groupInv);

     delete [] gv;
     delete [] gedg;
     delete [] levelPtr;
     delete [] levelSet;
    }


    StatSpMat::StatSpMat(CSR *L, std::vector<std::vector<int> > DAG, SpKerType kerType, int num_threads, int blksize) {
     omp_set_num_threads(num_threads);

     Setup(L, kerType);

     //profile the serial code
     fs_csr_stat(L->n, L->p, L->i, this->nFlops, this->nnz_access, this->nnz_reuse);


     int *groupPtr = (int *) malloc(sizeof(int) * (L->n + 1));
     memset(groupPtr, 0, sizeof(int) * (1 + L->n));
     int *groupSet = (int *) malloc(sizeof(int) * L->n);
     memset(groupSet, 0, sizeof(int) * L->n);
     int *groupInv = (int *) malloc(sizeof(int) * L->n);
     memset(groupInv, 0, sizeof(int) * L->n);

     int ngroup;
     group g(L->n, L->p, L->i);
     g.inspection_sptrsvcsr_v1(groupPtr, groupSet, ngroup, groupInv);
//        g.NaiveGrouping(L->n,  groupPtr, groupSet, ngroup, groupInv, blksize);
     this->ngroup = ngroup;
     this->npart = ngroup;

     /**
* apply grouping information to the DAG and generated a smaller DAG.
*/
     auto gDAG=Group_DAG(DAG, groupPtr, groupSet, groupInv, ngroup);


     size_t count=0;
     for (int j = 0; j < gDAG.size(); ++j) {
//            DAG[j].erase(std::unique(DAG[j].begin(), DAG[j].end()), DAG[j].end());
      count+=gDAG[j].size();
     }

     int *gv, *gedg;
     gv = new int[L->n+1]();
     gedg = new int[count+L->n]();
     int *levelPtr = new int[L->n+1]();
     int *levelSet = new int[L->n]();

     long int cti,edges=0;
     for(cti = 0, edges = 0; cti < ngroup; cti++){
      gv[cti] = edges;
      gedg[edges++] = cti;
      for (int ctj = 0; ctj < gDAG[cti].size(); ctj++) {
       gedg[edges++] = gDAG[cti][ctj];
      }
     }
     gv[cti] = edges;


     this->nlevels = buildLevelSet_CSC_Queue(ngroup, 0, gv, gedg, levelPtr, levelSet);

     this->num_sys = this->nlevels;


     this->nnzPerLevels = this->nnz * 1.0 /this->nlevels;
     this->averParallelism = this->ngroup * 1.0 / this->nlevels;

     std::vector<int> lcost;
     lcost.resize(this->nlevels);

     fs_csr_levelset_stat(L->p, L->i, groupPtr, groupSet, this->nlevels, levelPtr, levelSet, lcost.data());

     this->SumMaxDiff=0;
     for (auto &cost: lcost) {
      this->SumMaxDiff +=cost;
     }

     this->AverageMaxDiff = this->SumMaxDiff * 1.0 / lcost.size();

     double accum=0.0;
     std::for_each (std::begin(lcost), std::end(lcost), [&](const double d) {
         accum += (d - this->AverageMaxDiff) * (d - this->AverageMaxDiff);
     });
     this->VarianceMaxDiff = std::sqrt(accum/(lcost.size()-1));

    }


    StatSpMat::StatSpMat(CSR *L, SpKerType kerType, int num_threads, int *levelPtr, int *partPtr, int *nodePtr, int levelNo, int partNo,
                         int levelSetNo)
    {
     omp_set_num_threads(num_threads);

     Setup(L, kerType);

     this->ngroup=partNo;

     fs_csr_stat(L->n, L->p, L->i, this->nFlops, this->nnz_access, this->nnz_reuse);

     this->nlevels = levelSetNo;
     this->num_sys = levelNo;


     this->nnzPerLevels = this->nnz * 1.0 /this->num_sys;
     this->averParallelism = partNo* 1.0 / this->num_sys;


     std::vector<int> lcost;
     lcost.resize(this->num_sys);

     sptrsv_csr_lbc_stat(L->n, L->p, L->i,
                         levelNo, levelPtr,
                         partPtr, nodePtr, lcost.data());
//        printf("done\n");

     this->SumMaxDiff=0;
     for (auto &cost: lcost) {
      this->SumMaxDiff +=cost;
     }

//        this->SumMaxDiff = std::accumulate(lcost.begin(), lcost.end(), 0.0);
     this->AverageMaxDiff = this->SumMaxDiff * 1.0 / lcost.size();

     double accum=0.0;
     std::for_each (std::begin(lcost), std::end(lcost), [&](const double d) {
         accum += (d - this->AverageMaxDiff) * (d - this->AverageMaxDiff);
     });
     this->VarianceMaxDiff = std::sqrt(accum/(lcost.size()-1));

     lcost.clear();
    }


    StatSpMat::StatSpMat(CSR *L, SpKerType kerType, int num_threads, int *levelPtr, int *partPtr, int *nodePtr, int *groupPtr, int *groupSet,
                         int levelNo, int partNo, int levelSetNo, int groupNo)
    {
     omp_set_num_threads(num_threads);

     Setup(L, kerType);

     this->ngroup=groupNo;
     this->npart=levelPtr[levelNo];


     fs_csr_stat(L->n, L->p, L->i, this->nFlops, this->nnz_access, this->nnz_reuse);

     this->nlevels = levelSetNo;
     this->num_sys = levelNo;


     this->nnzPerLevels = this->nnz * 1.0 /this->num_sys;
     this->averParallelism = this->npart* 1.0 / this->num_sys;


     std::vector<int> lcost;
     lcost.resize(this->num_sys);

     sptrsv_csr_group_lbc_stat(L->n, L->p, L->i,
                               levelNo, levelPtr,
                               partPtr, nodePtr, groupPtr, groupSet, lcost.data());
//        printf("done\n");

     this->SumMaxDiff=0;
     for (auto &cost: lcost) {
      this->SumMaxDiff +=cost;
     }

//        this->SumMaxDiff = std::accumulate(lcost.begin(), lcost.end(), 0.0);
     this->AverageMaxDiff = this->SumMaxDiff * 1.0 / lcost.size();

     double accum=0.0;
     std::for_each (std::begin(lcost), std::end(lcost), [&](const double d) {
         accum += (d - this->AverageMaxDiff) * (d - this->AverageMaxDiff);
     });
     this->VarianceMaxDiff = std::sqrt(accum/(lcost.size()-1));

     lcost.clear();
    }


    void StatSpMat::Setup(CSR *L, SpKerType kerType) {
     this->t_serial=0;
     this->t_level=0;
     this->t_lbc=0;
     this->n = L->n;
     this->nnz = L->nnz;
     this->spkernel = kerType;
     this->NnzPerRows = L->nnz*1.0/L->n;
     this->numofcores = omp_get_max_threads();
    }


    void StatSpMat::PrintData() {
     PRINT_CSV(this->n);
     PRINT_CSV(this->nnz);
     PRINT_CSV(this->NnzPerRows);
     PRINT_CSV(this->ngroup);
     PRINT_CSV(this->npart);
     PRINT_CSV(this->nnz_access);
     PRINT_CSV(this->nnz_reuse);
     PRINT_CSV(this->nFlops);
     PRINT_CSV(this->nlevels);
     PRINT_CSV(this->num_sys);
     PRINT_CSV(this->averParallelism);
     PRINT_CSV(this->nnzPerLevels);
     PRINT_CSV(this->AverageMaxDiff);
     PRINT_CSV(this->VarianceMaxDiff);
     PRINT_CSV(this->SumMaxDiff);
     PRINT_CSV(this->numofcores);
     PRINT_CSV(this->t_serial);
     PRINT_CSV(this->t_group_level);
     PRINT_CSV(this->t_lbc);
     PRINT_CSV(this->t_level);

    }

}