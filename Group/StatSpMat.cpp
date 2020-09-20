//
// Created by labuser (Bangtian Liu) on 9/15/20.
//

#include <StatSpMat.h>

#include "../example/sparse_blas_lib.h"

using namespace sym_lib;
namespace group_cols{
    StatSpMat::StatSpMat(CSR *L, SpKerType kerType, int num_threads, int blksize)
    {
        omp_set_num_threads(num_threads);

        this->t_serial=0;
        this->t_level=0;
        this->t_lbc=0;
        this->n = L->n;
        this->nnz = L->nnz;
        this->spkernel = kerType;
        this->NnzPerRows = L->nnz*1.0/L->n;
        this->numofcores = omp_get_max_threads();

        fs_csr_stat(L->n, L->p, L->i, this->nFlops, this->nnz_access, this->nnz_reuse);


        int *groupPtr = (int *)malloc(sizeof(int)*(L->n+1));
        memset(groupPtr, 0, sizeof(int)*(1+L->n));
        int *groupSet = (int *)malloc(sizeof(int)*L->n);
        memset(groupSet, 0, sizeof(int)*L->n);
        int *groupInv = (int *)malloc(sizeof(int)*L->n);
        memset(groupInv, 0, sizeof(int)*L->n);

//        double *x = (double *)malloc(sizeof(double)*L->n);
//        memset(x, 0, sizeof(double)*L->n);

        int ngroup;
        group g(L->n, L->p, L->i);
//        g.inspection_sptrsvcsr(groupPtr, groupSet, ngroup, groupInv);

        NaiveGrouping(L->n,  groupPtr, groupSet, ngroup, groupInv, blksize);

        this->ngroup = ngroup;


        std::vector<std::vector<int>> DAG;
        DAG.resize(ngroup);

        fs_csr_inspector_dep(ngroup, groupPtr, groupSet, groupInv, L->p, L->i, DAG);

        size_t count=0;
        for (int j = 0; j < DAG.size(); ++j) {
            DAG[j].erase(std::unique(DAG[j].begin(), DAG[j].end()), DAG[j].end());
            count+=DAG[j].size();
        }
//   detectDAGCircle(DAG);

        int *gv, *gedg;
        gv = new int[L->n+1]();
        gedg = new int[count+L->n]();
        int *levelPtr = new int[L->n+1]();
        int *levelSet = new int[L->n]();

        long int cti,edges=0;
        for(cti = 0, edges = 0; cti < ngroup; cti++){
            gv[cti] = edges;
            gedg[edges++] = cti;
            for (int ctj = 0; ctj < DAG[cti].size(); ctj++) {
                gedg[edges++] = DAG[cti][ctj];
//                if(DAG[cti][ctj]==0)printf("cti=%d, ctj=%d\n", cti, ctj);
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


//        rhsInit_csr(L->n, L->p, L->i, L->x, x);

//        sptrsv_csr_levelset(L->n, L->p, L->i, L->x)


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
    }


    StatSpMat::StatSpMat(CSR *L, SpKerType kerType, int num_threads, int lparm, int divrate)
    {
        omp_set_num_threads(num_threads);

        this->t_serial=0;
        this->t_level=0;
        this->t_lbc=0;
        this->n = L->n;
        this->nnz = L->nnz;
        this->spkernel = kerType;
        this->NnzPerRows = L->nnz*1.0/L->n;
        this->numofcores = omp_get_max_threads();

        fs_csr_stat(L->n, L->p, L->i, this->nFlops, this->nnz_access, this->nnz_reuse);


        int *groupPtr = (int *)malloc(sizeof(int)*(L->n+1));
        memset(groupPtr, 0, sizeof(int)*(1+L->n));
        int *groupSet = (int *)malloc(sizeof(int)*L->n);
        memset(groupSet, 0, sizeof(int)*L->n);
        int *groupInv = (int *)malloc(sizeof(int)*L->n);
        memset(groupInv, 0, sizeof(int)*L->n);

//        double *x = (double *)malloc(sizeof(double)*L->n);
//        memset(x, 0, sizeof(double)*L->n);

        int ngroup;
        group g(L->n, L->p, L->i);
//        g.inspection_sptrsvcsr(groupPtr, groupSet, ngroup, groupInv);

        NaiveGrouping(L->n,  groupPtr, groupSet, ngroup, groupInv, 1);

        this->ngroup = ngroup;


        std::vector<std::vector<int>> DAG;
        DAG.resize(ngroup);

        fs_csr_inspector_dep(ngroup, groupPtr, groupSet, groupInv, L->p, L->i, DAG);

        size_t count=0;
        for (int j = 0; j < DAG.size(); ++j) {
            DAG[j].erase(std::unique(DAG[j].begin(), DAG[j].end()), DAG[j].end());
            count+=DAG[j].size();
        }
//   detectDAGCircle(DAG);

        int *gv, *gedg;
        gv = new int[L->n+1]();
        gedg = new int[count+L->n]();
        int *levelPtr = new int[L->n+1]();
        int *levelSet = new int[L->n]();

        long int cti,edges=0;
        for(cti = 0, edges = 0; cti < ngroup; cti++){
            gv[cti] = edges;
            gedg[edges++] = cti;
            for (int ctj = 0; ctj < DAG[cti].size(); ctj++) {
                gedg[edges++] = DAG[cti][ctj];
//                if(DAG[cti][ctj]==0)printf("cti=%d, ctj=%d\n", cti, ctj);
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


//        rhsInit_csr(L->n, L->p, L->i, L->x, x);

//        sptrsv_csr_levelset(L->n, L->p, L->i, L->x)


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
    }





    void StatSpMat::PrintData() {

        PRINT_CSV(this->n);
        PRINT_CSV(this->nnz);
        PRINT_CSV(this->NnzPerRows);
        PRINT_CSV(this->ngroup);
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
        PRINT_CSV(this->t_level);

    }

}