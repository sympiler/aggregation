//
// Created by labuser (Bangtian Liu) on 12/22/20.
//

//
// Created by kazem on 10/12/19.
//

#define DBG_LOG
#define CSV_LOG

#include <iostream>
#include "aggregation/sparse_io.h"
#include "aggregation/test_utils.h"
#ifdef ENABLE_OPENMP
#include <omp.h>
#endif
#include "aggregation/metis_interface.h"
#include "sptrsvprofiler_demo_utils.h"
#include "sptrsv_demo_utils.h"

using namespace sym_lib;

/// Evaluate spmv-sptrsv based on random matrices
/// \return
int sptrsv_csr_profiler_demo02(int argc, char *argv[]);

int main(int argc, char *argv[]){
    int ret_val;
    ret_val = sptrsv_csr_profiler_demo02(argc,argv);
    return ret_val;
}



int sptrsv_csr_profiler_demo02(int argc, char *argv[]){
    CSC *L1_csc, *A = NULLPNTR;
    CSR *L2_csr;
    size_t n;
    int num_threads = 6;
    int p2 = -1, p3 = 4000; // LBC params
    int header = 1;
    int *perm;
    std::string matrix_name;
    std::vector<timing_measurement> time_array;

    if (argc < 2) {
        PRINT_LOG("Not enough input args, switching to random mode.\n");
        n = 16;
        double density = 0.2;
        matrix_name = "Random_" + std::to_string(n);
        A = random_square_sparse(n, density);
        if (A == NULLPNTR)
            return -1;
        L1_csc = make_half(A->n, A->p, A->i, A->x);
    } else {
        std::string f1 = argv[1];
        matrix_name = f1;
        L1_csc = read_mtx(f1);
        if (L1_csc == NULLPNTR)
            return -1;
        n = L1_csc->n;
    }

    if(argc >= 3)
        p2 = atoi(argv[2]);
    if(argc >= 4)
        p3 = atoi(argv[3]);
    if(argc >= 5)
        num_threads = atoi(argv[4]);

    omp_set_num_threads(num_threads);
    /// Re-ordering L matrix
#ifdef METIS
    //We only reorder L since dependency matters more in l-solve.
    //perm = new int[n]();
    CSC *L1_csc_full = make_full(L1_csc);
    delete L1_csc;
    metis_perm_general(L1_csc_full, perm);
    L1_csc = make_half(L1_csc_full->n, L1_csc_full->p, L1_csc_full->i,
                       L1_csc_full->x);
    CSC *Lt = transpose_symmetric(L1_csc, perm);
    CSC *L1_ord = transpose_symmetric(Lt, NULLPNTR);
    delete L1_csc;
    L1_csc = L1_ord;
    delete Lt;
    delete L1_csc_full;
    delete[]perm;
#endif

    L2_csr = csc_to_csr(L1_csc);

    double *y_serial, *y_correct = new double[n];

    timing_measurement t_ser, t_levelset, t_levelset_group;

    timing_measurement t_c_tp, t_c_pp, t_c_sp, t_g_c_tp, t_g_c_sp;

    SptrsvSerial *ss = new SptrsvSerial(L2_csr, L1_csc, NULLPNTR, "serial"); //seq
    t_ser = ss->evaluate();
    y_serial = ss->solution();
    copy_vector(0,n,y_serial,y_correct);
    //print_vec("x:\n", 0, n, y_correct);

    /// Profiling
    int event_limit =1, instance_per_run = 5;
    std::vector<int> event_list = {PAPI_L1_DCM, PAPI_L1_TCM, PAPI_L2_DCM,
                                   PAPI_L2_DCA,PAPI_L2_TCM,PAPI_L2_TCA,
                                   PAPI_L3_TCM, PAPI_L3_TCA,
                                   PAPI_TLB_DM, PAPI_TOT_CYC, PAPI_LST_INS,
                                   PAPI_TOT_INS, PAPI_REF_CYC, PAPI_RES_STL,
                                   PAPI_BR_INS, PAPI_BR_MSP};

    auto *plevelset = new SptrsvLevelSetProfiler(event_list, event_limit, instance_per_run);

    plevelset->profile(L2_csr, L1_csc, y_correct, "levelset csr", num_threads); // 12 is number of threads, can be

    auto *pgroup = new SpTrsvCSR_Grouping_Profiler(event_list, event_limit, instance_per_run);

    pgroup->profile(L2_csr, L1_csc, y_correct, "levelset with grouping", num_threads);


    auto *plbcdag = new SptrsvLBCDAGProfiler(event_list, event_limit, instance_per_run);
    plbcdag->profile(L2_csr, L1_csc, y_correct, "coarsening", num_threads, p2, p3);

    auto *plbcdag_sorted = new SptrsvLBC_W_SortingProfiler(event_list, event_limit, instance_per_run);
    plbcdag_sorted->profile(L2_csr, L1_csc, y_correct, "coarsening", num_threads, p2, p3, true);

    auto *pglbcdag = new SpTrsvCSR_Grouping_H2_Profiler(event_list, event_limit, instance_per_run);
    pglbcdag->profile(L2_csr, L1_csc, y_serial, "grouping + coarsening", num_threads, p2, p3, false);

    auto *pglbcdag_sort = new SpTrsvCSR_Grouping_H2_Profiler(event_list, event_limit, instance_per_run);
    pglbcdag_sort->profile(L2_csr, L1_csc, y_serial, "grouping + coarsening", num_threads, p2, p3, true);



    if(header){
        print_common_header();
        plevelset->print_headers();
        PRINT_CSV("NNZ Cost,Unit Cost,Redundant Nodes,Redundant NNZ,NNZ Parallelism");
        PRINT_CSV("Node Parallelism,Critical Path,Max-diff");
        std::cout<<"\n";
    }
    size_t pos = matrix_name.find_last_of("/\\");
    matrix_name = matrix_name.substr(pos+1);

    auto print_stat = [&](Profiler *p){
        for (int i = 0; i < instance_per_run; ++i) {
            print_common(matrix_name,"IC0-TRSV-PROF",p->Name(),L1_csc,L1_csc,num_threads);
            p->print_counters(i);
            PRINT_CSV(p->CostNNZ());
            PRINT_CSV(p->CostUnit());
            PRINT_CSV(p->redundantNodes());
            PRINT_CSV(p->redundantNNZ());
//   PRINT_CSV(2*B->nnz);// we have only one level here
            PRINT_CSV(p->avgParallelism());
            PRINT_CSV(p->avgIterParallelism());
            PRINT_CSV(p->criticalPath());
            PRINT_CSV(p->maxDiff());
            std::cout<<"\n";
        }
    };

    print_stat(plevelset);
    print_stat(pgroup);
    print_stat(plbcdag);
    print_stat(plbcdag_sorted);
    print_stat(pglbcdag);
    print_stat(pglbcdag_sort);


    delete []y_correct;
    delete A;
    delete L1_csc;
    delete L2_csr;

    delete ss;
    delete plevelset;
    delete pgroup;
    delete plbcdag;
    delete plbcdag_sorted;
    delete pglbcdag;
    delete pglbcdag_sort;

    return 0;
}
