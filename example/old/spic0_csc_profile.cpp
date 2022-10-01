//
// Created by labuser (Bangtian Liu) on 10/18/20.
//

//
// Created by labuser (Bangtian Liu) on 9/16/20.
//

#define DBG_LOG
#define CSV_LOG

#include <iostream>
#include "aggregation/sparse_io.h"
#include "aggregation/test_utils.h"
#include <omp.h>
#include "aggregation/metis_interface.h"
#include <StatSpMat_v1.h>

#include "spic0_demo_utils.h"

using namespace sym_lib;


int sptrsv_csr_profile_demo02(int argc, char *argv[]);

int main(int argc, char *argv[]){
    int ret_val;
    ret_val = sptrsv_csr_profile_demo02(argc,argv);
    return ret_val;
}


int sptrsv_csr_profile_demo02(int argc, char *argv[]){
    CSC *L1_csc, *A = NULLPNTR;
    CSR *L2_csr;
    size_t n;
    int num_threads = 6;
    int p2 = -1, p3 = 4000; // LBC params
    int header = 0;
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

    int option=0;
    if(argc >= 3)
        num_threads=atoi(argv[2]);
    omp_set_num_threads(num_threads);
    if(argc >= 4)
        option=atoi(argv[3]);
    if(argc >= 5)
        p2 = atoi(argv[4]);
    if(argc >= 6)
        p3 = atoi(argv[5]);
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
    CSR *B_csr;
    CSC *B;
    if(true){
        B = copy_sparse(L1_csc); // will use the lower matrix here
        B_csr = csc_to_csr(B);
        //B = copy_sparse(L1_csc);
    } else {
        B = diagonal(L1_csc->n,1.0);
        B_csr = csc_to_csr(B);
    }
    auto *factor_correct = new double[L2_csr->nnz]();

    timing_measurement t_serial, t_lbc, t_lbc_sort, t_levelset, t_level_group, t_lbc_group, t_lbc_group_sort;


    auto *itsnf = new Spic0CSCSerial(L2_csr, L1_csc,B_csr, B,
                                     NULLPNTR,
                                     "Serial");
    t_serial = itsnf->evaluate();
    copy_vector(0, L2_csr->nnz, itsnf->Factor(), factor_correct);


    auto *itlevel = new SpicoCSCLevelSet(L2_csr, L1_csc, B_csr, B,
                                         factor_correct, "levelset");
    t_levelset = itlevel->evaluate();
    auto ins_levelset=itlevel->analysisTime();

    auto *itgroup = new SpicoCSC_Grouping(L2_csr, L1_csc, B_csr, B,
                                          factor_correct, "levelset grouping");
    t_level_group = itgroup->evaluate();
    auto ins_g = itgroup->groupTime();
    auto ins_g_levelset = itgroup->levelsetTime();

    auto *lbc = new Spic0CSCParallelLBC(L2_csr, L1_csc, B_csr, B,
                                        factor_correct,
                                        "Parallel LBC",num_threads,
                                        p2, p3, false);
    t_lbc = lbc->evaluate();
    auto ins_coarsen = lbc->analysisTime();

    auto *lbc_sort = new Spic0CSCParallelLBC(L2_csr, L1_csc, B_csr, B,
                                             factor_correct,
                                             "Parallel LBC",num_threads,
                                             p2, p3, true);
    t_lbc_sort = lbc_sort->evaluate();

    auto ins_coarsen_psort = lbc->analysisTime();

    auto *lbc_group = new SpicoCSC_Grouping_H2(L2_csr, L1_csc, B_csr, B,
                                               factor_correct,
                                               "Parallel LBC grouping",num_threads,
                                               p2, p3, false);
    t_lbc_group = lbc_group->evaluate();

    auto ins_glbc_group = lbc_group->groupTime();
    auto ins_glbc_coarsen = lbc_group->coarsenTime();


    auto *lbc_group_sort = new SpicoCSC_Grouping_H2(L2_csr, L1_csc, B_csr, B,
                                               factor_correct,
                                               "Parallel LBC grouping",num_threads,
                                               p2, p3, true);

    t_lbc_group_sort = lbc_group_sort->evaluate();

    auto ins_glbc_sort = lbc_group_sort->sortTime();

    size_t pos = matrix_name.find_last_of("/\\");
    matrix_name = matrix_name.substr(pos+1);

    PRINT_CSV(matrix_name);
    PRINT_CSV(p2);
    PRINT_CSV(p3);
    PRINT_CSV(ins_levelset.elapsed_time);
    PRINT_CSV(ins_g.elapsed_time);
    PRINT_CSV(ins_g_levelset.elapsed_time);
    PRINT_CSV(ins_coarsen.elapsed_time);
    PRINT_CSV(ins_coarsen_psort.elapsed_time-ins_coarsen.elapsed_time);
    PRINT_CSV(ins_glbc_group.elapsed_time);
    PRINT_CSV(ins_glbc_coarsen.elapsed_time);
    PRINT_CSV(ins_glbc_sort.elapsed_time);

    PRINT_CSV(itgroup->groupWidth());
    PRINT_CSV(lbc_group->groupWidth());
    PRINT_CSV(itlevel->numLevels());
    PRINT_CSV(itlevel->averparallelism());
    PRINT_CSV(itgroup->numLevels());
    PRINT_CSV(itgroup->averparallelism());


    PRINT_CSV(lbc->numLevels());
    PRINT_CSV(lbc->averparallelism());
    PRINT_CSV(lbc->consecutiveRatio());
    PRINT_CSV(lbc_sort->consecutiveRatio());

    PRINT_CSV(lbc_group->numLevels());
    PRINT_CSV(lbc_group->averparallelism());
    PRINT_CSV(lbc_group->consecutiveRatio());
    PRINT_CSV(lbc_group_sort->consecutiveRatio());


    SpKerType ktype = SpInChol_CSC;
    StatSpMat_v1 *profiler = new StatSpMat_v1();
    profiler->Setup(B, ktype);
    profiler->PrintData();

    delete []factor_correct;
    delete A;
    delete L1_csc;
    delete L2_csr;

    delete itsnf;
    delete itlevel;
    delete itgroup;
    delete lbc;
    delete lbc_sort;
    delete lbc_group;
    delete lbc_group_sort;

    return 0;
}
