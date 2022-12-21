#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <omp.h>
#include <iostream>

#include "sptrsv.h"
#include "aggregation/def.h"
#include "aggregation/hdagg.h"
#include "aggregation/sparse_inspector.h"
#include "aggregation/sparse_io.h"
#include "aggregation/sparse_utilities.h"
#include "aggregation/test_utils.h"


TEST_CASE("Check hdagg sptrsv", "[hdagg]") {
    // prepare input matrix
    auto n = GENERATE(20, 50, 100, 200);
    auto density = GENERATE(0.05, 0.1, 0.2, 0.5);
    auto seed = GENERATE(1U, 2U, 3U);
    INFO("n = " << n << "; density = " << density << "; seed = " << seed);

    sym_lib::CSC *A = sym_lib::random_square_sparse(n, density, 1.0, seed);
    REQUIRE(A != NULLPNTR);
    sym_lib::CSC *Lower_A_CSC = sym_lib::make_half(A->n, A->p, A->i, A->x);
    delete A;

    bool use_metis = GENERATE(false, true);
    INFO("use_metis = " << use_metis);
    if (use_metis) {
        A = sym_lib::make_full(Lower_A_CSC);
        delete Lower_A_CSC;
        int *perm;
        sym_lib::metis_perm_general(A, perm);
        Lower_A_CSC = sym_lib::make_half(A->n, A->p, A->i, A->x);
        delete A;
        sym_lib::CSC *Lt = sym_lib::transpose_symmetric(Lower_A_CSC, perm);
        delete Lower_A_CSC;
        Lower_A_CSC = sym_lib::transpose_symmetric(Lt, NULLPNTR);
        delete Lt;
        delete[] perm;
    } 
    sym_lib::CSR *Lower_A_CSR = sym_lib::csc_to_csr(Lower_A_CSC);

    // one thread will fail Tree_HDagg https://github.com/sympiler/aggregation/issues/5#issuecomment-1357255720
    int nthreads = GENERATE(2, 4);
    INFO("nthreads = " << nthreads);
    #pragma omp parallel default(shared)
    {
        omp_set_num_threads(nthreads);
    }

    // construct Tree_HDAGG partitoning
    // copy from `SpTrSv_LL_Tree_HDAGG` in example/SpTRSV_runtime.h
    bool isLfactor = false;
    auto bin_pack = GENERATE(false, true);
    INFO("bin_pack = " << bin_pack);

    int ngroups;
    std::vector<int> final_level_ptr, final_part_ptr, final_node_ptr;
    int final_level_no;
    std::vector<int> DAG_ptr;
    std::vector<int> DAG_set;
    std::vector<int> group_set, group_ptr;
    std::vector<int> level_ptr, level_set;
    int nlevels;
    std::vector<double> cost;

    sym_lib::CSR *L1_csr_ = Lower_A_CSR;
    sym_lib::CSC *L1_csc_ = Lower_A_CSC;
    int n_ = Lower_A_CSC->n;
    int nnz_ = Lower_A_CSC->nnz;
    HDAGG::partialSparsification(n_, nnz_, L1_csc_->p, L1_csc_->i, DAG_ptr,
                                DAG_set);
    HDAGG::treeBasedGrouping(n_, DAG_ptr, DAG_set, ngroups, group_ptr,
                             group_set, isLfactor);
    std::vector<int> group_DAG_ptr, group_DAG_set;
    HDAGG::buildGroupDAG(n_, ngroups, group_ptr.data(), group_set.data(),
                         DAG_ptr.data(), DAG_set.data(), group_DAG_ptr,
                         group_DAG_set);

    cost.resize(ngroups, 0);
    auto CSC_Lp = L1_csc_->p;
    auto CSC_Li = L1_csc_->i;
    auto CSR_Lp = L1_csr_->p;
    auto CSR_Li = L1_csr_->i;
    HDAGG::costComputation(ngroups, CSC_Lp, CSC_Li, CSR_Lp, CSR_Li,
                           HDAGG::SpTrSv_LL, group_ptr.data(), group_set.data(),
                           true, cost);
    HDAGG::HDAGG(ngroups, group_DAG_ptr[ngroups], group_DAG_ptr, group_DAG_set,
                 cost, nthreads, final_level_no, final_level_ptr,
                 final_part_ptr, final_node_ptr, false, false, bin_pack);
    HDAGG::ungroupingScheduleAndApplyOrdering(
        n_, final_level_no, final_level_ptr, final_part_ptr, final_node_ptr,
        group_ptr, group_set);

    // apply parallel schedule
    double *b = new double[n_];
    double *x_parallel = new double[n_];
    std::fill_n(b, n_, 1.0);
    sym_lib::copy_vector(0, n_, b, x_parallel);
    sym_lib::sptrsv_csr_lbc(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_parallel,
                   final_level_no, final_level_ptr.data(),
                   final_part_ptr.data(), final_node_ptr.data());

    // compare with serial reference result
    double *x_serial = new double[n_];
    sym_lib::copy_vector(0, n_, b, x_serial);
    sym_lib::sptrsv_csr(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_serial);
    CHECK(sym_lib::is_equal(0, n_, x_serial, x_parallel));

    delete[] b;
    delete[] x_parallel;
    delete[] x_serial;
    delete L1_csc_;  // same as Lower_A_CSC
    delete L1_csr_;  // same as Lower_A_CSR
}
