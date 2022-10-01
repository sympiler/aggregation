//
// Created by kazem on 12/5/21.
//

//#define DBG_LOG
#define CSV_LOG
#include <iostream>
#include "aggregation/sparse_io.h"
#include "aggregation/test_utils.h"
#include "aggregation/sparse_utilities.h"

#ifdef ENABLE_OPENMP
#include "omp.h"
#endif

#ifdef METIS
#include "aggregation/metis_interface.h"
#endif
#include "sptrsv_profiler_utils.h"

using namespace sym_lib;

int profiler_demo(int argc, char *argv[]);

int main(int argc, char *argv[]){
 int ret_val;
 ret_val = profiler_demo(argc,argv);
 return ret_val;
}


int profiler_demo(int argc, char *argv[]){
 CSC *L1_csc, *A = NULLPNTR;
 CSR *L2_csr = NULLPNTR;
 size_t n;
 int num_threads = 20;
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

 if(argc >= 3)
  p2 = atoi(argv[2]);
 if(argc >= 4)
  p3 = atoi(argv[3]);
 if(argc >= 5)
  num_threads = atoi(argv[4]);

#ifdef ENABLE_OPENMP
 omp_set_num_threads(num_threads);
#endif
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

 timing_measurement t_c_tp, t_c_pp, t_c_sp, t_g_c_tp, t_g_c_sp, lt_t;

 /// Correct answer computation
 // Serial
 SptrsvSerial *ss = new SptrsvSerial(L2_csr, L1_csc, NULLPNTR, "serial"); //seq
 t_ser = ss->evaluate();
 y_serial = ss->solution();
 copy_vector(0,n,y_serial,y_correct);

 /// Profiling
 int event_limit =1, instance_per_run = 5;
 std::vector<int> event_list = {PAPI_L1_DCM, PAPI_L1_TCM, PAPI_L2_DCM,
                                PAPI_L2_DCA,PAPI_L2_TCM,PAPI_L2_TCA,
                                PAPI_L3_TCM, PAPI_L3_TCA,
                                PAPI_TLB_DM, PAPI_TOT_CYC, PAPI_LST_INS,
                                PAPI_TOT_INS, PAPI_REF_CYC, PAPI_RES_STL,
                                PAPI_BR_INS, PAPI_BR_MSP};

 auto *ssp = new SptrsvSerialProfiler(event_list,
                                                     event_limit, instance_per_run);
 ssp->profile(L2_csr, L1_csc, NULLPNTR, NULLPNTR, NULLPNTR,
                 y_correct,"Serial");
 int p2_nf = 4;

 if(header){
  print_common_header();
  ssp->print_headers();
  PRINT_CSV("NNZ Cost,Unit Cost,Redundant Nodes,Redundant NNZ,NNZ Parallelism");
  PRINT_CSV("Node Parallelism,Critical Path,Max-diff");
  std::cout<<"\n";
 }

 auto print_stat = [&](Profiler *p){
  for (int i = 0; i < instance_per_run; ++i) {
   print_common(matrix_name,"TRSV-PROF",p->Name(),L1_csc,L1_csc,num_threads);
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
 print_stat(ssp);

 delete []y_correct;
 delete A;
 delete L1_csc;
 delete L2_csr;

 delete ss;
 delete ssp;
 return 0;
}
