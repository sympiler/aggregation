//
// Created by kazem on 10/12/19.
//

#define DBG_LOG
#define CSV_LOG

#include <iostream>
#include <sparse_io.h>
#include <test_utils.h>
#include <omp.h>
#include <metis_interface.h>

#include "sptrsv_demo_utils.h"

using namespace sym_lib;

/// Evaluate spmv-sptrsv based on random matrices
/// \return
int sptrsv_csc_demo02(int argc, char *argv[]);

int main(int argc, char *argv[]){
 int ret_val;
 ret_val = sptrsv_csc_demo02(argc,argv);
 return ret_val;
}



int sptrsv_csc_demo02(int argc, char *argv[]){
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

 timing_measurement t_c_tp, t_c_sp, t_g_c_tp, t_g_c_sp, lt_t;

 SptrsvSerial *ss = new SptrsvSerial(L2_csr, L1_csc, NULLPNTR, "serial"); //seq
 t_ser = ss->evaluate();
 y_serial = ss->solution();
 copy_vector(0,n,y_serial,y_correct);
 //print_vec("x:\n", 0, n, y_correct);

 auto *sls = new SptrsvLevelSet(L2_csr, L1_csc, y_correct, "levelset csc"); // levelset
 t_levelset = sls->evaluate();

 auto *sg = new SpTrsvCSR_Grouping(L2_csr, L1_csc, y_correct, "levelset with grouping", num_threads);
 t_levelset_group = sg->evaluate();

 auto *lbc_tree = new SptrsvLBC(L2_csr, L1_csc, y_serial, "LBC Tree",num_threads, p2, p3); // ng + c + tp
 lt_t = lbc_tree->evaluate();

 auto *sld = new SptrsvLBCDAG(L2_csr, L1_csc, y_serial, "coarsening levels",num_threads, p2, p3); // ng + c + tp
 t_c_tp = sld->evaluate();
//    auto t_coarsen_4=sld->analysisTime();

 auto *sld_sort = new SptrsvLBC_W_Sorting(L2_csr, L1_csc, y_serial, "c4+sorting", num_threads, p2, p3, true);
 t_c_sp = sld_sort->evaluate();

 auto *sglbc = new SpTrsvCSR_Grouping_H2(L2_csr, L1_csc, y_correct, "g_c4", num_threads, p2, p3, false);
 t_g_c_tp = sglbc->evaluate(); // g  + c + tp

 auto *sglbc_sort = new SpTrsvCSR_Grouping_H2(L2_csr, L1_csc, y_correct, "g_c4_sorting", num_threads, p2, p3, true);
 t_g_c_sp = sglbc_sort->evaluate(); // g + c + sp;


 if(header)
  std::cout<<"Matrix Name,Metis Enabled,"
             "Number of Threads,"
             "Serial Non-fused,Parallel Levelset CSC,Parallel LBC CSR,";
 size_t pos = matrix_name.find_last_of("/\\");
 matrix_name = matrix_name.substr(pos+1);
 PRINT_CSV(matrix_name);
#ifdef METIS
 PRINT_CSV("METIS");
#else
 PRINT_CSV("No");
#endif
 PRINT_CSV(num_threads);
 PRINT_CSV(p2);
 PRINT_CSV(p3);

 PRINT_CSV(t_ser.elapsed_time);
 PRINT_CSV(t_levelset.elapsed_time);
 PRINT_CSV(t_levelset_group.elapsed_time);

 PRINT_CSV(lt_t.elapsed_time);
 PRINT_CSV(t_c_tp.elapsed_time);
 PRINT_CSV(t_c_sp.elapsed_time);
 PRINT_CSV(t_g_c_tp.elapsed_time);
 PRINT_CSV(t_g_c_sp.elapsed_time);

 std::cout<<"\n";

 delete []y_correct;
 delete A;
 delete L1_csc;
 delete L2_csr;

 delete ss;
 delete sls;
 delete sg;
 delete sld;
 delete lbc_tree;
 delete sld_sort;
 delete sglbc;
 delete sglbc_sort;


 return 0;
}
