//
// Created by kazem on 2020-04-23.
//
#define CSV_LOG
#include <def.h>
#include <test_utils.h>
#include <sparse_io.h>
#include <omp.h>
#include <metis_interface.h>
#include "spic0_demo_utils.h"

using namespace sym_lib;

int spic0_demo01(int argc, char **argv);


int main(int argc, char *argv[]){
 int ret_val = 0;
 ret_val = spic0_demo01(argc,argv);
 return ret_val;
}

int spic0_demo01(int argc, char *argv[]) {
 CSC *L1_csc, *A = NULLPNTR;
 CSR *L2_csr;
 size_t n;
 int num_threads = 6;
 int p2 = -3, p3 = 4000; // LBC params
 int header = 0;
 int *perm;
 std::string matrix_name;
 std::vector<timing_measurement> time_array;

 /// Generating inputs
 if (argc < 2) {
  PRINT_LOG("Not enough input args, switching to random mode.\n");
  n = 15;
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
 if(argc == 6)
  header = atoi(argv[5]);
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
 timing_measurement t_serial, t_level_set, t_lbc ;

 //print_csc(1,"CSR:\n",n,B_csr->p, B_csr->i, B_csr->x);

 auto *itsnf = new Spic0Serial(L2_csr, L1_csc,B_csr, B,
                                             NULLPNTR,
                                             "Serial");
 t_serial = itsnf->evaluate();
 copy_vector(0, L2_csr->nnz, itsnf->Factor(), factor_correct);
 //print_vec("X: \n", 0, n, factor_correct);

 auto *spls = new Spic0ParallelLevelset(L2_csr, L1_csc, B_csr, B,
                                     factor_correct,
                                     "Parallel Level set");
 t_level_set = spls->evaluate();

 auto *itplnf = new Spic0ParallelLBC(L2_csr, L1_csc, B_csr, B,
                                             factor_correct,
                                             "Parallel LBC",num_threads,
                                             p2, p3);
 t_lbc = itplnf->evaluate();
 //print_vec("Y: \n", 0, n, itplnf->solution());



 if(header)
  std::cout<<"Matrix Name,Code Type,Data Type,Metis Enabled,"
             "Number of Threads,LBC Param1,LBC Param2,"
             "Serial Non-fused,"
             "Parallel Level set CSR Analysis Time,"
             "Parallel Level set CSR,"
             "Parallel LBC CSR Analysis Time,"
             "Parallel LBC CSR,"
             "\n";

 PRINT_CSV(matrix_name);
 PRINT_CSV("IC0");
#ifdef METIS
 PRINT_CSV("Metis");
#else
 PRINT_CSV("No Metis");
#endif
 PRINT_CSV(num_threads);
 PRINT_CSV(p2);
 PRINT_CSV(p3);

 PRINT_CSV(t_serial.elapsed_time);

 PRINT_CSV(spls->analysisTime().elapsed_time);
 PRINT_CSV(t_level_set.elapsed_time);

 PRINT_CSV(itplnf->analysisTime().elapsed_time);
 PRINT_CSV(t_lbc.elapsed_time);


 delete itsnf;
 delete spls;
 delete itplnf;

 delete A;
 delete B;
 delete B_csr;
 delete L1_csc;
 delete L2_csr;

 delete []factor_correct;

 return 1;
}

