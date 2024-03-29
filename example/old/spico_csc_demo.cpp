//
// Created by labuser (Bangtian Liu) on 10/3/20.
//

//#define CSV_LOG
#include "aggregation/def.h"
#include "aggregation/test_utils.h"
#include "aggregation/sparse_io.h"
#ifdef ENABLE_OPENMP
#include <omp.h>
#endif
#include "aggregation/metis_interface.h"
#include "spic0_demo_utils.h"

using namespace sym_lib;

int spic0_demo01(int argc, char **argv);


int main(int argc, char *argv[]) {
 int ret_val = 0;
 ret_val = spic0_demo01(argc, argv);
 return ret_val;
}

int spic0_demo01(int argc, char *argv[]) {
 CSC *L1_csc, *A = NULLPNTR;
 CSR *L2_csr;
 size_t n;
 int num_threads = 6;
 int p2 = -3, p3 = 4; // LBC params
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
 if (argc >= 3)
  p2 = atoi(argv[2]);
 if (argc >= 4)
  p3 = atoi(argv[3]);
 if (argc >= 5)
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
 CSR *B_csr;
 CSC *B;
 if (true) {
  B = copy_sparse(L1_csc); // will use the lower matrix here
  B_csr = csc_to_csr(B);
  //B = copy_sparse(L1_csc);
 } else {
  B = diagonal(L1_csc->n, 1.0);
  B_csr = csc_to_csr(B);
 }
 auto *factor_correct = new double[L2_csr->nnz]();
 timing_measurement t_serial, t_lbc, t_lbc_sort, t_levelset, t_level_group, t_lbc_group, t_lbc_group_sort;

 //print_csc(1,"CSR:\n",n,B_csr->p, B_csr->i, B_csr->x);

 auto *itsnf = new Spic0CSCSerial(L2_csr, L1_csc, B_csr, B,
                                  NULLPNTR,
                                  "Serial");
 t_serial = itsnf->evaluate();
 copy_vector(0, L2_csr->nnz, itsnf->Factor(), factor_correct);


 auto *itlevel = new SpicoCSCLevelSet(L2_csr, L1_csc, B_csr, B,
                                      factor_correct, "levelset");
 t_levelset = itlevel->evaluate();

 auto *itgroup = new SpicoCSC_Grouping(L2_csr, L1_csc, B_csr, B,
                                       factor_correct, "levelset grouping");
 t_level_group = itgroup->evaluate();
//    print_vec("X: \n", 0, 100, factor_correct);
//    print_vec("Y: \n", 0, 100, itlevel->Factor());
//    print_vec("Z: \n", 0, L1_csc->nnz, itgroup->Factor());

 auto *lbc = new Spic0CSCParallelLBC(L2_csr, L1_csc, B_csr, B,
                                     factor_correct,
                                     "Parallel LBC", num_threads,
                                     p2, p3, false);
 t_lbc = lbc->evaluate();

 auto *lbc_sort = new Spic0CSCParallelLBC(L2_csr, L1_csc, B_csr, B,
                                          factor_correct,
                                          "Parallel LBC sort", num_threads,
                                          p2, p3, true);
 t_lbc_sort = lbc_sort->evaluate();


 auto *lbc_group = new SpicoCSC_Grouping_H2(L2_csr, L1_csc, B_csr, B,
                                            factor_correct,
                                            "Parallel LBC grouping", num_threads,
                                            p2, p3, false);
 t_lbc_group = lbc_group->evaluate();

 auto *lbc_group_sort = new SpicoCSC_Grouping_H2(L2_csr, L1_csc, B_csr, B,
                                                 factor_correct,
                                                 "Parallel LBC grouping sort", num_threads,
                                                 p2, p3, true);
 t_lbc_group_sort = lbc_group_sort->evaluate();

 if (header)
  std::cout << "Matrix Name,Code Type,Data Type,Metis Enabled,"
               "Number of Threads,LBC Param1,LBC Param2,"
               "Serial Non-fused,"
               "Parallel LBC Non-fused CSR Analysis Time,"
               "Parallel LBC CSR,"
               "\n";
 size_t pos = matrix_name.find_last_of("/\\");
 matrix_name = matrix_name.substr(pos + 1);
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
 PRINT_CSV(itgroup->groupWidth());
 PRINT_CSV(lbc_group->groupWidth());
 PRINT_CSV(t_serial.elapsed_time);
 PRINT_CSV(t_levelset.elapsed_time);
 PRINT_CSV(t_level_group.elapsed_time);
 PRINT_CSV(t_lbc.elapsed_time);
 PRINT_CSV(t_lbc_sort.elapsed_time);
 PRINT_CSV(t_lbc_group.elapsed_time);
 PRINT_CSV(t_lbc_group_sort.elapsed_time);

 std::cout << "\n";

//
//    delete itsnf;
//    delete itplnf;

 delete A;
 delete B;
 delete B_csr;
 delete L1_csc;
 delete L2_csr;

 delete[]factor_correct;

 delete itlevel;
 delete itgroup;
 delete lbc;
 delete lbc_sort;
 delete lbc_group;
 delete lbc_group_sort;

 return 0;
}

