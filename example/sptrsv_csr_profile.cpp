//
// Created by labuser (Bangtian Liu) on 9/16/20.
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

 double *y_serial, *y_correct = new double[n];


 timing_measurement t_ser, t_levelset, t_group, t_lbc;

 SptrsvSerial *ss = new SptrsvSerial(L2_csr, L1_csc, NULLPNTR, "serial");
 t_ser = ss->evaluate();
 y_serial = ss->solution();
 copy_vector(0,n,y_serial,y_correct);

 auto *sls = new SptrsvLevelSet(L2_csr, L1_csc, y_correct, "levelset csc");
 t_levelset = sls->evaluate();

 auto LevelSetNo = sls->getLevelNo();

 if (option==0){
  /**
   * profiling for levelset method, in which grouping can be enabled by setting p2 > 1
   * **/
  auto *sg = new SpTrsvCSR_Grouping(L2_csr, L1_csc, y_correct, "grouping code", num_threads);
  t_group = sg->evaluate();
  SpKerType ktype = SpTrsv_CSR;
  //        if(num_threads==1)num_threads=16;
  StatSpMat profiler(L2_csr, ktype, num_threads, p2);
  profiler.set_seq_time(t_ser.elapsed_time);
  profiler.set_level_time(t_levelset.elapsed_time);
  profiler.set_lbc_time(t_lbc.elapsed_time);
  profiler.set_glevel_time(t_group.elapsed_time);

  size_t pos = matrix_name.find_last_of("/\\");
  matrix_name = matrix_name.substr(pos+1);
  PRINT_CSV(matrix_name);
  PRINT_CSV(p2);
  profiler.PrintData();
  std::cout<<"\n";

  delete sg;
 }
 else if (option==1){
  /**
   * profiling for lbc method combined with grouping method
   * */

  auto *sglbc = new SpTrsvCSR_Grouping_H2(L2_csr, L1_csc, y_correct, "grouping_lbc", num_threads, p2, p3,false);
  t_lbc = sglbc->evaluate();
  auto levelPtr = sglbc->getLevelPtr();
  auto partPtr = sglbc->getPartPtr();
  auto nodePtr = sglbc->getNodePtr();
  auto levelNo = sglbc->getLevelNo();
  auto partNo = sglbc->getPartNo();
  auto ngroup = sglbc->getGroupNo();
  auto groupPtr = sglbc->getGroupPtr();
  auto groupSet = sglbc->getGroupSetPtr();

  SpKerType ktype = SpTrsv_CSR;
  StatSpMat profiler(L2_csr, ktype, num_threads, levelPtr, partPtr, nodePtr, groupPtr, groupSet,
                     levelNo, partNo, LevelSetNo, ngroup);


  profiler.set_seq_time(t_ser.elapsed_time);
  profiler.set_level_time(t_levelset.elapsed_time);
  profiler.set_lbc_time(t_lbc.elapsed_time);
  profiler.set_glevel_time(t_group.elapsed_time);

  size_t pos = matrix_name.find_last_of("/\\");
  matrix_name = matrix_name.substr(pos+1);
  PRINT_CSV(matrix_name);
  PRINT_CSV(p2);
  PRINT_CSV(p3);
  profiler.PrintData();
  std::cout<<"\n";
  delete  sglbc;
 }
 else{
  /**
   * profling for lbc method
   * */
  auto *slbc = new SptrsvLBC(L2_csr, L1_csc, y_serial, "lbc",num_threads, p2, p3);
  t_lbc = slbc->evaluate();

  auto levelPtr = slbc->getLevelPtr();
  auto partPtr = slbc->getPartPtr();
  auto nodePtr = slbc->getNodePtr();
  auto levelNo = slbc->getLevelNo();
  auto partNo = slbc->getPartNo();

  SpKerType ktype = SpTrsv_CSR;
  StatSpMat profiler(L2_csr, ktype, num_threads, levelPtr, partPtr, nodePtr, levelNo, partNo, LevelSetNo);
  profiler.set_seq_time(t_ser.elapsed_time);
  profiler.set_level_time(t_levelset.elapsed_time);
  profiler.set_lbc_time(t_lbc.elapsed_time);
  profiler.set_glevel_time(t_group.elapsed_time);
  size_t pos = matrix_name.find_last_of("/\\");
  matrix_name = matrix_name.substr(pos+1);
  PRINT_CSV(matrix_name);
  PRINT_CSV(p2);
  profiler.PrintData();
  std::cout<<"\n";
  delete  slbc;
 }


 delete []y_correct;
 delete A;
 delete L1_csc;
 delete L2_csr;

 delete ss;
//    delete sl;
 delete sls;
 return 0;
}
