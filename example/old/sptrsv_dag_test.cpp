#include "aggregation/def.h"
#include "aggregation/utils.h"
#define DBG_LOG
#define CSV_LOG

#include <iostream>
#include "aggregation/metis_interface.h"
#ifdef ENABLE_OPENMP
#include <omp.h>
#endif
#include <regex>
#include "aggregation/sparse_io.h"
#include "aggregation/test_utils.h"

#include "sptrsv_demo_utils.h"

using namespace sym_lib;

int main(int argc, char *argv[]) {
 int ret_val;
 CSC *L1_csc, *A = NULLPNTR;
 CSR *L2_csr;
 size_t n;
 int num_threads = 6;
 int p2 = -1, p3 = 4000; // LBC params
 int header = 0;
 int *perm;
 int mode = 0; // 0: parallel l, 1: serial, 2: level set
 int creation_threads_temp = -1;
 std::string matrix_name;

 if (argc < 2) {
  PRINT_LOG("Not enough input args\n");
  return 1;
 }

 std::string f1 = argv[1];
 matrix_name = f1;
 L1_csc = read_mtx(f1);

 if (L1_csc == NULLPNTR)
  return -1;
 n = L1_csc->n;

 if (argc >= 3)
  p2 = atoi(argv[2]);
 if (argc >= 4)
  p3 = atoi(argv[3]);
 if (argc >= 5)
  num_threads = atoi(argv[4]);
 if (argc >= 6)
  mode = atoi(argv[5]);
 if (argc >= 7)
  creation_threads_temp = atoi(argv[6]);

 omp_set_num_threads(num_threads);

 int creation_threads =
  creation_threads_temp == -1 ? num_threads : creation_threads_temp;

 /// Re-ordering L matrix
#ifdef METIS
 // We only reorder L since dependency matters more in l-solve.
 // perm = new int[n]();
 CSC *L1_csc_full = make_full(L1_csc);
 delete L1_csc;
 metis_perm_general(L1_csc_full, perm);
 L1_csc =
  make_half(L1_csc_full->n, L1_csc_full->p, L1_csc_full->i, L1_csc_full->x);
 CSC *Lt = transpose_symmetric(L1_csc, perm);
 CSC *L1_ord = transpose_symmetric(Lt, NULLPNTR);
 delete L1_csc;
 L1_csc = L1_ord;
 delete Lt;
 delete L1_csc_full;
 delete[] perm;
#endif

 L2_csr = csc_to_csr(L1_csc);

 int final_level_no, *fina_level_ptr, *final_part_ptr, *final_node_ptr;
 int *level_set, *level_ptr;
 int part_no;

 std::regex re("^.+\\/(.*?)\\..*$");
 std::smatch m;
 std::regex_search(matrix_name, m, re);

 auto *cost = new double[n]();
 for (int i = 0; i < n; ++i) {
  cost[i] = L2_csr->p[i + 1] - L2_csr->p[i];
 }

 std::vector<timing_measurement> time_array;

 for (int i = 0; i < 10; ++i) {
  timing_measurement time;
  if (mode == 2) {
   time.start_timer();
   build_levelSet_CSC(L1_csc->n, L1_csc->p, L1_csc->i, level_ptr, level_set);
   time.measure_elapsed_time();
   delete[] level_set;
   delete[] level_ptr;
  } else {
   if (mode == 1) {
    time.start_timer();
    get_coarse_Level_set_DAG_CSC03(n, L1_csc->p, L1_csc->i, final_level_no,
                                   fina_level_ptr, part_no, final_part_ptr,
                                   final_node_ptr, num_threads, p2, p3, cost);
    time.measure_elapsed_time();
   } else {
    time.start_timer();
    get_coarse_Level_set_DAG_CSC03_parallel(
     n, L1_csc->p, L1_csc->i, final_level_no, fina_level_ptr, part_no,
     final_part_ptr, final_node_ptr, num_threads, p2, p3, cost,
     creation_threads);
    time.measure_elapsed_time();
   }

   delete[] fina_level_ptr;
   delete[] final_part_ptr;
   delete[] final_node_ptr;
  }
  time_array.emplace_back(time);
 }

 timing_measurement median_t = time_median(time_array);
 PRINT_CSV(m[1]);
 PRINT_CSV(median_t.elapsed_time);
 if (creation_threads_temp != -1) {
  PRINT_CSV(creation_threads_temp);
 }
 std::cout << std::endl;
 delete[] cost;
}
