//
// Created by behrooz on 1/17/21.
//

#include "aggregation/csv_utils.h"
#include <iostream>
#include "aggregation/metis_interface.h"
#include <omp.h>
#include "aggregation/sparse_io.h"
#include "aggregation/test_utils.h"

#include "SpTRSV_runtime.h"

using namespace sym_lib;

int main(int argc, char *argv[]) {

  CSC *Lower_A_CSC, *A = NULLPNTR;
  CSR *Lower_A_CSR;
  size_t n, nnz;
  int nthreads = -1;
  int *perm;
  std::string matrix_name;
  std::vector<timing_measurement> time_array;

  if (argc < 2) {
    PRINT_LOG("Not enough input args, switching to random mode.\n");
    n = 50;
    double density = 0.3;
    matrix_name = "Random_" + std::to_string(n);
    A = random_square_sparse(n, density);
    if (A == NULLPNTR)
      return -1;
    Lower_A_CSC = make_half(A->n, A->p, A->i, A->x);
    delete A;

    //        std::vector<std::vector<bool>> DAG(n, std::vector<bool>(n,false));
    //        addEdge(DAG, 0, 1);
    //        addEdge(DAG, 0, 3);
    //        addEdge(DAG, 0, 4);
    //        addEdge(DAG, 1, 2);
    //        addEdge(DAG, 3, 4);
    //        addEdge(DAG, 3, 8);
    //        addEdge(DAG, 3, 11);
    ////        addEdge(DAG, 4, 11);
    //        addEdge(DAG, 5, 7);
    //        addEdge(DAG, 5, 8);
    //        addEdge(DAG, 6, 7);
    //        addEdge(DAG, 6, 8);
    //        addEdge(DAG, 7, 8);
    //        addEdge(DAG, 7, 9);
    //        addEdge(DAG, 8, 9);
    //        addEdge(DAG, 8, 10);
    //        addEdge(DAG, 8, 11);
    //        addEdge(DAG, 8, 12);
    //        addEdge(DAG, 11, 12);
    //
    //        Lower_A_CSC = generateDAGfromEdgeList(n, DAG);
    //        std::vector<int> DAG_ptr(Lower_A_CSC->p, Lower_A_CSC->p + n + 1);
    //        std::vector<int> DAG_set(Lower_A_CSC->i, Lower_A_CSC->i +
    //        Lower_A_CSC->nnz);
    nnz = Lower_A_CSC->nnz;
  } else {
    std::cout << "Reading matrix " << argv[1] << std::endl;
    std::string f1 = argv[1];
    matrix_name = f1;
    Lower_A_CSC = read_mtx(f1);
    if (Lower_A_CSC == NULLPNTR) {
      std::cout << "Reading failed" << std::endl;
      return -1;
    }
    n = Lower_A_CSC->n;
  }

  if (argc >= 3) {
    nthreads = atoi(argv[2]);
  }

  /// Re-ordering matrix A
  //    std::cout << "METIS IS NOT ACTIVATED" << std::endl;
#ifdef METIS
  std::cout << "METIS IS ACTIVATED" << std::endl;
  // We only reorder A since dependency matters more in l-solve.
  A = make_full(Lower_A_CSC);
  delete Lower_A_CSC;
  metis_perm_general(A, perm);
  Lower_A_CSC = make_half(A->n, A->p, A->i, A->x);
  CSC *Lt = transpose_symmetric(Lower_A_CSC, perm);
  CSC *L1_ord = transpose_symmetric(Lt, NULLPNTR);
  delete Lower_A_CSC;
  Lower_A_CSC = L1_ord;
  Lower_A_CSR = csc_to_csr(Lower_A_CSC);
  delete Lt;
  delete A;
  delete[] perm;
#else
  CSC *tmp =
      make_half(Lower_A_CSC->n, Lower_A_CSC->p, Lower_A_CSC->i, Lower_A_CSC->x);
  delete Lower_A_CSC;
  Lower_A_CSC = tmp;
  Lower_A_CSR = csc_to_csr(Lower_A_CSC);
#endif

  nnz = Lower_A_CSC->nnz;

  bool isLfactor = false;

  double *y_serial, *y_correct = new double[n];
  double *y_correct_perfect = new double[n];
  std::cout << "Starting SpTrSv Runtime analysis" << std::endl;

  std::vector<std::string> Runtime_headers;
  Runtime_headers.emplace_back("Matrix_Name");
  Runtime_headers.emplace_back("Algorithm");
  Runtime_headers.emplace_back("Kernel");
  Runtime_headers.emplace_back("Core");
  Runtime_headers.emplace_back("Scheduling_Time");
  Runtime_headers.emplace_back("Executor_Runtime");
  Runtime_headers.emplace_back("nlevel");
  Runtime_headers.emplace_back("Profitable");

  std::string delimiter = "/";
  size_t matrix_name_start_pos = matrix_name.find(delimiter) + 1;
  delimiter = ".";
  size_t matrix_name_end_pos = matrix_name.find(delimiter);
  if (argc < 2) {
    matrix_name_start_pos = 0;
    matrix_name_end_pos = matrix_name.size();
  }
  std::string Mat_name(matrix_name.begin() + matrix_name_start_pos,
                       matrix_name.begin() + matrix_name_end_pos);
  std::string Data_name = "SpTrSv_Final_20";
  profiling_utils::CSVManager runtime_csv(Data_name, "some address",
                                          Runtime_headers, false);
#ifdef ANALYSE_STATISTIC
  std::vector<std::string> Parameter_headers;
  Parameter_headers.emplace_back("Matrix_Name");
  Parameter_headers.emplace_back("Algorithm");
  Parameter_headers.emplace_back("Kernel");
  Parameter_headers.emplace_back("Core");
  Parameter_headers.emplace_back("AvgParallelism");
  Parameter_headers.emplace_back("WeightedAvgParallelism");
  Parameter_headers.emplace_back("AvgGroupSize");

  Parameter_headers.emplace_back("AvgCostPerNNZSpTrSv");
  Parameter_headers.emplace_back("CostPerVertexSpTrSv");

  Parameter_headers.emplace_back("AvgCostPerNNZSpIC0");
  Parameter_headers.emplace_back("CostPerVertexSpIC0");

  Parameter_headers.emplace_back("CostPerVertexSpILU0");
  Parameter_headers.emplace_back("AvgCostPerNNZSpILU0");

  Parameter_headers.emplace_back("PG");

  profiling_utils::CSVManager parameter_csv("SpTrSv_Analysis", "some address",
                                            Parameter_headers, false);
#endif

  auto tmp = make_full(Lower_A_CSC);
  auto CSR_A = csc_to_csr(tmp);

  std::vector<int> Cores;
  if(nthreads == -1){
    Cores.emplace_back(20);
  } else {
    Cores.emplace_back(nthreads);
  }

  std::vector<bool> choice{true};
  std::vector<bool> bin_pack{true};

  //"********************* LL Serial  *********************"
  Sptrsv_LL_Serial LL_serial(Lower_A_CSR, Lower_A_CSC, NULLPNTR,
                             "LL serial"); // seq
  auto LL_serial_runtime = LL_serial.evaluate();
  y_serial = LL_serial.solution();
  copy_vector(0, n, y_serial, y_correct);
  std::cout << "Running LL Serial Code - The runtime:"
            << LL_serial_runtime.elapsed_time << std::endl;
  runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
  runtime_csv.addElementToRecord("Serial", "Algorithm");
  runtime_csv.addElementToRecord("LL", "Kernel");
  runtime_csv.addElementToRecord(1, "Core");
  runtime_csv.addElementToRecord(0, "Scheduling_Time");
  runtime_csv.addElementToRecord(LL_serial_runtime.elapsed_time,
                                 "Executor_Runtime");
  runtime_csv.addElementToRecord(1, "nlevel");
  runtime_csv.addElementToRecord(0, "Profitable");
  runtime_csv.addRecord();

  //"********************* LL Levelset *********************"
  for (auto &core : Cores) {
    SpTrsv_LL_Wavefront LL_lvl_obj(Lower_A_CSR, Lower_A_CSC, y_correct,
                                  "LL Levelset ", core);
    timing_measurement LL_lvl_runtime;
    omp_set_num_threads(core);
    LL_lvl_runtime = LL_lvl_obj.evaluate();
    std::cout << "Running LL Levelset Code with #core: " << core
              << " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
    int part_no, nlevels;
    LL_lvl_obj.getWaveStatistic(nlevels, part_no);
    runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
    runtime_csv.addElementToRecord("Wavefront", "Algorithm");
    runtime_csv.addElementToRecord("LL", "Kernel");
    runtime_csv.addElementToRecord(core, "Core");
    runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(),
                                   "Scheduling_Time");
    runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time,
                                   "Executor_Runtime");
    runtime_csv.addElementToRecord(nlevels, "nlevel");
    double profitable =
        (LL_lvl_obj.getSchedulingTime()) /
        (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
    runtime_csv.addElementToRecord(profitable, "Profitable");

    runtime_csv.addRecord();
  }
#ifdef SPMP
  //"********************* LL Levelset *********************"
  for (auto &core : Cores) {
    SpTrsv_LL_Parallel_Wavefront LL_lvl_obj(CSR_A, Lower_A_CSR, Lower_A_CSC,
                                           y_correct, "LL Levelset ", core);
    timing_measurement LL_lvl_runtime;
    omp_set_num_threads(core);
    LL_lvl_runtime = LL_lvl_obj.evaluate();
    std::cout << "Running LL parallel Levelset Code with #core: " << core
              << " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
    int part_no, nlevels;
    LL_lvl_obj.getWaveStatistic(nlevels, part_no);
    runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
    runtime_csv.addElementToRecord("Parallel_Wavefront", "Algorithm");
    runtime_csv.addElementToRecord("LL", "Kernel");
    runtime_csv.addElementToRecord(core, "Core");
    runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(),
                                   "Scheduling_Time");
    runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time,
                                   "Executor_Runtime");
    runtime_csv.addElementToRecord(nlevels, "nlevel");
    double profitable =
        (LL_lvl_obj.getSchedulingTime()) /
        (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
    runtime_csv.addElementToRecord(profitable, "Profitable");
    runtime_csv.addRecord();
  }
#endif
  //"********************* LL Tree + Levelset *********************"
  for (auto &core : Cores) {
    SpTrSv_LL_Tree_Wavefront_No_unPACK LL_lvl_obj(
        Lower_A_CSR, Lower_A_CSC, y_correct, "LL Tree Levelset ", core);
    timing_measurement LL_lvl_runtime;
    omp_set_num_threads(core);
    LL_lvl_runtime = LL_lvl_obj.evaluate();
    std::cout << "Running LL Levelset + Tree Code with #core: " << core
              << " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
    int part_no, nlevels;
    LL_lvl_obj.getWaveStatistic(nlevels, part_no);
    runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
    runtime_csv.addElementToRecord("Tree_Wavefront", "Algorithm");
    runtime_csv.addElementToRecord("LL", "Kernel");
    runtime_csv.addElementToRecord(core, "Core");
    runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(),
                                   "Scheduling_Time");
    runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time,
                                   "Executor_Runtime");
    runtime_csv.addElementToRecord(nlevels, "nlevel");
    double profitable =
        (LL_lvl_obj.getSchedulingTime()) /
        (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
    runtime_csv.addElementToRecord(profitable, "Profitable");
    runtime_csv.addRecord();
  }
  //"********************* LL HDAGG *********************"
  for (auto &core : Cores) {
    for (auto &&bin : bin_pack) {
      SpTrSv_LL_HDAGG LL_lvl_obj(Lower_A_CSR, Lower_A_CSC, y_correct,
                                 "LL HDAGG ", core, bin);
      timing_measurement LL_lvl_runtime;
      omp_set_num_threads(core);
      LL_lvl_runtime = LL_lvl_obj.evaluate();
      if (bin) {
        std::cout << "Running LL HDAGG with BIN Code with #core: " << core
                  << " - The runtime:" << LL_lvl_runtime.elapsed_time
                  << std::endl;
      } else {
        std::cout << "Running LL HDAGG Code with #core: " << core
                  << " - The runtime:" << LL_lvl_runtime.elapsed_time
                  << std::endl;
      }
      int part_no, nlevels;
      LL_lvl_obj.getWaveStatistic(nlevels, part_no);
      double avg_par = part_no * 1.0 / nlevels;
      runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
      if (bin) {
        runtime_csv.addElementToRecord("HDAGG_BIN", "Algorithm");
      } else {
        runtime_csv.addElementToRecord("HDAGG_NO_BIN", "Algorithm");
      }
      runtime_csv.addElementToRecord("LL", "Kernel");
      runtime_csv.addElementToRecord(core, "Core");
      runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(),
                                     "Scheduling_Time");
      runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time,
                                     "Executor_Runtime");

      runtime_csv.addElementToRecord(nlevels, "nlevel");
      double profitable =
          (LL_lvl_obj.getSchedulingTime()) /
          (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
      runtime_csv.addElementToRecord(profitable, "Profitable");

      runtime_csv.addRecord();
    }
  }
  //"********************* LL Tree HDAGG *********************"
  for (auto &core : Cores) {
    for (auto &&bin : bin_pack) {
      SpTrSv_LL_Tree_HDAGG LL_lvl_obj(Lower_A_CSR, Lower_A_CSC, y_correct,
                                      "LL Tree HDAGG ", core, isLfactor, bin);
      timing_measurement LL_lvl_runtime;
      omp_set_num_threads(core);
      LL_lvl_runtime = LL_lvl_obj.evaluate();
      if (bin) {
        std::cout << "Running LL Tree HDAGG with BIN Code with #core: " << core
                  << " - The runtime:" << LL_lvl_runtime.elapsed_time
                  << std::endl;
      } else {
        std::cout << "Running LL Tree HDAGG Code with #core: " << core
                  << " - The runtime:" << LL_lvl_runtime.elapsed_time
                  << std::endl;
      }
      int part_no, nlevels;
      LL_lvl_obj.getWaveStatistic(nlevels, part_no);
      runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
      if (bin) {
        runtime_csv.addElementToRecord("Tree_HDAGG_BIN", "Algorithm");
      } else {
        runtime_csv.addElementToRecord("Tree_HDAGG_NO_BIN", "Algorithm");
      }
      runtime_csv.addElementToRecord("LL", "Kernel");
      runtime_csv.addElementToRecord(core, "Core");
      runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(),
                                     "Scheduling_Time");
      runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time,
                                     "Executor_Runtime");

      runtime_csv.addElementToRecord(nlevels, "nlevel");
      double profitable =
          (LL_lvl_obj.getSchedulingTime()) /
          (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
      runtime_csv.addElementToRecord(profitable, "Profitable");
      runtime_csv.addRecord();
    }
  }
  //"********************* LL Tree BFS HDAGG *********************"
  for (auto &core : Cores) {
    for (auto &&bin : {true}) {
      SpTrSv_LL_Tree_HDAGG_BFS LL_lvl_obj(Lower_A_CSR, Lower_A_CSC, y_correct,
                                          "LL Tree BFS HDAGG ", core, isLfactor,
                                          bin);
      timing_measurement LL_lvl_runtime;
      omp_set_num_threads(core);
      LL_lvl_runtime = LL_lvl_obj.evaluate();
      if (bin) {
        std::cout << "Running LL Tree BFS HDAGG with BIN Code with #core: "
                  << core << " - The runtime:" << LL_lvl_runtime.elapsed_time
                  << std::endl;
      } else {
        std::cout << "Running LL Tree BFS HDAGG Code with #core: " << core
                  << " - The runtime:" << LL_lvl_runtime.elapsed_time
                  << std::endl;
      }
      int part_no, nlevels;
      LL_lvl_obj.getWaveStatistic(nlevels, part_no);

      runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
      if (bin) {
        runtime_csv.addElementToRecord("Tree_HDAGG_BFS_BIN", "Algorithm");
      } else {
        runtime_csv.addElementToRecord("Tree_HDAGG_BFS_NO_BIN", "Algorithm");
      }
      runtime_csv.addElementToRecord("LL", "Kernel");
      runtime_csv.addElementToRecord(core, "Core");
      runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(),
                                     "Scheduling_Time");
      runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time,
                                     "Executor_Runtime");

      runtime_csv.addElementToRecord(nlevels, "nlevel");
      double profitable =
          (LL_lvl_obj.getSchedulingTime()) /
          (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
      runtime_csv.addElementToRecord(profitable, "Profitable");
      runtime_csv.addRecord();
    }
  }

#ifdef SPMP
  //"********************* LL Tree BFS HDAGG P2P *********************"
  for (auto &core : Cores) {
    for (auto &&bin : {false}) {
      omp_set_num_threads(core);
      SpTrSv_LL_Tree_HDAGG_BFS_P2P LL_lvl_obj(
          Lower_A_CSR, Lower_A_CSC, y_correct, "LL Tree BFS HDAGG P2P ", core,
          isLfactor, bin);
      timing_measurement LL_lvl_runtime;
      LL_lvl_runtime = LL_lvl_obj.evaluate();
      if (bin) {
        std::cout << "Running LL Tree BFS HDAGG P2P with BIN Code with #core: "
                  << core << " - The runtime:" << LL_lvl_runtime.elapsed_time
                  << std::endl;
      } else {
        std::cout << "Running LL Tree BFS HDAGG P2P Code with #core: " << core
                  << " - The runtime:" << LL_lvl_runtime.elapsed_time
                  << std::endl;
      }
      int part_no, nlevels;
      LL_lvl_obj.getWaveStatistic(nlevels, part_no);

      runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
      if (bin) {
        runtime_csv.addElementToRecord("Tree_HDAGG_BFS_P2P_BIN", "Algorithm");
      } else {
        runtime_csv.addElementToRecord("Tree_HDAGG_BFS_P2P_NO_BIN", "Algorithm");
      }
      runtime_csv.addElementToRecord("LL", "Kernel");
      runtime_csv.addElementToRecord(core, "Core");
      runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(),
                                     "Scheduling_Time");
      runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time,
                                     "Executor_Runtime");
      runtime_csv.addElementToRecord(nlevels, "nlevel");
      double profitable =
          (LL_lvl_obj.getSchedulingTime()) /
          (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
      runtime_csv.addElementToRecord(profitable, "Profitable");
      runtime_csv.addRecord();
    }
  }
#endif

  //"********************* LL LBC *********************"
  std::vector<int> P3 = {5000};
  for (auto core : Cores) {
    for (auto p3 : P3) {
      SpTrSv_LL_LBC LL_lvl_obj(Lower_A_CSR, Lower_A_CSC, y_correct,
                               "LL LBC_Tree", core);
      omp_set_num_threads(core);
      LL_lvl_obj.setP2_P3(-1, p3);
      auto LL_lvl_runtime = LL_lvl_obj.evaluate();
      std::cout << "Running LL LBC_Tree Code with #core: " << core
                << " and P3=" << p3
                << " - The runtime:" << LL_lvl_runtime.elapsed_time
                << std::endl;
      int part_no, nlevels;
      LL_lvl_obj.getWaveStatistic(nlevels, part_no);

      runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
      runtime_csv.addElementToRecord("LBC", "Algorithm");
      runtime_csv.addElementToRecord("LL", "Kernel");
      runtime_csv.addElementToRecord(core, "Core");
      runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(),
                                     "Scheduling_Time");
      runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time,
                                     "Executor_Runtime");
      runtime_csv.addElementToRecord(nlevels, "nlevel");
      double profitable =
          (LL_lvl_obj.getSchedulingTime()) /
          (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
      runtime_csv.addElementToRecord(profitable, "Profitable");
      runtime_csv.addRecord();
    }
  }

#ifdef MKL
  //"********************* MKL *********************"
  for (auto &core : Cores) {
    omp_set_num_threads(core);
    SpTrSv_LL_MKL LL_lvl_obj(Lower_A_CSR, Lower_A_CSC, y_correct, "LL MKL",
                             core);
    timing_measurement LL_lvl_runtime;
    LL_lvl_runtime = LL_lvl_obj.evaluate();
    std::cout << "Running LL MKL Code with #core: " << core
              << " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
    runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
    runtime_csv.addElementToRecord("MKL", "Algorithm");
    runtime_csv.addElementToRecord("LL", "Kernel");
    runtime_csv.addElementToRecord(core, "Core");
    runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(),
                                   "Scheduling_Time");
    runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time,
                                   "Executor_Runtime");
    double profitable =
        (LL_lvl_obj.getSchedulingTime()) /
        (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
    runtime_csv.addElementToRecord(profitable, "Profitable");
    runtime_csv.addElementToRecord(1, "nlevel");
    runtime_csv.addRecord();
  }
#endif

#ifdef DAGP
  //"********************* DAGP *********************"
  for (auto &core : Cores) {
    std::vector<int> Parts{1000};
    for (auto &part : Parts) {
      SpTrSv_LL_DAGP LL_lvl_obj(Lower_A_CSR, Lower_A_CSC, y_correct, "LL DAGP",
                                part, core);
      omp_set_num_threads(core);
      timing_measurement LL_lvl_runtime;
      LL_lvl_runtime = LL_lvl_obj.evaluate();
      std::cout << "Running LL DAGP Code with #core: " << core
                << " - The runtime:" << LL_lvl_runtime.elapsed_time
                << std::endl;
      int part_no, nlevels;
      LL_lvl_obj.getWaveStatistic(nlevels, part_no);

      runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
      runtime_csv.addElementToRecord("DAGP", "Algorithm");
      runtime_csv.addElementToRecord("LL", "Kernel");
      runtime_csv.addElementToRecord(core, "Core");
      runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(),
                                     "Scheduling_Time");
      runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time,
                                     "Executor_Runtime");
      runtime_csv.addElementToRecord(nlevels, "nlevel");

      runtime_csv.addRecord();
    }
  }
#endif
#ifdef SPMP
  //"********************* SpMP *********************"
  for (auto &core : Cores) {
    omp_set_num_threads(core);
    SpTrsv_LL_SpMP LL_lvl_obj(CSR_A, Lower_A_CSR, Lower_A_CSC, y_correct,
                              "LL SpMP", core); // seq
    auto LL_lvl_runtime = LL_lvl_obj.evaluate();
    LL_lvl_runtime = LL_lvl_obj.evaluate();
    std::cout << "Running LL SpMP Code with #core: " << core
              << " - The runtime:" << LL_lvl_runtime.elapsed_time << std::endl;
    runtime_csv.addElementToRecord(Mat_name, "Matrix_Name");
    runtime_csv.addElementToRecord("SpMP", "Algorithm");
    runtime_csv.addElementToRecord("LL", "Kernel");
    runtime_csv.addElementToRecord(core, "Core");
    runtime_csv.addElementToRecord(LL_lvl_obj.getSchedulingTime(),
                                   "Scheduling_Time");
    runtime_csv.addElementToRecord(LL_lvl_runtime.elapsed_time,
                                   "Executor_Runtime");
    runtime_csv.addElementToRecord(LL_lvl_obj.getNumberOfP2P(), "nlevel");
    double profitable =
        (LL_lvl_obj.getSchedulingTime()) /
        (LL_serial_runtime.elapsed_time - LL_lvl_runtime.elapsed_time);
    runtime_csv.addElementToRecord(profitable, "Profitable");
    runtime_csv.addRecord();
  }
#endif

  delete tmp;
  delete CSR_A;
  delete[] y_correct;
  delete Lower_A_CSR;
  delete Lower_A_CSC;

  delete[] y_correct_perfect;
  return 0;
}
