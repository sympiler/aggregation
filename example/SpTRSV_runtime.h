//
// Created by behrooz on 1/17/21.
//

#ifndef SPTRSV_FINAL_H
#define SPTRSV_FINAL_H


//STD libraries
#include <algorithm>
#include <cmath>
#include <list>
#include <numeric>
#include <vector>


//Lib Codes
//#include "aggregation/Group.h"
//#include "aggregation/Group_Utils.h"
#include "aggregation/hdagg.h"
#include "aggregation/lbc.h"
#include "aggregation/sparse_inspector.h"
#include <sptrsv.h>
#include "aggregation/BCSCMatrix.h"
#include "aggregation/FusionDemo.h"
#include "aggregation/SuperNodalTools.h"

//Related works
#ifdef MKL
  #include <mkl.h>
#endif

#ifdef DAGP

#include "dagp_utils.h"
#include <dgraphReader.h>
#include <info.h>
#include <rvcycle.h>

#endif

#ifdef SPMP

#include "SpMP_Interface.h"

#endif

namespace sym_lib {

//============================== SpTrSv Serial LL ============================
class Sptrsv_LL_Serial : public FusionDemo {
protected:
  timing_measurement fused_code() override {
    timing_measurement t1;
    t1.start_timer();
    sptrsv_csr(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_);
    t1.measure_elapsed_time();
    copy_vector(0, n_, x_in_, x_);
    return t1;
  }

public:
  Sptrsv_LL_Serial(CSR *L, CSC *L_csc, double *correct_x, std::string name)
      : FusionDemo(L->n, L->nnz, name) {
    L1_csr_ = L;
    L1_csc_ = L_csc;
    correct_x_ = correct_x;
  };

  ~Sptrsv_LL_Serial() = default;
};

///======================================================================
///============================== Wavefront =============================
///======================================================================
//========================= SpTrSv LevelSet LL =================================
class SpTrsv_LL_Wavefront : public Sptrsv_LL_Serial {
protected:
  int nlevels, nthreads;
  std::vector<int> level_ptr;
  std::vector<int> level_set;
  std::vector<int> DAG_ptr;
  std::vector<int> DAG_set;
  timing_measurement scheduling_time;
  void build_set() override {
    // Create the DAG for automation stuff
    // auto DAG_time = computingDAG_CSR(L1_csr_, DAG_ptr, DAG_set);
    scheduling_time.start_timer();
    // Create the Levelset
    auto levelset_time = computingLevelSet_CSC(n_, L1_csc_->p, L1_csc_->i,
                                               level_ptr, level_set, nlevels);
    scheduling_time.measure_elapsed_time();
  }

  timing_measurement fused_code() override {
    timing_measurement t1;
    t1.start_timer();
    sptrsv_csr_levelset(L1_csr_->n, L1_csr_->p, L1_csr_->i, L1_csr_->x, nlevels,
                        level_ptr.data(), level_set.data(), x_in_);
    t1.measure_elapsed_time();
    sym_lib::copy_vector(0, n_, x_in_, x_);
    return t1;
  }

public:
  /*
   * @brief Class constructor
   * @param L Sparse Matrix with CSR format
   * @param L_csc Sparse Matrix with CSC format
   * @param correct_x The right answer for x
   * @param name The name of the algorithm
   * @param nt number of threads
   */
  SpTrsv_LL_Wavefront(CSR *L, CSC *L_csc, double *correct_x, std::string name,
                     int nt)
      : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
    L1_csr_ = L;
    L1_csc_ = L_csc;
    correct_x_ = correct_x;
    this->nthreads = nt;
  };
  double getSchedulingTime() { return scheduling_time.elapsed_time; };
  void getWaveStatistic(int &nlevels, int &part_no) {
    part_no = n_;
    nlevels = this->nlevels;
  }
  ~SpTrsv_LL_Wavefront() = default;
};
//====================== SpTrSv Parallel LevelSet LL ===========================
#ifdef SPMP
class SpTrsv_LL_Parallel_Wavefront : public Sptrsv_LL_Serial {
protected:
  int nlevels, nthreads;
  std::vector<int> DAG_ptr;
  std::vector<int> DAG_set;
  std::vector<int> level_ptr, level_set;
  timing_measurement scheduling_time;
  CSR *LBC_A;
  void build_set() override {
    // Convert lbc CSR version to SpMP CSR version
    SpMP::CSR *A = new SpMP::CSR();
    Convert_LBCCSR_to_SpMP(LBC_A, A);
    scheduling_time.start_timer();
    nlevels =
        HDAGG::levelsetCSRParallel_SpMP(A, nthreads, level_ptr, level_set);

    scheduling_time.measure_elapsed_time();
    delete A;
  }

  timing_measurement fused_code() override {
    timing_measurement t1;
    t1.start_timer();
    sptrsv_csr_levelset(L1_csr_->n, L1_csr_->p, L1_csr_->i, L1_csr_->x, nlevels,
                        level_ptr.data(), level_set.data(), x_in_);
    t1.measure_elapsed_time();
    sym_lib::copy_vector(0, n_, x_in_, x_);
    return t1;
  }

public:
  /*
   * @brief Class constructor
   * @param L Sparse Matrix with CSR format
   * @param L_csc Sparse Matrix with CSC format
   * @param correct_x The right answer for x
   * @param name The name of the algorithm
   * @param nt number of threads
   */
  SpTrsv_LL_Parallel_Wavefront(CSR *A, CSR *L, CSC *L_csc, double *correct_x,
                              std::string name, int nt)
      : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
    this->LBC_A = A;
    L1_csr_ = L;
    L1_csc_ = L_csc;
    correct_x_ = correct_x;
    this->nthreads = nt;
   #ifdef ENABLE_OPENMP
 omp_set_num_threads(nt);
#endif
  };

  double getSchedulingTime() { return scheduling_time.elapsed_time; };
  void getWaveStatistic(int &nlevels, int &part_no) {
    part_no = n_;
    nlevels = this->nlevels;
  }

  ~SpTrsv_LL_Parallel_Wavefront(){};
};
#endif
//======================= SpTRSV Levelset + Tree LL ============================
class SpTrSv_LL_Tree_Wavefront_No_unPACK : public Sptrsv_LL_Serial {
protected:
  int nthreads;
  timing_measurement Scheduling_time;
  int ngroups;
  std::vector<int> DAG_ptr;
  std::vector<int> DAG_set;
  std::vector<int> group_set, group_ptr;
  std::vector<int> level_ptr, level_set;
  int nlevels;
  void build_set() override {
    Scheduling_time.start_timer();
    // Create the DAG
    HDAGG::partialSparsification(n_, nnz_, L1_csc_->p, L1_csc_->i, DAG_ptr,
                                 DAG_set);

    HDAGG::treeBasedGrouping(n_, DAG_ptr, DAG_set, ngroups, group_ptr,
                             group_set, false);
    std::vector<int> group_DAG_ptr, group_DAG_set;
    HDAGG::buildGroupDAG(n_, ngroups, group_ptr.data(), group_set.data(),
                         DAG_ptr.data(), DAG_set.data(), group_DAG_ptr,
                         group_DAG_set);

    computingLevelSet_CSC(ngroups, group_DAG_ptr.data(), group_DAG_set.data(),
                          level_ptr, level_set, nlevels);

    for (int g = 0; g < ngroups; g++) {
      std::sort(group_set.data() + group_ptr[g],
                group_set.data() + group_ptr[g + 1]);
    }
    Scheduling_time.measure_elapsed_time();
  }

  timing_measurement fused_code() override {
    timing_measurement t1;
    t1.start_timer();
    sptrsv_csr_group_levelset(L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                              nlevels, level_ptr.data(), level_set.data(),
                              group_ptr.data(), group_set.data());
    t1.measure_elapsed_time();
    sym_lib::copy_vector(0, n_, x_in_, x_);
    return t1;
  }

public:
  /*
   * @brief Class constructor
   * @param L Sparse Matrix with CSR format
   * @param L_csc Sparse Matrix with CSC format
   * @param correct_x The right answer for x
   * @param name The name of the algorithm
   * @param nt number of threads
   */
  SpTrSv_LL_Tree_Wavefront_No_unPACK(CSR *L, CSC *L_csc, double *correct_x,
                                    std::string name, int nt)
      : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
    L1_csr_ = L;
    L1_csc_ = L_csc;
    correct_x_ = correct_x;
    nthreads = nt;
  };

  int getNlevels() { return nlevels; }
  void getWaveStatistic(int &nlevels, int &part_no) {
    part_no = ngroups;
    nlevels = this->nlevels;
  }
  double getSchedulingTime() { return Scheduling_time.elapsed_time; }
  ~SpTrSv_LL_Tree_Wavefront_No_unPACK() = default;
};

///======================================================================
///================================= HDAGG ==============================
///======================================================================
//================================ SpTRSV HDAGG LL ================================
class SpTrSv_LL_HDAGG : public Sptrsv_LL_Serial {
protected:
  int nthreads;
  timing_measurement Scheduling_time;
  int ngroups;
  std::vector<int> final_level_ptr, final_part_ptr, final_node_ptr;
  int final_level_no;
  std::vector<int> level_ptr, level_set;
  int nlevels;
  double avg_par, max_diff, var;
  int part_no;
  bool bin_pack;
  std::vector<double> cost;
  void build_set() override {
    Scheduling_time.start_timer();
    // Computing node Cost
    cost.resize(n_, 0);
    auto CSC_Lp = L1_csc_->p;
    auto CSC_Li = L1_csc_->i;
    auto CSR_Lp = L1_csr_->p;
    auto CSR_Li = L1_csr_->i;
    HDAGG::costComputation(n_, CSC_Lp, CSC_Li, CSR_Lp, CSR_Li, HDAGG::SpTrSv_LL,
                           NULL, NULL, false, cost);

    std::vector<int> DAG_ptr(L1_csc_->p, L1_csc_->p + n_ + 1);
    std::vector<int> DAG_set(L1_csc_->i, L1_csc_->i + nnz_);
    HDAGG::HDAGG(n_, nnz_, DAG_ptr, DAG_set, cost, nthreads, final_level_no,
                 final_level_ptr, final_part_ptr, final_node_ptr, false, false,
                 bin_pack);
    Scheduling_time.measure_elapsed_time();
  }

  timing_measurement fused_code() override {
    timing_measurement t1;
    t1.start_timer();
    sptrsv_csr_lbc(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                   final_level_no, final_level_ptr.data(),
                   final_part_ptr.data(), final_node_ptr.data());
    t1.measure_elapsed_time();
    sym_lib::copy_vector(0, n_, x_in_, x_);
    return t1;
  }

public:
  /*
   * @brief Class constructor
   * @param L Sparse Matrix with CSR format
   * @param L_csc Sparse Matrix with CSC format
   * @param correct_x The right answer for x
   * @param name The name of the algorithm
   * @param nt number of threads
   */
  SpTrSv_LL_HDAGG(CSR *L, CSC *L_csc, double *correct_x, std::string name,
                  int nt, bool bin_pack)
      : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
    L1_csr_ = L;
    L1_csc_ = L_csc;
    correct_x_ = correct_x;
    nthreads = nt;
    this->bin_pack = bin_pack;
  };
  // I'm not in the mood to delete codes
  void computeCharacteristics(double &first, double &second, double &third,
                              double &avg) {

    std::vector<double> max_array(final_level_no, 0);
    std::vector<double> whole_array(final_level_no, 0);
    for (int lvl = 0; lvl < final_level_no; lvl++) {
      double total_nodes = 0;
      double max_size = 0;
      for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1];
           ++w_ptr) {
        int cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
        if (cc_size > max_size) {
          max_size = cc_size;
        }
        total_nodes += cc_size;
      }
      if (lvl == 0) {
        first = max_size / total_nodes;
      }
      if (lvl == 1) {
        second = max_size / total_nodes;
      }
      if (lvl == 2) {
        third = max_size / total_nodes;
      }
      max_array[lvl] = max_size;
      whole_array[lvl] = total_nodes;
    }

    double max_cc_sizes = 0;
    for (auto &iter : max_array) {
      max_cc_sizes += iter;
    }
    avg = max_cc_sizes / n_;

#ifndef NDEBUG
    double sum = 0;
    for (auto &iter : whole_array) {
      sum += iter;
    }
    assert(sum == n_);
#endif
  }

  void computeStatistic() {
    double part_no = 0;
    for (int i = 0; i < final_part_ptr.size() - 1; i++) {
      if (final_part_ptr[i + 1] != final_part_ptr[i]) {
        part_no++;
      }
    }
    this->part_no = part_no;
    avg_par = part_no / final_level_no;
    std::vector<double> max_part(final_level_no, 0);
    std::vector<double> min_part(final_level_no, 0);
    std::vector<double> level_weight(final_level_no, 0);
    for (int lvl = 0; lvl < final_level_no; lvl++) {
      double total_nodes = 0;
      double max_size = 0;
      double min_size = n_;
      for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1];
           ++w_ptr) {
        double cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
        double cc_cost = 0;
        for (int i = final_part_ptr[w_ptr]; i < final_part_ptr[w_ptr + 1];
             i++) {
          int node = final_node_ptr[i];
          cc_cost += cost[node];
        }
        cc_size = cc_cost;
        if (cc_size > max_size) {
          max_size = cc_size;
        }
        if (cc_size < min_size) {
          min_size = cc_size;
        }
        total_nodes += cc_size;
      }
      max_part[lvl] = max_size;
      min_part[lvl] = min_size;
      level_weight[lvl] = total_nodes;
    }

    max_diff = 0;
    double total_weight = 0;
    for (int l = 0; l < final_level_no; l++) {
      max_diff += level_weight[l] * (max_part[l] - min_part[l]);
      total_weight += level_weight[l];
    }
    max_diff = max_diff / total_weight;
    var = 0;
  }
  void getWaveStatistic(int &nlevels, int &part_no) {
    part_no = this->part_no;
    nlevels = this->final_level_no;
  }
  double getSchedulingTime() { return Scheduling_time.elapsed_time; }
  void getStatistic(double &avg_par, double &max_diff, double &var) {
    avg_par = this->avg_par;
    max_diff = this->max_diff;
    var = this->var;
  }
  ~SpTrSv_LL_HDAGG() = default;
};
//================================ SpTRSV HDAGG LL ================================
class SpTrSv_LL_Tree_HDAGG : public Sptrsv_LL_Serial {
protected:
  int nthreads;
  timing_measurement Scheduling_time;
  int ngroups;
  std::vector<int> final_level_ptr, final_part_ptr, final_node_ptr;
  int final_level_no;
  std::vector<int> DAG_ptr;
  std::vector<int> DAG_set;
  std::vector<int> group_set, group_ptr;
  std::vector<int> level_ptr, level_set;
  int nlevels;
  bool isLfactor;
  double avg_par, max_diff, var;
  int part_no;
  bool bin_pack;
  std::vector<double> cost;
  void build_set() override {
    Scheduling_time.start_timer();
    // Create the DAG
    if (isLfactor) {
      std::cout << "Building Tree" << std::endl;
      HDAGG::buildETree(L1_csr_, nnz_, DAG_ptr, DAG_set);
    } else {
#ifndef NDEBUG
      std::cout << "Compute DAG" << std::endl;
#endif
      //                computingDAG_CSR(L1_csr_, L1_csc_->p, L1_csc_->i);
      HDAGG::partialSparsification(n_, nnz_, L1_csc_->p, L1_csc_->i, DAG_ptr,
                                   DAG_set);
    }
    HDAGG::treeBasedGrouping(n_, DAG_ptr, DAG_set, ngroups, group_ptr,
                             group_set, isLfactor);
    std::vector<int> group_DAG_ptr, group_DAG_set;
    HDAGG::buildGroupDAG(n_, ngroups, group_ptr.data(), group_set.data(),
                         DAG_ptr.data(), DAG_set.data(), group_DAG_ptr,
                         group_DAG_set);

    // Computing node Cost
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

    //                std::cout << "Static is COMPUTING" << std::endl;
    //                computeStatistic();

    HDAGG::ungroupingScheduleAndApplyOrdering(
        n_, final_level_no, final_level_ptr, final_part_ptr, final_node_ptr,
        group_ptr, group_set);
    Scheduling_time.measure_elapsed_time();
  }

  timing_measurement fused_code() override {
    timing_measurement t1;
    t1.start_timer();
    sptrsv_csr_lbc(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                   final_level_no, final_level_ptr.data(),
                   final_part_ptr.data(), final_node_ptr.data());
    t1.measure_elapsed_time();
    sym_lib::copy_vector(0, n_, x_in_, x_);
    return t1;
  }

public:
  /*
   * @brief Class constructor
   * @param L Sparse Matrix with CSR format
   * @param L_csc Sparse Matrix with CSC format
   * @param correct_x The right answer for x
   * @param name The name of the algorithm
   * @param nt number of threads
   */
  SpTrSv_LL_Tree_HDAGG(CSR *L, CSC *L_csc, double *correct_x, std::string name,
                       int nt, bool isLfactor, bool bin_pack)
      : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
    L1_csr_ = L;
    L1_csc_ = L_csc;
    correct_x_ = correct_x;
    nthreads = nt;
    this->isLfactor = isLfactor;
    this->bin_pack = bin_pack;
  };
  // I'm not in the mood to delete codes
  void computeCharacteristics(double &first, double &second, double &third,
                              double &avg) {

    std::vector<double> max_array(final_level_no, 0);
    std::vector<double> whole_array(final_level_no, 0);
    for (int lvl = 0; lvl < final_level_no; lvl++) {
      double total_nodes = 0;
      double max_size = 0;
      for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1];
           ++w_ptr) {
        int cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
        if (cc_size > max_size) {
          max_size = cc_size;
        }
        total_nodes += cc_size;
      }
      if (lvl == 0) {
        first = max_size / total_nodes;
      }
      if (lvl == 1) {
        second = max_size / total_nodes;
      }
      if (lvl == 2) {
        third = max_size / total_nodes;
      }
      max_array[lvl] = max_size;
      whole_array[lvl] = total_nodes;
    }

    double max_cc_sizes = 0;
    for (auto &iter : max_array) {
      max_cc_sizes += iter;
    }
    avg = max_cc_sizes / n_;

#ifndef NDEBUG
    double sum = 0;
    for (auto &iter : whole_array) {
      sum += iter;
    }
    assert(sum == n_);
#endif
  }

  void computeStatistic() {
    double part_no = 0;
    for (int i = 0; i < final_part_ptr.size() - 1; i++) {
      if (final_part_ptr[i + 1] != final_part_ptr[i]) {
        part_no++;
      }
    }
    this->part_no = part_no;
    avg_par = part_no / final_level_no;
    std::vector<double> max_part(final_level_no, 0);
    std::vector<double> min_part(final_level_no, 0);
    std::vector<double> level_weight(final_level_no, 0);
    for (int lvl = 0; lvl < final_level_no; lvl++) {
      double total_nodes = 0;
      double max_size = 0;
      double min_size = n_;
      for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1];
           ++w_ptr) {
        double cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
        double cc_cost = 0;
        for (int i = final_part_ptr[w_ptr]; i < final_part_ptr[w_ptr + 1];
             i++) {
          int node = final_node_ptr[i];
          cc_cost += cost[node];
        }
        cc_size = cc_cost;
        if (cc_size > max_size) {
          max_size = cc_size;
        }
        if (cc_size < min_size) {
          min_size = cc_size;
        }
        total_nodes += cc_size;
      }
      max_part[lvl] = max_size;
      min_part[lvl] = min_size;
      level_weight[lvl] = total_nodes;
    }

    max_diff = 0;
    double total_weight = 0;
    for (int l = 0; l < final_level_no; l++) {
      max_diff += level_weight[l] * (max_part[l] - min_part[l]);
      total_weight += level_weight[l];
    }
    max_diff = max_diff / total_weight;
    var = 0;
  }
  void getWaveStatistic(int &nlevels, int &part_no) {
    part_no = this->part_no;
    nlevels = this->final_level_no;
  }
  double getSchedulingTime() { return Scheduling_time.elapsed_time; }
  void getStatistic(double &avg_par, double &max_diff, double &var) {
    avg_par = this->avg_par;
    max_diff = this->max_diff;
    var = this->var;
  }
  ~SpTrSv_LL_Tree_HDAGG() = default;
};
//================================ SpTRSV HDAGG LL ================================
class SpTrSv_LL_Tree_HDAGG_BFS : public Sptrsv_LL_Serial {
protected:
  int nthreads;
  timing_measurement Scheduling_time;
  int ngroups;
  std::vector<int> final_level_ptr, final_part_ptr, final_node_ptr;
  int final_level_no;
  std::vector<int> DAG_ptr;
  std::vector<int> DAG_set;
  std::vector<int> group_set, group_ptr;
  std::vector<int> level_ptr, level_set;
  std::vector<int> group_DAG_ptr, group_DAG_set;
  double avgGroupSize;
  int nlevels;
  double avg_par, max_diff, var;
  int part_no = 0;
  std::vector<double> cost;
  bool isLfactor;
  bool bin_pack;
  void build_set() override {
    Scheduling_time.start_timer();
    // Create the DAG
    if (isLfactor) {
      HDAGG::buildETree(L1_csr_, nnz_, DAG_ptr, DAG_set);
    } else {
      //                computingDAG_CSR(L1_csr_, DAG_ptr_not_prune,
      //                DAG_set_not_prune);
      HDAGG::partialSparsification(n_, nnz_, L1_csc_->p, L1_csc_->i, DAG_ptr,
                                   DAG_set);
    }

    HDAGG::treeBasedGroupingBFS(n_, DAG_ptr, DAG_set, ngroups, group_ptr,
                                group_set, isLfactor);
    avgGroupSize = 0;
    for (int i = 0; i < ngroups; i++) {
      avgGroupSize += (group_ptr[i + 1] - group_ptr[i]);
    }
    avgGroupSize = avgGroupSize / ngroups;
    HDAGG::buildGroupDAG(n_, ngroups, group_ptr.data(), group_set.data(),
                         DAG_ptr.data(), DAG_set.data(), group_DAG_ptr,
                         group_DAG_set);

    // Computing node Cost
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

    //      std::cout << "Static is COMPUTING" << std::endl;
    //      computeStatistic();

    HDAGG::ungroupingScheduleAndApplyOrdering(
        n_, final_level_no, final_level_ptr, final_part_ptr, final_node_ptr,
        group_ptr, group_set, DAG_ptr.data(), DAG_set.data(), true, false);

    Scheduling_time.measure_elapsed_time();
  }

  timing_measurement fused_code() override {
    timing_measurement t1;
    t1.start_timer();
    sptrsv_csr_lbc(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                   final_level_no, final_level_ptr.data(),
                   final_part_ptr.data(), final_node_ptr.data());
    t1.measure_elapsed_time();
    sym_lib::copy_vector(0, n_, x_in_, x_);
    return t1;
  }

public:
  /*
   * @brief Class constructor
   * @param L Sparse Matrix with CSR format
   * @param L_csc Sparse Matrix with CSC format
   * @param correct_x The right answer for x
   * @param name The name of the algorithm
   * @param nt number of threads
   */
  SpTrSv_LL_Tree_HDAGG_BFS(CSR *L, CSC *L_csc, double *correct_x,
                           std::string name, int nt, bool isLfactor,
                           bool bin_pack)
      : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
    L1_csr_ = L;
    L1_csc_ = L_csc;
    correct_x_ = correct_x;
    nthreads = nt;
    this->isLfactor = isLfactor;
    this->bin_pack = bin_pack;
  };

  void computeCharacteristics(double &first, double &second, double &third,
                              double &avg) {

    std::vector<double> max_array(final_level_no, 0);
    std::vector<double> whole_array(final_level_no, 0);
    for (int lvl = 0; lvl < final_level_no; lvl++) {
      double total_nodes = 0;
      double max_size = 0;
      for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1];
           ++w_ptr) {
        int cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
        if (cc_size > max_size) {
          max_size = cc_size;
        }
        total_nodes += cc_size;
      }
      if (lvl == 0) {
        first = max_size / total_nodes;
      }
      if (lvl == 1) {
        second = max_size / total_nodes;
      }
      if (lvl == 2) {
        third = max_size / total_nodes;
      }
      max_array[lvl] = max_size;
      whole_array[lvl] = total_nodes;
    }

    double max_cc_sizes = 0;
    for (auto &iter : max_array) {
      max_cc_sizes += iter;
    }
    avg = max_cc_sizes / n_;

#ifndef NDEBUG
    double sum = 0;
    for (auto &iter : whole_array) {
      sum += iter;
    }
    assert(sum == n_);
#endif
  }

  void computeStatistic() {
    double part_no = 0;
    for (int i = 0; i < final_part_ptr.size() - 1; i++) {
      if (final_part_ptr[i + 1] != final_part_ptr[i]) {
        part_no++;
      }
    }
    this->part_no = part_no;
    avg_par = part_no / final_level_no;
    std::vector<double> max_part(final_level_no, 0);
    std::vector<double> min_part(final_level_no, 0);
    std::vector<double> level_weight(final_level_no, 0);
    for (int lvl = 0; lvl < final_level_no; lvl++) {
      double total_nodes = 0;
      double max_size = 0;
      double min_size = n_;
      for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1];
           ++w_ptr) {
        double cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
        double cc_cost = 0;
        for (int i = final_part_ptr[w_ptr]; i < final_part_ptr[w_ptr + 1];
             i++) {
          int node = final_node_ptr[i];
          cc_cost += cost[node];
        }
        cc_size = cc_cost;
        if (cc_size > max_size) {
          max_size = cc_size;
        }
        if (cc_size < min_size) {
          min_size = cc_size;
        }
        total_nodes += cc_size;
      }
      max_part[lvl] = max_size;
      min_part[lvl] = min_size;
      level_weight[lvl] = total_nodes;
    }

    max_diff = 0;
    double total_weight = 0;
    for (int l = 0; l < final_level_no; l++) {
      max_diff += level_weight[l] * (max_part[l] - min_part[l]);
      total_weight += level_weight[l];
    }
    max_diff = max_diff / total_weight;
    var = 0;
  }
  double getSchedulingTime() { return Scheduling_time.elapsed_time; }
  void getStatistic(double &avg_par, double &max_diff, double &var) {
    avg_par = this->avg_par;
    max_diff = this->max_diff;
    var = this->var;
  }
  void getWaveStatistic(int &nlevels, int &part_no) {
    part_no = this->part_no;
    nlevels = this->final_level_no;
  }
  double getPotentialGain() {
    double critical_path = 0;
    std::vector<double> cost(n_, 0);
    std::vector<double> cost_per_level(final_level_no, 0);
    std::vector<double> part_per_level(final_level_no, 0);
    std::vector<double> node_per_level(final_level_no, 0);
    std::vector<double> max_cost(final_level_no, 0);
    auto CSC_Lp = L1_csc_->p;
    auto CSC_Li = L1_csc_->i;
    auto CSR_Lp = L1_csr_->p;
    auto CSR_Li = L1_csr_->i;
    HDAGG::costComputation(n_, CSC_Lp, CSC_Li, CSR_Lp, CSR_Li, HDAGG::SpTrSv_LL,
                           NULL, NULL, false, cost);
    double Total_cost = 0;
    double total_part_number = 0;

    for (int lvl = 0; lvl < final_level_no; lvl++) {
      for (int prt_ptr = final_level_ptr[lvl];
           prt_ptr < final_level_ptr[lvl + 1]; prt_ptr++) {
        if (final_part_ptr[prt_ptr + 1] != final_part_ptr[prt_ptr]) {
          total_part_number++;
          part_per_level[lvl]++;
        }
        double part_cost = 0;
        for (int node_ptr = final_part_ptr[prt_ptr];
             node_ptr < final_part_ptr[prt_ptr + 1]; node_ptr++) {
          int node = final_node_ptr[node_ptr];
          assert(!std::isnan(cost[node]));
          Total_cost += cost[node];
          cost_per_level[lvl] += cost[node];
          part_cost += cost[node];
          node_per_level[lvl]++;
        }
        if (part_cost > max_cost[lvl]) {
          max_cost[lvl] = part_cost;
        }
      }
      critical_path += max_cost[lvl];
    }

    return 1 - (Total_cost / nthreads) / critical_path;
  }

  void computeMetrics(double &avgParallelism, double &WeightedAvgParallelism,
                      double &CostAvgParallelism, double &AvgGroupSize,
                      double &AvgCostPerVertexSpTrSv,
                      double &AvgCostPerNNZSpTrSv,
                      double &AvgCostPerVertexSpILU0,
                      double &AvgCostPerNNZSpLU0, double &AvgCostPerVertexSpIC0,
                      double &AvgCostPerNNZSpIC0, double &PG) {
    double Total_cost = 0;
    double total_part_number = 0;
    std::vector<double> cost_per_level(final_level_no, 0);
    std::vector<double> part_per_level(final_level_no, 0);
    std::vector<double> node_per_level(final_level_no, 0);
    std::vector<double> max_cost(final_level_no, 0);
    double critical_path = 0;
    std::vector<double> cost;
    cost.resize(n_, 0);
    auto CSC_Lp = L1_csc_->p;
    auto CSC_Li = L1_csc_->i;
    auto CSR_Lp = L1_csr_->p;
    auto CSR_Li = L1_csr_->i;
    HDAGG::costComputation(n_, CSC_Lp, CSC_Li, CSR_Lp, CSR_Li, HDAGG::SpTrSv_LL,
                           NULL, NULL, false, cost);
    for (int lvl = 0; lvl < final_level_no; lvl++) {
      for (int prt_ptr = final_level_ptr[lvl];
           prt_ptr < final_level_ptr[lvl + 1]; prt_ptr++) {
        if (final_part_ptr[prt_ptr + 1] != final_part_ptr[prt_ptr]) {
          total_part_number++;
          part_per_level[lvl]++;
        }
        double part_cost = 0;
        for (int node_ptr = final_part_ptr[prt_ptr];
             node_ptr < final_part_ptr[prt_ptr + 1]; node_ptr++) {
          int node = final_node_ptr[node_ptr];
          assert(!std::isnan(cost[node]));
          Total_cost += cost[node];
          cost_per_level[lvl] += cost[node];
          part_cost += cost[node];
          node_per_level[lvl]++;
        }
        if (part_cost > max_cost[lvl]) {
          max_cost[lvl] = part_cost;
        }
      }
      critical_path += max_cost[lvl];
    }

    PG = 1 - (Total_cost / nthreads) / critical_path;

    avgParallelism = total_part_number / final_level_no;
    for (int lvl = 0; lvl < final_level_no; lvl++) {
      CostAvgParallelism += (part_per_level[lvl] * cost_per_level[lvl]);
      WeightedAvgParallelism += (part_per_level[lvl] * node_per_level[lvl]);
    }

    CostAvgParallelism = CostAvgParallelism / Total_cost;
    WeightedAvgParallelism = WeightedAvgParallelism / n_;
    AvgGroupSize = this->avgGroupSize;

    AvgCostPerVertexSpTrSv = Total_cost / n_;
    AvgCostPerNNZSpTrSv = Total_cost / nnz_;

    cost.clear();
    cost.resize(n_, 0);
    HDAGG::costComputation(n_, CSC_Lp, CSC_Li, CSR_Lp, CSR_Li, HDAGG::SpICh0_UL,
                           NULL, NULL, false, cost);

    Total_cost = 0;
    for (auto &iter : cost) {
      Total_cost += iter;
    }
    AvgCostPerVertexSpIC0 = Total_cost / n_;
    AvgCostPerNNZSpIC0 = Total_cost / nnz_;

    auto tmp = make_full(L1_csc_);
    auto CSR_A = csc_to_csr(tmp);
    delete tmp;
    cost.clear();
    cost.resize(n_, 0);
    HDAGG::costComputation(n_, NULL, NULL, CSR_A->p, CSR_A->i, HDAGG::SpICh0_UL,
                           NULL, NULL, false, cost);

    Total_cost = 0;
    for (auto &iter : cost) {
      Total_cost += iter;
    }
    AvgCostPerVertexSpILU0 = Total_cost / n_;
    AvgCostPerNNZSpLU0 = Total_cost / nnz_;
    delete CSR_A;
  }
  ~SpTrSv_LL_Tree_HDAGG_BFS() = default;
};

///======================================================================
///=========================== RELATED WORKS ============================
///======================================================================

#ifdef SPMP
//================================ SpTRSV Tree BFS HDAGG LL ================================
class SpTrSv_LL_Tree_HDAGG_BFS_P2P : public Sptrsv_LL_Serial {
protected:
  int nthreads;
  timing_measurement Scheduling_time;
  int ngroups;
  std::vector<int> final_level_ptr, final_part_ptr, final_node_ptr;
  int final_level_no;
  std::vector<int> DAG_ptr;
  std::vector<int> DAG_set;
  std::vector<int> group_set, group_ptr;
  std::vector<int> level_ptr, level_set;
  std::vector<int> group_DAG_ptr, group_DAG_set;
  int nlevels;
  double avg_par, max_diff, var;
  int part_no = 0;
  std::vector<double> cost;
  bool isLfactor;
  bool bin_pack;
  SpMP::LevelSchedule *p2pScheduleWithTransitiveReduction;
  const int *invPerm;
  bool isfirst_time;

  void build_set() override {
    Scheduling_time.start_timer();

    if (isLfactor) {
      std::cout << "Building Tree" << std::endl;
      HDAGG::buildETree(L1_csr_, nnz_, DAG_ptr, DAG_set);
    } else {
      //                computingDAG_CSR(L1_csr_, DAG_ptr_not_prune,
      //                DAG_set_not_prune);
      HDAGG::partialSparsification(n_, nnz_, L1_csc_->p, L1_csc_->i, DAG_ptr,
                                   DAG_set);
    }

    HDAGG::treeBasedGroupingBFS(n_, DAG_ptr, DAG_set, ngroups, group_ptr,
                                group_set, isLfactor);

    HDAGG::buildGroupDAG(n_, ngroups, group_ptr.data(), group_set.data(),
                         DAG_ptr.data(), DAG_set.data(), group_DAG_ptr,
                         group_DAG_set);

    // Computing node Cost
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

    //      std::cout << "Static is COMPUTING" << std::endl;
    //      computeStatistic();

    HDAGG::ungroupingScheduleAndApplyOrdering(
        n_, final_level_no, final_level_ptr, final_part_ptr, final_node_ptr,
        group_ptr, group_set, DAG_ptr.data(), DAG_set.data(), true, false);

    // Create Group set
    HDAGG::getFinalScheduleDAG(n_, final_level_no, final_level_ptr.data(),
                               final_part_ptr.data(), final_node_ptr.data(),
                               L1_csc_->p, L1_csc_->i, ngroups, group_ptr,
                               group_set, DAG_ptr, DAG_set);
    if (!isfirst_time) {
      delete p2pScheduleWithTransitiveReduction;
    }

    // Create DAG
    CSC *CSC_DAG =
        new CSC(DAG_ptr.size() - 1, DAG_ptr.size() - 1, DAG_set.size());
    std::copy(DAG_ptr.begin(), DAG_ptr.end(), CSC_DAG->p);
    std::copy(DAG_set.begin(), DAG_set.end(), CSC_DAG->i);
    std::fill(CSC_DAG->x, CSC_DAG->x + DAG_set.size(), 1);
    CSC_DAG->stype = -1;
    auto tmp = make_full(CSC_DAG);
    auto CSR_A = csc_to_csr(tmp);
    delete tmp;
    delete CSC_DAG;
    // Create SpMP schedule
    // Convert lbc CSR version to SpMP CSR version
    auto *A = new SpMP::CSR();
    Convert_LBCCSR_to_SpMP(CSR_A, A);

    /////////////////////////////////////////////////////////////////////////////
    // Construct schedules
    /////////////////////////////////////////////////////////////////////////////
    bool hasZeroDiag = A->hasZeroDiag();
    if (hasZeroDiag) {
      fprintf(stderr,
              "Matrix has a zero diagonal element. Can't do Gauss Seidel\n");
    }

    int *symRowPtr = NULL, *symColIdx = NULL, *symDiagPtr = NULL,
        *symExtPtr = NULL;
    bool wasSymmetric = getSymmetricNnzPattern(A, &symRowPtr, &symDiagPtr,
                                               &symExtPtr, &symColIdx);

    p2pScheduleWithTransitiveReduction = new SpMP::LevelSchedule;
    if (wasSymmetric) {
      p2pScheduleWithTransitiveReduction->constructTaskGraph(*A);
    } else {
      p2pScheduleWithTransitiveReduction->constructTaskGraph(
          A->m, symRowPtr, symDiagPtr, symExtPtr, symColIdx,
          SpMP::PrefixSumCostFunction(symRowPtr));

      FREE(symRowPtr);
      FREE(symColIdx);
      FREE(symDiagPtr);
      FREE(symExtPtr);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Reorder matrix
    /////////////////////////////////////////////////////////////////////////////
    const int *perm = p2pScheduleWithTransitiveReduction->origToThreadContPerm;
    const int *invPerm =
        p2pScheduleWithTransitiveReduction->threadContToOrigPerm;
    assert(SpMP::isPerm(perm, A->m));
    assert(SpMP::isPerm(invPerm, A->m));
    this->invPerm = invPerm;

    delete A;
    delete CSR_A;
    isfirst_time = false;

    Scheduling_time.measure_elapsed_time();
  }

  timing_measurement fused_code() override {
    timing_measurement t1;
    auto Lp = L1_csr_->p;
    auto Li = L1_csr_->i;
    auto Lx = L1_csr_->x;
    const SpMP::LevelSchedule &schedule = *p2pScheduleWithTransitiveReduction;
    const int *perm = invPerm;
    double *x = x_in_;
    t1.start_timer();
#pragma omp parallel
    {
      #ifdef ENABLE_OPENMP
       int nthreads = omp_get_num_threads();
       int tid = omp_get_thread_num();
      #else
       int nthreads =1;
       int tid=1;
      #endif

      const int ntasks = schedule.ntasks;
      const short *nparents = schedule.nparentsForward;
      const std::vector<int> &threadBoundaries = schedule.threadBoundaries;
      const std::vector<int> &taskBoundaries = schedule.taskBoundaries;

      int nPerThread = (ntasks + nthreads - 1) / nthreads;
      int nBegin = std::min(nPerThread * tid, ntasks);
      int nEnd = std::min(nBegin + nPerThread, ntasks);

      volatile int *taskFinished = schedule.taskFinished;
      int **parents = schedule.parentsForward;

      memset((char *)(taskFinished + nBegin), 0, (nEnd - nBegin) * sizeof(int));

      synk::Barrier::getInstance()->wait(tid);

      for (int task = threadBoundaries[tid]; task < threadBoundaries[tid + 1];
           ++task) {
        SPMP_LEVEL_SCHEDULE_WAIT;

        for (int c = taskBoundaries[task]; c < taskBoundaries[task + 1]; ++c) {
          int group = perm[c];
          for (int g_ptr = group_ptr[group]; g_ptr < group_ptr[group + 1];
               g_ptr++) {
            int i = group_set[g_ptr];
            for (int j = Lp[i]; j < Lp[i + 1] - 1; j++) {
              x[i] -= Lx[j] * x[Li[j]];
            }
            x[i] /= Lx[Lp[i + 1] - 1];
          }
        }
        SPMP_LEVEL_SCHEDULE_NOTIFY;
      } // for each task
    }   // omp parallel
    t1.measure_elapsed_time();
    sym_lib::copy_vector(0, n_, x_in_, x_);
    return t1;
  }

public:
  /*
   * @brief Class constructor
   * @param L Sparse Matrix with CSR format
   * @param L_csc Sparse Matrix with CSC format
   * @param correct_x The right answer for x
   * @param name The name of the algorithm
   * @param nt number of threads
   */
  SpTrSv_LL_Tree_HDAGG_BFS_P2P(CSR *L, CSC *L_csc, double *correct_x,
                               std::string name, int nt, bool isLfactor,
                               bool bin_pack)
      : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
    L1_csr_ = L;
    L1_csc_ = L_csc;
    correct_x_ = correct_x;
    nthreads = nt;
    this->isLfactor = isLfactor;
    isfirst_time = true;
    this->bin_pack = bin_pack;
  };

  void computeCharacteristics(double &first, double &second, double &third,
                              double &avg) {

    std::vector<double> max_array(final_level_no, 0);
    std::vector<double> whole_array(final_level_no, 0);
    for (int lvl = 0; lvl < final_level_no; lvl++) {
      double total_nodes = 0;
      double max_size = 0;
      for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1];
           ++w_ptr) {
        int cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
        if (cc_size > max_size) {
          max_size = cc_size;
        }
        total_nodes += cc_size;
      }
      if (lvl == 0) {
        first = max_size / total_nodes;
      }
      if (lvl == 1) {
        second = max_size / total_nodes;
      }
      if (lvl == 2) {
        third = max_size / total_nodes;
      }
      max_array[lvl] = max_size;
      whole_array[lvl] = total_nodes;
    }

    double max_cc_sizes = 0;
    for (auto &iter : max_array) {
      max_cc_sizes += iter;
    }
    avg = max_cc_sizes / n_;
  }

  void computeStatistic() {
    double part_no = 0;
    for (int i = 0; i < final_part_ptr.size() - 1; i++) {
      if (final_part_ptr[i + 1] != final_part_ptr[i]) {
        part_no++;
      }
    }
    this->part_no = part_no;
    avg_par = part_no / final_level_no;
    std::vector<double> max_part(final_level_no, 0);
    std::vector<double> min_part(final_level_no, 0);
    std::vector<double> level_weight(final_level_no, 0);
    for (int lvl = 0; lvl < final_level_no; lvl++) {
      double total_nodes = 0;
      double max_size = 0;
      double min_size = n_;
      for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1];
           ++w_ptr) {
        double cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
        double cc_cost = 0;
        for (int i = final_part_ptr[w_ptr]; i < final_part_ptr[w_ptr + 1];
             i++) {
          int node = final_node_ptr[i];
          cc_cost += cost[node];
        }
        cc_size = cc_cost;
        if (cc_size > max_size) {
          max_size = cc_size;
        }
        if (cc_size < min_size) {
          min_size = cc_size;
        }
        total_nodes += cc_size;
      }
      max_part[lvl] = max_size;
      min_part[lvl] = min_size;
      level_weight[lvl] = total_nodes;
    }

    max_diff = 0;
    double total_weight = 0;
    for (int l = 0; l < final_level_no; l++) {
      max_diff += level_weight[l] * (max_part[l] - min_part[l]);
      total_weight += level_weight[l];
    }
    max_diff = max_diff / total_weight;
    var = 0;
  }
  double getSchedulingTime() { return Scheduling_time.elapsed_time; }
  void getStatistic(double &avg_par, double &max_diff, double &var) {
    avg_par = this->avg_par;
    max_diff = this->max_diff;
    var = this->var;
  }
  void getWaveStatistic(int &nlevels, int &part_no) {
    part_no = this->part_no;
    nlevels = this->final_level_no;
  }

  ~SpTrSv_LL_Tree_HDAGG_BFS_P2P() {
    delete p2pScheduleWithTransitiveReduction;
  };
};
#endif

#ifdef DAGP
//================================ SpTrSv DAGP LL
//================================
class SpTrSv_LL_DAGP : public Sptrsv_LL_Serial {
protected:
  int nthreads;
  timing_measurement Scheduling_time;
  int ngroups;
  int *final_level_ptr, *final_part_ptr, *final_node_ptr;
  int final_level_no;
  double avg_par, max_diff, var;
  int nlevels;
  int num_parts;
  void build_set() override {
    Scheduling_time.start_timer();
    // num_parts = num_threads_;
    int n_graphs = 1;
    int *parts;
    auto *sptrsv_spmv_oneb = convert_to_one_based(L1_csc_);
    CSC *p_graph = dagp_partition(sptrsv_spmv_oneb, n_, num_parts, parts);
    std::vector<int> h_level_ptr, h_par_ptr, h_partition;
    int n_levels;
    get_h_level_set(n_, p_graph, num_parts, parts, final_level_no,
                    final_level_ptr, final_part_ptr, final_node_ptr);
    std::vector<int> level_ptr(final_level_ptr,
                               final_level_ptr + final_level_no + 1);
    std::vector<int> part_ptr(final_part_ptr, final_part_ptr + num_parts + 1);
    std::vector<int> node_ptr(final_node_ptr, final_node_ptr + n_);

    if (ANALYSE_STATISTIC) {
      std::cout << "UNPACKING IS DEACTIVATED" << std::endl;
      computeStatistic();
      //            std::vector<int> final_level_ptr_copy(final_level_ptr,
      //            final_level_ptr + final_level_no + 1); std::vector<int>
      //            final_part_ptr_copy(final_part_ptr, final_part_ptr +
      //            final_level_ptr[final_level_no] + 1); std::vector<int>
      //            final_node_ptr_copy(final_node_ptr, final_node_ptr + n_);
      assert(test_unique(n_, final_node_ptr));
    }

    // print_hlevel_set("hlevel: \n", n_h_levels, h_level_ptr, h_part_ptr,
    //   h_partition );

    free(parts);
    delete p_graph;
    delete sptrsv_spmv_oneb;
    Scheduling_time.measure_elapsed_time();
  }

  timing_measurement fused_code() override {
    timing_measurement t1;
    t1.start_timer();
    sptrsv_csr_lbc(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                   final_level_no, final_level_ptr, final_part_ptr,
                   final_node_ptr);
    t1.measure_elapsed_time();
    sym_lib::copy_vector(0, n_, x_in_, x_);
    return t1;
  }

public:
  /*
   * @brief Class constructor
   * @param L Sparse Matrix with CSR format
   * @param L_csc Sparse Matrix with CSC format
   * @param correct_x The right answer for x
   * @param name The name of the algorithm
   * @param nt number of threads
   */
  SpTrSv_LL_DAGP(CSR *L, CSC *L_csc, double *correct_x, std::string name,
                 int num_parts, int nt)
      : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
    L1_csr_ = L;
    L1_csc_ = L_csc;
    correct_x_ = correct_x;
    nthreads = nt;
    this->num_parts = num_parts;
  };

  void computeStatistic() {
    double part_no = final_level_ptr[final_level_no];
    ngroups = part_no;
    avg_par = part_no / final_level_no;
    std::vector<double> max_part(final_level_no, 0);
    std::vector<double> min_part(final_level_no, 0);
    std::vector<double> level_weight(final_level_no, 0);
    for (int lvl = 0; lvl < final_level_no; lvl++) {
      double total_nodes = 0;
      double max_size = 0;
      double min_size = n_;
      for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1];
           ++w_ptr) {
        int cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
        if (cc_size > max_size) {
          max_size = cc_size;
        }
        if (cc_size < min_size) {
          min_size = cc_size;
        }
        total_nodes += cc_size;
      }
      max_part[lvl] = max_size;
      min_part[lvl] = min_size;
      level_weight[lvl] = total_nodes;
    }

    max_diff = 0;
    double total_weight = 0;
    for (int l = 0; l < final_level_no; l++) {
      max_diff += level_weight[l] * (max_part[l] - min_part[l]);
      total_weight += level_weight[l];
    }
    max_diff = max_diff / total_weight;
    var = 0;
  }

  double getSchedulingTime() { return Scheduling_time.elapsed_time; }
  void getStatistic(double &avg_par, double &max_diff, double &var) {
    avg_par = this->avg_par;
    max_diff = this->max_diff;
    var = this->var;
  }
  void getWaveStatistic(int &nlevels, int &part_no) {
    part_no = ngroups;
    nlevels = this->final_level_no;
  }
  ~SpTrSv_LL_DAGP() {
    delete[] final_level_ptr;
    delete[] final_part_ptr;
    delete[] final_node_ptr;
  }
};
#endif

#ifdef SPMP
//================================ SpTrSv SpMP LL ================================
class SpTrsv_LL_SpMP : public Sptrsv_LL_Serial {
protected:
  int nlevels, nthreads;
  std::vector<int> DAG_ptr;
  std::vector<int> DAG_set;
  timing_measurement scheduling_time;
  CSR *LBC_A;
  SpMP::LevelSchedule *p2pScheduleWithTransitiveReduction;
  const int *invPerm;
  void build_set() override {
    // Convert lbc CSR version to SpMP CSR version
    SpMP::CSR *A = new SpMP::CSR();
    Convert_LBCCSR_to_SpMP(LBC_A, A);
    scheduling_time.start_timer();
    /////////////////////////////////////////////////////////////////////////////
    // Construct schedules
    /////////////////////////////////////////////////////////////////////////////
    bool hasZeroDiag = A->hasZeroDiag();
    if (hasZeroDiag) {
      fprintf(stderr,
              "Matrix has a zero diagonal element. Can't do Gauss Seidel\n");
    }

    int *symRowPtr = NULL, *symColIdx = NULL, *symDiagPtr = NULL,
        *symExtPtr = NULL;
    bool wasSymmetric = getSymmetricNnzPattern(A, &symRowPtr, &symDiagPtr,
                                               &symExtPtr, &symColIdx);

    p2pScheduleWithTransitiveReduction = new SpMP::LevelSchedule;
    if (wasSymmetric) {
      p2pScheduleWithTransitiveReduction->constructTaskGraph(*A);
    } else {
      p2pScheduleWithTransitiveReduction->constructTaskGraph(
          A->m, symRowPtr, symDiagPtr, symExtPtr, symColIdx,
          SpMP::PrefixSumCostFunction(symRowPtr));

      FREE(symRowPtr);
      FREE(symColIdx);
      FREE(symDiagPtr);
      FREE(symExtPtr);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Reorder matrix
    /////////////////////////////////////////////////////////////////////////////

    const int *perm = p2pScheduleWithTransitiveReduction->origToThreadContPerm;
    const int *invPerm =
        p2pScheduleWithTransitiveReduction->threadContToOrigPerm;
    assert(SpMP::isPerm(perm, A->m));
    assert(SpMP::isPerm(invPerm, A->m));
    this->invPerm = invPerm;

    delete A;
    scheduling_time.measure_elapsed_time();
  }

  //        timing_measurement fused_code() override {
  //            timing_measurement t1;
  //            auto A = *L;
  //            int base = A.getBase();
  //            const int *rowptr = A.rowptr - base;
  //            const int *colidx = A.colidx - base;
  //            const double *values = A.values - base;
  //            const double *idiag = A.idiag - base;
  //            const SpMP::LevelSchedule &schedule =
  //            *p2pScheduleWithTransitiveReduction; const int* perm = invPerm;
  //            double* b = x_in_;
  //            t1.start_timer();
  //            #pragma omp parallel
  //            {
  //                int nthreads = omp_get_num_threads();
  //                int tid = omp_get_thread_num();
  //
  //                const int ntasks = schedule.ntasks;
  //                const short *nparents = schedule.nparentsForward;
  //                const std::vector<int>& threadBoundaries =
  //                schedule.threadBoundaries; const std::vector<int>&
  //                taskBoundaries = schedule.taskBoundaries; int nPerThread =
  //                (ntasks + nthreads - 1)/nthreads; int nBegin =
  //                std::min(nPerThread*tid, ntasks); int nEnd = std::min(nBegin
  //                + nPerThread, ntasks);
  //
  //                volatile int *taskFinished = schedule.taskFinished;
  //                int **parents = schedule.parentsForward;
  //
  //                memset((char *)(taskFinished + nBegin), 0, (nEnd -
  //                nBegin)*sizeof(int));
  //
  //                synk::Barrier::getInstance()->wait(tid);
  //
  //                for (int task = threadBoundaries[tid]; task <
  //                threadBoundaries[tid + 1]; ++task) {
  //                    SPMP_LEVEL_SCHEDULE_WAIT;
  //                    for (int i = taskBoundaries[task]; i <
  //                    taskBoundaries[task + 1]; ++i) {
  //                        int row = perm[i] + base;
  //                        double sum = b[row];
  //                        for (int j = rowptr[row]; j < rowptr[row + 1]; ++j)
  //                        {
  //                            sum -= values[j]*b[colidx[j]];
  //                        }
  //                        b[row] = sum*idiag[row];
  //                    }
  //                    SPMP_LEVEL_SCHEDULE_NOTIFY;
  //                } // for each task
  //            } // omp parallel
  //
  //
  //            t1.measure_elapsed_time();
  //            sym_lib::copy_vector(0, n_, x_in_, x_);
  //            return t1;
  //        }
  timing_measurement fused_code() override {
    timing_measurement t1;
    auto Lp = L1_csr_->p;
    auto Li = L1_csr_->i;
    auto Lx = L1_csr_->x;
    const SpMP::LevelSchedule &schedule = *p2pScheduleWithTransitiveReduction;
    const int *perm = invPerm;
    double *x = x_in_;
    t1.start_timer();
#pragma omp parallel
    {
      #ifdef ENABLE_OPENMP
       int nthreads = omp_get_num_threads();
       int tid = omp_get_thread_num();
      #else
       int nthreads =1;
       int tid=1;
      #endif


      const int ntasks = schedule.ntasks;
      const short *nparents = schedule.nparentsForward;
      const std::vector<int> &threadBoundaries = schedule.threadBoundaries;
      const std::vector<int> &taskBoundaries = schedule.taskBoundaries;
      int nPerThread = (ntasks + nthreads - 1) / nthreads;
      int nBegin = std::min(nPerThread * tid, ntasks);
      int nEnd = std::min(nBegin + nPerThread, ntasks);

      volatile int *taskFinished = schedule.taskFinished;
      int **parents = schedule.parentsForward;

      memset((char *)(taskFinished + nBegin), 0, (nEnd - nBegin) * sizeof(int));

      synk::Barrier::getInstance()->wait(tid);

      for (int task = threadBoundaries[tid]; task < threadBoundaries[tid + 1];
           ++task) {
        SPMP_LEVEL_SCHEDULE_WAIT;
        for (int c = taskBoundaries[task]; c < taskBoundaries[task + 1]; ++c) {
          int row = perm[c];
          for (int j = Lp[row]; j < Lp[row + 1] - 1; j++) {
            x[row] -= Lx[j] * x[Li[j]];
          }
          x[row] /= Lx[Lp[row + 1] - 1];
        }
        SPMP_LEVEL_SCHEDULE_NOTIFY;
      } // for each task
    }   // omp parallel
    t1.measure_elapsed_time();
    sym_lib::copy_vector(0, n_, x_in_, x_);
    return t1;
  }

public:
  /*
   * @brief Class constructor
   * @param L Sparse Matrix with CSR format
   * @param L_csc Sparse Matrix with CSC format
   * @param correct_x The right answer for x
   * @param name The name of the algorithm
   * @param nt number of threads
   */
  SpTrsv_LL_SpMP(CSR *A, CSR *L, CSC *L_csc, double *correct_x,
                 std::string name, int nt)
      : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
    this->LBC_A = A;
    L1_csr_ = L;
    L1_csc_ = L_csc;
    correct_x_ = correct_x;
    this->nthreads = nt;
   #ifdef ENABLE_OPENMP
 omp_set_num_threads(nt);
#endif
  };

  double getNumberOfP2P() {
    auto nparents = p2pScheduleWithTransitiveReduction->nparentsForward;
    const std::vector<int> &threadBoundaries =
        p2pScheduleWithTransitiveReduction->threadBoundaries;
    double incomming_edges = 0;
    for (int tid = 0; tid < nthreads; tid++) {
      for (int task = threadBoundaries[tid]; task < threadBoundaries[tid + 1];
           ++task) {
        incomming_edges += nparents[task];
      }
    }
    return incomming_edges;
  }

  double getSchedulingTime() { return scheduling_time.elapsed_time; };
  ~SpTrsv_LL_SpMP() { delete p2pScheduleWithTransitiveReduction; };
};
#endif
//================================ SpTRSV LBC LL ================================
class SpTrSv_LL_LBC : public Sptrsv_LL_Serial {
protected:
  int nthreads;
  timing_measurement Scheduling_time;
  std::vector<int> DAG_ptr;
  std::vector<int> DAG_set;
  int levelNo;
  int *final_level_ptr, *final_part_ptr, *final_node_ptr;
  int final_level_no;
  int part_no;
  int lp_, cp_, ic_;
  bool first_time;
  double avg_par, max_diff, var;
  std::vector<double> cost;
  void build_set() override {
    Scheduling_time.start_timer();
    cost.resize(n_, 0);
    auto CSC_Lp = L1_csc_->p;
    auto CSC_Li = L1_csc_->i;
    auto CSR_Lp = L1_csr_->p;
    auto CSR_Li = L1_csr_->i;
    HDAGG::costComputation(n_, CSC_Lp, CSC_Li, CSR_Lp, CSR_Li, HDAGG::SpTrSv_LL,
                           NULL, NULL, false, cost);

    get_coarse_levelSet_DAG_CSC_tree(
        n_, L1_csc_->p, L1_csc_->i, -1, final_level_no, final_level_ptr,
        part_no, final_part_ptr, final_node_ptr, lp_, ic_, cp_, cost.data());

    //      std::cout << "Static is COMPUTING" << std::endl;
    //      computeStatistic();

    Scheduling_time.measure_elapsed_time();
  }

  timing_measurement fused_code() override {
    timing_measurement t1;
    t1.start_timer();
    sptrsv_csr_lbc(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                   final_level_no, final_level_ptr, final_part_ptr,
                   final_node_ptr);
    t1.measure_elapsed_time();
    sym_lib::copy_vector(0, n_, x_in_, x_);
    return t1;
  }

public:
  /*
   * @brief Class constructor
   * @param L Sparse Matrix with CSR format
   * @param L_csc Sparse Matrix with CSC format
   * @param correct_x The right answer for x
   * @param name The name of the algorithm
   * @param nt number of threads
   */
  SpTrSv_LL_LBC(CSR *L, CSC *L_csc, double *correct_x, std::string name, int lp)
      : Sptrsv_LL_Serial(L, L_csc, correct_x, name) {
    L1_csr_ = L;
    L1_csc_ = L_csc;
    correct_x_ = correct_x;
    lp_ = lp;
    nthreads = lp;
  };

  void setP2_P3(int p2, int p3) {
    //            this->cp_ = p3;
    //            this->ic_ = p2;

    this->cp_ = p3;
    this->ic_ = p2;
  }

  void computeStatistic() {
    this->part_no = part_no;
    avg_par = part_no / final_level_no;
    std::vector<double> max_part(final_level_no, 0);
    std::vector<double> min_part(final_level_no, 0);
    std::vector<double> level_weight(final_level_no, 0);
    for (int lvl = 0; lvl < final_level_no; lvl++) {
      double total_nodes = 0;
      double max_size = 0;
      double min_size = n_;
      for (int w_ptr = final_level_ptr[lvl]; w_ptr < final_level_ptr[lvl + 1];
           ++w_ptr) {
        double cc_size = final_part_ptr[w_ptr + 1] - final_part_ptr[w_ptr];
        double cc_cost = 0;
        for (int i = final_part_ptr[w_ptr]; i < final_part_ptr[w_ptr + 1];
             i++) {
          int node = final_node_ptr[i];
          cc_cost += cost[node];
        }
        cc_size = cc_cost;
        if (cc_size > max_size) {
          max_size = cc_size;
        }
        if (cc_size < min_size) {
          min_size = cc_size;
        }
        total_nodes += cc_size;
      }
      max_part[lvl] = max_size;
      min_part[lvl] = min_size;
      level_weight[lvl] = total_nodes;
    }

    max_diff = 0;
    double total_weight = 0;
    for (int l = 0; l < final_level_no; l++) {
      max_diff += level_weight[l] * (max_part[l] - min_part[l]);
      total_weight += level_weight[l];
    }
    max_diff = max_diff / total_weight;
    var = 0;
  }
  void getWaveStatistic(int &nlevels, int &part_no) {
    part_no = this->part_no;
    nlevels = this->final_level_no;
  }
  double getSchedulingTime() { return Scheduling_time.elapsed_time; }
  void getStatistic(double &avg_par, double &max_diff, double &var) {
    avg_par = this->avg_par;
    max_diff = this->max_diff;
    var = this->var;
  }

  ~SpTrSv_LL_LBC() {
    delete final_level_ptr;
    delete final_node_ptr;
    delete final_part_ptr;
  };
};

#ifdef MKL
//================================ SpTRSV MKL LL ================================
class SpTrSv_LL_MKL : public FusionDemo {
protected:
  int *diagonal_locs_;
  sparse_matrix_t csrA;
  CSR *A_oneb_;
  int ierr;
  double alpha = 1.0;
  const double beta = 0;
  struct matrix_descr descrB, descrL;
  timing_measurement Scheduling_time;

  void build_set() override {
    Scheduling_time.start_timer();
    MKL_INT job[7] = {1, 0, 0, 0, 0, 1, 1};
    MKL_INT info = 0;
    const MKL_INT size_csr = n_;
    char transa;
    char matdescra[6];
    transa = 'n';
    matdescra[0] = 't';
    matdescra[1] = 'l';
    matdescra[2] = 'n';
    matdescra[3] = 'c';

    // Create handle with matrix stored in CSR format
    auto statl = mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ZERO,
                                         size_csr, // number of rows
                                         size_csr, // number of cols
                                         L1_csr_->p, L1_csr_->p + 1, L1_csr_->i,
                                         L1_csr_->x);
    if (statl != SPARSE_STATUS_SUCCESS) {
      printf("analysis failed with %d;", statl);
    }

    // Create matrix descriptor
    descrL.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
    descrL.mode = SPARSE_FILL_MODE_LOWER;
    descrL.diag = SPARSE_DIAG_NON_UNIT;

    if (mkl_sparse_set_sv_hint(csrA, SPARSE_OPERATION_NON_TRANSPOSE, descrL,
                               10000) != SPARSE_STATUS_SUCCESS) {
      printf("optimization failed with %d;", statl);
    }

    if (mkl_sparse_optimize(csrA) != SPARSE_STATUS_SUCCESS) {
      printf("optimization failed with %d;", statl);
    }

    Scheduling_time.measure_elapsed_time();
  }

  void setting_up() override { std::fill_n(x_, n_, 1.0); }

  timing_measurement fused_code() override {
    timing_measurement t1;
    auto *y = new double[n_]();

    t1.start_timer();
    auto stat = mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, csrA,
                                  descrL, x_, y);

    t1.measure_elapsed_time();
    if (stat != SPARSE_STATUS_SUCCESS) {
      printf("analysis failed with %d;", stat);
    }
    std::copy(y, y + n_, x_);
    delete[] y;
    return t1;
  }

public:
  SpTrSv_LL_MKL(CSR *A, CSC *A_csc, double *correct_x, std::string name, int nt)
      : FusionDemo(A->n, A->nnz, name) {
    L1_csr_ = A;
    L1_csc_ = A_csc;
    correct_x_ = correct_x;
    MKL_Set_Num_Threads(nt);
  };

  double getSchedulingTime() { return Scheduling_time.elapsed_time; }
  ~SpTrSv_LL_MKL() override{};
};
#endif

} // namespace sym_lib

#endif
