//
// Created by kazem on 2020-04-23.
//

#ifndef FUSION_SPIC0_SPTRSV_DEMO_UTILS_H
#define FUSION_SPIC0_SPTRSV_DEMO_UTILS_H

#include <lbc.h>
#include <sparse_utilities.h>
#include <sparse_blas_lib.h>
#include "FusionDemo.h"

#ifdef PROFILE
#include "Profiler.h"
#endif
namespace sym_lib {


 class Spic0Serial : public FusionDemo {
 protected:
  CSR *factor_;
  void build_set() override {
   factor_ = csc_to_csr(A_csc_);
  }

  void setting_up() override {
   std::fill_n(x_,n_,1.0);
   copy_from_to(A_csr_,factor_);
  }

  timing_measurement fused_code() override {
   timing_measurement t1;

   t1.start_timer();
   spic0_csr(n_, factor_->x, factor_->p, factor_->i);
   t1.measure_elapsed_time();
   //print_vec("Lx:\n",0,n_,x_);
   return t1;
  }

  void testing() override {
   if(correct_x_)
    if (!is_equal(0, factor_->nnz, correct_x_, factor_->x,1e-6))
      PRINT_LOG(name_ + " code != reference solution.\n");
  }

 public:
  Spic0Serial(CSR *L, CSC *L_csc, CSR *A, CSC *A_csc,
                            double *correct_x, std::string name) :
    FusionDemo(L->n, name) {
   L1_csr_ = L;
   L1_csc_ = L_csc;
   A_csr_ = A;
   A_csc_ = A_csc;
   correct_x_ = correct_x;
   factor_=NULLPNTR;
  };

  ~Spic0Serial() override {
   delete factor_;
  };

  double *Factor(){ return factor_->x;}
 };

 class Spic0ParallelLBC : public Spic0Serial {
 protected:
  int lp_, cp_, ic_;

  int final_level_no, *fina_level_ptr, *final_part_ptr, *final_node_ptr;
  int part_no;
  void build_set() override {
   Spic0Serial::build_set();
   auto *cost = new double[n_]();
   for (int i = 0; i < n_; ++i) {
    cost[i] = L1_csr_->p[i+1] - L1_csr_->p[i];
   }
   get_coarse_levelSet_DAG_CSC_tree(n_, L1_csr_->p, L1_csr_->i,
                                    final_level_no,
                                    fina_level_ptr,part_no,
                                    final_part_ptr,final_node_ptr,
                                    lp_,cp_, ic_, cost);

   delete []cost;

  }

  timing_measurement fused_code() override {
   timing_measurement t1;

   t1.start_timer();
   spic0_csr_lbc(n_, factor_->x, factor_->p, factor_->i,
                 final_level_no, fina_level_ptr,
                 final_part_ptr, final_node_ptr);
   t1.measure_elapsed_time();
   //print_vec("Lx:\n",0,n_,x_);
   return t1;
  }

 public:
  Spic0ParallelLBC(CSR *L, CSC *L_csc, CSR *A, CSC *A_csc,
                            double *correct_x, std::string name,
                                 int l_param, int c_param, int i_param) :
    Spic0Serial(L,L_csc,A,A_csc,correct_x,name),
    lp_(l_param),cp_(c_param),ic_(i_param){};

#ifdef PROFILE
  Spic0SptrsvParallelLBCNonFused(CSR *L, CSC *L_csc, CSR *A, CSC *A_csc,
  double *correct_x, std::string name,
  int l_param, int c_param, int i_param, PAPIWrapper *pw):
    Spic0SptrsvSerialNonFused(L,L_csc,A,A_csc,correct_x,name),
    lp_(l_param),cp_(c_param),ic_(i_param){
   pw_ = pw;
   num_threads_ = lp_;
  }
#endif
  ~Spic0ParallelLBC(){
   delete []fina_level_ptr;
   delete []final_part_ptr;
   delete []final_node_ptr;
  }
 };


}
#endif //FUSION_SPIC0_SPTRSV_DEMO_UTILS_H
