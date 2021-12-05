//
// Created by kazem on 12/5/21.
//

#ifndef LBC_LIB_SPTRSV_PROFILER_UTILS_H
#define LBC_LIB_SPTRSV_PROFILER_UTILS_H

#include <sparse_blas_lib.h>
#include "Profiler.h"

namespace sym_lib{

 class SptrsvSerial : public FusionDemo {
 protected:
  timing_measurement fused_code() override {
   timing_measurement t1;
   t1.start_timer();
   sptrsv_csr(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_);
   t1.measure_elapsed_time();
   copy_vector(0,n_,x_in_,x_);
   return t1;
  }

 public:
  SptrsvSerial(CSR *L, CSC *L_csc,
               double *correct_x, std::string name) :
    FusionDemo(L->n, name) {
   L1_csr_ = L;
   L1_csc_ = L_csc;
   correct_x_ = correct_x;
  };
#ifdef PROFILE
  SptrsvSerial(CSR *L, CSC *L_csc, double *correct_x, std::string name,
               PAPIWrapper *pw):
    SptrsvSerial(L, L_csc, correct_x, name){
   pw_ = pw;
  }
#endif
  ~SptrsvSerial() override {};
 };

 class SptrsvSerialProfiler: public Profiler{
 protected:

  FusionDemo *construct_demo(CSR *L, CSC* L_csc, CSC* Lt_csc,CSR *A, CSC *A_csc,
                             double *correct_x, std::string name,
                             PAPIWrapper *pw, int p1, int p2, int p3) override {
   auto *fd = new SptrsvSerial(L, L_csc, correct_x,name,pw);
   return fd;
  }

 public:
  SptrsvSerialProfiler(std::vector<int>  event_codes, int event_thr, int inst_no):
    Profiler(event_codes, event_thr, inst_no){};
 };

}

#endif //LBC_LIB_SPTRSV_PROFILER_UTILS_H
