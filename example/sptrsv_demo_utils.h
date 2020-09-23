//
// Created by kazem on 3/16/20.
//

#ifndef FUSION_SPTRSV_DEMO_UTILS_H
#define FUSION_SPTRSV_DEMO_UTILS_H

#include <sparse_inspector.h>
#include <lbc.h>

#include "FusionDemo.h"
#include "sparse_blas_lib.h"
#include <Group.h>
#include <Utils.h>
#include <executor.h>
#include <StatSpMat.h>

namespace sym_lib {

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

  ~SptrsvSerial() override {};
 };

 class SptrsvLevelSet : public SptrsvSerial {
 protected:
  int *level_set, *level_ptr, level_no;
  void build_set() override {

   level_no = build_levelSet_CSC(L1_csc_->n, L1_csc_->p, L1_csc_->i,
                                 level_ptr, level_set);
  }

  timing_measurement fused_code() override {
   timing_measurement t1;

   t1.start_timer();

   sptrsv_csr_levelset(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                           level_no, level_ptr, level_set);

   t1.measure_elapsed_time();
   copy_vector(0,n_,x_in_,x_);
   return t1;
  }

 public:
  SptrsvLevelSet (CSR *L, CSC *L_csc,
                  double *correct_x, std::string name) :
    SptrsvSerial(L, L_csc, correct_x, name) {
   L1_csr_ = L;
   L1_csc_ = L_csc;
   correct_x_ = correct_x;
  };

  int *getLeveSet(){
      return level_set;
  }

  int *getLevelPtr(){
      return level_ptr;
  }

  int getLevelNo(){
      return level_no;
  }


  ~SptrsvLevelSet () override {
   delete []level_ptr;
   delete []level_set;
  };
 };

 class SptrsvLBC : public SptrsvSerial {
 protected:
  int final_level_no, *fina_level_ptr, *final_part_ptr, *final_node_ptr;
  int part_no;
  int lp_, cp_, ic_;
  void build_set() override {
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
   
    sptrsv_csr_lbc(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                  final_level_no, fina_level_ptr,
                  final_part_ptr, final_node_ptr);
   t1.measure_elapsed_time();
   copy_vector(0,n_,x_in_,x_);
   return t1;
  }

 public:
  SptrsvLBC (CSR *L, CSC *L_csc,
               double *correct_x, std::string name,
               int lp, int cp, int ic) :
    SptrsvSerial(L, L_csc, correct_x, name) {
   L1_csr_ = L;
   L1_csc_ = L_csc;
   correct_x_ = correct_x;
   lp_=lp; cp_=cp; ic_=ic;
  };

  int *getLevelPtr(){
      return fina_level_ptr;
  }

  int *getPartPtr(){
      return final_part_ptr;
  }

  int * getNodePtr(){
      return final_node_ptr;
  }

  int getLevelNo(){
      return final_level_no;
  }

  int getPartNo(){
      return part_no;
  }

  ~SptrsvLBC () override {
   delete []fina_level_ptr;
   delete []final_part_ptr;
   delete []final_node_ptr;
  };
 };



}

namespace group_cols{
    using namespace sym_lib;

    class SpTrsvCSR_Grouping : public sym_lib::SptrsvSerial{
    protected:
        int *groupSet, *groupPtr, *groupInv, ngroup, nlevels, nthreads;
        int *levelPtr, *levelSet;
        int blksize=1;
        void build_set() override {
            groupPtr = (int *)malloc(sizeof(int)*(L1_csc_->n+1));
            memset(groupPtr, 0, sizeof(int)*(1+L1_csc_->n));
            groupSet = (int *)malloc(sizeof(int)*L1_csc_->n);
            memset(groupSet, 0, sizeof(int)*L1_csc_->n);
            groupInv = (int *)malloc(sizeof(int)*L1_csc_->n);
            memset(groupInv, 0, sizeof(int)*L1_csc_->n);

            group g(L1_csr_->n, L1_csr_->p, L1_csr_->i);

//            g.inspection_sptrsvcsr(groupPtr, groupSet, ngroup, groupInv);
            NaiveGrouping(L1_csr_->n,  groupPtr, groupSet, ngroup, groupInv, blksize);
            std::vector<std::vector<int>> DAG;
            DAG.resize(ngroup);

            fs_csr_inspector_dep(ngroup, groupPtr, groupSet, groupInv, L1_csr_->p, L1_csr_->i, DAG);

            size_t count=0;
            for (int j = 0; j < DAG.size(); ++j) {
                DAG[j].erase(std::unique(DAG[j].begin(), DAG[j].end()), DAG[j].end());
                count+=DAG[j].size();
            }
//   detectDAGCircle(DAG);

            int *gv, *gedg;
            gv = new int[L1_csc_->n+1]();
            gedg = new int[count+L1_csc_->n]();
            levelPtr = new int[L1_csc_->n+1]();
            levelSet = new int[L1_csc_->n]();

            long int cti,edges=0;
            for(cti = 0, edges = 0; cti < ngroup; cti++){
                gv[cti] = edges;
                gedg[edges++] = cti;
                for (int ctj = 0; ctj < DAG[cti].size(); ctj++) {
                    gedg[edges++] = DAG[cti][ctj];
//                if(DAG[cti][ctj]==0)printf("cti=%d, ctj=%d\n", cti, ctj);
                }
            }
            gv[cti] = edges;

            nlevels = buildLevelSet_CSC_Queue(ngroup, 0, gv, gedg, levelPtr, levelSet);

//            std::cout<<nlevels<<","<<ngroup*1.0/nlevels<<std::endl;

        }

        timing_measurement fused_code() override {
            //sym_lib::rhs_init(L1_csc_->n, L1_csc_->p, L1_csc_->i, L1_csc_->x, x_); // x is b
            timing_measurement t1;
            t1.start_timer();
            fs_csr_executor_sgroup(L1_csr_->n, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_, x_in_, groupPtr, groupSet, ngroup, nlevels, levelPtr, levelSet);
            t1.measure_elapsed_time();
            sym_lib::copy_vector(0,n_,x_in_,x_);
            return t1;
        }

    public:
        SpTrsvCSR_Grouping(CSR *L, CSC *L_csc,
                           double *correct_x, std::string name, int nt, int blksize_):
                SptrsvSerial(L, L_csc, correct_x, name){
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            nthreads = nt;
            blksize = blksize_;
        };

        ~SpTrsvCSR_Grouping() override{};
    };




}

#endif //FUSION_SPTRSV_DEMO_UTILS_H
