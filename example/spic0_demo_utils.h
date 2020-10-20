//
// Created by kazem on 2020-04-23.
//

#ifndef FUSION_SPIC0_SPTRSV_DEMO_UTILS_H
#define FUSION_SPIC0_SPTRSV_DEMO_UTILS_H

#include <lbc.h>
#include <sparse_utilities.h>
#include <sparse_blas_lib.h>

#include "FusionDemo.h"
#include <sparse_inspector.h>
#include <Group.h>
#include <Utils.h>
//
//#ifdef PROFILE
//#include "Profiler.h"
//#endif


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

     CSR * getFactor(){ return  factor_; }
 };

class Spic0CSCSerial : public FusionDemo {
    protected:
        CSC *factor_;
        void build_set() override {
            factor_ = csr_to_csc(A_csr_);
//            factor_ = A_csc_;
        }

        void setting_up() override {
            std::fill_n(x_,n_,1.0);
            copy_from_to(A_csc_,factor_);
        }

        timing_measurement fused_code() override {
            timing_measurement t1;

            t1.start_timer();
            ic0_csc(n_, factor_->x, factor_->p, factor_->i);
            t1.measure_elapsed_time();
            //print_vec("Lx:\n",0,n_,x_);
            return t1;
        }

        void testing() override {
            if(correct_x_)
                if (!is_equal(0, factor_->nnz, correct_x_, factor_->x,1e-3))
                       PRINT_CSV(name_ + " code != reference solution.\n");
        }

    public:
        Spic0CSCSerial(CSR *L, CSC *L_csc, CSR *A, CSC *A_csc,
                    double *correct_x, std::string name) :
                FusionDemo(L->n, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            A_csr_ = A;
            A_csc_ = A_csc;
            correct_x_ = correct_x;
            factor_=NULLPNTR;
        };

        ~Spic0CSCSerial() override {
            delete factor_;
        };

        double *Factor(){ return factor_->x;}
//
//        CSR * getFactor(){ return  factor_; }
};

class SpicoCSCLevelSet: public Spic0CSCSerial{
protected:
    int *level_set, *level_ptr, level_no;
    void build_set() override {
        Spic0CSCSerial::build_set();

        level_no = build_levelSet_CSC(L1_csc_->n, L1_csc_->p, L1_csc_->i,
                                      level_ptr, level_set);
    }

    timing_measurement fused_code() override {
        timing_measurement t1;

        t1.start_timer();

        spico_csc_levelset(n_, factor_->p, factor_->i, factor_->x,
                            level_no, level_ptr, level_set);

        t1.measure_elapsed_time();
//        copy_vector(0,n_,x_in_,x_);
        return t1;
    }

public:
    SpicoCSCLevelSet(CSR *L, CSC *L_csc, CSR *A, CSC *A_csc,
                     double *correct_x, std::string name): Spic0CSCSerial(L,L_csc,A,A_csc,correct_x,name){

    }

    int *getLeveSet(){
        return level_set;
    }

    int *getLevelPtr(){
        return level_ptr;
    }


    int numLevels(){
        return level_no;
    }

    double averparallelism(){
        return L1_csc_->n * 1.0/level_no;
    }

    ~SpicoCSCLevelSet(){
        delete [] level_ptr;
        delete [] level_set;
    }
};

class SpicoCSC_Grouping : public Spic0CSCSerial{
protected:
    int *groupSet, *groupPtr, *groupInv, ngroup, nlevels, nthreads;
    int *levelPtr, *levelSet;
    timing_measurement t_group, t_levelset;

    void build_set() override {
        Spic0CSCSerial::build_set();

        groupPtr = (int *)malloc(sizeof(int)*(L1_csc_->n+1));
        memset(groupPtr, 0, sizeof(int)*(1+L1_csc_->n));
        groupSet = (int *)malloc(sizeof(int)*L1_csc_->n);
        memset(groupSet, 0, sizeof(int)*L1_csc_->n);
        groupInv = (int *)malloc(sizeof(int)*L1_csc_->n);
        memset(groupInv, 0, sizeof(int)*L1_csc_->n);

        group g(L1_csc_->n, L1_csc_->p, L1_csc_->i);
        t_group.start_timer();
        g.inspection_spicocsc_v1(groupPtr, groupSet, ngroup, groupInv);
        t_group.measure_elapsed_time();

        t_levelset.start_timer();
        std::vector<std::vector<int>> DAG;
        DAG.resize(ngroup);

        fs_ic0csc_inspector_group(ngroup, groupPtr, groupSet, groupInv, L1_csc_->p, L1_csc_->i, DAG);

        size_t count=0;
        for (int j = 0; j < DAG.size(); ++j) {
            DAG[j].erase(std::unique(DAG[j].begin(), DAG[j].end()), DAG[j].end());
            count+=DAG[j].size();
        }

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
        t_levelset.measure_elapsed_time();
    }


    timing_measurement fused_code() override {
        timing_measurement t1;
        t1.start_timer();
        spico_csc_group_levelset(L1_csc_->n, factor_->p, factor_->i, factor_->x, nlevels,
                levelPtr, levelSet, groupPtr, groupSet);

        t1.measure_elapsed_time();
        return t1;
    }

public:
    SpicoCSC_Grouping(CSR *L, CSC *L_csc, CSR *A, CSC *A_csc,
                      double *correct_x, std::string name): Spic0CSCSerial(L,L_csc,A,A_csc,correct_x,name)
    {


    }

    timing_measurement groupTime(){
        return t_group;
    }

    timing_measurement levelsetTime(){
        return t_levelset;
    }


    double groupWidth(){
        return 1.0 * L1_csc_->n/ngroup;
    }

    int numLevels(){
        return nlevels;
    }

    double averparallelism(){
        return 1.0 * L1_csc_->n / nlevels;
    }

    ~SpicoCSC_Grouping(){
        free(groupSet);
        free(groupPtr);
        free(groupInv);

        delete [] levelPtr;
        delete [] levelSet;
    }
};

class SpicoCSC_Grouping_H2 : public Spic0CSCSerial{
protected:
    int *groupSet, *groupPtr, *groupInv, ngroup, nlevels, nthreads;
    int lp_, cp_, ic_;

    int final_level_no, *fina_level_ptr, *final_part_ptr, *final_node_ptr;
    int part_no;

    bool f_sort;

    timing_measurement t_group, t_coarsen, t_sort;

    void build_set() override {
        Spic0CSCSerial::build_set();

        groupPtr = (int *)malloc(sizeof(int)*(L1_csc_->n+1));
        memset(groupPtr, 0, sizeof(int)*(1+L1_csc_->n));
        groupSet = (int *)malloc(sizeof(int)*L1_csc_->n);
        memset(groupSet, 0, sizeof(int)*L1_csc_->n);
        groupInv = (int *)malloc(sizeof(int)*L1_csc_->n);
        memset(groupInv, 0, sizeof(int)*L1_csc_->n);

        group g(L1_csc_->n, L1_csc_->p, L1_csc_->i);
        t_group.start_timer();
        g.inspection_spicocsc_v1(groupPtr, groupSet, ngroup, groupInv);
        t_group.start_timer();

        t_coarsen.start_timer();
        std::vector<std::vector<int>> DAG;
        DAG.resize(ngroup);

        fs_ic0csc_inspector_group(ngroup, groupPtr, groupSet, groupInv, L1_csc_->p, L1_csc_->i, DAG);

        size_t count=0;
        for (int j = 0; j < DAG.size(); ++j) {
            DAG[j].erase(std::unique(DAG[j].begin(), DAG[j].end()), DAG[j].end());
            count+=DAG[j].size();
        }

        int *gv, *gedg;
        gv = new int[L1_csc_->n+1]();
        gedg = new int[count+L1_csc_->n]();

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

        auto *cost = new double[ngroup]();
        for (int i = 0; i < ngroup; ++i) {
            for (int j = groupPtr[i]; j < groupPtr[i+1]; ++j) {
                int k = groupSet[j];
                cost[i] += L1_csc_->p[k+1] - L1_csc_->p[k];
            }
        }

        get_coarse_Level_set_DAG_CSC03(ngroup,
                                       gv,
                                       gedg,
                                       final_level_no,
                                       fina_level_ptr,
                                       part_no,
                                       final_part_ptr,final_node_ptr,
                                       lp_,cp_, ic_, cost
        );

        nlevels = final_level_no;
        t_coarsen.measure_elapsed_time();

        t_sort.start_timer();
        if(f_sort)
        {
            part_no=fina_level_ptr[final_level_no];
            // Sorting the w partitions
            for (int i = 0; i < part_no; ++i) {
                std::sort(final_node_ptr + final_part_ptr[i],
                          final_node_ptr + final_part_ptr[i + 1]);
            }
        }
        t_sort.measure_elapsed_time();
        delete [] cost;
    }

    timing_measurement fused_code() override {
        timing_measurement t1;
        t1.start_timer();
        spic0_csc_group_lbc(
                factor_->n, factor_->p, factor_->i, factor_->x,
                final_level_no, fina_level_ptr, final_part_ptr, final_node_ptr, groupPtr, groupSet
                );
        t1.measure_elapsed_time();
        return t1;
    }

public:
    SpicoCSC_Grouping_H2(CSR *L, CSC *L_csc, CSR *A, CSC *A_csc,
                         double *correct_x, std::string name,
                         int l_param, int c_param, int i_param, bool sort_flag):
            Spic0CSCSerial(L,L_csc,A,A_csc,correct_x,name),
            lp_(l_param),cp_(c_param),ic_(i_param){
                f_sort=sort_flag;
    };

    double groupWidth(){
        return 1.0 * L1_csc_->n/ngroup;
    }

    int numLevels(){
        return nlevels;
    }

    double averparallelism(){
        return 1.0 *  fina_level_ptr[nlevels]/ nlevels;
    }

    double averWsize(){
        int part_no=fina_level_ptr[final_level_no];
        int ncol=0;
        // Sorting the w partitions
        for (int i = 0; i < part_no; ++i) {
            for(int j=final_part_ptr[i]; j<final_part_ptr[i+1]; j++)
            {
                int k = final_node_ptr[j];
                for (int l = groupPtr[k]; l < groupPtr[k+1]; ++l) {
                    int k1=groupSet[l];
                    ncol +=(L1_csr_->p[k1+1] - L1_csr_->p[k1]);
                }
            }
        }
        return ncol*1.0/part_no;
    }

    double consecutiveRatio(){
        part_no=fina_level_ptr[final_level_no];
        // Sorting the w partitions
        double aver_ratio=0;
        int sum=0;
        int t_sum=0;
        for (int i = 0; i < part_no; ++i) {
//             t_sum += (final_part_ptr[i+1]-final_part_ptr[i]-1);
            for (int k1 = final_part_ptr[i]; k1 < final_part_ptr[i+1]-1; ++k1) {
                if((final_node_ptr[k1+1]-final_node_ptr[k1])==1){
                    sum+=(groupPtr[k1+1]-groupPtr[k1]);
                    t_sum+=(groupPtr[k1+1]-groupPtr[k1]);
                }
                else{
                    sum+=(groupPtr[k1+1]-groupPtr[k1]-1);
                    t_sum+=(groupPtr[k1+1]-groupPtr[k1]);
                }
            }

            sum +=groupPtr[final_part_ptr[i+1]]-groupPtr[final_part_ptr[i+1]-1];
            t_sum += groupPtr[final_part_ptr[i+1]]-groupPtr[final_part_ptr[i+1]];
            t_sum = t_sum-1;
        }
        return 1.0*sum/t_sum;
    }


    timing_measurement groupTime(){
        return t_group;
    }

    timing_measurement coarsenTime(){
        return t_coarsen;
    }

    timing_measurement sortTime(){
        return t_sort;
    }

    ~SpicoCSC_Grouping_H2() override{
        delete [] fina_level_ptr;
        delete [] final_part_ptr;
        delete [] final_node_ptr;

        free(groupPtr);
        free(groupSet);
    }


};



class Spic0CSCParallelLBC : public Spic0CSCSerial {
protected:
    int lp_, cp_, ic_;

    int final_level_no, *fina_level_ptr, *final_part_ptr, *final_node_ptr;
    int part_no;

    bool f_sort;

    void build_set() override {
        Spic0CSCSerial::build_set();
        auto *cost = new double[n_]();

        for (int i = 0; i < n_; ++i) {
            cost[i] = L1_csc_->p[i + 1] - L1_csc_->p[i];
        }

        get_coarse_Level_set_DAG_CSC03(n_, L1_csc_->p, L1_csc_->i,
                                       final_level_no,
                                       fina_level_ptr, part_no,
                                       final_part_ptr, final_node_ptr,
                                       lp_, cp_, ic_, cost);
        if(f_sort){
            if(f_sort){
                part_no=fina_level_ptr[final_level_no];
                // Sorting the w partitions
                for (int i = 0; i < part_no; ++i) {
                    std::sort(final_node_ptr + final_part_ptr[i],
                              final_node_ptr + final_part_ptr[i + 1]);
                }
            }
        }


        delete[] cost;
    }

    timing_measurement fused_code() override {
        timing_measurement t1;
        t1.start_timer();
        spico_csc_lbc(n_, factor_->x, factor_->p, factor_->i,
                      final_level_no, fina_level_ptr,
                      final_part_ptr, final_node_ptr);
        t1.measure_elapsed_time();
        return t1;
    }

public:
    Spic0CSCParallelLBC(CSR *L, CSC *L_csc, CSR *A, CSC *A_csc,
    double *correct_x, std::string name,
    int l_param, int c_param, int i_param, bool sort_flag) :
    Spic0CSCSerial(L,L_csc,A,A_csc,correct_x,name),
    lp_(l_param),cp_(c_param),ic_(i_param){
        f_sort = sort_flag;
    };

    int numLevels(){
        return final_level_no;
    }

    double averparallelism(){
        return 1.0 * fina_level_ptr[final_level_no] / final_level_no;
    }


    double averWsize(){
        part_no=fina_level_ptr[final_level_no];
        // Sorting the w partitions
        int ncols=0;
        for (int i = 0; i < part_no; ++i) {
            for (int k1 = final_part_ptr[i]; k1 < final_part_ptr[i+1]; ++k1) {
                int k=final_node_ptr[k1];
                ncols+=(L1_csr_->p[k+1]-L1_csr_->p[k]);
            }
        }
        return 1.0*ncols/part_no;
    }

    double consecutiveRatio(){
        part_no=fina_level_ptr[final_level_no];
        // Sorting the w partitions
        double aver_ratio=0;
        int sum=0;
        int t_sum=0;
        for (int i = 0; i < part_no; ++i) {

            t_sum += (final_part_ptr[i+1]-final_part_ptr[i]-1);
            for (int k1 = final_part_ptr[i]; k1 < final_part_ptr[i+1]-1; ++k1) {
                if((final_node_ptr[k1+1]-final_node_ptr[k1])==1)sum+=1;
            }
        }
        return 1.0*sum/t_sum;
    }

    ~Spic0CSCParallelLBC(){
        delete []fina_level_ptr;
        delete []final_part_ptr;
        delete []final_node_ptr;
    }
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
      get_coarse_Level_set_DAG_CSC03(n_, L1_csc_->p, L1_csc_->i,
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
