//
// Created by kazem on 3/16/20.
//

#ifndef FUSION_SPTRSV_DEMO_UTILS_H
#define FUSION_SPTRSV_DEMO_UTILS_H

#include <algorithm>
#include "aggregation/sparse_inspector.h"
#include "aggregation/lbc.h"
#include <cstring>
#include "aggregation/FusionDemo.h"
#include "sptrsv.h"

#ifdef GROUPING_ENABLED
#include "aggregation/group.h"
#include "aggregation/group_utils.h"
#endif

#include <executor.h>

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
//   std::cout<<"=>"<<level_no<<"\n";
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
/*
#ifdef PAPI
        SptrsvLevelSet (CSR *L, CSC *L_csc,
                        double *correct_x, std::string name, PAPIWrapper *pw) :
                SptrsvSerial(L, L_csc, correct_x, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            pw_ = pw;
        };
#endif
*/

        int *getLeveSet(){
         return level_set;
        }

        int *getLevelPtr(){
         return level_ptr;
        }

        int getLevelNo(){
         return level_no;
        }

        double averParallelism()
        {
         return 1.0*L1_csc_->n/level_no;
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
         get_coarse_levelSet_DAG_CSC_tree(n_, L1_csc_->p, L1_csc_->i,
                                          L1_csc_->stype,
                                          final_level_no,
                                          fina_level_ptr,part_no,
                                          final_part_ptr,final_node_ptr,
                                          lp_,cp_, ic_, cost);

         delete []cost;
        }

        timing_measurement fused_code() override {
         timing_measurement t1;

         t1.start_timer();
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

        double AverParalleism(){
         return fina_level_ptr[final_level_no]*1.0/final_level_no;
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


    class SptrsvLBCDAG : public SptrsvLBC {
    protected:
        void build_set() override {
         auto *cost = new double[n_]();
         for (int i = 0; i < n_; ++i) {
          cost[i] = L1_csr_->p[i + 1] - L1_csr_->p[i];
         }
         get_coarse_Level_set_DAG_CSC03(n_, L1_csc_->p, L1_csc_->i, final_level_no,
                                        fina_level_ptr, part_no, final_part_ptr,
                                        final_node_ptr, lp_, cp_, ic_, cost);
         delete[] cost;
        }

    public:
        SptrsvLBCDAG(CSR *L, CSC *L_csc, double *correct_x, std::string name, int lp,
                     int cp, int ic)
                : SptrsvLBC(L, L_csc, correct_x, name, lp, cp, ic) {}

/*
#ifdef PAPI
        SptrsvLBCDAG(CSR *L, CSC *L_csc, double *correct_x, std::string name, int lp,
                     int cp, int ic, PAPIWrapper *pw)
                : SptrsvLBC(L, L_csc, correct_x, name, lp, cp, ic) {
            pw_ = pw;
        }

#endif
*/
        ~SptrsvLBCDAG() {}
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

    };

    class SptrsvLBCDAGParallel : public SptrsvLBCDAG {
    protected:
        void build_set() override {
         auto *cost = new double[n_]();
         for (int i = 0; i < n_; ++i) {
          cost[i] = L1_csr_->p[i + 1] - L1_csr_->p[i];
         }
         get_coarse_Level_set_DAG_CSC03_parallel(
          n_, L1_csc_->p, L1_csc_->i, final_level_no, fina_level_ptr, part_no,
          final_part_ptr, final_node_ptr, lp_, cp_, ic_, cost, -1);
         delete[] cost;
        }

    public:
        SptrsvLBCDAGParallel(CSR *L, CSC *L_csc, double *correct_x, std::string name, int lp,
                     int cp, int ic)
                : SptrsvLBCDAG(L, L_csc, correct_x, name, lp, cp, ic) {}
        ~SptrsvLBCDAGParallel() {}

    };

    class SptrsvLBC_W_Sorting : public SptrsvSerial {
    protected:
        int final_level_no, *fina_level_ptr, *final_part_ptr, *final_node_ptr;
        int part_no;
        int lp_, cp_, ic_;

        bool f_sort;

        void build_set() override {

         auto *cost = new double[n_]();
         for (int i = 0; i < n_; ++i) {
          cost[i] = L1_csr_->p[i + 1] - L1_csr_->p[i];
         }

         get_coarse_Level_set_DAG_CSC03(n_, L1_csc_->p, L1_csc_->i, final_level_no,
                                        fina_level_ptr, part_no, final_part_ptr,
                                        final_node_ptr, lp_, cp_, ic_, cost);

         if(f_sort){
          part_no=fina_level_ptr[final_level_no];
          // Sorting the w partitions
          for (int i = 0; i < part_no; ++i) {
           std::sort(final_node_ptr + final_part_ptr[i],
                     final_node_ptr + final_part_ptr[i + 1]);
          }
         }
         delete[] cost;
        }

        timing_measurement fused_code() override {
         timing_measurement t1;

         t1.start_timer();
         sptrsv_csr_lbc(n_, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                        final_level_no, fina_level_ptr,
                        final_part_ptr, final_node_ptr);
         t1.measure_elapsed_time();
         copy_vector(0,n_,x_in_,x_);
         return t1;
        }

    public:
        SptrsvLBC_W_Sorting(CSR *L, CSC *L_csc, double *correct_x, std::string name,
                            int lp, int cp, int ic, bool flag)
                : SptrsvSerial(L, L_csc, correct_x, name){
         L1_csr_ = L;
         L1_csc_ = L_csc;
         correct_x_ = correct_x;
         lp_=lp; cp_=cp; ic_=ic;

         f_sort = flag;
        };
/*
#ifdef PAPI
        SptrsvLBC_W_Sorting(CSR *L, CSC *L_csc, double *correct_x, std::string name,
                            int lp, int cp, int ic, bool flag, PAPIWrapper *pw)
                : SptrsvSerial(L, L_csc, correct_x, name){
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            lp_=lp; cp_=cp; ic_=ic;

            f_sort = flag;
            pw_ = pw;
        };
#endif
*/

        int levels(){
         return final_level_no;
        }

        double averParallelism(){
         return fina_level_ptr[final_level_no]*1.0/final_level_no;
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


        ~SptrsvLBC_W_Sorting() {
         delete []fina_level_ptr;
         delete []final_part_ptr;
         delete []final_node_ptr;
        };
    };


#ifdef GROUPING_ENABLED
    class SpTrsvCSR_Grouping : public sym_lib::SptrsvSerial{
    protected:
        int *groupSet, *groupPtr, *groupInv, ngroup, nlevels, nthreads;
        int *levelPtr, *levelSet;
        int blksize=1;
        timing_measurement t_group, t_levelset;
        void build_set() override {
         t_group.start_timer();
         groupPtr = (int *)malloc(sizeof(int)*(L1_csc_->n+1));
         memset(groupPtr, 0, sizeof(int)*(1+L1_csc_->n));
         groupSet = (int *)malloc(sizeof(int)*L1_csc_->n);
         memset(groupSet, 0, sizeof(int)*L1_csc_->n);
         groupInv = (int *)malloc(sizeof(int)*L1_csc_->n);
         memset(groupInv, 0, sizeof(int)*L1_csc_->n);

         group g(L1_csr_->n, L1_csr_->p, L1_csr_->i);
         g.inspection_sptrsvcsr_v1(groupPtr, groupSet, ngroup, groupInv);
         t_group.measure_elapsed_time();

         t_levelset.start_timer();
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


         t_levelset.measure_elapsed_time();
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
                           double *correct_x, std::string name, int nt):
                SptrsvSerial(L, L_csc, correct_x, name){
         L1_csr_ = L;
         L1_csc_ = L_csc;
         correct_x_ = correct_x;
         nthreads = nt;
         blksize = 1;
        };
/*
#ifdef PAPI
        SpTrsvCSR_Grouping(CSR *L, CSC *L_csc,
                           double *correct_x, std::string name, int nt, PAPIWrapper *pw):
                SptrsvSerial(L, L_csc, correct_x, name){
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            nthreads = nt;
            blksize = 1;
            pw_=pw;
        };
#endif
*/

        timing_measurement groupTime(){
         return t_group;
        }
        timing_measurement levelsetTime(){
         return t_levelset;
        }

        double groupWith(){
         return 1.0*L1_csc_->n /ngroup;
        }

        int levels(){
         return nlevels;
        }

        double averParallelism(){
         return 1.0*ngroup/nlevels;
        }

        ~SpTrsvCSR_Grouping() override{};
    };


    /**
     * @brief: this class is for running code, which combines the grouping and lbc
     */
    class SpTrsvCSR_Grouping_H2 : public sym_lib::SptrsvSerial
    {
    protected:
        int *groupSet, *groupPtr, *groupInv, ngroup, nlevels, nthreads;

        int final_level_no, *fina_level_ptr, *final_part_ptr, *final_node_ptr;
        int part_no;
        int lp_, cp_, ic_;

        bool f_sort; // flag for sorting w;

        timing_measurement t_group, t_coarsen, t_sort;

        void build_set() override {
         groupPtr = (int *)malloc(sizeof(int)*(L1_csc_->n+1));
         memset(groupPtr, 0, sizeof(int)*(1+L1_csc_->n));
         groupSet = (int *)malloc(sizeof(int)*L1_csc_->n);
         memset(groupSet, 0, sizeof(int)*L1_csc_->n);
         groupInv = (int *)malloc(sizeof(int)*L1_csc_->n);
         memset(groupInv, 0, sizeof(int)*L1_csc_->n);


         group g(L1_csr_->n, L1_csr_->p, L1_csr_->i);

         t_group.start_timer();
         g.inspection_sptrsvcsr_v1(groupPtr, groupSet, ngroup, groupInv);
//            g.NaiveGrouping(L1_csr_->n,  groupPtr, groupSet, ngroup, groupInv, 1);
         t_group.measure_elapsed_time();

         t_coarsen.start_timer();
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
           cost[i] = L1_csr_->p[k+1] - L1_csr_->p[k];
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
         if(f_sort){
          part_no=fina_level_ptr[final_level_no];
          // Sorting the w partitions
          for (int i = 0; i < part_no; ++i) {
           std::sort(final_node_ptr + final_part_ptr[i],
                     final_node_ptr + final_part_ptr[i + 1]);
          }
         }
         t_sort.measure_elapsed_time();
         delete []cost;
        }

        timing_measurement fused_code() override {
         timing_measurement t1;
         t1.start_timer();
         sptrsv_csr_group_lbc(L1_csr_->n, L1_csr_->p, L1_csr_->i, L1_csr_->x, x_in_,
                              final_level_no, fina_level_ptr, final_part_ptr, final_node_ptr, groupPtr, groupSet);
         t1.measure_elapsed_time();
         sym_lib::copy_vector(0,n_,x_in_,x_);
         return t1;
        }
    public:
        SpTrsvCSR_Grouping_H2(CSR *L, CSC *L_csc,
                              double *correct_x, std::string name,
                              int lp, int cp, int ic, bool flag) :
                SptrsvSerial(L, L_csc, correct_x, name) {
         L1_csr_ = L;
         L1_csc_ = L_csc;
         correct_x_ = correct_x;
         lp_=lp; cp_=cp; ic_=ic;
         f_sort=flag;
        };
/*
#ifdef PAPI
        SpTrsvCSR_Grouping_H2(CSR *L, CSC *L_csc,
                              double *correct_x, std::string name,
                              int lp, int cp, int ic, bool flag, PAPIWrapper *pw) :
                SptrsvSerial(L, L_csc, correct_x, name) {
            L1_csr_ = L;
            L1_csc_ = L_csc;
            correct_x_ = correct_x;
            lp_=lp; cp_=cp; ic_=ic;
            f_sort=flag;
            pw_ = pw;
        };
#endif
*/
        int *getLevelPtr(){
         return fina_level_ptr;
        }

        int *getPartPtr(){
         return final_part_ptr;
        }

        int * getNodePtr(){
         return final_node_ptr;
        }

        int * getGroupPtr(){
         return groupPtr;
        }

        int * getGroupSetPtr(){
         return groupSet;
        }

        int getLevelNo(){
         return final_level_no;
        }

        int getPartNo(){
         return part_no;
        }

        int getGroupNo(){
         return ngroup;
        }

        double averParallelism(){
         return fina_level_ptr[final_level_no] * 1.0 / final_level_no;
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

        ~SpTrsvCSR_Grouping_H2 () override {
         delete []fina_level_ptr;
         delete []final_part_ptr;
         delete []final_node_ptr;
         free(groupPtr);
         free(groupSet);
         free(groupInv);
        };


    };

#endif // grouping enabled


} // namespace sym_lib




#endif //FUSION_SPTRSV_DEMO_UTILS_H
