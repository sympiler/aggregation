//
// Created by labuser (Bangtian Liu) on 12/22/20.
//

#ifndef LBC_LIB_SPTRSVPROFILER_DEMO_UTILS_H
#define LBC_LIB_SPTRSVPROFILER_DEMO_UTILS_H

#include <Profiler.h>
#include "../example/sptrsvprofiler_demo_utils.h"
#include "../example/sptrsv_demo_utils.h"


namespace sym_lib {
    class SptrsvLevelSetProfiler : public Profiler {
    protected:
        FusionDemo *construct_demo(CSR *L, CSC *L_csc,
                                   double *correct_x, std::string name,
                                   PAPIWrapper *pw, int num_threads) override {
            auto *fd = new SptrsvLevelSet(L, L_csc, correct_x, name, pw);
            return fd;
        }

    public:
        SptrsvLevelSetProfiler(std::vector<int>  event_codes, int event_thr, int inst_no):
            Profiler(event_codes, event_thr, inst_no){};
    };


    class SpTrsvCSR_Grouping_Profiler : public Profiler {
    protected:
        FusionDemo *construct_demo(CSR *L, CSC *L_csc,
                                   double *correct_x, std::string name,
                                   PAPIWrapper *pw, int num_threads) override {
            auto *fd = new SpTrsvCSR_Grouping(L, L_csc, correct_x, name, num_threads, pw);
            return fd;
        }

    public:
        SpTrsvCSR_Grouping_Profiler(std::vector<int>  event_codes, int event_thr, int inst_no):
                Profiler(event_codes, event_thr, inst_no){};
    };

    class SptrsvLBCDAGProfiler : public Profiler{
    protected:
        FusionDemo *construct_demo(CSR *L, CSC *L_csc,
                                    double *correct_x, std::string name,
                                    PAPIWrapper *pw, int num_threads, int p2, int p3) override{
            auto *fd = new SptrsvLBCDAG(L, L_csc, correct_x, "coarsen DAG",num_threads, p2, p3, pw);
            return fd;
        }

    public:
        SptrsvLBCDAGProfiler(std::vector<int>  event_codes, int event_thr, int inst_no):
                Profiler(event_codes, event_thr, inst_no){};
    };

    class SptrsvLBC_W_SortingProfiler : public Profiler{
    protected:
        FusionDemo *construct_demo(CSR *L, CSC *L_csc,
                                   double *correct_x, std::string name,
                                   PAPIWrapper *pw, int num_threads, int p2, int p3, bool flag) override{
            auto *fd = new SptrsvLBC_W_Sorting(L, L_csc, correct_x, "sorted coarsen DAG", num_threads, p2, p3, flag, pw);
            return fd;
        }
    public:
        SptrsvLBC_W_SortingProfiler(std::vector<int>  event_codes, int event_thr, int inst_no):
                Profiler(event_codes, event_thr, inst_no){};

    };

    class SpTrsvCSR_Grouping_H2_Profiler : public Profiler{
    protected:
        FusionDemo *construct_demo(CSR *L, CSC *L_csc,
                                   double *correct_x, std::string name,
                                   PAPIWrapper *pw, int num_threads, int p2, int p3, bool flag){
            auto *fd=new SpTrsvCSR_Grouping_H2(L, L_csc, correct_x, "grouping + lbc", num_threads, p2, p3, flag, pw);
            return fd;
        }
    public:
        SpTrsvCSR_Grouping_H2_Profiler(std::vector<int>  event_codes, int event_thr, int inst_no):
        Profiler(event_codes, event_thr, inst_no){};
    };

}

#endif //LBC_LIB_SPTRSVPROFILER_DEMO_UTILS_H
