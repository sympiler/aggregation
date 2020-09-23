//
// Created by labuser (Bangtian Liu) on 9/15/20.
//

#ifndef LBC_LIB_GROUP_H
#define LBC_LIB_GROUP_H
#include <cstring>
#include <vector>
#include <cstdlib>

namespace group_cols
{

    class group {
    private:
        int *status;
        int *next;
        int *pre;
        int *mP;
        int *mI;
        int ncol;
        int *top;
        bool *visited;

        int *child;
        bool *lflag;

        int *levels;

        std::vector<std::vector<int>> levelset;


    public:
        group(int n, int *p, int *i) : ncol(n), mP(p), mI(i) {
            status = (int *)malloc(sizeof(int)*n);
            memset(status, 0, sizeof(int)*n);

            next = (int *)malloc(sizeof(int)*n);
            memset(next, 0, sizeof(int)*n);

            pre = (int *)malloc(sizeof(int)*n);
            memset(pre, 0, sizeof(int)*n);

            visited = (bool *)malloc(sizeof(bool)*n);
            memset(visited, 0, sizeof(bool)*n);

            top = (int *)malloc(sizeof(int)*n);
            memset(top, 0, sizeof(int)*n);

            child =(int *)malloc(sizeof(int)*ncol);
            memset(child, 0, sizeof(int)*ncol);

            levels=(int *)malloc(sizeof(int)*ncol);
            memset(levels, 0, sizeof(int)*ncol);

            lflag = (bool *)malloc(sizeof(bool)*ncol);
            memset(lflag, 0, sizeof(bool)*ncol);

            top = (int *)malloc(sizeof(int)*ncol);
            memset(top, 0, sizeof(int)*ncol);

            levelset.resize(ncol);
        }


        void groupInspection(int *groupPtr, int *groupSet, int &ngroup, int *groupInv);

        void nogroup(int *groupPtr, int *groupSet, int &ngroup, int *groupInv);

        void groupInspection_v1(int *groupPtr, int *groupSet, int &ngroup, int *groupInv);


        void groupInspection_v2(int *groupPtr, int *groupSet, int &ngroup, int *groupInv);


        void inspection_sptrsvcsr(int *groupPtr, int *groupSet, int &ngroup, int *groupInv);


        void inspection_sptrsvcsr_v1(int *groupPtr, int *groupSet, int &ngroup, int *groupInv);


        void inspection_sptrsvcsr_v2(int *groupPtr, int *groupSet, int &ngroup, int *groupInv);



    };

    /**
     * @brief group every blksize rows/columns
     * @param n  number of rows/cols
     * @param groupPtr the pointer to the starting address of one group
     * @param groupSet the pointer to index array
     * @param ngroup number of columns
     * @param ginv mapping column idx to group idx
     * @param blksize  parameter for grouping
     */
    void NaiveGrouping(int n, int *groupPtr, int *groupSet,  int &ngroup, int *ginv, int blksize=1);

}





#endif //LBC_LIB_GROUP_H
