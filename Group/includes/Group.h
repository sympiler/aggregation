//
// Created by labuser (Bangtian Liu) on 9/15/20.
//

#ifndef LBC_LIB_GROUP_H
#define LBC_LIB_GROUP_H
#include <cstring>
#include <vector>
#include <cstdlib>

namespace sym_lib
{

    class group {
    private:
        int *status; // status for each column to help to do grouping
        int *next; // pointer to its next column of one column/row
        int *pre; // pointer to its previous column of one column/row
        int *mP; // Pointer array in CSR/CSC
        int *mI; // Index array in CSR/CSC
        int ncol; // number of columns
        bool *visited; // indicate whether one column/row is visited or not, used for grouping

        int *child; // indicates its first dependent column/row based on first off-diagonal element

    public:
        /**
         * @brief construction function
         * @param n number of rows/columns
         * @param p  row/column pointer in the CSR/CSC format
         * @param i index array in the CSR/CSC format
         */
        group(int n, int *p, int *i) : ncol(n), mP(p), mI(i) {
            status = (int *)malloc(sizeof(int)*n);
            memset(status, 0, sizeof(int)*n);

            next = (int *)malloc(sizeof(int)*n);
            memset(next, 0, sizeof(int)*n);

            pre = (int *)malloc(sizeof(int)*n);
            memset(pre, 0, sizeof(int)*n);

            visited = (bool *)malloc(sizeof(bool)*n);
            memset(visited, 0, sizeof(bool)*n);

            child =(int *)malloc(sizeof(int)*ncol);
            memset(child, 0, sizeof(int)*ncol);


        }


        /**
         * @brief Grouping consecutive columns with dependence for sptrsv_csr
         * @param groupPtr Pointer to the starting location of one group
         * @param groupSet Pointer to the column indices in one group
         * @param ngroup Number of groups
         * @param groupInv  mapping the column Id to group Id
         */
        void inspection_sptrsvcsr_v1(int *groupPtr, int *groupSet, int &ngroup, int *groupInv);

        /**
         * @brief Grouping consecutive columns with dependence and also includes empty columns/rows for sptrsv_csr
         * @param groupPtr Pointer to the starting location of one group
         * @param groupSet Pointer to the column indices in one group
         * @param ngroup Number of groups
         * @param groupInv mapping the column Id to group Id
         */
        void inspection_sptrsvcsr_v2(int *groupPtr, int *groupSet, int &ngroup, int *groupInv);


        /**
         * @brief group every blksize rows/columns
         * @param n  number of rows/cols
         * @param groupPtr the pointer to the starting address of one group
         * @param groupSet the pointer to index array
         * @param ngroup number of columns
         * @param ginv mapping column idx to group idx
         * @param blksize  parameter for grouping
           **/
        void NaiveGrouping(int n, int *groupPtr, int *groupSet,  int &ngroup, int *ginv, int blksize=1);
    };



}





#endif //LBC_LIB_GROUP_H
