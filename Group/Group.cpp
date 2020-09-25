//
// Created by labuser (Bangtian Liu) on 9/16/20.
//

#include <Group.h>

namespace sym_lib
{

    void group::inspection_sptrsvcsr_v1(int *groupPtr, int *groupSet, int &ngroup, int *groupInv) {
        for (int i = 0; i < ncol; ++i) {
            int len = mP[i + 1] - mP[i];
            child[i] = len > 1 ? mI[mP[i+1] - 2] : -1;
        }

        for (int i = 1; i < ncol; ++i) {
            int curIdx = i;
            int chIdx = child[curIdx];
            if(chIdx!=-1 && curIdx-chIdx==1){
                next[chIdx]=curIdx;
                pre[curIdx]=chIdx;
            }
        }

        memset(visited, 0, sizeof(bool)*ncol);

        int pos=0; int index=0;

        for (int i = 0; i < ncol; ++i) {
            int curIdx = i;
            if (!visited[curIdx]) {
                groupPtr[index]=pos;
                visited[curIdx]=true;
                groupSet[pos++]=curIdx;
                groupInv[curIdx] = index;
                int ch = next[curIdx];

                while (ch!=0 && !visited[ch])
                {
                    curIdx=ch;
                    groupSet[pos++] = curIdx;
                    groupInv[curIdx] = index;
                    visited[curIdx]= true;
                    ch = next[curIdx];
                }
                index++;
            }
        }

        groupPtr[index] = pos;
        ngroup = index;
    }


    void group::inspection_sptrsvcsr_v2(int *groupPtr, int *groupSet, int &ngroup, int *groupInv) {
        for (int i = 0; i < ncol; ++i) {
            child[i]=-1;
        }

        bool *rCnt=(bool *)malloc(sizeof(bool)*ncol);
        memset(rCnt, 0, sizeof(bool)*ncol);

        for (int i = 0; i < ncol; ++i) {
            int len = mP[i + 1] - mP[i];
            if(len==1)rCnt[i]=true;
            int tidx = len > 1 ? mI[mP[i+1] - 2] : -1;
            if(tidx!=-1 ){
                if(child[tidx]==-1)child[tidx]=i;
            }
        }

        for (int i = 0; i < ncol; ++i) {
            int curIdx = i;

            if(!visited[curIdx])
            {
                visited[curIdx] = true;
                int chIdx = child[curIdx];

                while ((chIdx-curIdx==1 || rCnt[chIdx]) && !visited[chIdx] )
                {
                    visited[chIdx]=true;
                    pre[chIdx] = curIdx;
                    next[curIdx] = chIdx;

                    curIdx=chIdx;
                    chIdx=child[curIdx];
                }
//
            }
        }

        memset(visited, 0, sizeof(bool)*ncol);

        int pos=0; int index=0;
        int len;
        for (int i = 0; i < ncol; ++i) {
            int curIdx = i;
            len=0;
            if (!visited[curIdx]) {
                groupPtr[index]=pos;
                visited[curIdx]=true;
                groupSet[pos++]=curIdx;
                groupInv[curIdx] = index;
                int ch = next[curIdx];
                ++len;
                while (ch!=0 && !visited[ch])
                {
                    curIdx=ch;
                    groupSet[pos++] = curIdx;
                    ++len;
                    groupInv[curIdx] = index;
                    visited[curIdx]= true;
                    ch = next[curIdx];
                }
                index++;
            }
        }

        groupPtr[index] = pos;
        ngroup = index;
    }


    void group::NaiveGrouping(int n, int *groupPtr, int *groupSet,  int &ngroup, int *ginv, int blksize)
    {
        int index=0;
        int pos=0;
        for (int i = 0; i < n; ++i) {
            if(i%blksize==0) groupPtr[index++]=pos;
            groupSet[pos++]=i;
            ginv[i]=index-1;
        }

        ngroup=index;
        groupPtr[index]=n;
    }





}