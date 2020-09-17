//
// Created by labuser (Bangtian Liu) on 9/16/20.
//

#include <Group.h>

namespace group_cols
{
    void group::groupInspection(int *groupPtr, int *groupSet, int &ngroup, int *groupInv) {
        for (int i = 0; i < ncol; ++i) {
            int curIdx = i;
            int len = mP[i+1] - mP[i];
            int chIdx = len>1 ? mI[mP[i]+1] : -1;

            if(chIdx!=-1) {
                if(status[curIdx]==0 && status[chIdx]==0){
                    next[curIdx]=chIdx;
                    pre[chIdx]=i;
                    status[chIdx]=1;
                }
                else if (status[chIdx]==1 && status[curIdx]==0){
                    next[curIdx] = chIdx;
                    int tp = pre[chIdx];
                    pre[chIdx]=curIdx;
                    pre[curIdx]=tp;
                    next[tp]=curIdx;
                    status[curIdx]==2;
                }
                else if (status[curIdx]==1 && status[chIdx]==0 && (chIdx-curIdx==1)) // first consecutive column
                {
                    next[curIdx] = chIdx;
                    pre[chIdx]=curIdx;
                    status[chIdx]=3;
                }
                else if (status[curIdx]=3 && status[chIdx]==0 && (chIdx-curIdx)==1) // consecutive grouping column
                {
                    next[curIdx] = chIdx;
                    pre[chIdx]=curIdx;
                    status[chIdx]=3;
                }
//                else if(status[chIdx]==3 && status[curIdx]==0){
//                    next[curIdx] = chIdx;
//                    int tp = pre[chIdx];
//                    pre[chIdx]=curIdx;
//                    pre[curIdx]=tp;
//                    next[tp]=curIdx;
//                    status[curIdx]==2;
//                }
            }
        }

//        for (int i = 0; i < ncol; ++i) {
//            printf("i=%d, %d\n", i, next[i]);
//        }



        int pos=0; int index=0;
        for (int i = 0; i < ncol; ++i) {
            int curIdx = i;
            if (!visited[curIdx]) {

                groupPtr[index]=pos;
                visited[curIdx]=true;
                groupSet[pos++]=curIdx;
                groupInv[curIdx] = index;
                int ch = next[curIdx];
                while (ch!=0)
                {
//                    printf("pos=%d  %d %d\n", pos, curIdx, ch);
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


    void group::groupInspection_v1(int *groupPtr, int *groupSet, int &ngroup, int *groupInv) {

//#pragma omp parallel for
        for (int i = 0; i < ncol; ++i) {
            int len = mP[i + 1] - mP[i];
            child[i] = len > 1 ? mI[mP[i] + 1] : -1;
        }

        for (int i = 0; i < ncol; ++i) {
            int curIdx=i;


            int chIdx=child[curIdx];

            if(chIdx!=-1 && !visited[curIdx] ){
                visited[curIdx]= true;
                if(status[curIdx]==0 && status[chIdx]==0){

                    next[curIdx]=chIdx;
                    pre[chIdx]=curIdx;
                    status[chIdx]=1;


                    curIdx = chIdx;
                    visited[curIdx]= true;
                    chIdx = child[curIdx];
                    while (chIdx!=-1 && (chIdx-curIdx)==1 && !visited[chIdx])
                    {
                        visited[chIdx]=true;
                        pre[chIdx] = curIdx;
                        next[curIdx] = chIdx;
                        status[chIdx] = 3;
                        curIdx=chIdx;
                        chIdx=child[curIdx];
                    }

                    int tcurIdx = curIdx+1;
                    int len=0;
                    while (child[tcurIdx]==-1 && !visited[tcurIdx] && len<ncol/12)
                    {
                        visited[tcurIdx] = true;
                        next[tcurIdx-1]=tcurIdx;
                        pre[tcurIdx]= tcurIdx-1;
                        status[tcurIdx]=3;
                        ++tcurIdx;
                        ++len;
                    }
                }
                else if (status[chIdx]==1 && status[curIdx]==0){
                    visited[chIdx]=true;
                    next[curIdx] = chIdx;
                    int tp = pre[chIdx];
                    pre[chIdx]=curIdx;
                    pre[curIdx]=tp;
                    next[tp]=curIdx;
                    status[curIdx]==2;
                }
                else if(status[chIdx]==3 && status[curIdx]==0){
                    visited[chIdx]=true;
                    next[curIdx] = chIdx;
                    int tp = pre[chIdx];
                    pre[chIdx]=curIdx;
                    pre[curIdx]=tp;
                    next[tp]=curIdx;
                    status[curIdx]==2;
                }
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


    void group::groupInspection_v2(int *groupPtr, int *groupSet, int &ngroup, int *groupInv) {
        int *child =(int *)malloc(sizeof(int)*ncol);
        memset(child, 0, sizeof(int)*ncol);

//#pragma omp parallel for
        for (int i = 0; i < ncol; ++i) {
            int len = mP[i + 1] - mP[i];
            child[i] = len > 1 ? mI[mP[i] + 1] : -1;
            top[i] = -1;
        }

        for (int i = 0; i < ncol; ++i) {
            int curIdx=i;


            int chIdx=child[curIdx];

            if(chIdx!=-1 && !visited[curIdx] ){
                visited[curIdx]= true;
//                ++len;
                if(status[curIdx]==0 && status[chIdx]==0){

                    next[curIdx]=chIdx;
                    pre[chIdx]=curIdx;
                    status[chIdx]=1;
                    top[chIdx]=curIdx;
                    int root = curIdx;

//                    ++len;

                    curIdx = chIdx;
                    visited[curIdx]= true;
                    chIdx = child[curIdx];
                    while (chIdx!=-1 && (chIdx-curIdx)==1 && !visited[chIdx] )
                    {
                        visited[chIdx]=true;
                        pre[chIdx] = curIdx;
                        next[curIdx] = chIdx;
                        status[chIdx] = 3;
                        top[chIdx]=root;

//                        ++len;
                        curIdx=chIdx;
                        chIdx=child[curIdx];
                    }

                    int tcurIdx = curIdx+1;
                    int tlen=0;
                    while (child[tcurIdx]==-1 && !visited[tcurIdx] && tlen<10)
                    {
                        visited[tcurIdx] = true;
                        next[tcurIdx-1]=tcurIdx;
                        pre[tcurIdx]= tcurIdx-1;
                        status[tcurIdx]=3;
                        top[tcurIdx] = root;

                        ++tcurIdx;
                        ++tlen;
                    }
                }
                else if (status[chIdx]==1 && status[curIdx]==0){
                    visited[chIdx]=true;

//                    if(pre[chIdx]==top[chIdx]){
//                        next[curIdx] = chIdx;
//                    }
//                    else{
//
//                    }


                    next[curIdx] = chIdx;
                    int tp = pre[chIdx];
                    pre[chIdx]=curIdx;
                    pre[curIdx]=tp;
                    next[tp]=curIdx;
                    status[curIdx]==2;
                }
                else if(status[chIdx]==3 && status[curIdx]==0){
                    visited[chIdx]=true;

                    int tidx = top[chIdx];
                    while (next[tidx]<curIdx){
                        tidx=next[tidx];
                    }

                    tidx = next[tidx];
                    next[curIdx] = tidx;
                    int tp = pre[tidx];
                    pre[tidx] = curIdx;
                    pre[curIdx] = tp;
                    next[tp]=curIdx;


//                    next[curIdx] = chIdx;
//                    int tp = pre[chIdx];
//                    pre[chIdx]=curIdx;
//                    pre[curIdx]=tp;
//                    next[tp]=curIdx;
                    status[curIdx]==2;
                }
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

    void group::nogroup(int *groupPtr, int *groupSet, int &ngroup, int *groupInv) {
        int pos=0; int index=0;

        for (int i = 0; i < ncol; ++i) {
            groupPtr[i]=pos;
            groupSet[pos++]=i;
            groupInv[i] = i;
        }
        groupPtr[ncol]=pos;
        ngroup=ncol;
    }

    void group::inspection_sptrsvcsr(int *groupPtr, int *groupSet, int &ngroup, int *groupInv) {
        for (int i = 0; i < ncol; ++i) {
            int len = mP[i + 1] - mP[i];
            child[i] = len > 1 ? mI[mP[i+1] - 2] : -1;
        }

        for (int i = ncol-1; i >=1 ; --i) {
            int curIdx = i;
            int chIdx = child[curIdx];


            if(chIdx!=-1 && !visited[curIdx]){
                visited[curIdx]= true;

                if(status[curIdx]==0 && status[chIdx]==0){
//                        next[chIdx]=curIdx;
//                        pre[curIdx]=chIdx;
//                        status[chIdx]=1;
//
//                        curIdx=chIdx;
//                        visited[curIdx]= true;
                    chIdx = child[curIdx];

                    while (chIdx!=-1 && (curIdx-chIdx)==1 && !visited[chIdx])
                    {
                        visited[chIdx]=true;

                        next[chIdx] = curIdx;
                        pre[curIdx]=chIdx;

                        status[chIdx] = 3;
                        curIdx=chIdx;
                        chIdx=child[curIdx];
                    }

                    int tcurIdx = curIdx-1;
                    int len=0;
                    while (child[tcurIdx]==-1 && !visited[tcurIdx] && len<ncol/12)
                    {
                        visited[tcurIdx] = true;
                        next[tcurIdx] = tcurIdx+1;
                        pre[tcurIdx+1] = tcurIdx;
                        status[tcurIdx]=3;
                        --tcurIdx;
                    }
                }
                else if (status[chIdx]==1 && status[curIdx]==0)
                {
                    visited[chIdx]=true;
//                        next[chIdx] = curIdx;
                    int tp = next[chIdx];
                    next[chIdx]=curIdx;
                    pre[curIdx]=chIdx;
                    next[curIdx]=tp;
                    pre[tp]=curIdx;
                    status[curIdx]=2;
                }
                else if(status[chIdx]==3 && status[curIdx]==0)
                {
                    visited[chIdx]=true;
                    int tp = next[chIdx];

                    next[chIdx]=curIdx;
                    pre[curIdx]=chIdx;
                    next[curIdx]=tp;
                    pre[tp]=curIdx;

                    status[curIdx]=2;
                }
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


    void NaiveGrouping(int n, int *groupPtr, int *groupSet,  int &ngroup, int *ginv, int blksize)
    {
//        int blksize=1;
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