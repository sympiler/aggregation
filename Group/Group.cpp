//
// Created by labuser (Bangtian Liu) on 9/16/20.
//

#include <Group.h>

namespace sym_lib
{

 void group::inspection_spicocsc_v1(int *groupPtr, int *groupSet, int &ngroup, int *groupInv) {
     for (int i = 0; i < ncol; ++i) {
         int len = mP[i + 1] - mP[i];
         child[i] = len > 1 ? mI[mP[i] + 1] : -1;
     }

     for (int i = 0; i < ncol; ++i) {
         int curIdx = i;
         int chIdx = child[curIdx];
         if(chIdx!=-1&&chIdx-curIdx==1){
             next[curIdx]=chIdx;
             pre[chIdx]=curIdx;
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

 void group::inspection_spicocsc_v2(int *groupPtr, int *groupSet, int &ngroup, int *groupInv) {
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
