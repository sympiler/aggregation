//
// Created by labuser (Bangtian Liu) on 9/16/20.
//

#include <Utils.h>
#include <algorithm>
#include <cstdio>

namespace sym_lib
{

     void fs_csr_inspector_dep(int ngroup, int *groupPtr, int *groupSet, int *gInv, int *Lp, int *Li, std::vector<std::vector<int>> &DAG)
     {
          //#pragma omp parallel for
          for (int i = 0; i < ngroup; ++i)
          {
               for (int j = groupPtr[i]; j < groupPtr[i + 1]; j++)
               {
                    int iidx = groupSet[j];
                    for (int k = Lp[iidx]; k < Lp[iidx + 1] - 1; ++k)
                    {
                         auto sid = gInv[Li[k]];
                         //                long int idx = i*(long int)ngroup+sid;
                         if (sid != i)
                         {
                              connect(sid, i, DAG);
                         }
                    }
               }
          }
     }

     void fs_ic0csc_inspector_group(int ngroup, int *groupPtr, int *groupSet, int *gInv, int *Lp, int *Li, std::vector<std::vector<int>> &DAG)
     {
          //Inspector
          #pragma omp parallel for
          for (int i = 0; i < ngroup; ++i)
          {
               for (int j = groupPtr[i]; j < groupPtr[i + 1]; j++)
               {
                    int iidx = groupSet[j];
                    for (int k = Lp[iidx] + 1; k < Lp[iidx + 1]; ++k)
                    {
                         auto sid = gInv[Li[k]];
                         if (sid != i)
                         {
                              connect(i, sid, DAG);
                         }
                    }
               }
          }
     }

     void fs_csr_inspector_dep(int n, int *Lp, int *Li, std::vector<std::vector<int>> &DAG)
     {
          //#pragma omp parallel for
          for (int i = 0; i < n; ++i)
          {
               for (int k = Lp[i]; k < Lp[i + 1] - 1; ++k)
               {
                    auto sid = Li[k];
                    if (sid != i)
                    {
                         connect(sid, i, DAG);
                    }
               }
          }
     }

     std::vector<std::vector<int>> Group_DAG(std::vector<std::vector<int>> DAG, int *groupPtr, int *groupSet, int *groupInv, int ngroup)
     {
          std::vector<std::vector<int>> gDAG;
          //        gDAG.resize(ngroup);
          for (int i = 0; i < ngroup; ++i)
          {
               std::vector<int> tarray;
               for (int j = groupPtr[i]; j < groupPtr[i + 1]; ++j)
               {
                    int idx = groupSet[j];
                    for (auto &idy : DAG[idx])
                    {
                         if (groupInv[idy] != i)
                              tarray.push_back(groupInv[idy]);
                    }
                    //                tarray.insert(tarray.end(), DAG[idx].begin(), DAG[idx].end());
               }
               std::sort(tarray.begin(), tarray.end());
               tarray.erase(std::unique(tarray.begin(), tarray.end()), tarray.end());
               gDAG.push_back(tarray);
          }
          return gDAG;
     }

     void detectDAGCircle(std::vector<std::vector<int>> DAG)
     {
          for (int i = 0; i < DAG.size(); ++i)
          {
               for (int j = 0; j < DAG[i].size(); ++j)
               {
                    int src = i;
                    int tar = DAG[i][j];
                    auto pos = std::find(DAG[tar].begin(), DAG[tar].end(), src);
                    if (pos != DAG[tar].end())
                    {
                         printf("circle between %d and %d\n", src, tar);
                    }
               }
          }
     }

} // namespace sym_lib
