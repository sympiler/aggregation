//
// Created by labuser (Bangtian Liu) on 9/16/20.
//

#include <Utils.h>
#include <algorithm>

namespace group_cols {

    int buildLevelSet_CSC_Queue(int n, int nnz, int *Lp, int *Li, int *&levelPtr,
                                int *&levelSet) {
        int begin = 0, end = n - 1;
        int curLevel = 0, curLevelCol = 0;
//    levelPtr = new int[n+1]();
//    levelSet = new int[n]();;
        int *inDegree = (int *) malloc(sizeof(int) * n);
        memset(inDegree, 0, sizeof(int) * n);
//            new int[n]();
//    bool *visited = (bool *)malloc(sizeof(bool)*n);
//    memset(visited, 0, sizeof(bool)*n);

//            new bool[n]();
        for (int i = 0; i < Lp[n]; ++i) {//O(nnz)
            inDegree[Li[i]]++;
        }

#if 0
        for (int k = 0; k < n; ++k) {
        std::cout<<k<<":"<<inDegree[k]<<",";
    }
    std::cout<<"\n";
#endif


        std::deque<int> dq;


        for (int i = 0; i < n; ++i) {
            if (inDegree[i] == 1) {
                dq.push_back(i);
            }
        }

        while (dq.size() != 0) {
            int len = dq.size();
            for (int i = 0; i < len; ++i) {
                int idx = dq.front();
                dq.pop_front();
                levelSet[curLevelCol] = idx; //add it to current level
                curLevelCol++;//Adding to level-set
            }

            curLevel++;//all nodes with zero indegree are processed.
            levelPtr[curLevel] = curLevelCol;

            for (int l = levelPtr[curLevel - 1]; l < levelPtr[curLevel]; ++l) {
                int cc = levelSet[l];
                for (int j = Lp[cc] + 1; j < Lp[cc + 1]; ++j) {
                    inDegree[Li[j]]--;//removing corresponding edges
                    if (inDegree[Li[j]] == 1) {
                        dq.push_back(Li[j]);
                    }
                }
            }
        }


        free(inDegree);
        return curLevel;
    }


    void fs_csr_inspector_dep(int ngroup, int *groupPtr, int *groupSet, int *gInv, int *Lp, int *Li, std::vector<std::vector<int>> &DAG)
    {
//#pragma omp parallel for
        for (int i = 0; i < ngroup; ++i) {
            for(int j=groupPtr[i]; j<groupPtr[i+1]; j++)
            {
                int iidx = groupSet[j];
                for (int k = Lp[iidx]; k < Lp[iidx+1]-1; ++k) {
                    auto sid = gInv[Li[k]];
//                long int idx = i*(long int)ngroup+sid;
                    if(sid!=i ){
                        connect(sid,i, DAG);
                    }
                }
            }
        }
    }


    void fs_csr_inspector_dep(int n, int *Lp, int *Li, std::vector<std::vector<int>> &DAG)
    {
//#pragma omp parallel for
        for (int i = 0; i < n; ++i) {
                for (int k = Lp[i]; k < Lp[i+1]-1; ++k) {
                    auto sid = Li[k];
                    if(sid!=i ){
                        connect(sid,i, DAG);
                    }
                }
        }
    }

    std::vector<std::vector<int>> Group_DAG(std::vector<std::vector<int>> DAG, int *groupPtr, int *groupSet, int *groupInv, int ngroup)
    {
        std::vector<std::vector<int>> gDAG;
//        gDAG.resize(ngroup);
        for (int i = 0; i < ngroup; ++i) {
            std::vector<int> tarray;
            for (int j = groupPtr[i]; j < groupPtr[i+1]; ++j) {
                int idx = groupSet[j];
                for(auto &idy: DAG[idx])
                {
                    if(groupInv[idy]!=i)
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

    void rhsInit_csr(int n, int *Ap, int *Ai, double *Ax, double *b)
    {
        /*generating a rhs that produces a result of all 1 vector*/
        for (int j = 0; j < n; ++j) {
            b[j]=0;
        }
        for (int c = 0; c < n ; ++c) {
            for (int cc = Ap[c]; cc < Ap[c + 1]; ++cc) {
                b[c]+=Ax[cc];
            }
        }
    }

    bool detectDAGCircle(std::vector<std::vector<int>> DAG)
    {
        for (int i = 0; i < DAG.size(); ++i) {
            for(int j=0; j < DAG[i].size(); ++j){
                int src = i;
                int tar = DAG[i][j];
                auto pos = std::find(DAG[tar].begin(), DAG[tar].end(), src);
                if(pos!=DAG[tar].end()){
                    printf("circle between %d and %d\n", src, tar);
                }

            }
        }
    }


}