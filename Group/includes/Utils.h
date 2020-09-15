//
// Created by labuser (Bangtian Liu) on 9/15/20.
//

#ifndef LBC_LIB_UTILS_H
#define LBC_LIB_UTILS_H

namespace group_cols
{
    // Makes an edge inside dependence graph
    inline void connect(int v, int w, std::vector<std::vector<int>> &DAG){
        DAG[v].push_back( w );
    }

    int buildLevelSet_CSC_Queue(int n, int nnz, int *Lp, int *Li, int *&levelPtr,
                                int *&levelSet)
    {
        int begin=0,end=n-1;
        int curLevel=0, curLevelCol=0;
//    levelPtr = new int[n+1]();
//    levelSet = new int[n]();;
        int *inDegree = (int *)malloc(sizeof(int)*n);
        memset(inDegree, 0, sizeof(int)*n);
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
            if (inDegree[i]==1){
                dq.push_back(i);
            }
        }

        while (dq.size()!=0)
        {
            int len = dq.size();
            for (int i = 0; i < len; ++i) {
                int idx = dq.front();
                dq.pop_front();
                levelSet[curLevelCol] = idx; //add it to current level
                curLevelCol++;//Adding to level-set
            }

            curLevel++;//all nodes with zero indegree are processed.
            levelPtr[curLevel]=curLevelCol;

            for (int l = levelPtr[curLevel-1]; l < levelPtr[curLevel]; ++l) {
                int cc=levelSet[l];
                for (int j = Lp[cc] + 1; j < Lp[cc + 1]; ++j) {
                    inDegree[Li[j]]--;//removing corresponding edges
                    if(inDegree[Li[j]]==1){
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


}
#endif //LBC_LIB_UTILS_H
