//
// Created by behrooz on 2021-07-11.
//

#include "aggregation/hdagg.h"
#include <unordered_set>
#include "aggregation/sparse_utilities.h"
#include "aggregation/lbc.h"
#include <list>
#include <cmath>
#include <set>


namespace HDAGG
{
    int build_levelSet_CSC(size_t n, const int *Lp, const int *Li,
                           int *levelPtr, int *levelSet, int* node_to_level)
   {
        int begin = 0, end = n - 1;
        int cur_level = 0, cur_levelCol = 0;
        int *inDegree = new int[n]();
        bool *visited = new bool[n]();
        for (int i = 0; i < Lp[n]; ++i)
        { //O(nnz) -> but not catch efficient. This code should work well
            // on millions of none zeros to enjoy a gain in the parformance. Maybe we can find another way. Behrooz
            inDegree[Li[i]]++; // Isn't it the nnz in each row? or the rowptr[x + 1] - rowptr[x] in CSR?
        }
        //print_vec("dd\n",0,n,inDegree);
        while (begin <= end)
        {
            for (int i = begin; i <= end; ++i)
            { //For level cur_level
                if (inDegree[i] == 1 && !visited[i])
                {   //if no incoming edge
                    visited[i] = true;
                    levelSet[cur_levelCol] = i; //add it to current level
                    node_to_level[i] = cur_level;
                    cur_levelCol++;				//Adding to level-set - This is a cnt for the current level. Behrooz
                }
            }
            cur_level++; //all nodes_ with zero indegree are processed.
            //assert(cur_level < n);
            if (cur_level >= n)
                return -1; // The input graph has a cycle
                levelPtr[cur_level] = cur_levelCol; // The levelPtr starts from level 1. Behrooz
                while (inDegree[begin] == 1) // Why? Behrooz
                    {
                    begin++;
                    if (begin >= n)
                        break;
                    }
                while (inDegree[end] == 1 && begin <= end) // The same why as above. Behrooz
                    end--;
                //Updating degrees after removing the nodes_
                for (int l = levelPtr[cur_level - 1]; l < levelPtr[cur_level]; ++l) // I don't get this part. Behrooz
                    {
                    int cc = levelSet[l];
                    for (int j = Lp[cc]; j < Lp[cc + 1]; ++j)
                    {
                        if (Li[j] != cc)	   //skip diagonals
                            inDegree[Li[j]]--; //removing corresponding edges
                    }
                    }
                //print_vec("dd\n",0,n,inDegree);
        }
        delete[] inDegree;
        delete[] visited;
        return cur_level; //return number of level
    }

   void computePostOrderAndGroup(int n, const std::vector<int>& DAG_ptr, const std::vector<int>& DAG_set,
                                 int& ngroups, std::vector<int>& groupPtr, std::vector<int>& post_order,
                                 std::vector<int>& group_DAG_ptr, std::vector<int>& group_DAG_set){
        //TODO: Find a way to merge these functions
        //Find the leafs
        std::vector<int> in_degree(n, 0);

        for(int i = 0; i < n; i++){
            for(int child_ptr = DAG_ptr[i] + 1; child_ptr < DAG_ptr[i + 1]; child_ptr++){
                in_degree[DAG_set[child_ptr]]++;
            }
        }

        std::vector<int> leaf_id;
        for(int i = 0; i < n; i++){
            if(in_degree[i] == 0){
                leaf_id.push_back(i);
            }
        }
        //Apply dfs on the leaf in reverse order (from big to small)
        post_order.resize(n);
        int order_cnt = 0;
        std::vector<bool> visited(n, false);
        std::vector<int> stack;
        for(int leaf_ptr = 0; leaf_ptr < leaf_id.size(); leaf_ptr++){
            int leaf = leaf_id[leaf_id.size() - leaf_ptr - 1];
            //Apply DFS on each leaf
            int head = leaf;
            visited[leaf] = true;
            stack.push_back(leaf);
            while(!stack.empty()){
                bool all_parent_visited = true;
                for(int child_ptr = DAG_ptr[head + 1] - 1; child_ptr > DAG_ptr[head] ; child_ptr--){
                    int child = DAG_set[child_ptr];
                    if(!visited[child]){
                        visited[child] = true;
                        stack.push_back(child);
                        head = child;
                        all_parent_visited = false;
                        break;
                    }
                }
                if(all_parent_visited){
                    post_order[n - 1 - order_cnt] = stack.back();
                    order_cnt++;
                    stack.pop_back();
                    head = stack.back();
                }
            }
        }

        //Apply Grouping on the order
        assert(std::count(visited.begin(), visited.end(), false) == 0);
        groupPtr.reserve(n);
        groupPtr.push_back(0);
        std::vector<int> groupInv(n, 0);
        ngroups = 0;
        int nodes_cnt = 0;
        for(int i = 0; i < n; i++){
            int curr = post_order[i];
            int in_edges = 0;
            int next = 0;
            int child = 0;
            if(i < n - 1){
                next = post_order[i + 1];
                in_edges = in_degree[next];
                child = DAG_set[DAG_ptr[curr] + 1];//TODO: Ugly coding
            } else {
                in_edges = 1;
                child = next;
            }
            int out_edges = DAG_ptr[curr + 1] - DAG_ptr[curr] - 1;
            //TODO: Understand why post order needs the condition "AG_set[DAG_ptr[curr] + 1] == next"
            if(in_edges == out_edges && in_edges == 1 && child == next){
                nodes_cnt++;
                groupInv[curr] = ngroups;
            } else {
                nodes_cnt++;
                groupPtr.push_back(nodes_cnt);
                groupInv[curr] = ngroups;
                ngroups++;
            }
        }
        assert(groupPtr.back() == n);
        #ifndef NDEBUG
        //Check groups and make sure that they won't violate chain characteristics
        for(int group = 0; group < ngroups; group++){
            for(int node_ptr = groupPtr[group]; node_ptr < groupPtr[group + 1] - 1; node_ptr++){
                int curr = post_order[node_ptr];
                int next = post_order[node_ptr + 1];
                int out_degree = DAG_ptr[curr + 1] - DAG_ptr[curr] - 1;
                assert(in_degree[next] == out_degree == 1);
            }

            for(int node_ptr = groupPtr[group]; node_ptr < groupPtr[group + 1]; node_ptr++){
                assert(groupInv[post_order[node_ptr]] == group);
            }
        }
        #endif
        //Create the group DAG
        //Note that because of the special chain grouping, the real DAG just need a small modification
        group_DAG_ptr.push_back(0);
        group_DAG_set.reserve(DAG_ptr[n]);
        for(int group = 0; group < ngroups; group++){
            int node = post_order[groupPtr[group + 1] - 1];
            int tmp = group_DAG_ptr.back();
            group_DAG_ptr.push_back(tmp);
            for(int child = DAG_ptr[node]; child < DAG_ptr[node + 1]; child++){
                group_DAG_set.push_back(groupInv[DAG_set[child]]);
                group_DAG_ptr.back()++;
            }
        }
        assert(group_DAG_set.size() == group_DAG_ptr.back());
        assert(group_DAG_ptr.size() == (ngroups + 1));
        #ifndef NDEBUG
            //Aggressive DAG creation
            std::vector<std::vector<int>> DAG;
            DAG.resize(ngroups);

            for (int i = 0; i < ngroups; ++i)
            {
                for (int j = groupPtr[i]; j < groupPtr[i + 1]; j++)
                {
                    int iidx = post_order[j];
                    for (int k = DAG_ptr[iidx]; k < DAG_ptr[iidx + 1]; ++k)
                    {
                        auto sid = groupInv[DAG_set[k]];
                        DAG[i].push_back( sid );
                    }
                }
            }

            size_t count=0;
            for(int j = 0; j < DAG.size(); ++j) {
                DAG[j].erase(std::unique(DAG[j].begin(), DAG[j].end()), DAG[j].end());
                for(int i = 0; i < DAG[j].size(); i++){
                    assert(group_DAG_set[group_DAG_ptr[j] + i] == DAG[j][i]);
                }
            }

        #endif

    }

    void computePostOrder(int n, const int* DAG_ptr, const int* DAG_set,
                          std::vector<int>& post_order){
        //TODO: Find a way to merge these functions
        //Find the leafs
        std::vector<int> in_degree(n, 0);

        for(int i = 0; i < n; i++){
            for(int child_ptr = DAG_ptr[i] + 1; child_ptr < DAG_ptr[i + 1]; child_ptr++){
                in_degree[DAG_set[child_ptr]]++;
            }
        }

        std::vector<int> leaf_id;
        for(int i = 0; i < n; i++){
            if(in_degree[i] == 0){
                leaf_id.push_back(i);
            }
        }
        //Apply dfs on the leaf in reverse order (from big to small)
        post_order.resize(n);
        int order_cnt = 0;
        std::vector<bool> visited(n, false);
        std::vector<int> stack;
        for(int leaf_ptr = 0; leaf_ptr < leaf_id.size(); leaf_ptr++){
            int leaf = leaf_id[leaf_id.size() - leaf_ptr - 1];
            //Apply DFS on each leaf
            int head = leaf;
            visited[leaf] = true;
            stack.push_back(leaf);
            while(!stack.empty()){
                bool all_parent_visited = true;
                for(int child_ptr = DAG_ptr[head + 1] - 1; child_ptr > DAG_ptr[head] ; child_ptr--){
                    int child = DAG_set[child_ptr];
                    if(!visited[child]){
                        visited[child] = true;
                        stack.push_back(child);
                        head = child;
                        all_parent_visited = false;
                        break;
                    }
                }
                if(all_parent_visited){
                    post_order[n - 1 - order_cnt] = stack.back();
                    order_cnt++;
                    stack.pop_back();
                    head = stack.back();
                }
            }
        }
    }

    void parallel_cc(int n, const int *lC, const int *lR, int lbLevel, int ubLevel,
                     const int *node2Level, int numNodes, const int *nodes,
                     int *node2partition) {

        //Creating Stars
    #pragma omp parallel for
    for (int i = 0; i < numNodes; ++i) {
        int node = nodes[i];
        node2partition[node] = i;
    }

    bool change = true;
    while (change) {
        change = false;
        #pragma omp parallel for
        for (int i = 0; i < numNodes; ++i) {
            int u = nodes[i];
            // Now we go over all neighbors of u
            for (int r = lC[u]; r < lC[u + 1]; ++r) {
                int v = lR[r];

                int level = node2Level[v];
                if (level < lbLevel || level >= ubLevel)
                    continue;

                int comp_u = node2partition[u];
                int comp_v = node2partition[v];
                if (comp_u == comp_v)
                    continue;
                int high_comp = comp_u > comp_v ? comp_u : comp_v;
                int low_comp = comp_u + (comp_v - high_comp);
                if (high_comp == node2partition[nodes[high_comp]]) {
                    change = true;
                    node2partition[nodes[high_comp]] = low_comp;
                }
            }
        }

        #pragma omp parallel for
        for (int i = 0; i < numNodes; ++i) {
            int node = nodes[i];
            while (node2partition[node] != node2partition[nodes[node2partition[node]]]) {
                node2partition[node] = node2partition[nodes[node2partition[node]]];
            }
        }
    }
    }


    sym_lib::timing_measurement computingNode2Level(const std::vector<int>& LevelPtr, const std::vector<int>& LevelSet, int nlevels,
                                           int num_nodes, std::vector<int>& Node2Level){
        sym_lib::timing_measurement Node2Level_time;
        Node2Level_time.start_timer();
        Node2Level.resize(num_nodes);
        for (int lvl = 0; lvl < nlevels; ++lvl) {
            for (int j = LevelPtr[lvl]; j < LevelPtr[lvl + 1]; ++j) {
                int node = LevelSet[j];
                Node2Level[node] = lvl;
            }
        }

        Node2Level_time.measure_elapsed_time();
        return Node2Level_time;
    }


    void computeSchedule(int num_threads,
                         int n,
                         const int* DAG_ptr,
                         const int* DAG_set,
                         const int* level_ptr,
                         const int* level_set,
                         const int* node_to_level,
                         std::vector<int>& WM,
                         std::vector<std::vector<int>> Chunk_Sizes,
                         int& coarse_level_no,
                         std::vector<int>& coarse_level_ptr,
                         std::vector<int>& coarse_part_ptr,
                         std::vector<int>& coarse_node_ptr,
                         bool apply_serial_schedule,
                         bool postOrder){
        //TODO: Delete this part
        if(!coarse_level_ptr.empty()){
            coarse_level_ptr.clear();
        }

        if(!coarse_part_ptr.empty()){
            coarse_part_ptr.clear();
        }

        if(!coarse_node_ptr.empty()){
            coarse_node_ptr.clear();
        }

        assert(coarse_level_ptr.empty());
        assert(coarse_part_ptr.empty());
        assert(coarse_node_ptr.empty());


        coarse_level_no = WM.size() - 1;
        coarse_level_ptr.push_back(0);
        coarse_part_ptr.push_back(0);
        std::vector<int> post_order(n);
        std::vector<int> post_orderInv(n);
        std::vector<int> node_ptr_order(n);
        if(postOrder){
            computePostOrder(n, DAG_ptr, DAG_set, post_order);
            for(int i = 0; i < post_order.size(); i++){
                post_orderInv[post_order[i]] = i;
            }
        }
        for(int wm = 0; wm < WM.size() - 1; wm++){
            //Finding CCs
            std::vector<int> nodes_at_cur_level(n, 0);
            int numNodesAtCurLevel = 0;
            std::vector<int> node2partition(n, -1);
            int lb_level = WM[wm];
            int ub_level = WM[wm + 1];
            for (int ii = lb_level; ii < ub_level; ++ii) {
                for (int j = level_ptr[ii]; j < level_ptr[ii + 1]; ++j) {
                    int x = level_set[j];
                    nodes_at_cur_level[numNodesAtCurLevel++] = x;
                }
            }

            parallel_cc(n, DAG_ptr, DAG_set,
                        WM[wm], WM[wm + 1],
                        node_to_level,
                        numNodesAtCurLevel,
                        nodes_at_cur_level.data(),
                        node2partition.data());

            //Sorting the nodes based on their w-partitions
            std::vector<int> w_part_size_cnt(n + 1, 0);
            for(auto & iter: node2partition){
                if(iter != -1){
                    w_part_size_cnt[iter + 1]++;
                }
            }
            //Scan the counter vector
            for(int i = 1; i < w_part_size_cnt.size() - 1; i++){
                w_part_size_cnt[i + 1] += w_part_size_cnt[i];
            }

            std::vector<int> part_ptr;
            std::vector<int> node_ptr;
            //create part_ptr
            part_ptr.push_back(0);
            for(int i = 0; i < w_part_size_cnt.size() - 1; i++){
                if(w_part_size_cnt[i + 1] != w_part_size_cnt[i]){
                    part_ptr.push_back(w_part_size_cnt[i + 1]);
                }
            }
            //create node_ptr
            node_ptr.resize(numNodesAtCurLevel);
            for(int i = 0; i < numNodesAtCurLevel; i++){
                node_ptr[w_part_size_cnt[node2partition[nodes_at_cur_level[i]]]++] = nodes_at_cur_level[i];
            }

            //Apply specific serial schedule
            if(apply_serial_schedule){
                //Compute schedule chunk size
                std::vector<int> & chunk_array = Chunk_Sizes[wm];
                //Merge core's parts and apply schedule inside the merged parts
                std::vector<int> part_ptr_merged(num_threads + 1, 0);
                assert(chunk_array.back() == part_ptr.size() - 1);
                for(int i = 0; i < num_threads; i++){
                    size_t tmp = chunk_array[i + 1];
                    int maximum_part_idx = std::min(tmp, part_ptr.size() - 1);
                    part_ptr_merged[i + 1] = part_ptr[maximum_part_idx];
                }
                //Place arrays inside the coarse arrays
                auto tmp = coarse_level_ptr.back();
                assert(!part_ptr.empty());
                coarse_level_ptr.push_back(tmp + part_ptr_merged.size() - 1);
                auto last_part = coarse_part_ptr.back();
                for(int p = 1; p < part_ptr_merged.size(); p++){
                    auto tmp1 = part_ptr_merged[p] + last_part;
                    coarse_part_ptr.push_back(tmp1);
                }
                //Add the new nodes to the coarse node array
                if(postOrder){
                    for(int i = 0; i < part_ptr_merged.size() - 1; i++){
                        //TODO: Fix this copy later for faster code, it is too slow
                        int cnt = 0;
                        for(int j = part_ptr_merged[i]; j < part_ptr_merged[i + 1]; j++){
                            node_ptr_order[cnt++] = post_orderInv[node_ptr[j]];
                        }
                        std::sort(node_ptr_order.data(), node_ptr_order.data() + cnt);
                        cnt = 0;
                        for(int j = part_ptr_merged[i]; j < part_ptr_merged[i + 1]; j++){
                            node_ptr[j] = post_order[node_ptr_order[cnt++]];
                        }
                    }
                } else {
                    //std::cout << "Sorting" << std::endl;
                    for(int i = 0; i < part_ptr_merged.size() - 1; i++){
                        std::sort(node_ptr.data() + part_ptr_merged[i], node_ptr.data() + part_ptr_merged[i + 1]);
                    }
                }
            } else {
//                std::cout << "No Bin packing" << std::endl;
                //Place arrays inside the coarse arrays
                //TODO: use the coarse array directly
                auto tmp = coarse_level_ptr.back();
                assert(!part_ptr.empty());
                coarse_level_ptr.push_back(tmp + part_ptr.size() - 1);
                auto last_part = coarse_part_ptr.back();
                for(int p = 1; p < part_ptr.size(); p++){
                    auto tmp1 = part_ptr[p] + last_part;
                    coarse_part_ptr.push_back(tmp1);
                }
                //Add the new nodes to the coarse node array
            }
            coarse_node_ptr.insert(coarse_node_ptr.end(), node_ptr.begin(), node_ptr.end());
        }
    }



    void computeCC(int lb_level, int computed_up_to,
                   int ub_level,
                   int n,
                   const int* DAG_ptr,
                   const int* DAG_set,
                   const int* level_ptr,
                   const int* level_set,
                   const int* node_to_level,
                   std::vector<int>& node2partition,
                   int& part_cnt,
                   std::vector<double>& partition_cost,
                   const std::vector<double>& node_cost,
                   double& total_cost){


        const int* sub_DAG_nodes = level_set + level_ptr[lb_level];
        int num_nodes_processed = level_ptr[computed_up_to] - level_ptr[lb_level];
        int num_new_nodes = level_ptr[ub_level] - level_ptr[computed_up_to];
        int num_nodes_in_sub_DAG = num_new_nodes + num_nodes_processed;
        part_cnt = 0;
        // If there is only one level
        if(ub_level - lb_level == 1){
            part_cnt = level_ptr[ub_level] - level_ptr[lb_level];
            #pragma omp parallel for
            for(int node_ptr = level_ptr[lb_level]; node_ptr < level_ptr[ub_level]; node_ptr++){
                int node = level_set[node_ptr];
                node2partition[node] = node_ptr - level_ptr[lb_level];
                partition_cost[node_ptr - level_ptr[lb_level]] = node_cost[node];
                #pragma omp atomic
                total_cost += node_cost[node];
            }
        } else {
            bool change;
            std::vector<int> node2partition_map(n, 0);
            std::vector<int> part_exist(n, 0);
            #pragma omp parallel
            {
                //If there are multiple levels - Note that we have a bit of reusing right now
                #pragma omp for
                for (int i = num_nodes_processed; i < num_nodes_in_sub_DAG; ++i) {
                    int node = sub_DAG_nodes[i];
                    if(node2partition[node] == -1){
                        node2partition[node] = i;
                    } else {
                        std::cerr << "Somthing is wrong with your boundaries" << std::endl;
                    }
                }

                //Creating Stars
                #pragma omp single
                change = true;
                while (change) {
                    //Hooking
                    // If the change = false, then some of the thread may miss the while loop
                    //TODO change the while loop mechanism
                    #pragma omp barrier
                    #pragma omp single
                    change = false;
                    #pragma omp for
                    for (int i = 0; i < num_nodes_in_sub_DAG; ++i) {
                        int u = sub_DAG_nodes[i];
                        // Now we go over all neighbors of u
                        for (int r = DAG_ptr[u]; r < DAG_ptr[u + 1]; ++r) {
                            int v = DAG_set[r];

                            int level = node_to_level[v];
                            if (level < lb_level || level >= ub_level)
                                continue;

                            int comp_u = node2partition[u];
                            int comp_v = node2partition[v];
                            if (comp_u == comp_v)
                                continue;
                            int high_comp = comp_u > comp_v ? comp_u : comp_v;
                            int low_comp = comp_u + (comp_v - high_comp);
                            //Change the root partition indicator to the low component
                            if (high_comp == node2partition[sub_DAG_nodes[high_comp]]) {
                                change = true;
                                node2partition[sub_DAG_nodes[high_comp]] = low_comp;
                            }
                        }
                    }

                    #pragma omp for
                    for (int i = 0; i < num_nodes_in_sub_DAG; ++i) {
                        int node = sub_DAG_nodes[i];
                        while (node2partition[node] != node2partition[sub_DAG_nodes[node2partition[node]]]) {
                            node2partition[node] = node2partition[sub_DAG_nodes[node2partition[node]]];
                        }
                    }
                }

                #pragma omp for
                for(int i = 0; i < node2partition.size(); i++){
                    int node_partition = node2partition[i];
                    if(node_partition != -1){
                        part_exist[node_partition] = 1;
                    }
                }

                //This implementation is serial
                //TODO: one can use scan and compact to make it parallel
                #pragma omp single
                for(int i = 0; i < part_exist.size(); i++){
                    if(part_exist[i]){
                        node2partition_map[i] = part_cnt++;
                    }
                }

                assert(num_nodes_in_sub_DAG != 0);
                #pragma omp for reduction(+: total_cost)
                for(int node_ptr = 0; node_ptr < num_nodes_in_sub_DAG; node_ptr++){
                    int node = sub_DAG_nodes[node_ptr];
                    #pragma omp atomic
                    partition_cost[node2partition_map[node2partition[node]]] += node_cost[node];
                    total_cost += node_cost[node];
                }
            }
        }
    }


    std::vector<int> HDAGG(int n, int nnz,
                           std::vector<int> & DAG_ptr_not_prune, std::vector<int> & DAG_set_not_prune,
                           const std::vector<double> & node_cost,
                           int cores,
                           int & coarse_level_no,
                           std::vector<int> & coarse_level_ptr,
                           std::vector<int> & coarse_part_ptr,
                           std::vector<int> & coarse_node_ptr,
                           bool parallelLevelset,
                           bool postOrder,
                           bool bin_pack){

        auto DAG_ptr = DAG_ptr_not_prune.data();
        auto DAG_set = DAG_set_not_prune.data();

        std::vector<int> level_ptr;
        std::vector<int> level_set;
        std::vector<int> node_to_level(n, 0);
        int nlevels;

        if(parallelLevelset){
#ifdef SpMP
            //Create the CSC version of the codeto use make_full and csc_to_csr codes
            CSC* Lower_A_CSC = new CSC(n, n, nnz);
            std::copy(DAG_ptr, DAG_ptr + n + 1, Lower_A_CSC->p);
            std::copy(DAG_set, DAG_set + nnz, Lower_A_CSC->i);
            std::fill(Lower_A_CSC->x, Lower_A_CSC->x + nnz, 1);
            Lower_A_CSC->stype=-1;
            auto tmp = sym_lib::make_full(Lower_A_CSC);
            auto CSR_A = sym_lib::csc_to_csr(tmp);
            auto A = new SpMP::CSR();
            HDAGG::Convert_LBC_CSR_to_SpMP(CSR_A, A);
            nlevels = HDAGG::levelsetCSRParallel_SpMP(A, cores, level_ptr, level_set);
            for(int l = 0; l < nlevels; l++){
                for(int node_ptr = level_ptr[l]; node_ptr < level_ptr[l + 1]; node_ptr++){
                    int node = level_set[node_ptr];
                    node_to_level[node] = l;
                }
            }
            delete Lower_A_CSC;
            delete tmp;
            delete A;
            delete CSR_A;
#endif
        } else {
            //Create the Levelset and compute node_to_level
            level_ptr.resize(n + 1);
            level_set.resize(n);

            nlevels = build_levelSet_CSC(n, DAG_ptr, DAG_set,
                                         level_ptr.data(), level_set.data(), node_to_level.data());
        }

        std::vector<int> WM;
        WM.push_back(0);
        std::vector<std::vector<int>> chunk_sizes(nlevels);
        //Compute the schedule
        for(int lvl = 0; lvl < nlevels; ) {
            int lb_level = lvl;
            int computed_up_to_lvl = lb_level;
            int prev_upper_bound = 0;
            int WM_condidate = lb_level + 1;
            int ub_level = std::min(lb_level + 1, nlevels);
            assert(ub_level <= nlevels && lb_level < nlevels);
            bool done_flag = false;
            //TODO: Clean this code later
            //Reusing material
            std::vector<int> node2partition(n, -1);
            double curr_DKL = 0;
            std::vector<int> chunk_array(cores + 1, 0);
            int iteration_counter = 0;
            std::vector<double> tmp;
            int effective_cores = cores;
            while(!done_flag){
                iteration_counter++;
                chunk_sizes[WM.size() - 1] = chunk_array;
                std::vector<double> partition_cost(n, 0);
                //===================================== Greedy Merge ===================================
                //Compute the cc's and w-partitions of these levels and compute the critical path
                int part_cnt = 0;
                double total_cost = 0;

                if( (ub_level - lb_level == 1) && (level_ptr[ub_level] - level_ptr[ub_level - 1] < cores)){
                    effective_cores = level_ptr[ub_level] - level_ptr[ub_level - 1];
                } else if(level_ptr[ub_level] - level_ptr[ub_level - 1] > effective_cores){
                    effective_cores = cores;
                }

                computeCC(lb_level, computed_up_to_lvl, ub_level,
                          n, DAG_ptr, DAG_set, level_ptr.data(), level_set.data(), node_to_level.data(),
                          node2partition, part_cnt, partition_cost, node_cost, total_cost);

                //Compute cost
                //TODO: This is one of the ugliest code that I have ever done, and it is the core of my work :)))))
                std::vector<double> cost(cores, 0);
                double cost_per_core = total_cost / cores;
                int thread_id = 1;
                double cost_to_this_core = 0;
                //TODO: Consider the cases where it has a narrow part_cnt and then expand
                if(part_cnt < cores){
                    for(int part = 0; part < part_cnt; part++){
                        chunk_array[part + 1] = part + 1;
                        cost[part] = partition_cost[part];
                    }
                    thread_id = part_cnt;
                } else {
                    for(int part = 0; part < part_cnt; ){
                        if(cost[thread_id - 1] < cost_per_core){
                            cost[thread_id - 1] += partition_cost[part++];
                            chunk_array[thread_id] = part;
                            assert(cores != thread_id - 1);//I had a small optimization for last core
                        } else {
                            //Cost to this point to create balance workload based on the current condition
                            //Recompute the cost_per_core. Perhaps one part was too big so the if statement won't work
                            cost_to_this_core += cost[thread_id - 1];
                            cost_per_core = (total_cost - cost_to_this_core) / (cores - thread_id);
                            //Let's do the next thread
                            thread_id++;
                            if(cost_per_core == 0){
                                if(part != part_cnt - 1){
                                    std::cerr << "Not all the parts are assigned" << std::endl;
                                }
                                break;
                            }

                            if(thread_id == cores){
                                chunk_array[cores] = part_cnt;
                                cost[cores - 1] = total_cost - cost_to_this_core;
                                break;
                            }
                        }
                    }
                }

                //Fix the last partitions of chunk_array
                for(int i = thread_id; i < cores; i++){
                    chunk_array[i + 1] = part_cnt;
                }

                double critical_cost = 0;
                for(int i = 0; i < effective_cores; i++){
                    auto normalize = cost[i];
                    normalize = normalize / total_cost;
                    curr_DKL += normalize * std::log2(normalize * effective_cores);
                }



                //            if(ub_level < 200){
                //                std::cout << "Effective Cores: " << effective_cores << " Number of parts: " << part_cnt << std::endl;
                //                std::cout << " Divergence:" << curr_DKL << " ub_Level: " << ub_level << " lb_level: " << lb_level << "\t";
                //                int cnt = 0;
                //                for(auto & c : normalize_cost){
                //                     std::cout << "Core" << cnt << ": " << c << "\t";
                //                     cnt++;
                //                }
                //                std::cout << std::endl;
                //            }

                //==================================== Decision making procedure for WM ============================
                //For a bit of reusing in next iteration
                computed_up_to_lvl = ub_level;
                if(curr_DKL < 0.003){
                    WM_condidate = ub_level;
                    ub_level = std::min(ub_level + 1, nlevels);
                    if(prev_upper_bound == ub_level){
                        if(iteration_counter == 1){
                            chunk_sizes[WM.size() - 1] = chunk_array;
                        }
                        done_flag = true;
                    }
                    prev_upper_bound = ub_level;
                } else {
                    if(iteration_counter == 1){
                        chunk_sizes[WM.size() - 1] = chunk_array;
                    }
                    done_flag = true;
                }

            }
            WM.push_back(WM_condidate);
            lvl = WM_condidate;
        }

        computeSchedule(cores, n, DAG_ptr, DAG_set, level_ptr.data(), level_set.data(),
                        node_to_level.data(), WM, chunk_sizes,
                        coarse_level_no, coarse_level_ptr, coarse_part_ptr, coarse_node_ptr,
                        bin_pack, postOrder);
        #ifndef NDEBUG
            std::cout << "The wm is: " << std::endl;
            for(auto wm: WM){
                std::cout << wm << "\t";
            }
            std::cout << std::endl;
        #endif
        return WM;

    }

    std::vector<int> GLC_V2(int n, int nnz,
                         std::vector<int> & DAG_ptr_not_prune, std::vector<int> & DAG_set_not_prune,
                         const std::vector<double> & node_cost,
                         int cores,
                         int & coarse_level_no,
                         std::vector<int> & coarse_level_ptr,
                         std::vector<int> & coarse_part_ptr,
                         std::vector<int> & coarse_node_ptr,
                         bool parallelLevelset,
                         bool postOrder,
                         bool bin_pack){

        if(!coarse_level_ptr.empty()){
            coarse_level_ptr.clear();
        }

        if(!coarse_part_ptr.empty()){
            coarse_part_ptr.clear();
        }

        if(!coarse_node_ptr.empty()){
            coarse_node_ptr.clear();
        }

        coarse_level_ptr.push_back(0);
        coarse_part_ptr.push_back(0);

        auto DAG_ptr = DAG_ptr_not_prune.data();
        auto DAG_set = DAG_set_not_prune.data();

        std::vector<int> level_ptr;
        std::vector<int> level_set;
        std::vector<int> node_to_level(n, 0);
        int nlevels;

        if(parallelLevelset){
#ifdef SPMP
            std::cerr << "Parallel levelset is still under construction" << std::endl;
            //Create the CSC version of the codeto use make_full and csc_to_csr codes
            sym_lib::CSC* Lower_A_CSC = new sym_lib::CSC(n, n, nnz);
            std::copy(DAG_ptr, DAG_ptr + n + 1, Lower_A_CSC->p);
            std::copy(DAG_set, DAG_set + nnz, Lower_A_CSC->i);
            std::fill(Lower_A_CSC->x, Lower_A_CSC->x + nnz, 1);
            Lower_A_CSC->stype=-1;
            auto tmp = sym_lib::make_full(Lower_A_CSC);
            auto CSR_A = sym_lib::csc_to_csr(tmp);
            auto A = new SpMP::CSR();
            HDAGG::Convert_LBC_CSR_to_SpMP(CSR_A, A);
            nlevels = HDAGG::levelsetCSRParallel_SpMP(A, cores, level_ptr, level_set);
            for(int l = 0; l < nlevels; l++){
                for(int node_ptr = level_ptr[l]; node_ptr < level_ptr[l + 1]; node_ptr++){
                    int node = level_set[node_ptr];
                    node_to_level[node] = l;
                }
            }
            delete Lower_A_CSC;
            delete tmp;
            delete A;
            delete CSR_A;
#endif
        } else {
            //Create the Levelset and compute node_to_level
            level_ptr.resize(n + 1);
            level_set.resize(n);

            nlevels = build_levelSet_CSC(n, DAG_ptr, DAG_set,
                                         level_ptr.data(), level_set.data(), node_to_level.data());
        }

        std::vector<int> WM;
        WM.push_back(0);
        std::vector<std::vector<int>> chunk_sizes(nlevels);
        //Compute the schedule
        for(int lvl = 0; lvl < nlevels; ) {
            int lb_level = lvl;
            int computed_up_to_lvl = lb_level;
            int prev_upper_bound = 0;
            int WM_condidate = lb_level + 1;
            int ub_level = std::min(lb_level + 1, nlevels);
            assert(ub_level <= nlevels && lb_level < nlevels);
            bool done_flag = false;
            //TODO: Clean this code later
            //Reusing material
            std::vector<int> node2partition(n, -1);
            std::vector<int> node2partition_accepted(n);

            double curr_DKL = 0;
            std::vector<int> chunk_array(cores + 1, 0);
            int iteration_counter = 0;
            std::vector<double> tmp;
            int effective_cores = cores;
            while(!done_flag){
                iteration_counter++;
                chunk_sizes[WM.size() - 1] = chunk_array;
                std::vector<double> partition_cost(n, 0);
                //===================================== Greedy Merge ===================================
                //Compute the cc's and w-partitions of these levels and compute the critical path
                int part_cnt = 0;
                double total_cost = 0;

                //Optimize for maximum parallelism TODO: Fix for huge serial part at the end
                if( (ub_level - lb_level == 1) && (level_ptr[ub_level] - level_ptr[ub_level - 1] < cores)){
                    effective_cores = level_ptr[ub_level] - level_ptr[ub_level - 1];
                } else if(level_ptr[ub_level] - level_ptr[ub_level - 1] > effective_cores){
                    effective_cores = cores;
                }

                computeCC(lb_level, computed_up_to_lvl, ub_level,
                          n, DAG_ptr, DAG_set, level_ptr.data(), level_set.data(), node_to_level.data(),
                          node2partition, part_cnt, partition_cost, node_cost, total_cost);

                //Compute cost
                //TODO: This is one of the ugliest code that I have ever done, and it is the core of my work :)))))
                std::vector<double> cost(cores, 0);
                double cost_per_core = total_cost / cores;
                int thread_id = 1;
                double cost_to_this_core = 0;
                //TODO: Consider the cases where it has a narrow part_cnt and then expand
                if(part_cnt < cores){
                    for(int part = 0; part < part_cnt; part++){
                        chunk_array[part + 1] = part + 1;
                        cost[part] = partition_cost[part];
                    }
                    thread_id = part_cnt;
                } else {
                    for(int part = 0; part < part_cnt; ){
                        if(cost[thread_id - 1] < cost_per_core){
                            cost[thread_id - 1] += partition_cost[part++];
                            chunk_array[thread_id] = part;
                            assert(cores != thread_id - 1);//I had a small optimization for last core
                        } else {
                            //Cost to this point to create balance workload based on the current condition
                            //Recompute the cost_per_core. Perhaps one part was too big so the if statement won't work
                            cost_to_this_core += cost[thread_id - 1];
                            cost_per_core = (total_cost - cost_to_this_core) / (cores - thread_id);
                            //Let's do the next thread
                            thread_id++;
                            if(cost_per_core == 0){
                                if(part != part_cnt - 1){
                                    std::cerr << "Not all the parts are assigned" << std::endl;
                                }
                                break;
                            }

                            if(thread_id == cores){
                                chunk_array[cores] = part_cnt;
                                cost[cores - 1] = total_cost - cost_to_this_core;
                                break;
                            }
                        }
                    }
                }

                //Fix the last partitions of chunk_array
                for(int i = thread_id; i < cores; i++){
                    chunk_array[i + 1] = part_cnt;
                }

                double critical_cost = 0;
                for(int i = 0; i < effective_cores; i++){
                    auto normalize = cost[i];
                    normalize = normalize / total_cost;
                    curr_DKL += normalize * std::log2(normalize * effective_cores);
                }



                //            if(ub_level < 200){
                //                std::cout << "Effective Cores: " << effective_cores << " Number of parts: " << part_cnt << std::endl;
                //                std::cout << " Divergence:" << curr_DKL << " ub_Level: " << ub_level << " lb_level: " << lb_level << "\t";
                //                int cnt = 0;
                //                for(auto & c : normalize_cost){
                //                     std::cout << "Core" << cnt << ": " << c << "\t";
                //                     cnt++;
                //                }
                //                std::cout << std::endl;
                //            }

                //==================================== Decision making procedure for WM ============================
                //For a bit of reusing in next iteration
                computed_up_to_lvl = ub_level;
                if(curr_DKL < 0.003){
                    WM_condidate = ub_level;
                    ub_level = std::min(ub_level + 1, nlevels);
                    if(prev_upper_bound == ub_level){
                        if(iteration_counter == 1){
                            chunk_sizes[WM.size() - 1] = chunk_array;
                            node2partition_accepted = node2partition;
                        }
                        done_flag = true;
                    }
                    node2partition_accepted = node2partition;
                    prev_upper_bound = ub_level;
                } else {
                    if(iteration_counter == 1){
                        chunk_sizes[WM.size() - 1] = chunk_array;
                        node2partition_accepted = node2partition;
                    }
                    done_flag = true;
                }
            }
            WM.push_back(WM_condidate);
            lvl = WM_condidate;
            //In this part we want to create the coarsened levels
            //Nodes merged together
            int ubl = WM[WM.size() - 1];
            int lbl = WM[WM.size() - 2];
            int numNodesAtCurLevel = level_ptr[ubl] - level_ptr[lbl];
            //Sorting the nodes based on their w-partitions
            std::vector<int> w_part_size_cnt(n + 1, 0);
            for(auto & iter: node2partition_accepted){
                if(iter != -1){
                    w_part_size_cnt[iter + 1]++;
                }
            }
            //Scan the counter vector
            for(int i = 1; i < w_part_size_cnt.size() - 1; i++){
                w_part_size_cnt[i + 1] += w_part_size_cnt[i];
            }

            std::vector<int> part_ptr;
            std::vector<int> node_ptr;
            //create part_ptr
            part_ptr.push_back(0);
            for(int i = 0; i < w_part_size_cnt.size() - 1; i++){
                if(w_part_size_cnt[i + 1] != w_part_size_cnt[i]){
                    part_ptr.push_back(w_part_size_cnt[i + 1]);
                }
            }
            //create node_ptr
            int* nodes_at_cur_level = &level_set[level_ptr[lbl]];
            node_ptr.resize(numNodesAtCurLevel);
            for(int i = 0; i < numNodesAtCurLevel; i++){
                assert(node2partition_accepted[nodes_at_cur_level[i]] != -1);
                node_ptr[w_part_size_cnt[node2partition_accepted[nodes_at_cur_level[i]]]++]
                = nodes_at_cur_level[i];
            }

            if(bin_pack){
                //Compute schedule chunk size
                std::vector<int> & chunk_array = chunk_sizes[WM.size() - 2];
                //Merge core's parts and apply schedule inside the merged parts
                std::vector<int> part_ptr_merged(cores + 1, 0);
                assert(chunk_array.back() == part_ptr.size() - 1);
                for(int i = 0; i < cores; i++){
                    size_t tmp = chunk_array[i + 1];
                    int maximum_part_idx = std::min(tmp, part_ptr.size() - 1);
                    part_ptr_merged[i + 1] = part_ptr[maximum_part_idx];
                }
                //Place arrays inside the coarse arrays
                auto tmp = coarse_level_ptr.back();
                assert(!part_ptr.empty());
                coarse_level_ptr.push_back(tmp + part_ptr_merged.size() - 1);
                auto last_part = coarse_part_ptr.back();
                for(int p = 1; p < part_ptr_merged.size(); p++){
                    auto tmp1 = part_ptr_merged[p] + last_part;
                    coarse_part_ptr.push_back(tmp1);
                }
            } else {
                //TODO: use the coarse array directly
                auto tmp = coarse_level_ptr.back();
                assert(!part_ptr.empty());
                coarse_level_ptr.push_back(tmp + part_ptr.size() - 1);
                auto last_part = coarse_part_ptr.back();
                for(int p = 1; p < part_ptr.size(); p++){
                    auto tmp1 = part_ptr[p] + last_part;
                    coarse_part_ptr.push_back(tmp1);
                }
            }
            coarse_node_ptr.insert(coarse_node_ptr.end(), node_ptr.begin(), node_ptr.end());
        }
        coarse_level_no = WM.size() - 1;

        #ifndef NDEBUG
        std::cout << "The wm is: " << std::endl;
        for(auto wm: WM){
            std::cout << wm << "\t";
        }
        std::cout << std::endl;
        #endif
        return WM;

    }

    void convertCNScheduleToSuperNode(int n, int num_supernode,
                                      int* group_ptr,
                                      int* group_set,
                                      int* supernode_ptr,
                                      int & coarse_level_no,
                                      std::vector<int> & coarse_level_ptr,
                                      std::vector<int> & coarse_part_ptr,
                                      std::vector<int> & coarse_node_ptr){
        //Create the supernode inverse
        std::vector<int> supernode_inv(n, 0);
        #pragma omp parallel for
        for(int super = 0; super < num_supernode; super++){
            for(int start_node = supernode_ptr[super]; start_node < supernode_ptr[super + 1]; start_node++){
                supernode_inv[start_node] = super;
            }
        }


        //Unpacking and blocking
        auto current_part_ptr = coarse_part_ptr;
        auto current_node_ptr = coarse_node_ptr;
        coarse_node_ptr.clear();
        for (int lvl = 0; lvl < coarse_level_no; ++lvl) {
            for (int lvl_ptr = coarse_level_ptr[lvl]; lvl_ptr < coarse_level_ptr[lvl + 1]; ++lvl_ptr) {
                int new_num_parts = 0;
                for (int part_ptr = current_part_ptr[lvl_ptr]; part_ptr < current_part_ptr[lvl_ptr + 1]; ++part_ptr) {
                    int group = current_node_ptr[part_ptr];
                    //Convert CN to supernode
                    int prev_node = -1;
                    for(int node_ptr = group_ptr[group]; node_ptr < group_ptr[group + 1]; node_ptr++){
                        int node = group_set[node_ptr];//TODO: Delete group set because it is unnecessary in CN grouping
                        if(supernode_inv[node] != prev_node){
                            prev_node = supernode_inv[node];
                            coarse_node_ptr.push_back(prev_node);
                            new_num_parts++;
                        }
                    }
                }
                coarse_part_ptr[lvl_ptr + 1] = new_num_parts + coarse_part_ptr[lvl_ptr];
            }
        }

        //TODO: Assert that every node in supernode exists

    }


    double aggressiveSupernodeGrouping(int num_nodes,
                                       int* group_ptr,
                                       int* group_set,
                                       int* DAG_ptr,
                                       int* DAG_set,
                                       int & coarse_level_no,
                                       std::vector<int> & coarse_level_ptr,
                                       std::vector<int> & coarse_part_ptr,
                                       std::vector<int> & coarse_node_ptr) {


        //Assign nodes to a core and to a level and search for supernodes inside that core and level
        std::vector<std::pair<int, int>> node_assigned(num_nodes, std::pair<int, int>());
        for (int lvl = 0; lvl < coarse_level_no; ++lvl) {
            int core = 0;
            for (int lvl_ptr = coarse_level_ptr[lvl]; lvl_ptr < coarse_level_ptr[lvl + 1]; ++lvl_ptr) {
                for (int part_ptr = coarse_part_ptr[lvl_ptr]; part_ptr < coarse_part_ptr[lvl_ptr + 1]; ++part_ptr) {
                    int group = coarse_node_ptr[part_ptr];
                    for (int node_ptr = group_ptr[group]; node_ptr < group_ptr[group + 1]; node_ptr++) {
                        int node = group_set[node_ptr];
                        node_assigned[node].first = core;
                        node_assigned[node].second = lvl;
                    }
                }
                core++;
            }
        }


        //Create supernode Group set and group ptr
        std::vector<bool> grouped(num_nodes, false);
        //It is a CSC/CSR like format
        std::vector<int> supernode_ptr;
        supernode_ptr.push_back(0);
        std::vector<int> supernode_set;
        for (int col = 0; col < num_nodes;) {
            //Nodes that are children of the current node
            std::vector<int> candidates(DAG_set + DAG_ptr[col] + 1, DAG_set + DAG_ptr[col + 1]);
            grouped[col] = true;
            auto tmp = supernode_ptr.back();
            supernode_ptr.push_back(tmp + 1);
            supernode_set.push_back(col);
            for(auto & iter: candidates){
                if(node_assigned[col] != node_assigned[iter]) {
                    continue;
                }
                if(grouped[iter]) {
                    std::cout << "Interesting" << std::endl;
                    continue;
                }
                int diff = iter - col;
                // compare number of off diagonals
                int off_last = DAG_ptr[iter + 1] - DAG_ptr[iter];
                int off_first = DAG_ptr[col + 1] - DAG_ptr[col] - diff;
                // If you don't understand this line
                // get back and read about supernode .. don't waste your time
                if (off_first != off_last)
                    break;
                assert(off_first > 0 && off_last > 0);
                // now examine column pattern
                bool can_form_supernode = true;
                for (int row = 0; row < off_first; row++) {
                    if (DAG_set[DAG_ptr[col] + row + diff] != DAG_set[DAG_ptr[iter] + row]) {
                        can_form_supernode = false;
                        break;
                    }
                }

                if (can_form_supernode) {
                    supernode_ptr.back()++;
                    supernode_set.push_back(iter);
                    grouped[iter] = true;
                }
            }
            //The starting point of the next supernode
            while (grouped[col]) {
                col++;
                if (col == num_nodes) {
                    break;
                }
            }
        }

        double avg_supernode_size = 0;
        double cnt = 0;
        for(int i = 0; i < supernode_ptr.size() - 1; i++) {
            assert( (supernode_ptr[i + 1] - supernode_ptr[i]) >= 1 );
            avg_supernode_size += (supernode_ptr[i + 1] - supernode_ptr[i]);
            cnt++;
        }
        std::cout << "WTF? " << avg_supernode_size << "\t" << cnt << std::endl;
        avg_supernode_size = avg_supernode_size / cnt;
        std::cout << "Average supernode size after HDAGG is: " << avg_supernode_size << std::endl;

        return avg_supernode_size;
    }


    void buildGroupDAG(const int& n, const int& ngroups, const int* group_ptr,
                         const int* group_set,
                         const int* DAG_ptr, const int* DAG_set,
                         std::vector<int>& group_DAG_ptr, std::vector<int>& group_DAG_set){
        //Computing group inverse
        std::vector<int> group_inv(n);
        #pragma omp parallel for
        for(int g = 0; g < ngroups; g++){
            for(int i = group_ptr[g]; i < group_ptr[g + 1]; i++){
                int node = group_set[i];
                group_inv[node] = g;
            }
        }

        std::vector<std::vector<int>> DAG(ngroups);
        int nnz = 0;
        #pragma omp parallel for reduction(+: nnz)
        for(int g = 0; g < ngroups; g++){
            for(int node_ptr = group_ptr[g]; node_ptr < group_ptr[g + 1]; node_ptr++){
                int node = group_set[node_ptr];
                for(int child_ptr = DAG_ptr[node]; child_ptr < DAG_ptr[node + 1]; child_ptr++){
                    int child = group_inv[DAG_set[child_ptr]];
                    DAG[g].push_back(child);
                }
            }
            std::sort(DAG[g].begin(), DAG[g].end());
            DAG[g].erase(std::unique(DAG[g].begin(), DAG[g].end()), DAG[g].end());
            nnz += DAG[g].size();
        }

        group_DAG_ptr.resize(ngroups + 1, 0);
        group_DAG_set.resize(nnz);

        long int cti,edges=0;
        for(cti = 0, edges = 0; cti < ngroups; cti++){
            group_DAG_ptr[cti] = edges;
            for(int ctj = 0; ctj < DAG[cti].size(); ctj++) {
                group_DAG_set[edges++] = DAG[cti][ctj];
            }
        }
        group_DAG_ptr[cti] = edges;
    }

    void computeParents(int& n, const int* Ap, const int* Ai, std::vector<int>& parents){
        int inext;
        parents.resize(n);
        std::vector<int> ancestor(n);
        for(int k = 0; k < n; k++){
            parents[k] = -1;
            ancestor[k] = -1;
            for(int p = Ap[k]; p < Ap[k + 1]; p++){
                for(int i = Ai[p]; i != -1 && i < k; i = inext){
                    assert(i < n);
                    assert(i >= 0);
                    inext = ancestor[i];
                    ancestor[i] = k;
                    if(inext == -1){
                        parents[i] = k;
                    }
                }
            }
        }
    }

    void buildETree(sym_lib::CSR* A, int& nedges, std::vector<int>& tree_ptr, std::vector<int>& tree_set){
        auto Ap = A->p;
        auto Ai = A->i;
        int n = A->n;

        std::vector<int> parents;
        computeParents(n, Ap, Ai, parents);

        tree_ptr.clear();
        tree_ptr.reserve(n);
        tree_set.clear();
        tree_set.reserve(2 * n);
        nedges = 0;
        for(int i = 0; i < n; i++){
            tree_ptr.push_back(nedges);
            tree_set.push_back(i);
            nedges++;
            if(parents[i] != -1){
                assert(parents[i] > i);
                tree_set.push_back(parents[i]);
                nedges++;
            }
        }
        tree_ptr.push_back(nedges);
        assert(tree_ptr.size() == n + 1);
    }

    int computeNoneZeroPatternOfRowK(const int &n, const int* Ap, const int* Ai,
                                     const int k, const int* parent, int* s, std::vector<bool>& w){
        int top = n - 1;
        assert(n != 1);
        int len;
        w[k] = true;
        for(int p = Ap[k]; p < Ap[k + 1]; p++){
            int i = Ai[p];
            for(len = 0; !w[i]; i = parent[i]){
                s[len++] = i;
                w[i] = true;
            }
            while(len > 0) s[--top] = s[--len];
        }
        for(int p = top; p < n; p++){
            w[s[p]] = false;
        }
        w[k] = false;
        s[n - 1] = k;
        return top;
    }

    bool treeBasedGrouping(int n, std::vector<int>& DAG_ptr, std::vector<int>& DAG_set,
                                      int& ngroups, std::vector<int>& group_ptr, std::vector<int>& group_set,
                                      bool Lfactor){

        int size_restriction = 1;
        if(Lfactor){
            size_restriction = 0.001 * n;
            if(size_restriction == 0){
                size_restriction = 1;
            }
        } else {
            size_restriction = n;
        }

        auto edge_counter = DAG_ptr[n];
        //Create tree grouping
        std::vector<int> group_DAG_ptr;
        std::vector<int> group_DAG_set;
        //Create the inverse DAG
        sym_lib::CSC DAG(n,n, edge_counter);
        std::copy(DAG_ptr.begin(), DAG_ptr.end(), DAG.p);
        std::copy(DAG_set.begin(), DAG_set.end(), DAG.i);
        auto DAG_inv = sym_lib::csc_to_csr(&DAG);
        std::vector<int> DAG_inv_ptr(DAG_inv->p, DAG_inv->p + n + 1);
        std::vector<int> DAG_inv_set(DAG_inv->i, DAG_inv->i + edge_counter);
        delete DAG_inv;
        //Find the leaves of the inverse DAG (roots of the DAG)
        std::vector<int> root_id;
        for(int i = 0; i < n; i++){
            if(DAG_ptr[i + 1] - DAG_ptr[i] == 1){
                root_id.push_back(i);
            }
        }
        //Apply a traversal to this DAG and group sub-trees
        //All the nodes in each group should have outgoing edge of one except the root
        std::vector<bool> visited(n, false);
        std::list<int> group_stack;
        ngroups = 1;
        if(!group_ptr.empty()){
            group_ptr.clear();
        }
        if(!group_set.empty()){
            group_set.clear();
        }
        group_set.resize(n);
        group_ptr.reserve(n);
        //TODO: Clean this code .. for a single group restriction size you mess the whole code up!
        int set_cnt = 0;
        group_ptr.push_back(0);
        while(!root_id.empty()){
            auto head = root_id.back();
            root_id.pop_back();
            visited[head] = true;
            group_stack.push_back(head);
            if(group_ptr.back() != set_cnt){
                group_ptr.push_back(set_cnt);
                ngroups++;
            }
            int last_size = set_cnt;
            while(!group_stack.empty()){
                head = group_stack.front();
                group_stack.pop_front();
                for(int child_ptr = DAG_inv_ptr[head]; child_ptr < DAG_inv_ptr[head + 1] - 1; child_ptr++){
                    int child = DAG_inv_set[child_ptr];
                    if(!visited[child]){
                        visited[child] = true;
                        if(DAG_ptr[child + 1] - DAG_ptr[child] > 2){ //The node can only become a root
                            root_id.push_back(child);
                        } else if (DAG_ptr[child + 1] - DAG_ptr[child] == 2){
                            group_stack.push_back(child);
                        } else {
                            std::cerr << "Invalid DAG in CSC format" << std::endl;
                        }
                    }
                }
                group_set[set_cnt++] = head;
                if(set_cnt - last_size == size_restriction){
                    group_ptr.push_back(set_cnt);
                    ngroups++;
                    last_size = set_cnt;
                }
            }
        }
        if(group_ptr.back() != set_cnt){
            group_ptr.push_back(set_cnt);
        }

        if(group_ptr.size() != ngroups + 1){
            if(group_ptr.size() != ngroups){
                std::cerr << "The group_ptr size is: " << group_ptr.size() << " The number of groups is: " << ngroups << std::endl;
                throw std::logic_error("Something is wrong in the Tree Grouping Function");
            } else {
                ngroups--;
            }
        }

        #ifndef NDEBUG
        std::vector<bool> marks(n, false);
        for(auto & iter: group_set){
            marks[iter] = true;
        }

        for(auto && mark : marks){
            assert(mark == true);
        }
        #endif
        std::vector<int> group_id(n, 0);
        std::vector<int> group_size(n + 1, 0);
        //NOTE: We have to do these stuff so that the resulted DAG will be a lower triangular one
        //Sort Groups and give each group a name based on the minimum of the ids inside a group
        #pragma omp parallel for
        for(int g = 0; g < ngroups; g++){
            //Don't sort right now, it will be ordered at the end of coarsening.
            //std::sort(group_set.data() + group_ptr[g], group_set.data() + group_ptr[g + 1]);
            //Smallest id
            int id = n;
            for(int i = group_ptr[g]; i < group_ptr[g + 1]; i++){
                int node = group_set[i];
                if(node < id){
                    id = node;
                }
            }
            group_size[id + 1] = group_ptr[g + 1] - group_ptr[g];
            for(int i = group_ptr[g]; i < group_ptr[g + 1]; i++){
                int node = group_set[i];
                group_id[node] = id;
            }
        }

        //Sort the nodes based on their ids using counting sort
        for(int i = 1; i < n; i++){
            group_size[i + 1] += group_size[i];
        }

        group_ptr.clear();
        group_ptr.push_back(0);
        for(int i = 0; i < n; i++){
            if(group_size[i + 1] != group_size[i]){
                group_ptr.push_back(group_size[i + 1]);
            }
        }

        std::vector<int> unordered_group_set(n);
        for(int i = 0; i < n; i++){
            unordered_group_set[group_size[group_id[group_set[i]]]++] = group_set[i];
        }

        group_set = unordered_group_set;
        return true;
    }

    bool treeBasedGrouping_Aggressive(int n, std::vector<int>& DAG_ptr, std::vector<int>& DAG_set,
                           int& ngroups, std::vector<int>& group_ptr, std::vector<int>& group_set,
                           bool apply_size_restriction){

        int size_restriction = 1;
        if(apply_size_restriction){
            size_restriction = 0.002 * n;
            if(size_restriction == 0){
                size_restriction = 1;
            }
        } else {
            size_restriction = n;
        }

        auto edge_counter = DAG_ptr[n];
        //Create tree grouping
        std::vector<int> group_DAG_ptr;
        std::vector<int> group_DAG_set;
        //Create the inverse DAG
        sym_lib::CSC DAG(n,n, edge_counter);
        std::copy(DAG_ptr.begin(), DAG_ptr.end(), DAG.p);
        std::copy(DAG_set.begin(), DAG_set.end(), DAG.i);
        auto DAG_inv = sym_lib::csc_to_csr(&DAG);
        std::vector<int> DAG_inv_ptr(DAG_inv->p, DAG_inv->p + n + 1);
        std::vector<int> DAG_inv_set(DAG_inv->i, DAG_inv->i + edge_counter);
        delete DAG_inv;
        //Find the leaves of the inverse DAG (roots of the DAG)
        std::vector<int> root_id;
        for(int i = 0; i < n; i++){
            if(DAG_ptr[i + 1] - DAG_ptr[i] == 1){
                root_id.push_back(i);
            }
        }
        //Apply a traversal to this DAG and group sub-trees
        //All the nodes in each group should have outgoing edge of one except the root
        std::vector<bool> visited(n, false);
        std::list<int> group_stack;
        ngroups = 1;
        if(!group_ptr.empty()){
            group_ptr.clear();
        }
        if(!group_set.empty()){
            group_set.clear();
        }
        group_set.resize(n);
        group_ptr.reserve(n);
        //TODO: Clean this code .. for a single group restriction size you mess the whole code up!
        int set_cnt = 0;
        group_ptr.push_back(0);
        while(!root_id.empty()){
            auto head = root_id.back();
            root_id.pop_back();
            visited[head] = true;
            group_stack.push_back(head);
            if(group_ptr.back() != set_cnt){
                group_ptr.push_back(set_cnt);
                ngroups++;
            }
            int last_size = set_cnt;
            while(!group_stack.empty()){
                head = group_stack.front();
                group_stack.pop_front();
                for(int child_ptr = DAG_inv_ptr[head]; child_ptr < DAG_inv_ptr[head + 1] - 1; child_ptr++){
                    int child = DAG_inv_set[child_ptr];
                    if(!visited[child]){
                        visited[child] = true;
                        if(DAG_ptr[child + 1] - DAG_ptr[child] > 2){ //The node can only become a root
                            root_id.push_back(child);
                        } else if (DAG_ptr[child + 1] - DAG_ptr[child] == 2){
                            group_stack.push_back(child);
                        } else {
                            std::cerr << "Invalid DAG in CSC format" << std::endl;
                        }
                    }
                }
                group_set[set_cnt++] = head;
                if(set_cnt - last_size == size_restriction){
                    group_ptr.push_back(set_cnt);
                    ngroups++;
                    last_size = set_cnt;
                }
            }
        }
        if(group_ptr.back() != set_cnt){
            group_ptr.push_back(set_cnt);
        }

        if(group_ptr.size() != ngroups + 1){
            if(group_ptr.size() != ngroups){
                std::cerr << "The group_ptr size is: " << group_ptr.size() << " The number of groups is: " << ngroups << std::endl;
                throw std::logic_error("Something is wrong in the Tree Grouping Function");
            } else {
                ngroups--;
            }
        }

        #ifndef NDEBUG
        std::vector<bool> marks(n, false);
        for(auto & iter: group_set){
            marks[iter] = true;
        }

        for(auto && mark : marks){
            assert(mark == true);
        }
        #endif
        std::vector<int> group_id(n, 0);
        std::vector<int> group_size(n + 1, 0);
        //NOTE: We have to do these stuff so that the resulted DAG will be a lower triangular one
        //Sort Groups and give each group a name based on the minimum of the ids inside a group
        for(int g = 0; g < ngroups; g++){
            //Don't sort right now, it will be ordered at the end of coarsening.
//            std::sort(group_set.data() + group_ptr[g], group_set.data() + group_ptr[g + 1]);
            //Smallest id
            int id = n;
            for(int i = group_ptr[g]; i < group_ptr[g + 1]; i++){
                int node = group_set[i];
                if(node < id){
                    id = node;
                }
            }
            group_size[id + 1] = group_ptr[g + 1] - group_ptr[g];
            for(int i = group_ptr[g]; i < group_ptr[g + 1]; i++){
                int node = group_set[i];
                group_id[node] = id;
            }
        }

        //Sort the nodes based on their ids using counting sort
        for(int i = 1; i < n; i++){
            group_size[i + 1] += group_size[i];
        }

        group_ptr.clear();
        group_ptr.push_back(0);
        for(int i = 0; i < n; i++){
            if(group_size[i + 1] != group_size[i]){
                group_ptr.push_back(group_size[i + 1]);
            }
        }

        std::vector<int> unordered_group_set(n);
        for(int i = 0; i < n; i++){
            unordered_group_set[group_size[group_id[group_set[i]]]++] = group_set[i];
        }

        group_set = unordered_group_set;
        return true;
    }

    bool chainGrouping(int n, std::vector<int>& DAG_ptr, std::vector<int>& DAG_set,
                              int& ngroups, std::vector<int>& group_ptr, std::vector<int>& group_set){
        auto edge_counter = DAG_ptr[n];

        //Find the leaves of the inverse DAG (roots of the DAG)
        std::vector<int> root_id;
        std::vector<int> in_degree(n, 0);
        for(int i = 0; i < n; i++){
            for(int child_ptr = DAG_ptr[i] + 1; child_ptr < DAG_ptr[i + 1]; child_ptr++){
                int child = DAG_set[child_ptr];
                in_degree[child]++;
            }
        }

        //Find the roots of the tree
        for(int j = 0; j < n; j++){
            if(in_degree[j] == 0){
                root_id.push_back(j);
            }
        }

        //Apply a traversal to this DAG and group sub-trees
        //All the nodes in each group should have outgoing edge of one except the root
        std::vector<bool> visited(n, false);
        std::vector<int> group_stack;

        if(!group_ptr.empty()){
            group_ptr.clear();
        }
        if(!group_set.empty()){
            group_set.clear();
        }
        group_set.resize(n);
        group_ptr.reserve(n);
        //TODO: Clean this code .. for a single group restriction size you mess the whole code up!
        int set_cnt = 0;
        ngroups = 0;
        while(!root_id.empty()){
            auto head = root_id.back();
            root_id.pop_back();
            visited[head] = true;
            group_stack.push_back(head);
            group_ptr.push_back(set_cnt);
            ngroups++;
            int finish_group = false;
            while(!group_stack.empty()){
                head = group_stack.back();
                group_stack.pop_back();
                visited[head] = true;
                if(DAG_ptr[head + 1] - DAG_ptr[head] == 2){
                    int child = DAG_set[DAG_ptr[head] + 1];
                    if(in_degree[child] == 1){
                        assert(visited[child] == false);
                        group_stack.push_back(child);
                    } else if(in_degree[child] > 1){
                        if(!visited[child]){
                            root_id.push_back(child);
                        }
                        finish_group = true;
                    } else {
                        std::cerr << "Invalid DAG in CSC format" << std::endl;
                    }
                } else {
                    for(int i = DAG_ptr[head] + 1; i < DAG_ptr[head + 1]; i++){
                        int child = DAG_set[i];
                        if( !visited[child] ){
                            root_id.push_back(child);
                        }
                    }
                    finish_group = true;
                }
                group_set[set_cnt++] = head;
                if(finish_group == true){
                    assert(group_stack.empty());
                }
            }
        }
        if(group_ptr.back() != set_cnt){
            group_ptr.push_back(set_cnt);
        }

        #ifndef NDEBUG
        std::vector<bool> marks(n, false);
        for(auto & iter: group_set){
            marks[iter] = true;
        }

        for(auto && mark : marks){
            assert(mark == true);
        }
        #endif
        std::vector<int> group_id(n, 0);
        std::vector<int> group_size(n + 1, 0);
        //NOTE: We have to do these stuff so that the resulted DAG will be a lower triangular one
        //Sort Groups and give each group a name based on the minimum of the ids inside a group
        for(int g = 0; g < ngroups; g++){
            //Don't sort right now, it will be ordered at the end of coarsening.
            //            std::sort(group_set.data() + group_ptr[g], group_set.data() + group_ptr[g + 1]);
            //Smallest id
            int id = n;
            for(int i = group_ptr[g]; i < group_ptr[g + 1]; i++){
                int node = group_set[i];
                if(node < id){
                    id = node;
                }
            }
            group_size[id + 1] = group_ptr[g + 1] - group_ptr[g];
            for(int i = group_ptr[g]; i < group_ptr[g + 1]; i++){
                int node = group_set[i];
                group_id[node] = id;
            }
        }

        //Sort the nodes based on their ids using counting sort
        for(int i = 1; i < n; i++){
            group_size[i + 1] += group_size[i];
        }

        group_ptr.clear();
        group_ptr.push_back(0);
        for(int i = 0; i < n; i++){
            if(group_size[i + 1] != group_size[i]){
                group_ptr.push_back(group_size[i + 1]);
            }
        }

        std::vector<int> unordered_group_set(n);
        for(int i = 0; i < n; i++){
            unordered_group_set[group_size[group_id[group_set[i]]]++] = group_set[i];
        }

        group_set = unordered_group_set;
        return true;
    }

    bool treeBasedGroupingBFS(int n, std::vector<int>& DAG_ptr, std::vector<int>& DAG_set,
                              int& ngroups, std::vector<int>& group_ptr, std::vector<int>& group_set, bool Lfactor){

        int size_restriction = 1;
        if(Lfactor){
            size_restriction = 0.0001 * n;
            if(size_restriction == 0){
                size_restriction = 1;
            }
        } else {
            size_restriction = n;
        }
//        timing_measurement t1;
//        t1.start_timer();
        auto edge_counter = DAG_ptr[n];
        //Create tree grouping
        std::vector<int> group_DAG_ptr;
        std::vector<int> group_DAG_set;
        //Create the inverse DAG
        sym_lib::CSC DAG(n,n, edge_counter);
        std::copy(DAG_ptr.begin(), DAG_ptr.end(), DAG.p);
        std::copy(DAG_set.begin(), DAG_set.end(), DAG.i);
        auto DAG_inv = sym_lib::csc_to_csr(&DAG);
        std::vector<int> DAG_inv_ptr(DAG_inv->p, DAG_inv->p + n + 1);
        std::vector<int> DAG_inv_set(DAG_inv->i, DAG_inv->i + edge_counter);
        delete DAG_inv;
//        t1.measure_elapsed_time();
//        std::cout << "Inverse time " << t1.measure_elapsed_time() << std::endl;
        //Find the leaves of the inverse DAG (roots of the DAG)
        std::list<int> root_id;
        for(int i = 0; i < n; i++){
            if(DAG_ptr[i + 1] - DAG_ptr[i] == 1){
                root_id.push_back(i);
            }
        }

//        timing_measurement t2;
//        t2.start_timer();
        //Apply a traversal to this DAG and group sub-trees
        //All the nodes in each group should have outgoing edge of one except the root
        std::vector<bool> visited(n, false);
        std::list<int> group_stack;

        if(!group_ptr.empty()){
            group_ptr.clear();
        }
        if(!group_set.empty()){
            group_set.clear();
        }
        group_set.resize(n);
        group_ptr.reserve(n);
        //TODO: Clean this code .. for a single group restriction size you mess the whole code up!
        int set_cnt = 0;
        ngroups = 0;
        //The BFS Traversal
        while(!root_id.empty()){
            //First visit, first serve (BFS style)
            auto head = root_id.front();
            root_id.pop_front();
            visited[head] = true;
            group_stack.push_back(head);
            group_ptr.push_back(set_cnt);
            int last_size = set_cnt;
            ngroups++;
            while(!group_stack.empty()){
                head = group_stack.front();
                group_stack.pop_front();
                int finish_group = false;
                group_set[set_cnt++] = head;
                for(int child_ptr = DAG_inv_ptr[head]; child_ptr < DAG_inv_ptr[head + 1] - 1; child_ptr++){
                    int child = DAG_inv_set[child_ptr];
                    if(!visited[child]){
                        visited[child] = true;
                        if(DAG_ptr[child + 1] - DAG_ptr[child] > 2){ //The node can only become a root
                            root_id.push_back(child);
                            finish_group = true;
                        } else if (DAG_ptr[child + 1] - DAG_ptr[child] == 2){
                            group_stack.push_back(child);
                        } else {
                            std::cerr << "Invalid DAG in CSC format" << std::endl;
                        }
                    }
                }
                if(finish_group == true || set_cnt - last_size == size_restriction){
                    while(!group_stack.empty()){
                        auto tmp = group_stack.front();
                        group_stack.pop_front();
                        root_id.push_back(tmp);
                    }
                }
            }
        }
        if(group_ptr.back() != set_cnt){
            group_ptr.push_back(set_cnt);
        }

//        t2.measure_elapsed_time();
//        std::cout << "Grouping time " << t2.elapsed_time << std::endl;

//        timing_measurement t3;
//        t3.start_timer();
        #ifndef NDEBUG
        std::vector<bool> marks(n, false);
        for(auto & iter: group_set){
            marks[iter] = true;
        }

        for(auto && mark : marks){
            assert(mark == true);
        }
        #endif
        std::vector<int> group_id(n, 0);
        std::vector<int> group_size(n + 1, 0);
        //NOTE: We have to do these stuff so that the resulted DAG will be a lower triangular one
        //Sort Groups and give each group a name based on the minimum of the ids inside a group
        for(int g = 0; g < ngroups; g++){
            //Don't sort right now, it will be ordered at the end of coarsening.
            //            std::sort(group_set.data() + group_ptr[g], group_set.data() + group_ptr[g + 1]);
            //Smallest id
            int id = n;
            for(int i = group_ptr[g]; i < group_ptr[g + 1]; i++){
                int node = group_set[i];
                if(node < id){
                    id = node;
                }
            }
            group_size[id + 1] = group_ptr[g + 1] - group_ptr[g];
            for(int i = group_ptr[g]; i < group_ptr[g + 1]; i++){
                int node = group_set[i];
                group_id[node] = id;
            }
        }

        //Sort the nodes based on their ids using counting sort
        for(int i = 1; i < n; i++){
            group_size[i + 1] += group_size[i];
        }

        group_ptr.clear();
        group_ptr.push_back(0);
        for(int i = 0; i < n; i++){
            if(group_size[i + 1] != group_size[i]){
                group_ptr.push_back(group_size[i + 1]);
            }
        }

        std::vector<int> unordered_group_set(n);
        for(int i = 0; i < n; i++){
            unordered_group_set[group_size[group_id[group_set[i]]]++] = group_set[i];
        }

        group_set = unordered_group_set;
//        t3.measure_elapsed_time();
//        std::cout << "The sorting part " << t3.elapsed_time << std::endl;
        return true;
    }

    void ungroupingScheduleAndApplyOrdering(int n, int final_level_no, std::vector<int>& final_level_ptr,
                                            std::vector<int>& final_part_ptr, std::vector<int>& final_node_ptr,
                                            std::vector<int>& group_ptr, std::vector<int>& group_set,
                                            const int* DAG_ptr, const int* DAG_set,
                                            bool apply_sort_order, bool apply_post_order){

        if(apply_post_order){
            if(DAG_ptr == NULLPNTR || DAG_set == NULLPNTR){
                throw std::logic_error("For post ordering the prune DAG (with no grouping) is required");
            }
        }

        //Unpack nodes from groups
        std::vector<int> unpack_part_ptr;
        std::vector<int> unpack_node_ptr;
        int node_cnt = 0;
        for (int i1 = 0; i1 < final_level_no; ++i1){
            for (int j1 = final_level_ptr[i1]; j1 < final_level_ptr[i1 + 1]; ++j1){
                unpack_part_ptr.push_back(node_cnt);
                for (int k1 = final_part_ptr[j1]; k1 < final_part_ptr[j1 + 1]; ++k1){
                    int g = final_node_ptr[k1];
                    for (int k = group_ptr[g]; k < group_ptr[g + 1]; ++k){
                        int node = group_set[k];
                        unpack_node_ptr.push_back(node);
                        node_cnt++;
                    }
                }
            }
        }
        unpack_part_ptr.push_back(node_cnt);

        if(apply_post_order && apply_sort_order){
            throw std::logic_error("The code cannot apply both ordering to a single sequence");
        }

        //Recreate final part and node pointers and sort
        assert(final_level_ptr[final_level_no] + 1 == unpack_part_ptr.size());
        assert(unpack_node_ptr.size() == n);
        final_part_ptr = unpack_part_ptr;
        final_node_ptr = unpack_node_ptr;
        if(apply_sort_order){
            for (int i1 = 0; i1 < final_level_no; ++i1){
                for (int j1 = final_level_ptr[i1]; j1 < final_level_ptr[i1 + 1]; ++j1){
                    if(final_part_ptr[j1] != final_part_ptr[j1 + 1]){
                        std::sort(final_node_ptr.data() + final_part_ptr[j1],
                                  final_node_ptr.data() + final_part_ptr[j1 + 1]);
                    }
                }
            }
        } else if(apply_post_order){
            #ifndef NDEBUG
            std::cout << "Applying post ordering after unpacking the groups" << std::endl;
            #endif
            std::vector<int> post_order;
            std::vector<int> post_orderInv(n);
            std::vector<int> node_ptr_order(n);
            computePostOrder(n,  DAG_ptr, DAG_set, post_order);
            for(int i = 0; i < post_order.size(); i++){
                post_orderInv[post_order[i]] = i;
            }
            for (int lvl = 0; lvl < final_level_no; ++lvl){
                for (int p = final_level_ptr[lvl]; p < final_level_ptr[lvl + 1]; ++p){
                    int cnt = 0;
                    for(int np = final_part_ptr[p]; np < final_part_ptr[p + 1]; ++np){
                        int node = final_node_ptr[np];
                        node_ptr_order[cnt++] = post_orderInv[node];
                    }
                    std::sort(node_ptr_order.data(), node_ptr_order.data() + cnt);
                    cnt = 0;
                    for(int np = final_part_ptr[p]; np < final_part_ptr[p + 1]; ++np){
                        final_node_ptr[np] = post_order[node_ptr_order[cnt++]];
                    }
                }
            }
        }
    }

    void partialSparsification(int n, int nnz, const int* DAG_ptr_not_prune, const int* DAG_set_not_prune,
                               std::vector<int>& DAG_ptr, std::vector<int>& DAG_set, bool cliqueSimplification){

        if(cliqueSimplification){
            //find cliques -> after simplification they are just chains
            //For now it only support lower matrices
            auto Lp = DAG_ptr_not_prune;
            auto Li = DAG_set_not_prune;
            int limit = n;
            //It is a CSC/CSR like format
            std::vector<int> clique_ptr;
            clique_ptr.reserve(n);
            clique_ptr.push_back(0);
            std::vector<int> node_to_clique(n);
            #pragma omp parallel
            {

#ifdef ENABLE_OPENMP
          int bins = omp_get_num_threads();
          int tid = omp_get_thread_num();
#else
          int bins =1;
       int tid=1;
#endif
                int start_col = (n / bins) * tid;
                int end_col = (n / bins) * (tid + 1);
                std::vector<int> clique_ptr_per_thread;
                if(tid == bins - 1){
                    end_col = n;
                }
                clique_ptr_per_thread.reserve(end_col - start_col);
                clique_ptr_per_thread.push_back(start_col);
                for (int col = start_col; col < end_col; ) {
                    int width = 1;
                    int first = col;
                    int last = col + 1;
                    while (last < n && width < limit) {
                        int diff = last - first;
                        // compare number of off diagonals
                        int off_last  = Lp[last + 1] - Lp[last];
                        int off_first = Lp[first + 1] - Lp[first] - diff;
                        // If you don't understand this line
                        // get back and read about clique .. don't waste your time
                        if (off_first != off_last)
                            break;
                        assert(off_first > 0 && off_last > 0);
                        // now examine column pattern
                        bool can_form_clique = true;
                        for (int row = 0; row < off_first; row++) {
                            if (Li[Lp[first] + row + diff] != Li[Lp[last] + row]) {
                                can_form_clique = false;
                                break;
                            }
                        }
                        if (can_form_clique) {
                            //Increase the width of the clique
                            width++;
                            //Lets check the next row
                            last++;
                        } else {
                            //Lets go to the next row, this clique is done
                            break;
                        }
                    }
                    //The starting point of the next clique
                    col = first + width;
                    assert(col <= n);
                    clique_ptr_per_thread.push_back(col);
                    assert(Lp[first + 1] - Lp[first] - width >= 0);
                }

                for(int i = 0; i < clique_ptr_per_thread.size() - 1; i++){
                    int start_node = clique_ptr_per_thread[i];
                    int end_node = clique_ptr_per_thread[i + 1];
                    for(int j = start_node; j < end_node; j++){
                        node_to_clique[j] = i + start_col;
                    }
                }
            }

            int back = 0;
            int node_cnt = 0;
            for(int i : node_to_clique){
                if(i != back){
                    clique_ptr.push_back(node_cnt);
                }
                node_cnt++;
                back = i;
            }
            clique_ptr.push_back(node_cnt);
            assert(clique_ptr.back() == n);

//            for (int col = 0; col < n; ) {
//                int width = 1;
//                int first = col;
//                int last = col + 1;
//                while (last < n && width < limit) {
//                    int diff = last - first;
//                    // compare number of off diagonals
//                    int off_last  = Lp[last + 1] - Lp[last];
//                    int off_first = Lp[first + 1] - Lp[first] - diff;
//                    // If you don't understand this line
//                    // get back and read about clique .. don't waste your time
//                    if (off_first != off_last)
//                        break;
//                    assert(off_first > 0 && off_last > 0);
//                    // now examine column pattern
//                    bool can_form_clique = true;
//                    for (int row = 0; row < off_first; row++) {
//                        if (Li[Lp[first] + row + diff] != Li[Lp[last] + row]) {
//                            can_form_clique = false;
//                            break;
//                        }
//                    }
//                    if (can_form_clique) {
//                        //Increase the width of the clique
//                        width++;
//                        //Lets check the next row
//                        last++;
//                    } else {
//                        //Lets go to the next row, this clique is done
//                        break;
//                    }
//                }
//                //The starting point of the next clique
//                col = first + width;
//                assert(col <= n);
//                clique_ptr.push_back(col);
//                assert(Lp[first + 1] - Lp[first] - width >= 0);
//            }
//
//            assert(clique_ptr.back() == n);

            //Creating the prune DAG based on the clique pointers
            std::vector<int> DAG_ptr_no_clique(n + 1, 0);
            std::vector<int> DAG_set_no_clique;
            DAG_set_no_clique.reserve(nnz);
            int edge_cnt = 0;
            for(int i = 0; i < clique_ptr.size() - 1; i++){
                for(int j = clique_ptr[i]; j < clique_ptr[i + 1] - 1; j++){
                    DAG_ptr_no_clique[j] = edge_cnt;
                    DAG_set_no_clique.push_back(j);
                    DAG_set_no_clique.push_back(j + 1);
                    edge_cnt+=2;
                }
                int last_node = clique_ptr[i + 1] - 1;
                DAG_ptr_no_clique[last_node] = edge_cnt;
                for(int j = DAG_ptr_not_prune[last_node]; j < DAG_ptr_not_prune[last_node + 1]; j++){
                    DAG_set_no_clique.push_back(DAG_set_not_prune[j]);
                    edge_cnt++;
                }
            }
            DAG_ptr_no_clique[n] = edge_cnt;

            std::vector<int> deleted_edge(edge_cnt, 0);
            #pragma omp parallel for
            for(int i = 0; i < n; i++){
                int start_i_child_ptr = DAG_ptr_no_clique[i] + 1;
                int end_i_child_ptr = DAG_ptr_no_clique[i + 1];
                for(int k_ptr = DAG_ptr_no_clique[i] + 1; k_ptr < DAG_ptr_no_clique[i + 1]; k_ptr++){
                    int k = DAG_set_no_clique[k_ptr];//Child of i
                    int l2 = DAG_ptr_no_clique[k] + 1;
                    int u2 = DAG_ptr_no_clique[k + 1];
                    int l1 = start_i_child_ptr;
                    int u1 = end_i_child_ptr;
                    while ( l1 < u1 && l2 < u2) {
                        if (DAG_set_no_clique[l1] == DAG_set_no_clique[l2]){
                            deleted_edge[l1] = 1;
                            l1++;
                            l2++;
                        } else if (DAG_set_no_clique[l1] < DAG_set_no_clique[l2]){
                            l1++;
                        } else {
                            l2++;
                        }
                    }
                }
            }

            DAG_ptr.resize(n + 1, 0);
            DAG_set.reserve(edge_cnt);
            int edge_counter = 0;
            for(int i = 0; i < n; i++) {
                for (int k_ptr = DAG_ptr_no_clique[i] ; k_ptr < DAG_ptr_no_clique[i + 1]; k_ptr++) {
                    if(!deleted_edge[k_ptr]){
                        int k = DAG_set_no_clique[k_ptr];//Child of i
                        DAG_set.push_back(k);
                        edge_counter++;
                    }
                }
                DAG_ptr[i + 1] = edge_counter;
            }
        } else {
            std::cout << "SLOW Sparsification" << std::endl;
            std::vector<int> deleted_edge(nnz, 0);
            #pragma omp parallel for
            for(int i = 0; i < n; i++){
                std::vector<int> child(DAG_ptr_not_prune[i + 1] - (DAG_ptr_not_prune[i] + 1));
                int cnt = 0;
                int start_k_ptr = DAG_ptr_not_prune[i] + 1;
                for(int k_ptr = DAG_ptr_not_prune[i] + 1; k_ptr < DAG_ptr_not_prune[i + 1]; k_ptr++){
                    int k = DAG_set_not_prune[k_ptr];//Child of i
                    child[cnt++] = k;
                }

                for(int k_ptr = DAG_ptr_not_prune[i] + 1; k_ptr < DAG_ptr_not_prune[i + 1]; k_ptr++){
                    int k = DAG_set_not_prune[k_ptr];//Child of i
                    for(int j_ptr = DAG_ptr_not_prune[k] + 1; j_ptr < DAG_ptr_not_prune[k + 1]; j_ptr++){
                        int j = DAG_set_not_prune[j_ptr];//Child of k
                        auto iter = std::lower_bound(child.begin(),child.end(), j);
                        if(iter != child.end() && *iter <= j){
                            auto id = std::distance(child.begin(), iter);
                            deleted_edge[start_k_ptr + id] = 1;
                        }
                    }
                }
            }
//            int num_edges = 0;
//            for(auto && i : deleted_edge){
//                if(i == 0){
//                    num_edges++;
//                }
//            }

            DAG_ptr.resize(n + 1, 0);
            DAG_set.reserve(nnz);
            int edge_counter = 0;
            for(int i = 0; i < n; i++) {
                for (int k_ptr = DAG_ptr_not_prune[i] ; k_ptr < DAG_ptr_not_prune[i + 1]; k_ptr++) {
                    if(!deleted_edge[k_ptr]){
                        int k = DAG_set_not_prune[k_ptr];//Child of i
                        DAG_set.push_back(k);
                        edge_counter++;
                    }
                }
                DAG_ptr[i + 1] = edge_counter;
            }
        }
    }


    void costComputation(int nodes, const int* CSC_Lp, const int* CSC_Li, const int* CSR_Lp, const int* CSR_Li,
                         Kernel kernel, const int* group_ptr, const int* group_set, bool grouped,
                         std::vector<double>& cost){
        if(!cost.empty()){
            cost.clear();
            cost.resize(nodes, 0);
        }
        if(kernel == SpTrSv_LL || kernel == SpTrSv_RL){
            const int* Lp;
            if(kernel == SpTrSv_LL){
                Lp = CSR_Lp;
            } else {
                Lp = CSC_Lp;
            }
            if(grouped){
                #pragma omp parallel for
                for(int group = 0; group < nodes; group++){
                    for(int node_ptr = group_ptr[group]; node_ptr < group_ptr[group + 1]; node_ptr++){
                        int node = group_set[node_ptr];
                        cost[group] += 1 * (Lp[node + 1] - Lp[node]);
                    }
                }
            } else {
                #pragma omp parallel for
                for(int row = 0; row < nodes; row++){
                    cost[row] += Lp[row + 1] - Lp[row];
                }
            }
        } else if(kernel == SpICh0_LL){
            if(grouped){
                #pragma omp parallel for
                for(int group = 0; group < nodes; group++){
                    for(int node_ptr = group_ptr[group]; node_ptr < group_ptr[group + 1]; node_ptr++){
                        int node = group_set[node_ptr];
                        for (int j = CSR_Lp[node]; j < CSR_Lp[node + 1]; j++) {
                            auto dep = CSR_Li[j];
                            cost[group] += (CSC_Lp[dep + 1] - CSC_Lp[dep]);
                        }
                    }
                }
            } else {
                #pragma omp parallel for
                for(int node = 0; node < nodes; node++){
                    for (int j = CSR_Lp[node]; j < CSR_Lp[node + 1]; j++) {
                        auto dep = CSR_Li[j];
                        cost[node] += (CSC_Lp[dep + 1] - CSC_Lp[dep]);
                    }
                }
            }
        } else if (kernel == SpICh0_RL){
            if(grouped){
                #pragma omp parallel for
                for(int group = 0; group < nodes; group++){
                    for(int node_ptr = group_ptr[group]; node_ptr < group_ptr[group + 1]; node_ptr++){
                        int node = group_set[node_ptr];
                        for (int j = CSC_Lp[node]; j < CSC_Lp[node + 1]; j++) {
                            auto dep = CSC_Li[j];
                            cost[group] += (CSC_Lp[dep + 1] - CSC_Lp[dep]);
                        }
                    }
                }
            } else {
                #pragma omp parallel for
                for(int node = 0; node < nodes; node++){
                    for (int j = CSC_Lp[node]; j < CSC_Lp[node + 1]; j++) {
                        auto dep = CSC_Li[j];
                        cost[node] += (CSC_Lp[dep + 1] - CSC_Lp[dep]);
                    }
                }
            }
        } else if (kernel == SpICh0_UL) {
            if (grouped) {
                #pragma omp parallel for
                for (int group = 0; group < nodes; group++) {
                    for (int node_ptr = group_ptr[group]; node_ptr < group_ptr[group + 1]; node_ptr++) {
                        int node = group_set[node_ptr];
                        for (int j = CSR_Lp[node]; j < CSR_Lp[node + 1]; j++) {
                            auto dep = CSR_Li[j];
                            cost[group] += (CSR_Lp[dep + 1] - CSR_Lp[dep]);
                        }
                    }
                }
            } else {
                #pragma omp parallel for
                for (int node = 0; node < nodes; node++) {
                    for (int j = CSR_Lp[node]; j < CSR_Lp[node + 1]; j++) {
                        auto dep = CSR_Li[j];
                        cost[node] += (CSR_Lp[dep + 1] - CSR_Lp[dep]);
                    }
                }
            }
        } else if (kernel == SpILU0_UL){
            if(grouped){
                #pragma omp parallel for
                for(int group = 0; group < nodes; group++){
                    for(int node_ptr = group_ptr[group]; node_ptr < group_ptr[group + 1]; node_ptr++){
                        int node = group_set[node_ptr];
                        for (int j = CSR_Lp[node]; j < CSR_Lp[node + 1]; j++) {
                            if(node > j){
                                break;
                            }
                            auto dep = CSR_Li[j];
                            cost[group] += (CSR_Lp[dep + 1] - CSR_Lp[dep]);
                        }
                    }
                }
            } else {
                    #pragma omp parallel for
                for(int node = 0; node < nodes; node++){
                    for (int j = CSR_Lp[node]; j < CSR_Lp[node + 1]; j++) {
                        auto dep = CSR_Li[j];
                        cost[node] += (CSR_Lp[dep + 1] - CSR_Lp[dep]);
                    }
                }
            }
        } else if (kernel == General){
                if(grouped){
                    #pragma omp parallel for
                    for(int group = 0; group < nodes; group++){
                        for(int node_ptr = group_ptr[group]; node_ptr < group_ptr[group + 1]; node_ptr++){
                            int node = group_set[node_ptr];
                            for (int j = CSR_Lp[node]; j < CSR_Lp[node + 1]; j++) {
                                auto dep = CSR_Li[j];
                                cost[group] += (CSR_Lp[dep + 1] - CSR_Lp[dep]);
                            }
                        }
                    }
                } else {
                    #pragma omp parallel for
                    for(int node = 0; node < nodes; node++){
                        for (int j = CSR_Lp[node]; j < CSR_Lp[node + 1]; j++) {
                            auto dep = CSR_Li[j];
                            cost[node] += (CSR_Lp[dep + 1] - CSR_Lp[dep]);
                        }
                    }
                }
        } else{
            std::wcerr << "The kernel is not supported, using cost of 1 for each node" << std::endl;
        }
    }


    void twoDReordering(int n,
                        int nnz,
                        int* Lp,
                        int* Li,
                        int final_level_no,
                        std::vector<int>& final_level_ptr,
                        std::vector<int>& final_part_ptr,
                        std::vector<int>& final_node_ptr,
                        int T1_size,
                        std::vector<int>& final_tuple_ptr,
                        std::vector<twoDIteration>& final_iter_ptr){

        //Unpack the loops iterations and define consecutive groups
        final_iter_ptr.clear();
        final_iter_ptr.reserve(nnz);
        int cnt = 0;
        std::vector<int> cnc_group_ptr;
        cnc_group_ptr.push_back(0);
        final_tuple_ptr.push_back(cnt);
        for (int i1 = 0; i1 < final_level_no; ++i1){
            for (int j1 = final_level_ptr[i1]; j1 < final_level_ptr[i1 + 1]; ++j1){
                int prev_i = -2;//Because we have i = 0 so -2 + 1 = -1 makes sense for initial condition
                for (int k1 = final_part_ptr[j1]; k1 < final_part_ptr[j1 + 1]; ++k1){
                    int i = final_node_ptr[k1];
                    //Body of the Code
                    if( i != (prev_i + 1)){//Current group
                        if(cnc_group_ptr.back() != cnt){
                            cnc_group_ptr.push_back(cnt);
                        }
                    }
                    prev_i = i;
                    for (int k = Lp[i]; k < Lp[i+1]; ++k) {
                        twoDIteration iteration;
                        iteration.i = i;
                        iteration.k = k;
                        iteration.j = Li[k];
                        final_iter_ptr.push_back(iteration);
                        cnt++;
                    }
                }
                if(cnc_group_ptr.back() != cnt){//|| cnt - cnc_group_ptr.back() > T1_size
                    cnc_group_ptr.push_back(cnt);
                }
                final_tuple_ptr.push_back(cnt);
            }
        }

//        int break_cnt = 0;
//        for(auto & iter: final_iter_ptr){
//            std::cout << iter.i << "\t" << iter.j << std::endl;
//            break_cnt++;
//            if(break_cnt == 20){
//                break;
//            }
//        }

        for(int g = 0; g < cnc_group_ptr.size() - 1; g++){
            std::stable_sort(final_iter_ptr.data() + cnc_group_ptr[g], final_iter_ptr.data() + cnc_group_ptr[g + 1],
                      sortbysec);
        }

//        break_cnt = 0;
//        std::cout << "New era" << std::endl;
//        for(auto & iter: final_iter_ptr){
//            std::cout << iter.i << "\t" << iter.j << std::endl;
//            break_cnt++;
//            if(break_cnt == 20){
//                break;
//            }
//        }

    }

    bool sortbyfirst(const twoDIteration& a, const twoDIteration& b){
        return (a.i < b.i);
    }

    bool sortbysec(const twoDIteration& a, const twoDIteration& b){
        return (a.j < b.j);
    }


    void getFinalScheduleDAG(int n, int final_level_no, const int* final_level_ptr,
                                const int* final_part_ptr, const int* final_node_ptr,
                                const int* orig_DAG_ptr, const int* orig_DAG_set,
                                int& ngroups, std::vector<int>& group_ptr, std::vector<int>& group_set,
                                std::vector<int>& DAG_ptr, std::vector<int>& DAG_set){

        if(!group_ptr.empty()){
            group_ptr.clear();
        }

        if(!group_set.empty()){
            group_set.clear();
        }

        if(!DAG_ptr.empty()){
            DAG_ptr.clear();
        }

        if(!DAG_set.empty()){
            DAG_set.clear();
        }
        //Create Group set
        int node_cnt = 0;
        group_ptr.push_back(node_cnt);
        group_set.resize(n, 0);
        for (int i1 = 0; i1 < final_level_no; ++i1)
        {
            // Iterate over all the w-partitions of a l-partition
            #pragma omp for schedule(auto)
            for (int j1 = final_level_ptr[i1]; j1 < final_level_ptr[i1 + 1]; ++j1)
            {
                // Iterate over all the node of a w-partition
                for (int k1 = final_part_ptr[j1]; k1 < final_part_ptr[j1 + 1]; ++k1)
                {
                    //Detect the node
                    int i = final_node_ptr[k1];
                    group_set[node_cnt++] = i;
                }
                if(group_ptr.back() != node_cnt){
                    group_ptr.push_back(node_cnt);
                }
            }
        }

        ngroups = group_ptr.size() - 1;
        assert(node_cnt == n);
        assert(group_ptr.back() == n);

        //Create DAG
        HDAGG::buildGroupDAG(n, ngroups, group_ptr.data(), group_set.data(),
                             orig_DAG_ptr, orig_DAG_set, DAG_ptr, DAG_set);
    }


 int build_coarsened_level_parallel(const int n, const int* Lp, const int* Li,
                                    const int cores,
                                    int& coarse_level_no,
                                    std::vector<int>& coarse_level_ptr,
                                    std::vector<int>& coarse_part_ptr,
                                    std::vector<int>& coarse_node_ptr
 ){
  int coarsen_factor, num_wparts=cores, i_cut=0;
  bool bin_pack=false;
  bool postOrder=false;//set to false manually here
  sym_lib::lbc_config(n, Lp[n], cores, num_wparts, coarsen_factor, i_cut, bin_pack);
  if(coarsen_factor > 0){
   std::vector<int> level_ptr;
   std::vector<int> level_set, node_to_level, WM;
   level_ptr.resize(n + 1);
   level_set.resize(n);
   node_to_level.resize(n);
   int n_levels = HDAGG::build_levelSet_CSC(n, Lp, Li, level_ptr.data(), level_set.data(),
                                            node_to_level.data() );
   for (int i = 0; i < n_levels;) {
    WM.push_back(i);
    i+=coarsen_factor;
    if(i >= n_levels){
     WM.push_back(n_levels);
     break;
    }
   }
   HDAGG::computeSchedule(cores, n, Lp, Li, level_ptr.data(), level_set.data(),
                          node_to_level.data(), WM, {{0}},
                          coarse_level_no, coarse_level_ptr, coarse_part_ptr, coarse_node_ptr,
                          bin_pack, postOrder);
  } else{ // fall back to basic wavefront for very sparse DAGs
   std::vector<int> node_to_level;
   coarse_part_ptr.resize(n + 1);
   coarse_node_ptr.resize(n);
   node_to_level.resize(n);
   coarse_level_no = HDAGG::build_levelSet_CSC(n, Lp, Li, coarse_part_ptr.data(),
                                               coarse_node_ptr.data(), node_to_level.data());
   coarse_level_ptr.resize( coarse_level_no + 1);
   for (int i = 0; i < coarse_level_no+1; ++i) {
    coarse_level_ptr[i] = i;
   }
  }
  return 1;
 }


#ifdef SPMP
    int levelsetCSRParallel_SpMP(SpMP::CSR* A, int nthreads,
                                 std::vector<int>& level_ptr, std::vector<int>& level_set){
        int nlevels = 0;
        int *symRowPtr = NULL, *symColIdx = NULL, *symDiagPtr = NULL, *symExtPtr = NULL;
        bool wasSymmetric = getSymmetricNnzPattern(A, &symRowPtr, &symDiagPtr, &symExtPtr, &symColIdx);
       #ifdef ENABLE_OPENMP
 omp_set_num_threads(nthreads);
#endif
        SpMP::LevelSchedule *barrierSchedule = new SpMP::LevelSchedule;
        barrierSchedule->useBarrier = true;
        barrierSchedule->transitiveReduction = false;
        if (wasSymmetric) {
            FREE(symRowPtr);
            FREE(symColIdx);
            FREE(symDiagPtr);
            FREE(symExtPtr);

            barrierSchedule->constructTaskGraph(*A);
        }
        else {
            barrierSchedule->constructTaskGraph(
                    A->m, symRowPtr, symDiagPtr, symExtPtr, symColIdx,
                    SpMP::PrefixSumCostFunction(symRowPtr));
        }

        const std::vector<int>& threadBoundaries = barrierSchedule->threadBoundaries;
        const std::vector<int>& taskBoundaries = barrierSchedule->taskBoundaries;
        const std::vector<int>& levIndices = barrierSchedule->levIndices;
        auto lvl_set = barrierSchedule->threadContToOrigPerm;
        nlevels = barrierSchedule->levIndices.size() - 1;
        level_set.resize(A->n);
        int cnt = 0;
        for(int lvl = 0; lvl < nlevels; lvl++){
            for(int tid = 0; tid < nthreads; tid++){
                auto task = threadBoundaries[tid] + lvl;
                assert(task < threadBoundaries[tid + 1]);
                for (int i = taskBoundaries[task]; i < taskBoundaries[task + 1]; ++i) {
                    int row = lvl_set[i];
                    level_set[cnt++] = row;
                } // for each row
            }
        }
        level_ptr.resize(levIndices.size());
        std::copy(levIndices.begin(), levIndices.end(), level_ptr.begin());
        assert(level_ptr.back() == A->m);
        delete barrierSchedule;
        return nlevels;
    }

    void Convert_LBC_CSR_to_SpMP(sym_lib::CSR* LBC_A, SpMP::CSR* SpMP_A){
        auto n = LBC_A->n;
        auto nnz = LBC_A->nnz;
        // Allocate space and copy data
        SpMP_A->alloc(LBC_A->n, LBC_A->nnz, true);
        SpMP_A->n = n;
        std::copy(LBC_A->p, LBC_A->p + n + 1, SpMP_A->rowptr);
        std::copy(LBC_A->i, LBC_A->i + nnz, SpMP_A->colidx);
        std::copy(LBC_A->x, LBC_A->x + nnz, SpMP_A->values);

        // Initialize meta data for SpMP_A
        int base = 0;
        #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            for (int j = SpMP_A->rowptr[i] - base; j < SpMP_A->rowptr[i + 1] - base; ++j) {
                if (SpMP_A->colidx[j] - base == i) {
                    SpMP_A->diagptr[i] = j + base;
                    SpMP_A->idiag[i] = 1 / SpMP_A->values[j];
                    SpMP_A->diag[i] = SpMP_A->values[j];
                }
            }
        }
    }
#endif
}
