//
// Created by labuser (Bangtian Liu) on 9/16/20.
//

#include <Group.h>


namespace sym_lib
{

	void group::inspection_spicocsc_v1(int *groupPtr, int *groupSet, int &ngroup, int *groupInv)
	{
		for (int i = 0; i < ncol; ++i)
		{
			int len = mP[i + 1] - mP[i];
			child[i] = len > 1 ? mI[mP[i] + 1] : -1;
		}

		for (int i = 0; i < ncol; ++i)
		{
			int curIdx = i;
			int chIdx = child[curIdx];
			if (chIdx != -1 && chIdx - curIdx == 1)
			{
				next[curIdx] = chIdx;
				pre[chIdx] = curIdx;
			}
		}

		memset(visited, 0, sizeof(bool) * ncol);

		int pos = 0;
		int index = 0;

		for (int i = 0; i < ncol; ++i)
		{
			int curIdx = i;
			if (!visited[curIdx])
			{
				groupPtr[index] = pos;
				visited[curIdx] = true;
				groupSet[pos++] = curIdx;
				groupInv[curIdx] = index;
				int ch = next[curIdx];

				while (ch != 0 && !visited[ch])
				{
					curIdx = ch;
					groupSet[pos++] = curIdx;
					groupInv[curIdx] = index;
					visited[curIdx] = true;
					ch = next[curIdx];
				}
				index++;
			}
		}

		groupPtr[index] = pos;
		ngroup = index;
	}

	void group::inspection_spicocsc_v2(int *groupPtr, int *groupSet, int &ngroup, int *groupInv)
	{
		for (int i = 0; i < ncol; ++i)
		{
			int len = mP[i + 1] - mP[i];
			child[i] = len > 1 ? mI[mP[i] + 1] : -1;
		}

		for (int i = 0; i < ncol; ++i)
		{
			int curIdx = i;

			int chIdx = child[curIdx];

			if (chIdx != -1 && !visited[curIdx])
			{
				visited[curIdx] = true;
				if (status[curIdx] == 0 && status[chIdx] == 0)
				{

					next[curIdx] = chIdx;
					pre[chIdx] = curIdx;
					status[chIdx] = 1;

					curIdx = chIdx;
					visited[curIdx] = true;
					chIdx = child[curIdx];
					while (chIdx != -1 && (chIdx - curIdx) == 1 && !visited[chIdx])
					{
						visited[chIdx] = true;
						pre[chIdx] = curIdx;
						next[curIdx] = chIdx;
						status[chIdx] = 3;
						curIdx = chIdx;
						chIdx = child[curIdx];
					}

					int tcurIdx = curIdx + 1;
					int len = 0;
					while (child[tcurIdx] == -1 && !visited[tcurIdx] && len < ncol / 12)
					{
						visited[tcurIdx] = true;
						next[tcurIdx - 1] = tcurIdx;
						pre[tcurIdx] = tcurIdx - 1;
						status[tcurIdx] = 3;
						++tcurIdx;
						++len;
					}
				}
				else if (status[chIdx] == 1 && status[curIdx] == 0)
				{
					visited[chIdx] = true;
					next[curIdx] = chIdx;
					int tp = pre[chIdx];
					pre[chIdx] = curIdx;
					pre[curIdx] = tp;
					next[tp] = curIdx;
					status[curIdx] == 2;
				}
				else if (status[chIdx] == 3 && status[curIdx] == 0)
				{
					visited[chIdx] = true;
					next[curIdx] = chIdx;
					int tp = pre[chIdx];
					pre[chIdx] = curIdx;
					pre[curIdx] = tp;
					next[tp] = curIdx;
					status[curIdx] == 2;
				}
			}
		}

		memset(visited, 0, sizeof(bool) * ncol);

		int pos = 0;
		int index = 0;

		for (int i = 0; i < ncol; ++i)
		{
			int curIdx = i;
			if (!visited[curIdx])
			{
				groupPtr[index] = pos;
				visited[curIdx] = true;
				groupSet[pos++] = curIdx;
				groupInv[curIdx] = index;
				int ch = next[curIdx];

				while (ch != 0 && !visited[ch])
				{
					curIdx = ch;
					groupSet[pos++] = curIdx;
					groupInv[curIdx] = index;
					visited[curIdx] = true;
					ch = next[curIdx];
				}
				index++;
			}
		}

		groupPtr[index] = pos;
		ngroup = index;
	}

	void group::inspection_sptrsvcsr_v1(int *groupPtr, int *groupSet, int &ngroup, int *groupInv)
	{
		for (int i = 0; i < ncol; ++i)
		{
			int len = mP[i + 1] - mP[i];
			child[i] = len > 1 ? mI[mP[i + 1] - 2] : -1; //There is always a diagonal element, so len > 1
		}
        // It starts from 1 because the first row always has 1 element (the diagonal  one)
		for (int i = 1; i < ncol; ++i)
		{
			int curIdx = i;
			int chIdx = child[curIdx];
			if (chIdx != -1 && curIdx - chIdx == 1)
			{
				next[chIdx] = curIdx;
				pre[curIdx] = chIdx;
			}
		}

		memset(visited, 0, sizeof(bool) * ncol);

		int pos = 0;
		int index = 0;

		for (int i = 0; i < ncol; ++i)
		{
			int curIdx = i;
			if (!visited[curIdx])
			{
				groupPtr[index] = pos;
				visited[curIdx] = true;
				groupSet[pos++] = curIdx;
				groupInv[curIdx] = index;
				int ch = next[curIdx];

				while (ch != 0 && !visited[ch])
				{
					curIdx = ch;
					groupSet[pos++] = curIdx;
					groupInv[curIdx] = index;
					visited[curIdx] = true;
					ch = next[curIdx];
				}
				index++;
			}
		}

		groupPtr[index] = pos;
		ngroup = index;
	}



	void group::inspection_sptrsvcsr_v2(int *groupPtr, int *groupSet, int &ngroup, int *groupInv)
	{
		for (int i = 0; i < ncol; ++i)
		{
			child[i] = -1;
		}

		bool *rCnt = (bool *)malloc(sizeof(bool) * ncol);
		memset(rCnt, 0, sizeof(bool) * ncol);

		for (int i = 0; i < ncol; ++i)
		{
			int len = mP[i + 1] - mP[i];
			if (len == 1)
				rCnt[i] = true;
			int tidx = len > 1 ? mI[mP[i + 1] - 2] : -1;
			if (tidx != -1)
			{
				if (child[tidx] == -1)
					child[tidx] = i;
			}
		}

		for (int i = 0; i < ncol; ++i)
		{
			int curIdx = i;

			if (!visited[curIdx])
			{
				visited[curIdx] = true;
				int chIdx = child[curIdx];

				while ((chIdx - curIdx == 1 || rCnt[chIdx]) && !visited[chIdx])
				{
					visited[chIdx] = true;
					pre[chIdx] = curIdx;
					next[curIdx] = chIdx;

					curIdx = chIdx;
					chIdx = child[curIdx];
				}
				//
			}
		}

		memset(visited, 0, sizeof(bool) * ncol);

		int pos = 0;
		int index = 0;
		int len;
		for (int i = 0; i < ncol; ++i)
		{
			int curIdx = i;
			len = 0;
			if (!visited[curIdx])
			{
				groupPtr[index] = pos;
				visited[curIdx] = true;
				groupSet[pos++] = curIdx;
				groupInv[curIdx] = index;
				int ch = next[curIdx];
				++len;
				while (ch != 0 && !visited[ch])
				{
					curIdx = ch;
					groupSet[pos++] = curIdx;
					++len;
					groupInv[curIdx] = index;
					visited[curIdx] = true;
					ch = next[curIdx];
				}
				index++;
			}
		}

		groupPtr[index] = pos;
		ngroup = index;
	}

	void group::NaiveGrouping(int n, int *groupPtr, int *groupSet, int &ngroup, int *ginv, int blksize)
	{
		int index = 0;
		int pos = 0;
		for (int i = 0; i < n; ++i)
		{
			if (i % blksize == 0)
				groupPtr[index++] = pos;
			groupSet[pos++] = i;
			ginv[i] = index - 1;
		}

		ngroup = index;
		groupPtr[index] = n;
	}

	// ============================ Greedy Grouping =============================
    GreedyGrouping::GreedyGrouping(int cols_rows, int* csc_ptr, int* csc_idx,
                                   int max_search_space_size, double minimum_group_cost_threshold,
                                   bool merge_less_th){
	    this-> merge_less_th = merge_less_th;
	    this->ptr_size=cols_rows;
	    this->csc_ptr=csc_ptr;
	    this->csc_idx=csc_idx;
	    this->current_candidates.reserve(cols_rows);
	    this->group_set.reserve(cols_rows + 1);
	    this->node_group_idx.resize(cols_rows);
	    this->total_num_child = 0;
	    this->max_search_space_size = max_search_space_size;
	    this->minimum_group_cost_threshold = minimum_group_cost_threshold;
	    this->group_size = 0;
	}

    GreedyGrouping::~GreedyGrouping() {}

    void GreedyGrouping::startGrouping(){
	    int current_group_idx = 0;
	    int group_head = 0; // Columns 0
	    group_set.push_back(group_head);
	    int next_candidate;
        while( group_head < ptr_size){//There is sill a possibility to advance
            addChildesToCandidate(group_head);
            node_group_idx[group_head] = current_group_idx;
            next_candidate = getBestCandidate(group_head);
            if(next_candidate - group_head == 1){
                group_head = next_candidate;
            } else { //Start a new group
                current_group_idx++;
                group_head++;
                group_set.push_back(group_head);
                //Start fresh
                candidate_pool.clear();
                current_candidates.clear();
            }
        }
        // group_size is needed in the following functions
        group_set.shrink_to_fit();
        group_size = group_set.size() - 1;
        // This function will recompute the group_size
//        if(merge_less_th){
//            mergeSmallGroups();
//        }
        finalizeGroupingProcedure();
	}

	int GreedyGrouping::getBestCandidate(int group_head) {
        if(!current_candidates.empty()) {
            //delete candidates less than group_head
            while(!current_candidates.empty()){
                auto delete_iter = current_candidates.end();
                if(current_candidates.back()<= group_head){
                    delete_iter--;
                } else {
                    break;
                }
                current_candidates.erase(delete_iter, current_candidates.end());
            }
            // If nothing left, return
            if(current_candidates.empty()) return -1;

            if (current_candidates.back() - group_head == 1) {//if the next best candidate is the next columns
                int best_candidate = current_candidates.back();
                current_candidates.pop_back();
                return best_candidate;
            } else {
                for (auto iter: current_candidates) {
                    candidate_pool.insert(iter);
                }
                current_candidates.clear();
                current_candidates.reserve(candidate_pool.size());
                current_candidates.assign(candidate_pool.begin(), candidate_pool.end());
                //Sort in descending order
                std::sort(current_candidates.begin(),
                          current_candidates.end(), std::greater<int>());//We can use radix sort too

                //delete candidates less than group_head
                while(!current_candidates.empty()){
                    auto delete_iter = current_candidates.end();
                    if(current_candidates.back() <= group_head){
                        delete_iter--;
                    } else {
                        break;
                    }
                    current_candidates.erase(delete_iter, current_candidates.end());
                }
                if(current_candidates.empty()) return -1;
                int best_candidate = current_candidates.back();
                current_candidates.pop_back();
                return best_candidate;
            }
        }else{
            if(!candidate_pool.empty()){
                current_candidates.assign(candidate_pool.begin(), candidate_pool.end());
                //Sort in descending order
                std::sort(current_candidates.begin(),
                          current_candidates.end(), std::greater<int>());//We can use radix sort too

                //delete candidates less than group_head
                while(!current_candidates.empty()){
                    auto delete_iter = current_candidates.end();
                    if(current_candidates.back()<= group_head){
                        delete_iter--;
                    } else {
                        break;
                    }
                    current_candidates.erase(delete_iter, current_candidates.end());
                }
                if(current_candidates.empty()) return -1;

                int best_candidate = current_candidates.back();
                current_candidates.pop_back();
                return best_candidate;
            } else {
                return -1;
            }
        }
	}

	void GreedyGrouping::addChildesToCandidate(int node) {
//	    int current_group_size = node - group_set.back() + 1;
	    int acceptable_nodes = max_search_space_size + node;
        for(int rows_iter = csc_ptr[node]; rows_iter < csc_ptr[node + 1]; rows_iter++){
            if(node < csc_idx[rows_iter] && csc_idx[rows_iter] <= acceptable_nodes){
                candidate_pool.insert(csc_idx[rows_iter]);
            }
        }
	}

	void GreedyGrouping::finalizeGroupingProcedure(){
        //the last column in group_set is for showing the group_set(last_group + 1) - group_set(last_group)
        group_children.reserve(group_size);
        for(int group_idx = 0; group_idx < group_size; group_idx++){
            int group_begin = group_set[group_idx];
            int group_end = group_set[group_idx + 1];
            std::unordered_set<int> child_buff;
            child_buff.reserve(csc_ptr[group_end] - csc_ptr[group_begin]);
            //add children of the new aggregated node group_idx
            for(int rows_iter = csc_ptr[group_begin]; rows_iter < csc_ptr[group_end]; rows_iter++){
                if(csc_idx[rows_iter] >= group_end) {
                    child_buff.insert(node_group_idx[csc_idx[rows_iter]]);
                }
            }
            std::vector<int> children;
            children.assign(child_buff.begin(), child_buff.end());
            children.push_back(group_idx);//We basically need to fill the diagonal elements too
            std::sort(children.begin(), children.end());
            group_children.push_back(children);
            total_num_child += children.size();
        }
        assert(group_children.size()==group_set.size() - 1);
        group_nnz = total_num_child;
	}

    void GreedyGrouping::getGroupedCSCFormat(int* group_csc_ptr, int* group_csc_idx,
                                             int group_csc_ptr_size, int group_csc_idx_size){
        assert(group_csc_ptr_size == this->group_size + 1);
        assert(group_csc_idx_size == this->group_nnz );
        group_csc_ptr[0] = 0;
        for(int group_idx = 0; group_idx < group_size; group_idx++){
            group_csc_ptr[group_idx + 1] = group_children[group_idx].size() + group_csc_ptr[group_idx];
            std::copy(group_children[group_idx].begin(), group_children[group_idx].end(),
                      group_csc_idx + group_csc_ptr[group_idx]);
        }
	}

	void GreedyGrouping::getGroupedRangeArray(int *group_range, int group_range_size) {
	    assert(group_range_size==group_size + 1);
	    std::copy(group_set.begin(), group_set.end(), group_range);
	}

	double GreedyGrouping::getGroupCost(int group_idx) {
	    if(group_idx >= group_size || group_idx < 0){
	        return INFINITY;
	    }
        int group_begin = group_set[group_idx];
        int group_end = group_set[group_idx + 1];
        return csc_ptr[group_end] - csc_ptr[group_begin];
	}

	void GreedyGrouping::mergeSmallGroups() {
        // No merging is possible
	    if(group_size == 1){return;}
	    // We are merging nodes like this: if we merge node to its previous column
	    // Then for the next group, the previous full group will be the previous column of current group.
	    // The next group for a current group is always valid.
	    int pre_group = -1;
	    int next_group;
        for(int current_group = 0; current_group < group_size; current_group++){
            next_group = current_group + 1;
            double current_group_cost = getGroupCost(current_group);
            double next_group_cost = getGroupCost(next_group);
            double pre_group_cost = getGroupCost(pre_group);
            //Check whether the current group should merge into another group
            if(current_group_cost < minimum_group_cost_threshold){
                // Merge into the next or previous column (whichever has lower cost)
                if(pre_group_cost < next_group_cost){// Merge into the previous group
                    //Extend the previous group boundary
                    // Make the current group empty an re-adjust the boundary of all the
                    // empty spaces between the current group and the previous group
                    for(int gap = pre_group + 1; gap <= current_group; gap++){
                        group_set[gap] = group_set[next_group];
                    }
                } else {
                    // Merge into the next group
                    group_set[next_group] = group_set[current_group];
                }
            } else {
                // no merge happened
                pre_group = current_group;
            }
        }
        // TODO: Make this function a bit better, there is no need for this group_set_tmp
        // TODO: we can change the function above to have the same result with one group_set
        // TODO: Also we can reuse the pre_cost in the above chunk of the code
        std::vector<int> group_set_tmp (group_set.begin(), group_set.end());
        group_set.clear();
        group_set.reserve(group_size + 1);
        for(int i = 0; i < group_size; i++){
            if(group_set_tmp[i] < group_set_tmp[i + 1]){
                group_set.push_back(group_set_tmp[i]);
            }
        }
        group_set.push_back(group_set_tmp.back());
        group_set.shrink_to_fit();
        group_size = group_set.size() - 1;

        //Update the node_group_index
        for(int group_idx = 0; group_idx < group_size; group_idx++){
            for(int node = group_set[group_idx]; node < group_set[group_idx + 1]; node++){
                node_group_idx[node] = group_idx;
            }
        }


    }
    // ============================ Collection Grouping =============================
    CollectionGrouping::CollectionGrouping(int num_nodes, int *DAG_nodes, int *DAG_edges) {
        this->input_DAG_num_nodes = num_nodes;
        this->input_DAG_num_edges = DAG_nodes[num_nodes]; //The nnz/edges is in the last element
        this->input_DAG_nodes = DAG_nodes;
        this->input_DAG_edges = DAG_edges;

        grouped_DAG_num_edges = 0;
        grouped_DAG_num_nodes = 0;
        this->group_set.resize(num_nodes);
        this->group_set_ptr.resize(num_nodes + 1);
        this->grouped_DAG_edges.resize(input_DAG_num_edges);
        this->grouped_DAG_nodes.resize(num_nodes + 1);
        this->node_to_group.resize(num_nodes);
        grouping_happened = false;
	}

    //*NOTE:The DAG should be like CSC format
	void CollectionGrouping::startCNGrouping(){
	    if(grouping_happened){
	        std::cout << "Grouping Happened before!" << std::endl;
	        return;
	    }
	    grouping_time.start_timer();
        grouping_happened = true;
	    group_set_ptr[0] = 0;
	    this->grouped_DAG_num_nodes = 0;
        for (int node = 0; node < input_DAG_num_nodes; ++node)
        {
            group_set[node] = node;
            int num_childs = input_DAG_nodes[node + 1] - input_DAG_nodes[node];
            int child_idx = num_childs > 1 ? input_DAG_edges[input_DAG_nodes[node] + 1] : -1;
            node_to_group[node] = grouped_DAG_num_nodes;
            if (child_idx != -1 && child_idx - node == 1)
            {
                group_set_ptr[grouped_DAG_num_nodes + 1]++;
            } else {
                group_set_ptr[grouped_DAG_num_nodes + 1] += group_set_ptr[grouped_DAG_num_nodes];
                group_set_ptr[grouped_DAG_num_nodes + 1]++;//We add the effect of the current node at the end
                grouped_DAG_num_nodes++;
            }
        }
        grouping_time.measure_elapsed_time();
        cmp_DAG_creation_time.start_timer();

        //Building the Compressed DAG
        grouped_DAG_nodes[0] = 0;
        for(int cmp_node = 0; cmp_node < grouped_DAG_num_nodes; cmp_node++){
            int group_begin = group_set[group_set_ptr[cmp_node]];
            int group_end = group_set[group_set_ptr[cmp_node + 1]];
            std::unordered_set<int> child_buff;
            child_buff.reserve(input_DAG_nodes[group_end] - input_DAG_nodes[group_begin]);
            //add children of the new aggregated node group_idx
            for(int child_idx = input_DAG_nodes[group_begin]; child_idx < input_DAG_nodes[group_end]; child_idx++){
                if(input_DAG_edges[child_idx] >= group_end) {
                    child_buff.insert(node_to_group[input_DAG_edges[child_idx]]);
                }
            }
            std::vector<int> children;
            child_buff.insert(cmp_node);
            grouped_DAG_nodes[cmp_node + 1] = grouped_DAG_nodes[cmp_node] + child_buff.size();
            int start_child_idx = grouped_DAG_nodes[cmp_node];
            int end_child_idx = grouped_DAG_nodes[cmp_node + 1];
            std::copy(child_buff.begin(), child_buff.end(), &grouped_DAG_edges[start_child_idx]);
            std::sort(&grouped_DAG_edges[start_child_idx], &grouped_DAG_edges[end_child_idx]);
        }
        grouped_DAG_num_edges = grouped_DAG_nodes[grouped_DAG_num_nodes];

        //Adjusting the size
        group_set_ptr.resize(grouped_DAG_num_nodes + 1);
        grouped_DAG_nodes.resize(grouped_DAG_num_nodes + 1);
        grouped_DAG_edges.resize(grouped_DAG_num_edges);
        cmp_DAG_creation_time.measure_elapsed_time();
    }
    //*NOTE:The DAG should be like CSC format
    void CollectionGrouping::startSingleChildGrouping(){
        if(grouping_happened){
            std::cout << "Grouping Happened before!" << std::endl;
            return;
        }
        grouping_time.start_timer();
        grouping_happened = true;
        //0 -> 1, 0 is the parent and 1 is child
        //Calculating number of incoming edges into a node
        std::vector<int> Input_Edges(input_DAG_num_nodes, 0);
        for(int node = 0; node < input_DAG_num_nodes; node++){
            for(int child_idx = input_DAG_nodes[node]; child_idx < input_DAG_nodes[node + 1]; child_idx++){
                Input_Edges[input_DAG_edges[child_idx]]++;
            }
        }
        //Start Grouping
        std::vector<int> grouped(input_DAG_num_nodes, false);
        group_set_ptr[0] = 0;
        this->grouped_DAG_num_nodes = 0;
        int group_head = 0;
        int group_set_cnt = 0;
        for (int node = 0; node < input_DAG_num_nodes; ++node) {
            if (grouped[node] == true) {
                continue;
            }
            group_head = node;
            grouped[group_head] = true;
            node_to_group[group_head] = grouped_DAG_num_nodes;
            group_set[group_set_cnt] = group_head;
            group_set_cnt++;
            group_set_ptr[grouped_DAG_num_nodes + 1]++;
            while (true) {
                int num_childs = input_DAG_nodes[group_head + 1] - input_DAG_nodes[group_head];
                int child = num_childs > 1 ? input_DAG_edges[input_DAG_nodes[node] + 1] : -1;
                int tmp_child = 0;
                if(child == -1){//Stop the bad memory access
                    tmp_child = 0;
                } else {
                    tmp_child = child;
                }
                if (child != -1 && num_childs == 1 && Input_Edges[tmp_child] == 1) {
                    group_head = child;
                    group_set[group_set_cnt] = group_head;
                    group_set_cnt++;
                    node_to_group[group_head] = grouped_DAG_num_nodes;
                    group_set_ptr[grouped_DAG_num_nodes + 1]++;
                    grouped[node] = true;
                } else {
                    group_set_ptr[grouped_DAG_num_nodes + 1] += group_set_ptr[grouped_DAG_num_nodes];
                    grouped_DAG_num_nodes++;
                    break;
                }
            }
        }
        assert(group_set_cnt == input_DAG_num_nodes);
        grouping_time.measure_elapsed_time();
        cmp_DAG_creation_time.start_timer();
        //Building the Compressed DAG
        grouped_DAG_nodes[0] = 0;
        for(int cmp_node = 0; cmp_node < grouped_DAG_num_nodes; cmp_node++){
            int group_begin = group_set[group_set_ptr[cmp_node]];
            int group_end = group_set[group_set_ptr[cmp_node + 1]];
            std::unordered_set<int> child_buff;
            child_buff.reserve(input_DAG_nodes[group_end] - input_DAG_nodes[group_begin]);
            //add children of the new aggregated node group_idx
            for(int child_idx = input_DAG_nodes[group_begin]; child_idx < input_DAG_nodes[group_end]; child_idx++){
                if(input_DAG_edges[child_idx] >= group_end) {
                    child_buff.insert(node_to_group[input_DAG_edges[child_idx]]);
                }
            }
            std::vector<int> children;
            child_buff.insert(cmp_node);
            grouped_DAG_nodes[cmp_node + 1] = grouped_DAG_nodes[cmp_node] + child_buff.size();
            int start_child_idx = grouped_DAG_nodes[cmp_node];
            int end_child_idx = grouped_DAG_nodes[cmp_node + 1];
            std::copy(child_buff.begin(), child_buff.end(), &grouped_DAG_edges[start_child_idx]);
            std::sort(&grouped_DAG_edges[start_child_idx], &grouped_DAG_edges[end_child_idx]);
        }
        grouped_DAG_num_edges = grouped_DAG_nodes[grouped_DAG_num_nodes];

        //Adjusting the size
        group_set_ptr.resize(grouped_DAG_num_nodes + 1);
        grouped_DAG_nodes.resize(grouped_DAG_num_nodes + 1);
        grouped_DAG_edges.resize(grouped_DAG_num_edges);
        cmp_DAG_creation_time.measure_elapsed_time();

    }

    void CollectionGrouping::getGroupedDAG(std::vector<int>& DAG_nodes, std::vector<int>& DAG_edges) {
	    assert(DAG_edges.size() >= grouped_DAG_edges.size());
        assert(DAG_nodes.size() >= grouped_DAG_nodes.size());
        std::copy(grouped_DAG_nodes.begin(), grouped_DAG_nodes.end(), DAG_nodes.begin());
        std::copy(grouped_DAG_edges.begin(), grouped_DAG_edges.end(), DAG_edges.begin());
	}

	void CollectionGrouping::getGroups(std::vector<int>& groupSet, std::vector<int>& groupPtr) {
        assert(groupSet.size() >= group_set.size());
        assert(groupPtr.size() >= group_set_ptr.size());
        std::copy(group_set.begin(), group_set.end(), groupSet.begin());
        std::copy(group_set_ptr.begin(), group_set_ptr.end(), groupPtr.begin());
	}



}



