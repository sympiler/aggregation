//
// Created by kazem on 3/1/21.
//

#ifndef FUSION_DAGP_UTILS_H
#define FUSION_DAGP_UTILS_H


//#include <dagP.h>
//#include <dgraphReader.h>
//#include <info.h>
//#include <rvcycle.h>
#include "aggregation/sparse_inspector.h"
#include <set>


namespace sym_lib {


    CSC *convert_to_pgraph(const dgraph *G, const int *part) {
        int maxParts = -1, i = 1, j = 0;
        for (i = 1; i <= G->nVrtx; ++i) {
            if (maxParts < part[i])
                maxParts = part[i];
        }
        double **adj = new double *[maxParts + 1]();
        for (i = 0; i <= maxParts; ++i)
            adj[i] = new double[maxParts + 1]();
        for (i = 1; i <= G->nVrtx; ++i)
            for (j = G->outStart[i]; j <= G->outEnd[i]; ++j) {
                if (part[i] != part[G->out[j]])
                    //adj[part[i]][part[G->out[j]]] += G->ecOut[j];
                    // adj[part[G->out[j]]][part[i]] += G->ecOut[j];
                    adj[part[i]][part[G->out[j]]] += G->ecOut[j];
            }
        int ttt = 0;
#if 0
        for (i=0; i<=maxParts; i++) {
         for (j=0; j<=maxParts; j++) {
          if(adj[i][j] > 0) {
           std::cout << i << " -> " << j << "[" << adj[i][j] << ", " << adj[i][j] << "];\n";
           ttt++;
          }
         }
        }
#endif
        CSC *p_graph = dense_to_csc(maxParts + 1, maxParts + 1, adj);
        for (i = 0; i <= maxParts; ++i)
            delete[]adj[i];
        delete[]adj;
        return p_graph;
    }


/*
 * argv is required if pragma inside is enabled, by default will be ignored
 * L1_csc is the input graph in CSC form
 * n: number of columns
 * n_parts: specifies the number of partitions
 * parts: allocated inside, represent vertex to partition
 */
    CSC *dagp_partition(const CSC *L1_csc,
                        size_t n, int num_parts, int *&parts) {
        dgraph G;
        MLGP_option opt;
        idxType nbParts = num_parts;
        opt.runs = 1;
        dagP_init_parameters(&opt, nbParts); // initializes default parameters
        //dagP_init_filename(&opt, argv[1]); // initialize the input file name and then the  output file name
        //dagP_opt_reallocUBLB(&opt, nbParts);

        idxType N = n, NNZL = L1_csc->nnz;
        allocateDGraphData(&G, N, NNZL, 3);
        if (G.vw == NULL) {
            G.vw = (vwType *) malloc(sizeof(vwType) * (N + 1));
        }
        for (int i = 0; i <= N; i++)
            G.vw[i] = 1;

        loadFromCSC(&G, N, NNZL, L1_csc->p, L1_csc->i, NULL);

/*
  delete []tmp_lp;
  delete []tmp_li;
*/
        set_dgraph_info(&G);
        int maxindegree, minindegree, maxoutdegree, minoutdegree;
        double aveindegree, aveoutdegree;
        dgraph_info(&G, &maxindegree, &minindegree, &aveindegree, &maxoutdegree, &minoutdegree, &aveoutdegree);
        G.maxindegree = maxindegree;
        G.maxoutdegree = maxoutdegree;
        parts = (idxType *) calloc((G.nVrtx + 1), sizeof(idxType));
        if (parts == NULL)
            printf("Could not allocate `parts` array.\n");

        ecType x = dagP_partition_from_dgraph(&G, &opt, parts);

        printf("edge cut: %d\n", (int) x);

        CSC *p_graph = convert_to_pgraph(&G, parts);
        // print_csc(0,"P_graph \n",p_graph);
        dagP_free_option(&opt);
        dagP_free_graph(&G);
        return p_graph;
    }


    CSC *dagp_partition_from_file(char *const *argv, int num_parts, int *&parts) {
        dgraph G;
        MLGP_option opt;
        idxType nbParts = num_parts;
        dagP_init_parameters(&opt, nbParts); // initializes default parameters
        //dagP_init_filename(&opt, argv[1]); // initialize the input file name and then the  output file name
        //dagP_opt_reallocUBLB(&opt, nbParts);

        dagP_read_graph(argv[1], &G, &opt);
        parts = (idxType *) calloc((G.nVrtx + 1), sizeof(idxType));
        if (parts == NULL)
            printf("Could not allocate `parts` array.\n");

        ecType x = dagP_partition_from_dgraph(&G, &opt, parts);

        //printf("edge cut: %d\n", (int) x);

        CSC *p_graph = convert_to_pgraph(&G, parts);
        dagP_free_option(&opt);
        dagP_free_graph(&G);
        return p_graph;
    }


    int convert_partition_to_p_graph(CSC *G, int n_parts, int *parts) {
        //print_csc(1, "G:\n", G);
        std::vector<std::set<int>> partitioned_dag;
        partitioned_dag.resize(n_parts);
        for (int i = 0; i < G->n; ++i) {
            int src_p = parts[i + 1];
            assert(src_p < n_parts);
            for (int j = G->p[i] + 1; j < G->p[i + 1]; ++j) {
                int cur_v = G->i[j] + 1;
                int cur_p = parts[cur_v];
                assert(cur_p < n_parts);
                if (cur_p != src_p)
                    partitioned_dag[src_p].insert(cur_p);
            }
        }
        std::cout << "\n";
        for (auto k : partitioned_dag) {
            for (auto l : k) {
                std::cout << l << ", ";
            }
            std::cout << "\n";
        }
        return 0;
    }


    int build_levelSet_CSC_pgraph(size_t n, int *Lp, int *Li,
                                  int *&levelPtr, int *&levelSet) {
        int begin = 0, end = n - 1;
        int cur_level = 0, cur_levelCol = 0;
        levelPtr = new int[n + 1]();
        levelSet = new int[n]();
        int *inDegree = new int[n]();
        bool *visited = new bool[n]();
        for (int i = 0; i < n; ++i) {
            for (int j = Lp[i]; j < Lp[i + 1]; ++j) {
                int cn = Li[j];
                if (cn != i) {
                    inDegree[cn]++;
                }
            }
        }
        // for (int i = 0; i < Lp[n]; ++i) {//O(nnz)
        //  inDegree[Li[i]]++;
        // }
        //print_vec("dd\n",0,n,inDegree);
        while (begin <= end) {
            for (int i = begin; i <= end; ++i) {//For level cur_level
                if (inDegree[i] == 0 && !visited[i]) {//if no incoming edge
                    visited[i] = true;
                    levelSet[cur_levelCol] = i; //add it to current level
                    cur_levelCol++;//Adding to level-set
                }
            }
            cur_level++;//all nodes_ with zero indegree are processed.
            //assert(cur_level < n);
            if (cur_level > n) {
                std::cout << cur_level << ", " << n << std::endl;
                return -1; // The input graph has a cycle
            }
            levelPtr[cur_level] = cur_levelCol;
            while (inDegree[begin] == 0) {
                begin++;
                if (begin >= n)
                    break;
            }
            while (inDegree[end] == 0 && begin <= end)
                end--;
            //Updating degrees after removing the nodes_
            for (int l = levelPtr[cur_level - 1]; l < levelPtr[cur_level]; ++l) {
                int cc = levelSet[l];
                for (int j = Lp[cc]; j < Lp[cc + 1]; ++j) {
                    if (Li[j] != cc) //skip diagonals
                        inDegree[Li[j]]--;//removing corresponding edges
                }
            }
            //print_vec("dd\n",0,n,inDegree);
        }
        delete[]inDegree;
        delete[]visited;
        return cur_level;//return number of levels
    }

    int get_h_level_set(int n_vertices, CSC *p_graph, int n_parts,
                        const int *parts, int &n_levels, int *&h_level_ptr,
                        int *&h_par_ptr, int *&h_partition) {
        int *level_ptr, *level_set;
        n_levels = build_levelSet_CSC_pgraph(p_graph->n, p_graph->p, p_graph->i, level_ptr, level_set);
        h_level_ptr = new int[n_levels + 1]();
        h_par_ptr = new int[n_parts + 1]();
        h_partition = new int[n_vertices]();
        std::vector<std::vector<int>> partitions(n_parts);
        for (int i = 1; i < n_vertices + 1; ++i) {
            auto c_p = parts[i];
            partitions[c_p].push_back(i - 1);
        }
        int c_part = 1, c_vert = 0;
        for (int j = 0; j < n_levels; ++j) {
            h_level_ptr[j + 1] = level_ptr[j + 1];
            for (int i = level_ptr[j]; i < level_ptr[j + 1]; ++i) {
                auto tmp_c = level_set[i];
                h_par_ptr[c_part] += h_par_ptr[c_part - 1] + partitions[tmp_c].size();
                for (int k = 0; k < partitions[tmp_c].size(); ++k) {
                    h_partition[c_vert] = partitions[tmp_c][k];
                    c_vert++;
                }
                c_part++;
            }
        }
        std::cout << c_part << ", " << n_parts << std::endl;
//  assert(c_part-1 == n_parts);
        assert(c_vert == n_vertices);
        delete[]level_ptr;
        delete[]level_set;
        return 1;
    }


    void get_h_level_set_V2(int n, CSC *org_graph, int n_parts,
                            const int *parts, int &n_levels, int *&h_level_ptr,
                            int *&h_par_ptr, int *&h_partition) {


        std::vector<std::vector<int>> partitions(n_parts, std::vector<int>());
        for (int i = 1; i < n + 1; ++i) {
            auto c_p = parts[i];
            assert(i - 1 >= 0 && i - 1 < n);
            assert(c_p < n_parts);
            partitions[c_p].push_back(i - 1);
        }

        std::vector<int> group_ptr;
        std::vector<int> group_set(n, 0);
        int cnt = 0;
        for (int i = 0; i < n_parts; i++) {
            group_ptr.push_back(cnt);
            for (auto &iter: partitions[i]) {
                group_set[cnt++] = iter;
            }
        }
        group_ptr.push_back(cnt);
        std::vector<int> DAG_ptr, DAG_set;
        int ngroups = group_ptr.size() - 1;
        GLC::buildGroupDAG(n, ngroups, group_ptr.data(), group_set.data(),
                           org_graph->p, org_graph->i, DAG_ptr, DAG_set);
        std::vector<int> level_ptr, level_set;

        level_ptr.resize(ngroups + 1);
        level_set.resize(ngroups);
        n_levels = build_levelSet_CSC_V2(ngroups, DAG_ptr.data(), DAG_set.data(),
                                         level_ptr.data(), level_set.data());

        h_level_ptr = new int[n_levels + 1]();
        h_par_ptr = new int[n_parts + 1]();
        h_partition = new int[n]();

        int c_part = 1, c_vert = 0;
        for (int j = 0; j < n_levels; ++j) {
            h_level_ptr[j + 1] = level_ptr[j + 1];
            for (int i = level_ptr[j]; i < level_ptr[j + 1]; ++i) {
                auto tmp_c = level_set[i];
                h_par_ptr[c_part] += h_par_ptr[c_part - 1] + partitions[tmp_c].size();
                for (int k = 0; k < partitions[tmp_c].size(); ++k) {
                    h_partition[c_vert] = partitions[tmp_c][k];
                    c_vert++;
                }
                c_part++;
            }
        }
        std::cout << c_part << ", " << n_parts << std::endl;


    }

}

#endif //FUSION_DAGP_UTILS_H
