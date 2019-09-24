//
// Created by george on 2019-09-24.
//

#ifndef LBC_MERGEGRAPH_H
#define LBC_MERGEGRAPH_H

/**
 * Assumes sequential dependency of the graphs
 * Incorporate with graph class later
 *
 * @param ngraphs
 * @param n
 * @param Gp
 * @param Gi
 */
void merge_graph(int ngraphs, int n, int **Gps, int **Gis, int *&nGp, int *&nGi);

#endif //LBC_MERGEGRAPH_H
