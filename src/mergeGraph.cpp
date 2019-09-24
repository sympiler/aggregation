//
// Created by george on 2019-09-24.
//

#include <assert.h>
#include <iostream>

void merge_graph(int ngraphs, int n, int **Gps, int **Gis, int *&nGp, int *&nGi) {
    const int *Gp, *Gi;

    int nnz = 0;
    for(int i = 0; i < ngraphs; i++)
        nnz += Gps[i][n];
    nnz += (ngraphs-1) * n;
    /** allocate new graph space **/
    nGp = new int[ngraphs * n + 1]();
    nGi = new int[nnz]();

    int p_counter = 1;
    int i_counter = 0;
    nGp[0] = 0;
    /** first <ngraphs>-1 graphs**/
    for(int i = 0; i < ngraphs-1; i++) {
       Gp = Gps[i];
       Gi = Gis[i];

       for(int j = 0; j < n; j++) {
           int diff = Gp[j+1] - Gp[j] + 1;
           nGp[p_counter] = nGp[p_counter-1] + diff;
           p_counter++;
           for(int p = Gp[j]; p < Gp[j+1]; p++) {
               int row = Gi[p] + (i * n);
               nGi[i_counter] = row;
               i_counter++;
           }
           nGi[i_counter] = j + (i+1) * n;
           i_counter++;
       }
    }
    /** last graph **/
    Gp = Gps[ngraphs-1];
    Gi = Gis[ngraphs-1];
    for(int j = 0; j < n; j++) {
        int diff = Gp[j+1] - Gp[j];
        nGp[p_counter] = nGp[p_counter-1] + diff;
        p_counter++;
        for(int p = Gp[j]; p < Gp[j+1]; p++) {
            int row = Gi[p] + (ngraphs-1) * n;
            nGi[i_counter] = row;
            i_counter++;
        }
    }
    assert(p_counter == ngraphs * n+1);
    assert(nGp[ngraphs*n] == nnz);
    assert(i_counter == nnz);
}