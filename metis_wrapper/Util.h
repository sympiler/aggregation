//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_UTIL_H
#define CHOLOPENMP_UTIL_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iomanip>
#include "def.h"

/*
 * reading a CSC matrix from a coordinate file, stored col-ordered
 */
bool readMatrix(std::string fName, size_t &n, size_t &NNZ, int *&col,//FIXME change col type to size_t
                int *&row, double *&val) {
    /*This function reads the input matrix from "fName" file and
     * allocate memory for matrix A, L and U.
     * - The input file is a coordinate version and e
     * ach row of the file shows (col, row, nnz)
     * - The matrices are zero-indexed
     */

    std::ifstream inFile;
    inFile.open(fName);
    std::string line, banner, mtx, crd, arith, sym;
    /*  File format:
     *    %%MatrixMarket matrix coordinate real general/symmetric/...
     *    % ...
     *    % (optional comments)
     *    % ...
     *    #rows    #non-zero
     *    Triplet in the rest of lines: row    col    value
     */
    std::getline(inFile, line);
    for (unsigned i = 0; i < line.length(); line[i] = tolower(line[i]), i++);
    std::istringstream iss(line);
    if (!(iss >> banner >> mtx >> crd >> arith >> sym)) {
        std::cout << "Invalid header (first line does not contain 5 tokens)\n";
        return false;
    }

    if (banner.compare("%%matrixmarket")) {
        std::cout << "Invalid header (first token is not \"%%%%MatrixMarket\")\n";
        return false;
    }
    if (mtx.compare("matrix")) {
        std::cout << "Not a matrix; this driver cannot handle that.\"\n";
        return false;
    }
    if (crd.compare("coordinate")) {
        std::cout << "Not in coordinate format; this driver cannot handle that.\"\n";
        return false;
    }
    if (arith.compare("real") and arith.compare("integer")) {
        if (!arith.compare("complex")) {
            std::cout << "Complex matrix; use zreadMM instead!\n";
            return false;
        } else if (!arith.compare("pattern")) {
            std::cout << "Pattern matrix; values are needed!\n";
            return false;
        } else {
            std::cout << "Unknown arithmetic\n";
            return false;
        }
    }
    while (!line.compare(0, 1, "%")) {
        std::getline(inFile, line);
    }
    std::istringstream issDim(line);
    if (!(issDim >> n >> n >> NNZ)) {
        std::cout << "The matrix dimension is missing\n";
        return false;
    }
    if (n <= 0 || NNZ <= 0)
        return false;
    col = new int[n + 1]();
    // colL = new int[n + 1]; colU = new int[n + 1];
    row = new int[NNZ];
    // rowL = new int[factorSize]; rowU = new int[factorSize];
    val = new double[NNZ];
    // valL = new double[factorSize]; valU = new double[factorSize];
    if (!val || !col || !row)
        return false;
    //Initializing the result vector
    int y, x, colCnt = 0, nnzCnt = 0;
    double value;

    col[0] = 0;
    int sw = 1;
    for (int i = 0; nnzCnt < NNZ;) {//Reading from file row by row
        inFile >> x;
        x--;
        inFile >> y;
        y--;//zero indexing
        inFile >> value;
        if (i == 0 && sw) {//matrix starts with empty cols
            if (y != 0) {
                i = y;
            }
            sw = 0;
        }
        if (y > n)
            return false;
        if (y == i) {
            val[nnzCnt] = value;
            row[nnzCnt] = x;
            colCnt++;
            nnzCnt++;
        } else {//New col
            int nnz_c = col[i] + colCnt;
            col[i + 1] = nnz_c;
            i++;//next iteration
            for (int j = i; j < y; ++j) {//empty cols
                col[j + 1] = nnz_c;
                i++;
            }
            colCnt = 1;
            val[nnzCnt] = value;
            row[nnzCnt] = x;
            nnzCnt++;
        }
    }
    //col[y+1] = col[y] + colCnt;
    if (y == n - 1) {//each column has something in it
        col[n] = col[n - 1] + colCnt;//last col
    } else { // we have some empty columns
        int nnzLast = col[y] + colCnt;
        col[y + 1] = nnzLast;
        y++;
        assert(nnzLast == NNZ);
        for (int i = y + 1; i < n + 1; ++i) {
            col[i] = nnzLast;
        }
    }

    return true;
}

#endif //CHOLOPENMP_UTIL_H
