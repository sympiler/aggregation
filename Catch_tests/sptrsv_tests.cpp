//
// Created by george on 2019-10-09.
//
#include <cstring>
#include <utility>
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "sptrsv.h"
#include <cstdio>
#include <dirent.h>
#include <iostream>
#ifdef ENABLE_OPENMP
#include <omp.h>
#endif
#include "aggregation/sparse_inspector.h"
#include "aggregation/sparse_io.h"
#include "aggregation/sparse_utilities.h"
#include <sys/types.h>
#include "aggregation/test_utils.h"

using namespace sym_lib;
using namespace std;

static string dataPath;

bool testSptrsvCsr(CSC *L) {
 int *level_set, *level_ptr, level_no;

 level_no = build_levelSet_CSC(L->n, L->p, L->i, level_ptr, level_set);
 CSR *L_csr = csc_to_csr(L);

 auto *y_serial = new double[L->n]();
 auto *y_parallel = new double[L->n]();
 std::fill_n(y_serial, L->n, 1.0);
 sptrsv_csr(L_csr->n, L_csr->p, L_csr->i, L_csr->x, y_serial);

 std::fill_n(y_parallel, L->n, 1.0);
 sptrsv_csr_levelset(L_csr->n, L_csr->p, L_csr->i, L_csr->x, y_parallel,
                     level_no, level_ptr, level_set);

 bool is_correct = is_equal(0, L->n, y_serial, y_parallel);

 delete[] y_serial;
 delete[] y_parallel;
 delete[] level_ptr;
 delete[] level_set;
 delete L_csr;

 return is_correct;
}

TEST_CASE("Check lower triangular cases", "[sptrsvCorrectnessChecks]") {
 SECTION("Sptrsv CSR, parallel vs. serial (random matrices)") {
  vector<pair<int, double>> configs = {
   make_pair(10, 0.2), make_pair(100, 0.03), make_pair(500, 0.004),
   make_pair(1000, 0.0005), make_pair(10000, .0001)};
  // omp_set_num_threads(1);
  for (auto i : configs) {
   size_t n = i.first;
   double density = i.second;

   CSC *A = random_square_sparse(n, density);
   CSC *L = make_half(A->n, A->p, A->i, A->x);

   CHECK(testSptrsvCsr(L));

   delete A;
   delete L;
  }
 }

 // For running this test, you need to set the data path
 SECTION("Sptrsv CSR, parallel vs serial ( using datapath )") {
  DIR *dataDir = opendir(dataPath.c_str());
  if (dataDir) {
   struct dirent *dp;
   while ((dp = readdir(dataDir)) != NULL) {
    // Don't check . and ..
    if (!strcmp(dp->d_name, ".") || !strcmp(dp->d_name, "..")) {
     continue;
    } else if (strcmp(dp->d_name + (strlen(dp->d_name) - 4), ".mtx")) {
     continue;
    }
    char matrixFile[1024];
    const char *slash = dataPath.at(dataPath.length() - 1) == '/' ? "" : "/";
    sprintf(matrixFile, "%s%s%s", dataPath.c_str(), slash, dp->d_name);

    cout << "Testing: " << matrixFile << endl;

    CSC *L = read_mtx(matrixFile);

    if (!L) {
     cout << "Error in reading matrix file: " << matrixFile << endl;
    } else {
     CHECK(testSptrsvCsr(L));
     delete L;
    }
   }
  } else {
   cout << "Data path not set, so [datapath] tests won't run." << endl;
  }
 }
}

int main(int argc, char *argv[]) {
 Catch::Session session;

 using namespace Catch::clara;
 auto cli = session.cli() | Opt(dataPath, "datapath")["--datapath"](
                             "path where matrix data resides");

 session.cli(cli);

 int returnCode = session.applyCommandLine(argc, argv);
 if (returnCode != 0) // Indicates a command line error
  return returnCode;

 return session.run();
}
