//
// Created by Kazem on 10/10/19.
//
#include <iostream>
#include "includes/sparse_io.h"
#include "def.h"
#include "external/includes/mmio.h"

namespace sym_lib {

 CSC *read_mtx(std::string fname) {
  FILE *mf = fopen(fname.c_str(), "r");
  if (!mf) exit(1);

  MM_typecode mcode;
  if (mm_read_banner(mf, &mcode) != 0) {
   std::cerr << "Error processing matrix banner\n";
   fclose(mf);
   return nullptr;
  }

  int m, n, nnz;
  if (mm_read_mtx_crd_size(mf, &m, &n, &nnz) != 0) exit(1);
  CSC *A = new CSC(m, n, nnz);
  int *J = new int[nnz]();

  // Copy matrix data into COO format
  for (int i = 0; i < nnz; i++) {
   if (fscanf(mf, "%d %d %lg\n", (A->i) + i, &J[i], (A->x) + i) == EOF) {
    std::cerr << "Failed to load matrix at " << i + 1 << "th element\n";
    fclose(mf);
    delete (A);
    return nullptr;
   }
   (A->i)[i]--;
   J[i]--;
  }
  int i = 0, j;
  int index = 0, cur;
  for (; i < nnz; i++) {
   A->p[index] = i;
   cur = J[i];
   for (j = i + 1; j < nnz; j++) {
    if (J[j] != cur)
     break;
    else
     i++;
   }
   index += 1;
  }
  A->p[n] = nnz;
  delete[]J;

  A->m = m;
  A->n = n;
  A->nnz = nnz;
  A->stype = -1;

  return A;
 }

 void CSC_to_mtx(std::string fname, CSC *A) {
  FILE *fp = fopen(fname.c_str(), "w");

  MM_typecode matcode;
  mm_initialize_typecode(&matcode);
  mm_set_matrix(&matcode);
  mm_set_coordinate(&matcode);
  mm_set_real(&matcode);

  mm_write_banner(fp, matcode);
  mm_write_mtx_crd_size(fp, A->m, A->n, A->nnz);

  for (int i = 0; i < A->n; i++)
   for (int j = A->p[i]; j < A->p[i + 1]; j++)
    fprintf(fp, "%d %d %10.3g\n", A->i[j] + 1, i + 1, A->x[j]);
  fclose(fp);
 }

 void BCSC_to_mtx(std::string fname, BCSC *A) {
  FILE *f = fopen(fname.c_str(), "w");

  MM_typecode matcode;
  mm_initialize_typecode(&matcode);
  mm_set_matrix(&matcode);
  mm_set_coordinate(&matcode);
  mm_set_real(&matcode);

  mm_write_banner(f, matcode);
  mm_write_mtx_crd_size(f, A->m, A->n, A->nnz);

  for (int i = 0; i < A->nodes; i++) {
   int index = A->p[i];
   int width = A->supernodes[i + 1] - A->supernodes[i];
   int nrows = (A->p[i + 1] - A->p[i]) / width;

   for (int j = 0; j < width; j++) {
    for (int k = 0; k < nrows; k++) {
     int pos = index + j * nrows + k;
     fprintf(f, "%d %d %10.3g\n", A->i[pos] + 1, A->supernodes[i] + j + 1,
             A->x[pos]);
    }
   }
  }
  fclose(f);
 }

 void
 print_csc(int fd, std::string beg, size_t n, int *Ap, int *Ai, double *Ax) {
  dprintf(fd, "%s\n", beg.c_str());
  int nnz = n > 0 ? Ap[n] : 0;
  dprintf(fd, "%zu %zu %d\n", n, n, nnz);
  for (int i = 0; i < n; ++i) {
   for (int j = Ap[i]; j < Ap[i + 1]; ++j) {
    double x = Ax != NULLPNTR ? Ax[j] : 0;
    dprintf(fd, "%d %d %.12f", Ai[j] + 1, i + 1, x);
//    std::cout<<Ai[j]+1<<" "<<i+1<<" "<<std::setprecision(12)<< x;
    if (j + 1 != Ap[n])
     dprintf(fd, "\n");
   }
  }
 }

    void
    print_csr(int fd, std::string beg, size_t n, int *Ap, int *Ai, double *Ax) {
        dprintf(fd, "%s\n", beg.c_str());
        int nnz = n > 0 ? Ap[n] : 0;
        dprintf(fd, "%zu %zu %d\n", n, n, nnz);
        for (int i = 0; i < n; ++i) {
            for (int j = Ap[i]; j < Ap[i + 1]; ++j) {
                double x = Ax != NULLPNTR ? Ax[j] : 0;
                dprintf(fd, "%d %d %.12f", i + 1, Ai[j] + 1,  x);
//    std::cout<<Ai[j]+1<<" "<<i+1<<" "<<std::setprecision(12)<< x;
                if (j + 1 != Ap[n])
                    dprintf(fd, "\n");
            }
        }
    }




 void print_csc(int fd, std::string beg, CSC *A) {
  print_csc(fd, beg, A->n, A->p, A->i, A->x);
 }

 void print_dense(int fd, Dense *A) {
  for (int i = 0; i < A->row; i++) {
   for (int j = 0; j < A->col; j++) {
    dprintf(fd, "%f ", A->a[i + j * A->row]);
   }
   dprintf(fd, "\n");
  }
 }


 void print_level_set(std::string beg, int n, int *level_ptr, int *level_set) {
  std::cout << beg;
  for (int i = 0; i < n; ++i) {
   for (int j = level_ptr[i]; j < level_ptr[i + 1]; ++j) {
    std::cout << level_set[j] << ",";
   }
   std::cout << "\n";
  }
 }


 void print_hlevel_set(std::string beg, int n,
   const int *level_ptr, const int *level_part_ptr, int *level_set) {
  std::cout << beg;
  for (int i = 0; i < n; ++i) {
   for (int j = level_ptr[i]; j < level_ptr[i + 1]; ++j) {
    for (int k = level_part_ptr[j]; k < level_part_ptr[j + 1]; ++k) {
     std::cout << level_set[k] << ",";
    }
    std::cout << "; \n";
   }
   std::cout << "\n\n";
  }
 }

}