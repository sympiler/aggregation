//
// Created by Kazem on 10/10/19.
//
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include "aggregation/sparse_io.h"
#include "aggregation/def.h"
#include "aggregation/mmio.h"
#include "aggregation/utils.h"
#include "aggregation/exceptions.h"

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

  A->stype = mm_is_symmetric(mcode) ? -1 : 0;

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

  CSC * convert_to_one_based(const CSC *A){
  CSC *A_one = new CSC(A->m+1, A->n+1, A->nnz, A->is_pattern);
  if(A->is_pattern)
   for (int j = 0; j < A->nnz; ++j) {
    A_one->i[j] = A->i[j]+1;
   }
  else
   for (int j = 0; j < A->nnz; ++j) {
    A_one->i[j] = A->i[j]+1;
    A_one->x[j] = A_one->x[j];
   }

  for (int k = 0; k < A->n + 1; ++k) {
   A_one->p[k+1] = A->p[k];
  }
  A_one->m--; A_one->n--;
  return A_one;
 }

 CSC *dense_to_csc(int rows, int cols, double **val){
  int n_nz = 0; int *nz_col = new int[cols]();
  for (int j = 0; j < cols; ++j) {
   for (int i = 0; i < rows; ++i) {
    if(val[i][j] != 0){
     n_nz++;
     nz_col[j]++;
    }
   }
  }
  CSC *out = new CSC(rows, cols, n_nz);
  for (int k = 0; k < cols; ++k) {
   out->p[k+1] = out->p[k] + nz_col[k];
  }
  assert(out->p[cols] == n_nz);
  n_nz=0;
  for (int j = 0; j < cols; ++j) {
   for (int i = 0; i < rows; ++i) {
    if(val[i][j] != 0){
     out->i[n_nz] = i;
     out->x[n_nz] = val[i][j];
     n_nz++;
    }
   }
  }
  delete []nz_col;
  return out;
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


 /// New IO

 std::string type_str(int type) {
  switch (type) {
   case REAL:
    return "REAL";
   case INT:
    return "INT";
   case COMPLEX:
    return "COMPLEX";
   case PATTERN:
    return "PATTERN";
   default:
    return "UNKNOWN";
  }
 }


 std::string format_str(int fmt) {
  switch(fmt) {
   case COORDINATE:
    return "COORDINATE";
   case ARRAY:
    return "ARRAY";
   default:
    return "UNKNOWN";
  }
 }

 void print_string(std::string context, std::streambuf* out = std::cout.rdbuf(),
                   std::string indent = "  "){
  std::streambuf* sb_cout_backup = std::cout.rdbuf();
  std::cout.rdbuf(out);
  std::stringstream ss(context);
  std::string line;
  while (std::getline(ss,line)){
   std::cout<<indent<<line<<"\n";
  }
  std::cout.rdbuf(sb_cout_backup);
 }

 void read_header(std::ifstream &inFile, int &n_row, int &n_col,
                  size_t &n_nnz, int &type, int &shape, int &mtx_format){
  std::string line,banner, mtx, crd, arith, sym;
  std::getline(inFile,line);
  trim(line);
  for (unsigned i=0; i<line.length(); line[i]=tolower(line[i]),i++);
  std::istringstream iss(line);
  if (!(iss >> banner >> mtx >> crd >> arith >> sym)){
   throw mtx_header_error("Unknown", "First line does not contain 5 tokens");
  }
  if(!(banner =="%%matrixmarket")) {
   throw mtx_header_error("Unknown", "first token is not \"%%%%MatrixMarket\"");
  }
  if(!(mtx =="matrix")) {
   throw mtx_header_error("Unknown", "Not a matrix, unable to handle");
  }
  if(crd == "coordinate") {
   mtx_format = COORDINATE;
  } else if(crd == "array") {
   mtx_format = ARRAY;
  } else{
   throw mtx_header_error("Unknown", "Unknown matrix format, unable to handle");
  }
  if(arith == "real")
   type = REAL;
  else if(arith == "integer")
   type = INT;
  else if (arith == "complex")
   type = COMPLEX;
  else if(arith == "pattern")
   type = PATTERN;
  else{
   throw mtx_header_error("Unknown",
                          "Unknown arithmetic, unable to handle");
  }
  if(sym == "symmetric")
   shape = LOWER;
  else if(sym == "general")
   shape = GENERAL;
  else{
   throw mtx_header_error("Unknown", "Unknown shape, unable to handle");
  }
  while (!line.compare(0,1,"%"))
  {
   std::getline(inFile, line);
   trim(line);
  }
  std::istringstream issDim(line);
  if(mtx_format != ARRAY){
   if (!(issDim >> n_row >> n_col >> n_nnz)){
    throw mtx_header_error("Unknown", "The matrix dimension is missing");
   }
  } else{
   if (!(issDim >> n_row >> n_col)){
    throw mtx_header_error("Unknown", "The matrix dimension is missing");
   }
   n_nnz = n_row*n_col;
  }
 }


 int shape2int(int shape){
  int st = 0;
  switch (shape) {
   case LOWER:
    st=-1;
    break;
   case UPPER:
    st =1;
    break;
   case GENERAL:
    st = 0;
    break;
   default:
    st=0;
  }
  return st;
 }


 void read_triplets_real(std::ifstream &inFile, int nnz,
                         std::vector<triplet>& triplet_vec, bool zero_indexing=false){
  for (int i = 0; i < nnz; ++i) {
   triplet tmp;
   inFile >> tmp.row;
   inFile >> tmp.col;
   inFile >> tmp.val;
   if(!zero_indexing){
    tmp.col--; tmp.row--;
   }
   triplet_vec.push_back(tmp);
  }
 }


 void compress_triplets_to_csc(std::vector<triplet>& triplet_vec, CSC *A,
                               bool add_diags= true){
  assert(A->nnz == triplet_vec.size());
  std::sort(triplet_vec.begin(), triplet_vec.end(),
            [](const triplet& a, const triplet& b){return (a.col<b.col) || (a.col==b.col && a.row<b.row);});
  auto *count = new int[A->n]();
  for (auto i = 0; i < A->nnz; ++i) {
   count[triplet_vec[i].col]++;
  }
  A->p[0] = 0;
  for (auto j = 0; j < A->n; ++j) {
   if(count[j] == 0 && add_diags){ // insert zero diag for empty cols
    triplet tmp; tmp.col = tmp.row = j; tmp.val=0;
    triplet_vec.insert(triplet_vec.begin()+A->p[j], tmp);
    A->p[j+1] = A->p[j] + 1;
   }else{
    A->p[j+1] = A->p[j] + count[j];
   }
  }
  delete []count;
  for (auto k = 0; k < A->nnz; ++k) {
   A->i[k] = triplet_vec[k].row;
   A->x[k] = triplet_vec[k].val;
  }
 }

 void read_mtx_csc_real(std::ifstream &in_file, CSC *&A, bool insert_diag){
  int n, m;
  int shape, arith, mtx_format;
  size_t nnz;
  std::vector<triplet> triplet_vec;

  read_header(in_file, m, n, nnz, arith, shape, mtx_format);
  if(arith != REAL)
   throw mtx_arith_error("REAL", type_str(arith));
  else if (mtx_format != COORDINATE)
   throw mtx_format_error("COORDINATE", format_str(mtx_format));
  A = new CSC(m,n,nnz,false, shape2int(shape));
  read_triplets_real(in_file, nnz, triplet_vec);
  compress_triplets_to_csc(triplet_vec, A, insert_diag);
  A->nnz = A->p[n]; // if insert diag is true, it will be different.
  //print_csc(A->n, A->p, A->i, A->x);
 }


 void read_mtx_array_real(std::ifstream &in_file, Dense *&A) {
  int n, m;
  int shape, arith, mtx_format;
  size_t nnz;
  read_header(in_file, m, n, nnz, arith, shape, mtx_format);
  if(arith != REAL)
   throw mtx_arith_error("REAL", type_str(arith));
  else if (mtx_format != ARRAY)
   throw mtx_format_error("ARRAY", format_str(mtx_format));
  A = new Dense(m, n, 1);//
  for (int i = 0; i < m * n; i++) {//writing from file row by row
   in_file >> A->a[i];
  }
  std::ofstream file;
  //print_dense(A->row, A->col, A->lda, A->a);
 }


 void load_mtx_array_real(std::string &filename, Dense *&A) {
  std::ifstream fin(filename);
  if(fin.is_open()){
   try {
    read_mtx_array_real(fin, A);
   } catch (const mtx_format_error& e) {
    throw mtx_format_error(e.expected_format(),
                           e.got_format(),
                           filename,
                           e.what());
   } catch (const mtx_arith_error& e) {
    throw mtx_arith_error(e.expected_arith(),
                          e.got_arith(),
                          filename,
                          e.what());
   } catch (const mtx_header_error& e) {
    throw mtx_header_error(filename, e.what());
   }
  } else {
   fin.close();
   throw read_file_error(filename);
  }
  fin.close();
 }


/// Reads a real constant from input
/// \param in_file
/// \param val
 void read_real_constant(std::ifstream &in_file, double &val){
  in_file >> val;
 }


 void read_string(std::ifstream &in_file, std::string &context, std::string indent = "  "){
  std::string line;
  while (!in_file.eof()){
   std::streampos oldpos = in_file.tellg();  // stores the position
   std::getline(in_file,line);
   if(line.compare(0, indent.size(), indent)){
    if(line[0] == '-') //TODO not sure if this is always true
     continue;
    in_file.seekg (oldpos);   // get back to the position
    break;
   } else{
    ltrim(line);
    line+="\n";
    context.append(line);
   }
  }
 }


}
