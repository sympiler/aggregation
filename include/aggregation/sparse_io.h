//
// Created by Kazem on 10/10/19.
//

#ifndef PROJECT_SPARSE_IO_H
#define PROJECT_SPARSE_IO_H

#include <iostream>
#include <string>
#include <iomanip>
#include "def.h"
namespace sym_lib{


///
/// \param fName
/// \return
 CSC* read_mtx(std::string fName);

 ///
 /// \param fname
 /// \param A
 void CSC_to_mtx(std::string fname, CSC *A);

 ///
 /// \param fname
 /// \param A
 void BCSC_to_mtx(std::string fname, BCSC *A);
 /// Converts dense to CSC
 /// \param rows in
 /// \param cols in
 /// \param val in
 /// \return
 CSC *dense_to_csc(int rows, int cols, double **val);

 /// Converts indices to one-based
 /// \param A input CSC
 /// \return
 CSC * convert_to_one_based(const CSC *A);

 ///
/// \tparam type
/// \param header the string at the beginning of output
/// \param beg the beginning index to print
/// \param end the last index to print
/// \param vec the vector of values
 template<class type> void print_vec(std::string header,
                                     int beg,int end, type *vec){

  std::cout<<header;
  for (int i = beg; i < end; ++i) {
   std::cout<<std::setprecision(15)<<vec[i]<<", ";
  }
  std::cout<<"\n";
 }


 /// copying vec_in to vec_out from beg to end
 /// \tparam type
 /// \param beg
 /// \param end
 /// \param vec_in
 /// \param vec_out
 template<class type> void copy_vector(int beg, int end, type *vec_in, type *vec_out){
  for (int i = beg; i < end; ++i) {
   vec_out[i] = vec_in[i];
  }
 }

/// Print CSC matrix into output
/// \param beg
/// \param n
/// \param Ap
/// \param Ai
/// \param Ax
 void print_csc(int fd, std::string beg, size_t n, int *Ap, int *Ai, double *Ax);
 void print_csc(int fd, std::string beg, CSC *A);
 /**
  * @brief Print CSR Matrix into CSR
  * @param fd
  * @param beg
  * @param n
  * @param Ap
  * @param Ai
  * @param Ax
  */
 void print_csr(int fd, std::string beg, size_t n, int *Ap, int *Ai, double *Ax);
    ///
 /// \param fd
 /// \param A
 void print_dense(int fd, Dense *A);

 /// Printing a level-set
 /// \param beg
 /// \param n
 /// \param level_ptr
 /// \param level_set
 void print_level_set(std::string beg, int n, int *level_ptr, int *level_set);

/// Prints an H-level set
/// \param beg
/// \param n
/// \param level_ptr
/// \param level_part_ptr
/// \param level_set
 void print_hlevel_set(std::string beg, int n,
                       const int *level_ptr, const int *level_part_ptr,
                       int *level_set);

//// NEW IO section. This will replace old read matrix gradually since this is more general.

#define print_precision  48
 enum TYPE{
  REAL,INT,COMPLEX,PATTERN
 };

 enum SHAPE{// LOWER and UPPER both are symmetric matrices.
  LOWER,UPPER,GENERAL
 };
 enum FORMAT{
  COORDINATE,ARRAY
 };

 struct triplet{
  int row{}; int col{}; double val{};
 };


 std::string type_str(int type);

 std::string format_str(int fmt);


 /*
  * Printing a compressed matrix, CSC to output
  * output can be stdout or a file (The out will be file.rdbuf())
  */
 template<class type> void print_csc_templated( size_t m, size_t n, int *Ap,
                                                int *Ai, type *Ax,
                                                std::streambuf* out = std::cout.rdbuf(),
                                                const std::string indent="  ",
                                                const std::string &beg="%%MatrixMarket matrix coordinate real symmetric"){
  std::streambuf* sb_cout_backup = std::cout.rdbuf();
  std::cout.rdbuf(out);
  std::cout<<indent<<beg << "\n";
  size_t nnz = n>0?Ap[n]:0;
  std::cout<<indent<<m<<" "<<n<<" "<<nnz<<"\n";
  for (auto i = 0; i < n; ++i) {
   for (auto j = Ap[i]; j < Ap[i+1]; ++j) {
    assert(j<nnz);
    std::cout<<indent<<Ai[j]+1<<" "<<i+1<<" "<<std::setprecision(print_precision)<<Ax[j];
    std::cout<<"\n";
   }
  }
  std::cout.rdbuf(sb_cout_backup);
 }

 /*
 * Print dense matrix, stored col-wise
 * The out can be file.rdbuf()
 */
 template<class type> void print_dense( int row_no,
                                        int col_no, int lda,
                                        type *mat,
                                        std::streambuf* out = std::cout.rdbuf(),
                                        const std::string indent="  ",
                                        const std::string &header="%%MatrixMarket matrix array real general"){
  std::streambuf* sb_cout_backup = std::cout.rdbuf();
  std::cout.rdbuf(out);
  std::cout<<indent<<header<<"\n";
  std::cout<<indent<<row_no<<" "<<col_no<<"\n";
  for (int i = 0; i < col_no*row_no; i+=lda) {
   for (int j = 0; j < lda; ++j) {
    std::cout<<indent<<std::setprecision(print_precision)<<mat[i+j];
   }
   std::cout<<"\n";
  }
  std::cout.rdbuf(sb_cout_backup);
 }

 /// prints a constant value with indention
 /// \tparam type
 /// \param val
 /// \param out
 /// \param indent
 template<class type> void print_constant( type val,
                                           std::streambuf* out,
                                           const std::string indent){
  std::streambuf* sb_cout_backup = std::cout.rdbuf();
  std::cout.rdbuf(out);
  std::cout<<indent<<std::setprecision(print_precision)<<val<<"\n";
  std::cout.rdbuf(sb_cout_backup);
 }



 /// print a string to output with an indent
 /// \param context
 /// \param out
 /// \param indent
 void print_string(std::string context, std::streambuf* out,
                   std::string indent);


/// Reads the header information of a mtx format.
/// \param inFile input stream
/// \param n_row
/// \param n_col
/// \param n_nnz
/// \param type
/// \param shape
/// \param mtx_format
/// \return
 void read_header(std::ifstream &inFile, int &n_row, int &n_col,
                  size_t &n_nnz, int &type, int &shape, int &mtx_format);


 int shape2int(int shape);


/// Builds triplet from coordinate mtx file
/// \param inFile
/// \param nnz
/// \param triplet_vec
/// \param zero_indexing
 void read_triplets_real(std::ifstream &inFile, int nnz,
                         std::vector<triplet>& triplet_vec,
                         bool read_val,
                         bool zero_indexing);



 void compress_triplets_to_csc(std::vector<triplet>& triplet_vec, CSC *A,
                               bool add_diags);



 void read_mtx_csc_real(std::ifstream &in_file, CSC *&A, bool insert_diag=false);


 /*
 * Reads an array stored in matrix market format.
 */
 void read_mtx_array_real(std::ifstream &in_file, Dense *&A);


 void load_mtx_array_real(std::string &filename, Dense *&A);


/// Reads a real constant from input
/// \param in_file
/// \param val
 void read_real_constant(std::ifstream &in_file, double &val);


 /// Reading sting from a file with indention
 /// \param in_file
 /// \param context
 /// \param indent
 void read_string(std::ifstream &in_file, std::string &context, std::string indent);

}


#endif //PROJECT_SPARSE_IO_H
