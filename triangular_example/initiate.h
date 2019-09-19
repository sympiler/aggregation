#include "mmio.h"
#include "def.h"

/**
 * Process a matrix and vector from files in Matrix Market Format
 * and stores in Ap, Ai, Ax
 *
 * @param mf file pointer to matrix
 * @param Ap pointer to column pointer of matrix A
 * @param Ai pointer to row indices of matrix A
 * @param Ax pointer to values of matrix A
 * @param n pointer to column rank of matrix A
 * @param nz pointer to number non-zero of matrix A
 */
void process_matrix(FILE *mf, CSC *&A);

/**
 * Reads in a matrix from mtx file and convert
 * from COO format to CSC format
 *
 * @param f file pointer
 * @param n column rank
 * @param nz number of non-zeros
 * @param Ai row indices of A
 * @param Ap column pointer of A
 * @param Ax explicit values of A
 */
void load_matrix(FILE *f, int n, int nz, int *&Ai, int *&Ap, double *&Ax);


/**
 * Converts COO format to CSC format
 *
 * @param n column rank of matrix
 * @param nz number of non-zeros
 * @param J column indices
 * @param Lp resulting column pointer
 */
void to_csc(int n, int nz, int *J, int *Lp);


/**
 * Initialize RHS to produce a
 * resulting vector of Lx = b
 * of all 1's
 *
 * @param n column rank of A
 * @param Ap column pointer of A
 * @param Ai row indices of A
 * @param Ax explicit values of A
 * @param b RHS of Lx = b
 */
void rhsInit(int n, int *Lp, int *Li, double *Lx, double *b);


/**
 * Load a sparse vector from file
 *
 * @param f file pointer
 * @param n vector dimension
 * @param vz number of non-zero in vector
 * @param b vector to be loaded into
 * @param index indices of values
 */
void load_vector_sparse(FILE *f, int n, int vz, double *b, int **index);


/**
 * Load a dense vector from file
 *
 * @param f file pointer
 * @param n vector dimension
 * @param b vector to be loaded to
 */
void load_vector_dense(FILE *f, int n, double *b);


int transpose(CSC *L, CSC *U);


int allocateLC(BCSC *L, int sw);


int allocateAC(CSC *A, int nrow, int nnz, int sytpe, int sw);