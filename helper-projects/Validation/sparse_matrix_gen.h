/*
 * sparse_matrix_gen.h
 *
 *  Created on: Aug, 2024
 *      Author: yan
 */

#ifndef SPARSE_MAT_GEN_H_
#define SPARSE_MAT_GEN_H_

/**
 * Generate a sparse matrix in compressed row format which has sparsity percentage of values in each 
 * row filled with zeros. The data type of the sparse matrix can be 1) int, 2) float, 3) double. CSR
 * format has three different linear arrays: value array, column index array, and the row range array.
 * The final boolean argument specifies whether the generated files for the sparse matrix should be
 * saved as binary files or text files.
 * */
void genCSRMatrix(int sparsity, int rowCount, int colCount, int dataType, bool binary);

#endif /* SPARSE_MAT_GEN_H_ */
