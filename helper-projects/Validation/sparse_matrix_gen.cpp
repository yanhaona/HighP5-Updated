#include <cstdlib>
#include <cctype>
#include <stdlib.h>
#include <string.h>
#include <deque>
#include <sys/time.h>
#include <time.h>
#include <string>

#include "sparse_matrix_gen.h"
#include "fileUtility.h"
#include "structures.h"

namespace sparse_matrix {

	const int INT_TYPE = 1;
	const int FLOAT_TYPE = 2;
	const int DOUBLE_TYPE = 3;
}

using namespace std;
using namespace sparse_matrix;

void genCSRMatrix(int sparsity, int rowCount, int colCount, int dataType, const char* dirPath) {


        srand(time(NULL));

	int matrixSize = rowCount * colCount;
	int nonZeroEntries = matrixSize * sparsity / 100;

	// determine the file paths of the three output files
	string directory = "";
	if (dirPath != NULL) {
		directory.assign(dirPath);
	}	
	string valueFile = directory + "/values.txt";
	string colIndexFile = directory + "/colIndexes.txt";
	string rowRangeFile = directory + "/rowRanges.txt";

	// generate and write the value array
	Dimension valueDim;
	valueDim.length = nonZeroEntries;
	valueDim.range.min = 0;
	valueDim.range.max = nonZeroEntries - 1;
	if (dataType == INT_TYPE) {
		int *valueArray = allocate<int>(1, &valueDim);
		for (int i = 0; i < nonZeroEntries; i++) {
			valueArray[i] = rand() % 100;
		}	
		writeArrayToFile<int>(valueFile.c_str(), valueArray, 1, &valueDim);	
	} else if (dataType == FLOAT_TYPE) {
		float *valueArray = allocate<float>(1, &valueDim);
		for (int i = 0; i < nonZeroEntries; i++) {
			valueArray[i] = (rand() % 100) * 1.0 / (rand() % 100);
		}
		writeArrayToFile<float>(valueFile.c_str(), valueArray, 1, &valueDim);	
	} else if (dataType == DOUBLE_TYPE) {
		double *valueArray = allocate<double>(1, &valueDim);
		for (int i = 0; i < nonZeroEntries; i++) {
			valueArray[i] = (rand() % 100) * 1.0 / (rand() % 100);
		}	
		writeArrayToFile<double>(valueFile.c_str(), valueArray, 1, &valueDim);	
	}
}

