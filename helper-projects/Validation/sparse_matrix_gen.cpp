#include <cstdlib>
#include <cctype>
#include <stdlib.h>
#include <string.h>
#include <deque>
#include <sys/time.h>
#include <time.h>
#include <string>
#include <random>

#include "sparse_matrix_gen.h"
#include "fileUtility.h"
#include "structures.h"
#include "stream.h"
#include "list.h"

namespace sparse_matrix {

	const int INT_TYPE = 1;
	const int FLOAT_TYPE = 2;
	const int DOUBLE_TYPE = 3;
}

using namespace std;
using namespace sparse_matrix;

void saveDoubleArrayInBinaryFile(const char* fileName, Dimension dimension, double *data) {
	
	List<Dimension*> *dimensionList = new List<Dimension*>();
	dimensionList->Append(&dimension);
	TypedOutputStream<double> stream = TypedOutputStream<double>(fileName, dimensionList, true);
	stream.open();
	for (int i = 0; i < dimension.length; i++) {
		stream.writeNextElement(data[i]);
	}
	stream.close();
}

void saveFloatArrayInBinaryFile(const char* fileName, Dimension dimension, float *data) {

	List<Dimension*> *dimensionList = new List<Dimension*>();
	dimensionList->Append(&dimension);
	TypedOutputStream<float> stream = TypedOutputStream<float>(fileName, dimensionList, true);
	stream.open();
	for (int i = 0; i < dimension.length; i++) {
		stream.writeNextElement(data[i]);
	}
	stream.close();
}

void saveIntArrayInBinaryFile(const char* fileName, Dimension dimension, int *data) {

	List<Dimension*> *dimensionList = new List<Dimension*>();
	dimensionList->Append(&dimension);
	TypedOutputStream<int> stream = TypedOutputStream<int>(fileName, dimensionList, true);
	stream.open();
	for (int i = 0; i < dimension.length; i++) {
		stream.writeNextElement(data[i]);
	}
	stream.close();
}

void genCSRMatrix(int sparsity, int rowCount, int colCount, int dataType, bool binary) {


        srand(time(NULL));

	int nonZerosPerColumn = colCount * (100 - sparsity) / 100;
	int nonZeroEntries = nonZerosPerColumn;

	// ------------------------------------------------------------- generate and write the value array
	Dimension valueDim;
	valueDim.length = nonZeroEntries;
	valueDim.range.min = 0;
	valueDim.range.max = nonZeroEntries - 1;
	if (dataType == INT_TYPE) {
		int *valueArray = allocate<int>(1, &valueDim);
		for (int i = 0; i < nonZeroEntries; i++) {
			valueArray[i] = rand() % 100 + 1;	
		}
		if (binary == true) {
			cout << "saving the value array in file 'values'\n";
			saveIntArrayInBinaryFile("values", valueDim, valueArray);
		} else {
			writeArrayToFile<int>("values", valueArray, 1, &valueDim);
		}	
		delete valueArray;
	} else if (dataType == FLOAT_TYPE) {
		float *valueArray = allocate<float>(1, &valueDim);
		for (int i = 0; i < nonZeroEntries; i++) {
			valueArray[i] = (rand() % 100) * 1.0 / (rand() % 100 + 1);
		}
		if (binary == true) {
			cout << "saving the value array in file 'values'\n";
			saveFloatArrayInBinaryFile("values", valueDim, valueArray);
		} else {
			writeArrayToFile<float>("values", valueArray, 1, &valueDim);
		}	
		delete valueArray;
	} else if (dataType == DOUBLE_TYPE) {
		double *valueArray = allocate<double>(1, &valueDim);
		for (int i = 0; i < nonZeroEntries; i++) {
			valueArray[i] = (rand() % 100) * 1.0 / (rand() % 100 + 1);
		}	
		if (binary == true) {
			cout << "saving the value array in file 'values'\n";
			saveDoubleArrayInBinaryFile("values", valueDim, valueArray);
		} else {
			writeArrayToFile<double>("values", valueArray, 1, &valueDim);
		}	
		delete valueArray;
	}

	// --------------------------------------------------------------- Generate the Column Index Array
	int *colIndexArray = allocate<int>(1, &valueDim);

	// initialize a random number generator
	random_device seed;
	mt19937 gen{seed()};
	uniform_int_distribution<> dist{0, colCount - 1};

	int filled = 0;
	for (int row = 0; row < rowCount; row++) {

		// generate array of random indexes for each row
		for (int i = 0; i < nonZerosPerColumn && filled < nonZeroEntries; i++) {
			
			// generate a random guess
			int guess = dist(gen);

			// check if this number is already entered in the sorted array of column indices
			// for the current row
			int currRowStart = row * nonZerosPerColumn;
			bool duplicate = false;
			for (int j = currRowStart; j < filled; j++) {
				if (colIndexArray[j] == guess) {
					duplicate = true;
					break;
				}
			}

			// if the random guess is a duplicate then generate another index
			if (duplicate == true) {
				i--;
				continue;
			}


			// otherwise add the index in sorted order in the part of the array reserved for
			// the current row
			if (filled == currRowStart) {
				colIndexArray[currRowStart] = guess;
			} else {
				for (int j = filled - 1; j >= currRowStart; j--) {
					if (colIndexArray[j] > guess) {
						colIndexArray[j + 1] = colIndexArray[j];
						if (j == currRowStart) {
							colIndexArray[j] = guess;
						}	
					} else {
						colIndexArray[j + 1] = guess;
						break;
					}
				}
			}
			
			// increment the array filler tracker
			filled++;
		}
	}
	if (binary == true) {
		cout << "saving the column index array in file 'columns'\n";
		saveIntArrayInBinaryFile("columns", valueDim, colIndexArray);
	} else {
		writeArrayToFile<int>("colIndexes", colIndexArray, 1, &valueDim);
	}	
	delete colIndexArray;
	
	// ------------------------------------------------------------------ Generate the Row Range Array
	
	Dimension rangeDim;
	rangeDim.length = rowCount;
	rangeDim.range.min = 0;
	rangeDim.range.max = rowCount - 1;
	int *rowRangeArray = allocate<int>(1, &rangeDim);
	int rowEnd = nonZerosPerColumn;
	for (int row = 0; row < rowCount; row++) {
		rowRangeArray[row] = nonZerosPerColumn * (row + 1) - 1;
	}
	rowRangeArray[rowCount - 1] = nonZeroEntries - 1;
	if (binary == true) {
		cout << "saving the row range array in file 'rows'\n";
		saveIntArrayInBinaryFile("rows", rangeDim, rowRangeArray);
	} else {
		writeArrayToFile<int>("rowRanges", rowRangeArray, 1, &rangeDim);
	}	
	delete rowRangeArray;
}

int main(int argc, const char* argv[]) {

	if (argc < 6) {
                std::cout << "provide the following information\n";
                std::cout << "1. the number of rows in the matrix,\n";
                std::cout << "2. the number of columns in the matrix\n";
                std::cout << "3. the sparsity of matrix as a percentage from 0 to 100 (a higher percentage means more sparse)\n";
		std::cout << "4. data type of the matrix\n";
		std::cout << "\t1 = int\n";
		std::cout << "\t2 = float\n";
		std::cout << "\t3 = double\n";
		std::cout << "5. File type binary=1, text=0\n";
                std::exit(EXIT_FAILURE);
        }

	int rowCount = atoi(argv[1]);
	int colCount = atoi(argv[2]);
	int sparsity = atoi(argv[3]);
	int dataType = atoi(argv[4]);
	bool binaryFile = (atoi(argv[5]) == 1);

	genCSRMatrix(sparsity, rowCount, colCount, dataType, binaryFile);
	return 0;
}

