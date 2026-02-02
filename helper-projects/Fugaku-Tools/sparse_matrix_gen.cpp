#include <cstdlib>
#include <cctype>
#include <stdlib.h>
#include <string.h>
#include <deque>
#include <sys/time.h>
#include <time.h>
#include <string>
#include <random>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

void genCSRMatrix(int sparsity, int rowCount, int colCount) {

        srand(time(NULL));

	int nonZerosPerColumn = colCount * (100 - sparsity) / 100;
	int nonZeroEntries = nonZerosPerColumn * rowCount;

	// ------------------------------------------------------------- generate and write the value array

	cout << "saving the value array in file 'values'\n";
	ofstream vstream("values", ios_base::binary);
	ostringstream str;
        str << nonZeroEntries;
        int strLength = str.str().length();
        vstream.write(str.str().c_str(), sizeof(char) * strLength);
        vstream.flush();
	char lineEnd = '\n';
        vstream.write(&lineEnd, sizeof(char));
        vstream.flush();
	for (int i = 0; i < nonZeroEntries; i++) {
		double item = (rand() % 100 * 1.0) / (rand() % 100 + 1.0);
                vstream.write(reinterpret_cast<char*>(&item), sizeof(item));
	}	
	vstream.close();

	// --------------------------------------------------------------- Generate the Column Index Array
	int *colIndexArray = new int[nonZeroEntries];

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
	cout << "saving the column index array in file 'columns'\n";
	ofstream cstream("columns", ios_base::binary);
        cstream.write(str.str().c_str(), sizeof(char) * strLength);
        cstream.flush();
        cstream.write(&lineEnd, sizeof(char));
        cstream.flush();

	for (int i = 0; i < nonZeroEntries; i++) {
		int item = colIndexArray[i];
                cstream.write(reinterpret_cast<char*>(&item), sizeof(item));
	}	
	delete[] colIndexArray;
	cstream.close();
	
	// ------------------------------------------------------------------ Generate the Row Range Array
	
	int *rowRangeArray = new int[rowCount];
	for (int row = 0; row < rowCount; row++) {
		rowRangeArray[row] = nonZerosPerColumn * (row + 1) - 1;
	}
	rowRangeArray[rowCount - 1] = nonZeroEntries - 1;
	cout << "saving the row range array in file 'rows'\n";
	ofstream rstream("rows", ios_base::binary);
	ostringstream str2;
        str2 << rowCount;
        strLength = str2.str().length();
        rstream.write(str2.str().c_str(), sizeof(char) * strLength);
        rstream.flush();
        rstream.write(&lineEnd, sizeof(char));
        rstream.flush();
	for (int i = 0; i < rowCount; i++) {
		int item = rowRangeArray[i];
                rstream.write(reinterpret_cast<char*>(&item), sizeof(item));
	}	
	delete[] rowRangeArray;
	rstream.close();
}

int main(int argc, const char* argv[]) {

	if (argc < 4) {
                std::cout << "provide the following information\n";
                std::cout << "1. the number of rows in the matrix,\n";
                std::cout << "2. the number of columns in the matrix\n";
                std::cout << "3. the sparsity of matrix as a percentage from 0 to 100 (a higher percentage means more sparse)\n";
                std::exit(EXIT_FAILURE);
        }

	int rowCount = atoi(argv[1]);
	int colCount = atoi(argv[2]);
	int sparsity = atoi(argv[3]);

	genCSRMatrix(sparsity, rowCount, colCount);
	return 0;
}

