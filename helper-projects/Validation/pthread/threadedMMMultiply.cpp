#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cctype>
#include <stdlib.h>
#include <string.h>
#include <deque>
#include <sys/time.h>
#include <time.h>
#include "../utils.h"
#include "../structures.h"
#include "../fileUtility.h"


// ---------------------------------------------------------------------------------- Matrix Information
double *m1;
double *m2;
double *c;
Dimension a_Dims[2];
Dimension c_Dims[2];
Dimension b_Dims[2];

//------------------------------------------------------------------------------ Thread Interaction Data
int threadCountMMM;
int blockSizeMMM;

//-------------------------------------------------------------------------------------- Thread Function

void *computeBMMM(void *arg) {

        int threadId = *((int*) arg);

	// different threads get different chunks of rows from the result matrix to process
	int totalRows = c_Dims[0].length;
	int rowsPerThread = (totalRows + threadCountMMM - 1) / threadCountMMM;
	int rowStart = rowsPerThread * threadId;
	int rowEnd = rowStart + rowsPerThread - 1;
	if (rowEnd > totalRows - 1) {
		rowEnd = totalRows - 1;
	}
	
	// run the block matrix-matrix multiplication algorithm for the rows allocated to the thread
	for (int iB = rowStart; iB <= rowEnd; iB += blockSizeMMM) {
		int rStart = iB;
		int rEnd = rStart + blockSizeMMM - 1;
		if (rEnd >= a_Dims[0].length) rEnd = a_Dims[0].length - 1;
		for (int jB = 0; jB < b_Dims[1].length; jB += blockSizeMMM) {
			int cStart = jB;
			int cEnd = cStart + blockSizeMMM - 1;
			if (cEnd >= b_Dims[1].length) cEnd = b_Dims[1].length - 1;
			for (int kB = 0; kB < a_Dims[1].length; kB += blockSizeMMM) {
				int startIndex = kB;
				int endIndex = startIndex + blockSizeMMM - 1;
				if (endIndex >= a_Dims[1].length) endIndex = a_Dims[1].length - 1;
				for (int i = rStart; i <= rEnd; i++) {
					int aRowIndex = i * a_Dims[1].length;
					int cRowIndex = i * c_Dims[1].length;
					for (int j = cStart; j <= cEnd; j++) {
						for (int k = startIndex; k <= endIndex; k++) {
							int bRowIndex = k * b_Dims[1].length;
							c[cRowIndex + j] += m1[aRowIndex + k] * m2[bRowIndex + j];
						}
					}
				}
			}
		}
	}

        // exit thread
        pthread_exit(NULL);
}

//----------------------------------------------------------------------------------------------- main function
int mainTMMM(int argc, const char *argv[]) {

	if (argc < 5) {
                std::cout << "provide input file 1, input file 2, and blocking size\n";
                std::cout << "then specify the number of threads to be used as the last argument\n";
                std::exit(EXIT_FAILURE);
        }

        // interpret command line parameters
        const char* matrix1File = argv[1];
        const char* matrix2File= argv[2];
        blockSizeMMM = atoi(argv[3]);
        threadCountMMM = atoi(argv[4]);


	// load original inputs and output from generated files
	m1 = readArrayFromFile <double> (matrix1File, 2, a_Dims);
	m2 = readArrayFromFile <double> (matrix2File, 2, b_Dims);

	// starting execution timer clock
	struct timeval start;
	gettimeofday(&start, NULL);

	// declare and initialize c for current computation
	c_Dims[0] = a_Dims[0]; c_Dims[1] = b_Dims[1];
	int cSize = a_Dims[0].length * b_Dims[1].length;
	c = new double[cSize];
	for (int i = 0; i < cSize; i++) c[i] = 0;

	// start the threads
        int threadIds[threadCountMMM];
        pthread_t threads[threadCountMMM];
        for (int i = 0; i < threadCountMMM; i++) {
                threadIds[i] = i;
                int status = pthread_create(&threads[i], NULL, computeBMMM, (void*) &threadIds[i]);
                if (status != 0) {
                        std::cout << "Could not create some pthreads\n";
                        std::exit(EXIT_FAILURE);
                }
        }

        // join threads
        for (int i = 0; i < threadCountMMM; i++) {
                pthread_join(threads[i], NULL);
        }


	//-------------------------------- calculate running time
	struct timeval end;
	gettimeofday(&end, NULL);
	double runningTime = ((end.tv_sec + end.tv_usec / 1000000.0)
			- (start.tv_sec + start.tv_usec / 1000000.0));
	std::cout << "PThreaded Execution Time: " << runningTime << " Seconds" << std::endl;

	return 0;
}



