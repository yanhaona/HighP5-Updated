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
#include <array>
#include <mpi.h>
#include <pthread.h>

#include "../utils.h"
#include "../structures.h"
#include "../fileUtility.h"
#include "../cpuCores.h"
#include "../stream.h"


// ------------------------------------------------------------------------------------- MPI Information
int processId;
int processCount;

// ---------------------------------------------------------------------------------- Matrix Information
double *a;
double *b;
double *c;
Dimension aDims[2];
Dimension cDims[2];
Dimension bDims[2];

//------------------------------------------------------------------------------ Thread Interaction Data
int threadCountMMM;
int blockSizeMMM;

//------------------------------------------------------------------------------ Input Files Reader Code
void readAFromFile(const char *filePath) {

	int rowsPerProcess = (aDims[0].length + processCount - 1) / processCount;
	int rowStart = processId * rowsPerProcess;
	int rowEnd = rowStart + rowsPerProcess - 1;
	if (rowEnd >= aDims[0].length) {
		rowEnd = aDims[0].length - 1;
	}
	int rowCount = rowEnd - rowStart + 1;
	int localEntries = rowCount * aDims[1].length;
	a = new double[localEntries];

	TypedInputStream<double> *stream = new TypedInputStream<double>(filePath);
	int storeIndex = 0;
	stream->open();
	List<int> *indexList = new List<int>();
        for (int i = rowStart; i <= rowEnd; i++) {
		indexList->clear();
		indexList->Append(i);
		indexList->Append(0);
		a[storeIndex] = stream->readElement(indexList);
		storeIndex++;
		int count = 1;
                while (count < aDims[1].length) {
			a[storeIndex] = stream->readNextElement();
			storeIndex++;
			count++;
		}
	}
	stream->close();

	delete indexList;
	delete stream;
}

void readBFromFile(const char *filePath) {

	int bSize = bDims[0].length * bDims[1].length;
	b = new double[bSize];
	TypedInputStream<double> *stream = new TypedInputStream<double>(filePath);
	int storeIndex = 0;
	stream->open();
	for (int i = 0; i < bSize; i++) {
		b[storeIndex] = stream->readNextElement();
	}
	stream->close();
	delete stream;
}


//-------------------------------------------------------------------------------------- Thread Function

void *computeBMMM(void *arg) {

        int threadId = *((int*) arg);

	// different threads get different chunks of rows from the result matrix to process
	int totalRows = cDims[0].length;
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
		if (rEnd >= aDims[0].length) rEnd = aDims[0].length - 1;
		for (int jB = 0; jB < bDims[1].length; jB += blockSizeMMM) {
			int cStart = jB;
			int cEnd = cStart + blockSizeMMM - 1;
			if (cEnd >= bDims[1].length) cEnd = bDims[1].length - 1;
			for (int kB = 0; kB < aDims[1].length; kB += blockSizeMMM) {
				int startIndex = kB;
				int endIndex = startIndex + blockSizeMMM - 1;
				if (endIndex >= aDims[1].length) endIndex = aDims[1].length - 1;
				for (int i = rStart; i <= rEnd; i++) {
					int aRowIndex = i * aDims[1].length;
					int cRowIndex = i * cDims[1].length;
					for (int j = cStart; j <= cEnd; j++) {
						for (int k = startIndex; k <= endIndex; k++) {
							int bRowIndex = k * bDims[1].length;
							c[cRowIndex + j] += a[aRowIndex + k] * b[bRowIndex + j];
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
int main(int argc, char *argv[]) {

        MPI_Init(&argc, &argv);
	
	// this code only works for a single MPI process. So checking if that restriction is met
        MPI_Comm_rank(MPI_COMM_WORLD, &processId);
        MPI_Comm_size(MPI_COMM_WORLD, &processCount);
	if (processCount > 1) {
                std::cout << "you cannot have more than one MPI process for this program. \n";
		MPI_Finalize();
                std::exit(EXIT_FAILURE);
	}

	if (argc < 5) {
                std::cout << "provide input file 1, input file 2, and blocking size\n";
                std::cout << "then specify the number of threads to be used \n";
                std::exit(EXIT_FAILURE);
        }
	
	// starting execution timer clock
	struct timeval start;
	gettimeofday(&start, NULL);

	// read input matrices
	const char *filePathA = argv[1];
	std::ifstream fileA(filePathA);
        if (!fileA.is_open()) {
                std::cout << "could not open the specified file\n";
                std::exit(EXIT_FAILURE);
        }
	readArrayDimensionInfoFromFile(fileA, 2, aDims);
	fileA.close();
	const char *filePathB = argv[2];
	std::ifstream fileB(filePathB);
        if (!fileB.is_open()) {
                std::cout << "could not open the specified file\n";
                std::exit(EXIT_FAILURE);
        }
	readArrayDimensionInfoFromFile(fileB, 2, bDims);
	fileB.close();

	readAFromFile(filePathA);
	readBFromFile(filePathB);

	struct timeval memEnd;
        gettimeofday(&memEnd, NULL);

        // interpret command line parameters
        blockSizeMMM = atoi(argv[3]);
        threadCountMMM = atoi(argv[4]);


	// declare and initialize c for current computation
	cDims[0] = aDims[0]; cDims[1] = bDims[1];
	int cSize = aDims[0].length * bDims[1].length;
	c = new double[cSize];
	for (int i = 0; i < cSize; i++) c[i] = 0;

	std::cout << "executing for  " << threadCountMMM << " number of threads" << std::endl;


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
	std::cout << "MPI + PThread Execution Time: " << runningTime << " Seconds" << std::endl;

	MPI_Finalize();
	return 0;
}
