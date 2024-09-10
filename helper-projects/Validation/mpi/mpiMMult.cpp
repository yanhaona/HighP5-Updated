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
#include <pthread.h>
#include <mpi.h>

#include "../utils.h"
#include "../structures.h"
#include "../fileUtility.h"
#include "../stream.h"

namespace mpi_mmult {

//-------------------------------------------------------------- Data Dimension Information

Dimension aDims[2];
Dimension bDims[2];
Dimension cDims[2];

//--------------------------------------------------------------- Data Structure References

double *a;
double *b;
double *c;

//-------------------------------------------------------------- Partition Config Variables

int blockSize;
int processId;
int processCount;


//-------------------------------------------------------------------- Supporting functions

void allocateC() {

	int rowsPerProcess = (aDims[0].length + processCount - 1) / processCount;
	int localEntries = rowsPerProcess * bDims[1].length;
	c = new double[localEntries];
}

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
                for (int j = 0; j < aDims[1].length; j++) {
			indexList->clear();
			indexList->Append(i);
			indexList->Append(j);
			a[storeIndex] = stream->readElement(indexList);
			storeIndex++;
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

void writeCToFile() {

	// wait for your turn
        if (processId != 0) {
                int predecessorDone = 0;
                MPI_Status status;
                MPI_Recv(&predecessorDone, 1, MPI_INT, processId - 1, 0, MPI_COMM_WORLD, &status);
        }

        List<Dimension*> *dimLengths = new List<Dimension*>;
        dimLengths->Append(&cDims[0]);
        dimLengths->Append(&cDims[1]);
        TypedOutputStream<double> *stream = new TypedOutputStream<double>("c.bin", dimLengths, processId == 0);
	
	int rowsPerProcess = (aDims[0].length + processCount - 1) / processCount;
	int rowStart = processId * rowsPerProcess;
	int rowEnd = rowStart + rowsPerProcess - 1;
	if (rowEnd >= aDims[0].length) {
		rowEnd = aDims[0].length - 1;
	}

        List<int> *indexList = new List<int>;
        int storeIndex = 0;
        stream->open();
        for (int i = rowStart; i <= rowEnd; i++) {
        	for (int j = 0; j < cDims[1].length; j++) {
                	indexList->clear();
                        indexList->Append(i);
                        indexList->Append(j);
                        stream->writeElement(c[storeIndex], indexList);
                        storeIndex++;
                }
        }
        stream->close();

        delete indexList;
        delete stream;

        // notify the next in line
        if (processId < processCount - 1) {
                int writingDone = 1;
                MPI_Send(&writingDone, 1, MPI_INT, processId + 1, 0, MPI_COMM_WORLD);
        }
}

//--------------------------------------------------------------------- Block MMULT function

void multiplyMatrices() {
	
	// determine the starting and ending row for the current process
	int rowsPerProcess = (aDims[0].length + processCount - 1) / processCount;
	int rowStart = processId * rowsPerProcess;
	int rowEnd = rowStart + rowsPerProcess - 1;
	if (rowEnd >= aDims[0].length) {
		rowEnd = aDims[0].length - 1;
	}


	// run the block matrix-matrix multiplication algorithm for the rows allocated to the thread
        for (int iB = rowStart; iB <= rowEnd; iB += blockSize) {
                int rStart = iB;
                int rEnd = rStart + blockSize - 1;
                if (rEnd >= aDims[0].length) rEnd = aDims[0].length - 1;
                for (int jB = 0; jB < bDims[1].length; jB += blockSize) {
                        int cStart = jB;
                        int cEnd = cStart + blockSize - 1;
                        if (cEnd >= bDims[1].length) cEnd = bDims[1].length - 1;
                        for (int kB = 0; kB < aDims[1].length; kB += blockSize) {
                                int startIndex = kB;
                                int endIndex = startIndex + blockSize - 1;
                                if (endIndex >= aDims[1].length) endIndex = aDims[1].length - 1;

				// regular matrix matrix multiplication within a block
                                for (int i = rStart; i <= rEnd; i++) {
                                        int aRowIndex = (i - rowStart) * aDims[1].length;
                                        int cRowIndex = (i - rowStart) * cDims[1].length;
                                        for (int j = cStart; j <= cEnd; j++) {
                                                for (int k = startIndex; k <= endIndex; k++) {
                                                        int bRowIndex = k * bDims[1].length;
                                                        c[cRowIndex + j] += a[aRowIndex + k] * b[bRowIndex + j];
                                                }
                                        }
                                } // end of regular matix matrix multiplication
                        } // end of result row blocking
                } // end of column blocking
        } // end of row blocking
}

} // end of namespace

using namespace mpi_mmult;


//--------------------------------------------------------------------------- Main Function

int mainMMMult(int argc, char *argv[]) {
	
	// do MPI intialization
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &processId);
	MPI_Comm_size(MPI_COMM_WORLD, &processCount);

	if (argc < 5) {
		if (processId == 0) {
			std::cout << "Pass arguments as follows\n";
			std::cout << "\t1. block size for partitioning\n";
			std::cout << "\t2. path of the binary input file for matrix 1\n";
			std::cout << "\t3. path of the binary input file for matrix 2\n";
			std::cout << "\t4. a 0 to disable file writing or a 1 to write the result matrix to a file\n";
		}
		MPI_Finalize();
		exit(EXIT_FAILURE);		
	}

	if (processId == 0) {
		std::cout << "Running MPI Block MMult experiment with " << processCount << " MPI processes\n";
	}

	// start timer
	struct timeval start;
        gettimeofday(&start, NULL);

	// parse command line arguments
	blockSize = atoi(argv[1]);
	bool fileWriteMode (atoi(argv[4]) == 1);	
	const char *filePathA = argv[2];
	std::ifstream fileA(filePathA);
        if (!fileA.is_open()) {
                std::cout << "could not open the specified file\n";
                std::exit(EXIT_FAILURE);
        }
	readArrayDimensionInfoFromFile(fileA, 2, aDims);
	fileA.close();
	const char *filePathB = argv[3];
	std::ifstream fileB(filePathB);
        if (!fileB.is_open()) {
                std::cout << "could not open the specified file\n";
                std::exit(EXIT_FAILURE);
        }
	readArrayDimensionInfoFromFile(fileB, 2, bDims);
	fileB.close();

	// setup the result matrix
	cDims[0] = aDims[0];
	cDims[1] = bDims[1];
	allocateC();

	// read input matrices
	readAFromFile(filePathA);
	readBFromFile(filePathB);

	//------------------------------------------------------------ Computation Starts 

	multiplyMatrices();

	//-------------------------------------------------------------- Computation Ends


	MPI_Barrier(MPI_COMM_WORLD);
	
	// end timer
	struct timeval end;
        gettimeofday(&end, NULL);
	if (processId == 0) {
		double executionTime = ((end.tv_sec + end.tv_usec / 1000000.0)
				- (start.tv_sec + start.tv_usec / 1000000.0));
		std::cout << "Execution time: " << executionTime << " Seconds\n";
		std::cout << "Matrix dimension M1: " << aDims[0].length << " by " << aDims[1].length << "\n";
		std::cout << "Matrix dimension M2: " << bDims[0].length << " by " << bDims[1].length << "\n";
		std::cout << "Block size for partition: " << blockSize << "\n";
	}
	
	// write outputs to files
	if (fileWriteMode) {
		if (processId == 0) {
			std::cout << "Writing the result in file c.bin\n";
		}
		writeCToFile();
	}

	// do MPI teardown
	MPI_Finalize();
	return 0;
}



