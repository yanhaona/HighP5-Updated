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

namespace mpi_cg {

//-------------------------------------------------------------- Data Dimension Information

// array dimension objects
Dimension valDims[1];
Dimension colDims[1];
Dimension rowDims[1];
Dimension bDims[1];
Dimension xDims[1];

//--------------------------------------------------------------- Data Structure References

// arrays for the sparse matrix 
double *values = NULL;
int *columns = NULL;
int *rows = NULL;

// arrays for the known and prediction vectors
double *b = NULL;
double *x_i = NULL;

// intermediate vectors needed by the algorithm
double *r_i = NULL;
double *a_r_i = NULL;

int maxIterations;

// descent progress tracker
double norm;
double denorm;
double alpha_i;

//-------------------------------------------------------------- Partition Config Variables

int processId;
int processCount;

//-------------------------------------------------------------------- Supporting functions


double *readDoubleArrayFromFile(const char *filePath, Dimension *dim) {

        std::ifstream file(filePath);
        if (!file.is_open()) {
                std::cout << "could not open file" << filePath << "\n";
                std::exit(EXIT_FAILURE);
        }
        readArrayDimensionInfoFromFile(file, 1, dim);
        file.close();
	
        double *array = new double[dim[0].length];
        TypedInputStream<double> *stream = new TypedInputStream<double>(filePath);
        int storeIndex = 0;
        stream->open();
        for (int i = 0; i < dim[0].length; i++) {
                array[storeIndex] = stream->readNextElement();
                storeIndex++;
        }
        stream->close();
        delete stream;

        return array;
}

int *readIntArrayFromFile(const char *filePath, Dimension *dim) {

        std::ifstream file(filePath);
        if (!file.is_open()) {
                std::cout << "could not open file" << filePath << "\n";
                std::exit(EXIT_FAILURE);
        }
        readArrayDimensionInfoFromFile(file, 1, dim);
        file.close();
        int *array = new int[dim[0].length];
        TypedInputStream<int> *stream = new TypedInputStream<int>(filePath);
        int storeIndex = 0;
        stream->open();
        for (int i = 0; i < dim[0].length; i++) {
                array[storeIndex] = stream->readNextElement();
                storeIndex++;
        }
        stream->close();
        delete stream;

        return array;
}


void writeResultToFile() {

        if (processId != 0) {
		return;
        }

        List<Dimension*> *dimLengths = new List<Dimension*>;
        dimLengths->Append(&xDims[0]);
        TypedOutputStream<double> *stream = new TypedOutputStream<double>("last-pred.bin", dimLengths, processId == 0);
	
        int storeIndex = 0;
        stream->open();
        for (int j = 0; j < xDims[0].length; j++) {
        	stream->writeNextElement(x_i[storeIndex]);
                storeIndex++;
        }

        stream->close();
        delete stream;
}

//--------------------------------------------------------------------- Conjugate Gradient Function

void performConjugateGradientIterations() {

	int rowsPerProcess = (xDims[0].length + processCount - 1) / processCount;
	int rowStart = rowsPerProcess * processId;
	int rowEnd = rowStart + rowsPerProcess - 1;
	if (rowEnd >= xDims[0].length) {
		rowEnd = xDims[0].length - 1;
	}
	int selfRowCount = rowEnd - rowStart + 1;

	// this two buffers are here to facilitate MPI allgather communication
	double *sendBuffer = new double[rowsPerProcess];
	double *receiveBuffer = new double[rowsPerProcess * processCount];

	int iteration = 0;
	do {
                //--------------------------------------------------------- r_i = b - A * x_i computation
                //
                for (int i = rowStart; i <= rowEnd; i++) r_i[i] = 0;
                
		for (int i = rowStart; i <= rowEnd; i++) {
                        int start = (i > 0) ? rows[i - 1] + 1 : 0;
                        int end = rows[i];
                        for (int j = start; j <= end; j++) {
                                r_i[i] = r_i[i] + values[j] * x_i[columns[j]];
                        }
                }

                for (int i = rowStart; i <= rowEnd; i++) {
                        r_i[i] = b[i] - r_i[i];
                }

                // do an all to all scatter-getter communication so that everyone get the updated full copy
		// of r_i and also know that a_r_i can be recomputed now
		double *sendStart = &r_i[rowStart];
		memcpy(sendBuffer, sendStart, selfRowCount * sizeof(double));
                MPI_Allgather(sendBuffer, rowsPerProcess, MPI_DOUBLE, receiveBuffer, rowsPerProcess, MPI_DOUBLE, MPI_COMM_WORLD);
		memcpy(r_i, receiveBuffer, xDims[0].length * sizeof(double));

		//----------------------------------------------------------- a_r_i = A * r_i computation
		//
	
		for (int i = rowStart; i <= rowEnd; i++) a_r_i[i] = 0;
                
		for (int i = rowStart; i <= rowEnd; i++) {
                        int start = (i > 0) ? rows[i - 1] + 1 : 0;
                        int end = rows[i];
                        for (int j = start; j <= end; j++) {
                                a_r_i[i] = a_r_i[i] + values[j] * r_i[columns[j]];
                        }
                }
                
		// do an all to all scatter-getter communication so that everyone get the updated full copy
		// of a_r_i
		sendStart = &a_r_i[rowStart];
		memcpy(sendBuffer, sendStart, selfRowCount * sizeof(double));
                MPI_Allgather(sendBuffer, rowsPerProcess, MPI_DOUBLE, receiveBuffer, rowsPerProcess, MPI_DOUBLE, MPI_COMM_WORLD);
		memcpy(a_r_i, receiveBuffer, xDims[0].length * sizeof(double));

		//------------------------------------- alpha_i = (r_i * r_i) / (r_i * a_r_i) computation
                //

		// do partial r_i * r_i and r_i * a_r_i computations in all processes
		double partNorm = 0.0;
                double partDenorm = 0.0;
                for (int i = rowStart; i <= rowEnd; i++) {
                        partNorm += r_i[i] * r_i[i];
                        partDenorm += r_i[i] * a_r_i[i];
                }

		// accumulate the partial dot products in all processes
		MPI_Allreduce(&partNorm, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&partDenorm, &denorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		alpha_i = norm / denorm;

		//------------------------------------------------ x_i = x_i + alpha_i * r_i computation
                //
                for (int i = rowStart; i <= rowEnd; i++) {
                        x_i[i] = x_i[i] + alpha_i * r_i[i];
                }

		// do an all to all scatter-getter communication so that everyone get the updated full copy
		// of updated x_i
		sendStart = &x_i[rowStart];
		memcpy(sendBuffer, sendStart, selfRowCount * sizeof(double));
                MPI_Allgather(sendBuffer, rowsPerProcess, MPI_DOUBLE, receiveBuffer, rowsPerProcess, MPI_DOUBLE, MPI_COMM_WORLD);
		memcpy(x_i, receiveBuffer, xDims[0].length * sizeof(double));

		iteration++;
        } while (iteration < maxIterations);
}

} // end of namespace

using namespace mpi_cg;


//--------------------------------------------------------------------------- Main Function

int mainMConjGrad(int argc, char *argv[]) {
	
	// do MPI intialization
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &processId);
	MPI_Comm_size(MPI_COMM_WORLD, &processCount);

	if (argc < 8) {
		if (processId == 0) {
			std::cout << "Pass arguments as follows\n";
			std::cout << "\t1. path of the binary input file for the value array of the sparse matrix\n";
			std::cout << "\t2. path of the binary input file for the column array of the sparse matrix\n";
			std::cout << "\t3. path of the binary input file for the row array of the sparse matrix\n";
			std::cout << "\t4. path of the binary input file for the known vector\n";
			std::cout << "\t5. path of the binary input file for the prediction vector\n";
			std::cout << "\t6. the number of conjugate gradient refinement iterations\n";
			std::cout << "\t7. a 0 to disable file writing or a 1 to write the final version of the plate to a file\n";
		}
		MPI_Finalize();
		exit(EXIT_FAILURE);		
	}

	if (processId == 0) {
		std::cout << "Running MPI Stencil experiment with " << processCount << " MPI processes\n";
	}

	// start timer
	struct timeval start;
        gettimeofday(&start, NULL);

	// parse command line arguments
	maxIterations = atoi(argv[6]);
	bool fileWriteMode (atoi(argv[7]) == 1);

	const char *filePathV= argv[1];
	const char *filePathC= argv[2];
	const char *filePathR= argv[3];
	const char *filePathK= argv[4];
	const char *filePathP= argv[5];
	std::ifstream fileV(filePathV);
        if (!fileV.is_open()) {
                std::cout << "could not open the file " << filePathV << "\n";
                std::exit(EXIT_FAILURE);
        }
	readArrayDimensionInfoFromFile(fileV, 1, valDims);
	fileV.close();
	std::ifstream fileC(filePathC);
        if (!fileC.is_open()) {
                std::cout << "could not open the file " << filePathC << "\n";
                std::exit(EXIT_FAILURE);
        }
	readArrayDimensionInfoFromFile(fileC, 1, colDims);
	fileC.close();
	std::ifstream fileR(filePathR);
        if (!fileR.is_open()) {
                std::cout << "could not open the file " << filePathR << "\n";
                std::exit(EXIT_FAILURE);
        }
	readArrayDimensionInfoFromFile(fileR, 1, rowDims);
	fileR.close();
	std::ifstream fileK(filePathK);
        if (!fileK.is_open()) {
                std::cout << "could not open the file " << filePathK << "\n";
                std::exit(EXIT_FAILURE);
        }
	readArrayDimensionInfoFromFile(fileK, 1, bDims);
	fileK.close();
	std::ifstream fileP(filePathP);
        if (!fileP.is_open()) {
                std::cout << "could not open the file " << filePathP << "\n";
                std::exit(EXIT_FAILURE);
        }
	readArrayDimensionInfoFromFile(fileP, 1, xDims);
	fileP.close();

	// read all arrays from file
        // reading elements of the CSR sparse matrix
        values = readDoubleArrayFromFile(filePathV, valDims);
        columns = readIntArrayFromFile(filePathC, colDims);
        rows = readIntArrayFromFile(filePathR, rowDims);
        // reading other vectors
        b = readDoubleArrayFromFile(filePathK, bDims);
        x_i = readDoubleArrayFromFile(filePathP, xDims);

        // create other variables used in the algorithm
        r_i = new double[xDims[0].length];
        a_r_i = new double[xDims[0].length];
	
	struct timeval memEnd;
        gettimeofday(&memEnd, NULL);

	//------------------------------------------------------------ Computation Starts 

	performConjugateGradientIterations();

	//-------------------------------------------------------------- Computation Ends


	MPI_Barrier(MPI_COMM_WORLD);
	
	// end timer
	struct timeval end;
        gettimeofday(&end, NULL);
	if (processId == 0) {
		double dataReadingTime = ((memEnd.tv_sec + memEnd.tv_usec / 1000000.0)
				- (start.tv_sec + start.tv_usec / 1000000.0));
		double executionTime = ((end.tv_sec + end.tv_usec / 1000000.0)
				- (start.tv_sec + start.tv_usec / 1000000.0));
		std::cout << "Memory initialization time: " << dataReadingTime << " Seconds\n";
		std::cout << "Execution time: " << executionTime << " Seconds\n";
		std::cout << "Matrix dimension: " << rowDims[0].length << " by " << rowDims[0].length << "\n";
		std::cout << "Refinement iterations: " << maxIterations << "\n";
	}
	
	// write outputs to files
	if (fileWriteMode) {
		if (processId == 0) {
			std::cout << "Writing the result in file plate-out.bin\n";
		}
		writeResultToFile();
	}

	// do MPI teardown
	MPI_Finalize();
	return 0;
}



