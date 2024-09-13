#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cctype>
#include <stdlib.h>
#include <string.h>
#include <deque>
#include "../utils.h"
#include "../structures.h"
#include "../fileUtility.h"
#include <sys/time.h>
#include <time.h>
#include <array>
#include "../cpuCores.h"
#include "../stream.h"

namespace threaded_cg {

// ------------------------------------------------------------------------------------ Shared Input Data

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

// configuration parameters
int sparsity;
int dimLength;
int maxIterations;

// array dimension objects
Dimension valDims[1];
Dimension colDims[1];
Dimension rowDims[1];
Dimension bDims[1];
Dimension xDims[1];

// descent progress tracker
double norm;
double alpha_i;

// ------------------------------------------------------------------------------- Thread Interaction Data

int threadCount;
pthread_barrier_t barrier;

// a variable for participating in dot product
double *partialSums1;
double *partialSums2;

// -------------------------------------------------------------------------------------- Helper Functions

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


// --------------------------------------------------------------------------------------- Thread Function

void *computeConjugateGradient(void *arg) {

        int threadId = *((int*) arg);
	
	// start the iteration counter;
	int iteration = 0;

	do {
		//--------------------------------------------------------- r_i = b - A * x_i computation
		//
		for (int i = threadId; i < xDims[0].length; i = i + threadCount) r_i[i] = 0;
		
		for (int i = threadId; i < rowDims[0].length; i = i + threadCount) {
			int start = (i > 0) ? rows[i - 1] + 1 : 0;
			int end = rows[i];
			for (int j = start; j <= end; j++) {
				r_i[i] = r_i[i] + values[j] * x_i[columns[j]];
			}
		}
		
		for (int i = threadId; i < bDims[0].length; i = i + threadCount) {
			r_i[i] = b[i] - r_i[i];
		}

		// synchronize the threads to ensure r_i computation is over in all threads and thus
		// a_r_i can be recomputed now
		pthread_barrier_wait(&barrier);

		//----------------------------------------------------------- a_r_i = A * r_i computation
		//
		for (int i = threadId; i < xDims[0].length; i = i + threadCount) a_r_i[i] = 0;
		for (int i = threadId; i < rowDims[0].length; i = i + threadCount) {
			int start = (i > 0) ? rows[i - 1] + 1 : 0;
			int end = rows[i];
			for (int j = start; j <= end; j++) {
				a_r_i[i] = a_r_i[i] + values[j] * r_i[columns[j]];
			}
		}

		// barrier synchronize to ensure all threads completed computing a_r_i
		pthread_barrier_wait(&barrier);

		//------------------------------------- alpha_i = (r_i * r_i) / (r_i * a_r_i) computation
		//
		double partNorm = 0.0;
		double partDenom = 0.0;
		for (int i = threadId; i < xDims[0].length; i = i + threadCount) {
			partNorm += r_i[i] * r_i[i];
			partDenom += r_i[i] * a_r_i[i];
		}
		partialSums1[threadId] = partNorm;
		partialSums2[threadId] = partDenom;
		
		// barrier synchronize to ensure all threads completed computing partial dot products
		pthread_barrier_wait(&barrier);

		// let the first thread accumulate the partial sums
		if (threadId == 0) {

			norm = 0;
			double denominator = 0;
			for (int i = 0; i < threadCount; i++) {
				norm += partialSums1[i];
				denominator += partialSums2[i];
			}	
			alpha_i = norm / denominator;
		}

		// barrier synchronize again to ensure that the computation of alpha_i is available for
		// other threads
		pthread_barrier_wait(&barrier);

		//------------------------------------------------ x_i = x_i + alpha_i * r_i computation
		//
		for (int i = threadId; i < xDims[0].length; i = i + threadCount) {
			x_i[i] = x_i[i] + alpha_i * r_i[i];
		}

		// barrier synchronize again to ensure that the new x_i is available for all  threads
		pthread_barrier_wait(&barrier);

		iteration++;
	} while (iteration < maxIterations);
	
	// exit pthread
        pthread_exit(NULL);
}


} // end of the namespace


// ----------------------------------------------------------------------------------------- Main Function

using namespace threaded_cg;
using namespace std;

int mainTConjGrad(int argc, const char *argv[]) {

	struct timeval start;
        gettimeofday(&start, NULL);

         if (argc < 9) {
                cout << "Provide the following command line arguments \n";
                cout << "\t1. file value array of CSR matrix\n";
                cout << "\t2. file for column array of CSR matrix\n";
                cout << "\t3. file for row array of CSR matrix\n";
                cout << "\t4. file for the known vector array\n";
                cout << "\t5. file for the pred vector array\n";
                cout << "\t6. max gradient descent iterations\n";
                cout << "\t7. The number of threads to use\n";
                cout << "\t8. Provide the machine name tozammel|brac\n";
                exit(EXIT_FAILURE);
        }

        // parse command line arguments
        maxIterations = atoi(argv[6]);
        threadCount = atoi(argv[7]);
        const char *machine = argv[8];

	// read all arrays from file
	// reading elements of the CSR sparse matrix
	values = readDoubleArrayFromFile(argv[1], valDims);
	columns = readIntArrayFromFile(argv[2], colDims);
	rows = readIntArrayFromFile(argv[3], rowDims);
	// reading other vectors
	b = readDoubleArrayFromFile(argv[4], bDims);
	x_i = readDoubleArrayFromFile(argv[5], xDims);

	// create other variables used in the algorithm
	r_i = new double[xDims[0].length];
	a_r_i = new double[xDims[0].length];

	// allocate variables for storing partial dot products
	partialSums1 = new double[threadCount];
	partialSums2 = new double[threadCount];

	// -----------------------------------------------------------------------------parallel programming starts
        
	struct timeval compute;
	gettimeofday(&compute, NULL);

        // initialize thread barrier
        pthread_barrier_init(&barrier, NULL, threadCount);

        // retrieve the processor core numbering for the machine and thread count configuration
        int* coreOrder = getCoreNumberingArray(machine, threadCount);
        int coreCount = getCoreCount(machine, threadCount);
        int coreJump = coreCount / threadCount;
        std::cout << "executing for " << machine << " machine with " << coreCount << " core model" << std::endl;

        // prepare attribute for thread pinning
        pthread_attr_t attr;
        cpu_set_t cpus;
        pthread_attr_init(&attr);

        // start the threads
        int threadIds[threadCount];
        pthread_t threads[threadCount];
        for (int i = 0; i < threadCount; i++) {
                threadIds[i] = i;

                // get the CPU core's physical ID for pinning
                int cpuId = i * coreJump;
                int physicalId = coreOrder[cpuId];
                CPU_ZERO(&cpus);
                CPU_SET(physicalId, &cpus);
                pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpus);

                // run thread
                int status = pthread_create(&threads[i], &attr, computeConjugateGradient, (void*) &threadIds[i]);
                if (status != 0) {
                        cout << "Could not create some pthreads\n";
                        exit(EXIT_FAILURE);
                }
        }

	// -------------------------------------------------------------------------------- parallel programming ends

        // join threads
        for (int i = 0; i < threadCount; i++) {
                pthread_join(threads[i], NULL);
        }
        cout << "number of threads: " << threadCount << "\n";

	struct timeval end;
        gettimeofday(&end, NULL);
        double executionTime = ((end.tv_sec + end.tv_usec / 1000000.0)
                        - (start.tv_sec + start.tv_usec / 1000000.0));
        double computationTime = ((end.tv_sec + end.tv_usec / 1000000.0)
                        - (compute.tv_sec + compute.tv_usec / 1000000.0));
        cout << "Execution Time: " << executionTime << " Seconds\n";
        cout << "Computation Time: " << computationTime << " Seconds\n";

	return 0;
}

