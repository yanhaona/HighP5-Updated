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

namespace threaded_cg {

// ------------------------------------------------------------------------------------ Shared Input Data

// arrays for the sparse matrix	
double *values;
int *columns;
int *rows;

// arrays for the known and prediction vectors
double *b;
double *x_i;

// intermediate vectors needed by the algorithm
double *r_i;
double *a_r_i;

// configuration parameters
int sparsity;
int dimLength;
int maxIterations;
double precision;

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

// --------------------------------------------------------------------------------------- Thread Function

void *computeConjugateGradient(void *arg) {

        int threadId = *((int*) arg);

	// start the iteration counter;
	int iteration = 0;

	do {
		// calculating r_i = b - A * x_i using strided partition of rows among threads
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

		// calculating alpha_i = (r_i * r_i) / (r_i * (A * r_i)) again using strided partition
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

		// let the first thread compute the vector dot product and the norm of update
		if (threadId == 0) {
			norm = 0;
			double denominator = 0;
			for (int i = 0; i < xDims[0].length; i++) {
				norm += r_i[i] * r_i[i];
				denominator += r_i[i] * a_r_i[i];
			}
			alpha_i = norm / denominator;
		}

		// barrier synchronize again to ensure that the computation of alpha_i is available for
		// other threads
		pthread_barrier_wait(&barrier);

		// calculating x_i = x_i + alpha_i * r_i using strided partitioning of columns among threads
		for (int i = threadId; i < xDims[0].length; i = i + threadCount) {
			x_i[i] = x_i[i] + alpha_i * r_i[i];
		}

		iteration++;
	} while (iteration < maxIterations && norm > precision);

	// exit pthread
        pthread_exit(NULL);
}


} // end of the namespace


using namespace threaded_cg;
using namespace std;

int mainTConjGrad(int argc, const char *argv[]) {

	struct timeval start;
        gettimeofday(&start, NULL);

         if (argc < 5) {
                cout << "Provide the following command line arguments \n";
                cout << "\t1. max gradient descent iterations\n";
                cout << "\t2. the fractional precision to stop iterations\n";
                cout << "\t3. The number of threads to use\n";
                cout << "\t4. Provide the machine name tozammel|brac\n";
                exit(EXIT_FAILURE);
        }

        // parse command line arguments
        maxIterations = atoi(argv[1]);
	precision = atof(argv[2]);
        threadCount = atoi(argv[3]);
        const char *machine = argv[4];

	// read all arrays from file
	// reading elements of the CSR sparse matrix
	columns = readArrayFromFile <int> ("columns", 1, colDims);
	rows = readArrayFromFile <int> ("rows", 1, rowDims);
	values = readArrayFromFile <double> ("values", 1, valDims);
	// reading other vectors
	b = readArrayFromFile <double> ("b", 1, bDims);
	x_i = readArrayFromFile <double> ("x_0", 1, xDims);

	// declare iteration limiter variables
	const int maxIterations = 10;
	const double precision = 1;

	// declare other variables used in the algorithm
	r_i = new double[xDims[0].length];
	a_r_i = new double[xDims[0].length];

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
        cout << "execution time: " << executionTime << " Seconds\n";
        cout << "computation time: " << computationTime << " Seconds\n";

	return 0;
}

