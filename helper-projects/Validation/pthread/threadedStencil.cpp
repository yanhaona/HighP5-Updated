#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cctype>
#include <string.h>
#include <deque>
#include <sys/time.h>
#include <time.h>

#include "../utils.h"
#include "../structures.h"
#include "../fileUtility.h"
#include "../cpuCores.h"

using namespace std;

namespace pthread_stencil {

// ------------------------------------------------------------------------------------ Shared Input Data
double *plate;
int plateDimLength;
double **plate0Parts;
double **plate1Parts;
int iterations;
int padding;

// ------------------------------------------------------------------------------- Thread Interaction Data

int threadCount;
pthread_barrier_t barrier;

// --------------------------------------------------------------------------------------- Thread Function

void *computeStencil(void *arg) {

        int threadId = *((int*) arg);

	// determine the row index range of plate belonging to the current thread;
	int rowsPerThread = (plateDimLength + threadCount - 1) / threadCount;
	int rowStart = rowsPerThread * threadId;
	int rowEnd = rowStart + rowsPerThread - 1;
	if (rowEnd >= plateDimLength) {
		rowEnd = plateDimLength - 1;
	}	
	
	// determine if there is upper and lower row padding to be synchronized
	int padRowStart = rowStart;
	int padRowEnd = rowEnd;
	if (threadId > 0) {
		padRowStart = rowStart - padding;
	}
	if (threadId < threadCount - 1) {
		padRowEnd = rowEnd + padding;
	}
	int rowCount = padRowEnd - padRowStart + 1;
	int partSize = rowCount * plateDimLength;

	// initialize plate part belonging to the current thread
	plate0Parts[threadId] = new double[partSize];
	plate1Parts[threadId] = new double[partSize];

	// copy data from the original plate variable to the two parts of the thread
	for (int i = padRowStart; i <= padRowEnd; i++) {
		for (int j = 0; j < plateDimLength; j++) {
			int plateIndex = i * plateDimLength + j;
			double data = plate[plateIndex];
			int storeIndex = (i - padRowStart) * plateDimLength + j;
			plate0Parts[threadId][storeIndex] = data;
			plate1Parts[threadId][storeIndex] = data;
		}
	}

	// synchronize all threads after parts loading
	pthread_barrier_wait(&barrier);
	
	// run Jacobi iterations
	int totalIteration = 0;
	for (int i = 0; i < iterations; i++) {

		// setup the plate part pointers for the current and next iterations
		double *input, *output;
		if (totalIteration % 2 == 0) {
			input = plate0Parts[threadId];
			output = plate1Parts[threadId];
		} else {
			input = plate1Parts[threadId];
			output = plate0Parts[threadId];
		}

		// perform padding number of internal iterations before a synchronization attempt
		for (int p = 0; p < padding; p++) {
				
			// update the output cells based on Jacobi iteration logic
			for (int y = padRowStart + 1; y <  padRowEnd - 1; y++) {

				int yIndex0 = (y - padRowStart) * plateDimLength;
				int yIndex1 = yIndex0 - plateDimLength;
				int yIndex2 = yIndex0 + plateDimLength;

				for (int x = 1; x < plateDimLength - 1; x++) {

					int index0 = yIndex0 + x;
					int index1 = yIndex0 + (x - 1);
					int index2 = yIndex0 + (x + 1);
					int index3 = yIndex1 + x;
					int index4 = yIndex2 + x;
					output[index0] = 0.25 * (input[index1] +
							input[index2] +
							input[index3] +
							input[index4]);
				}
			}

			// swap the input and output parts
			double *temp = input;
			input = output;
			output = input;

			// update total iteration counter
			totalIteration++;
		}

		// ensure that all threads have completed padding number of internal iterations
		pthread_barrier_wait(&barrier);

		// copy from current thread's part to previous thread's padding region
		if (threadId > 0) {
			double *upperThreadsPart = plate0Parts[threadId - 1];
			if (totalIteration % 2 == 0) {
				upperThreadsPart = plate1Parts[threadId - 1];
			}
			for (int x = 1; x < plateDimLength - 1; x++) {
				for (int p = 1; p <= padding; p++) {
					int myBoundaryY = padding - 1 + p;
					int receiverPaddingY = rowsPerThread - 1 + p;
					if (threadId - 1 > 0) {
						receiverPaddingY += padding;
					}
					int cellIndex = myBoundaryY * plateDimLength + x;
					double data = output[cellIndex];
					int storeIndex = receiverPaddingY * plateDimLength + x;
					upperThreadsPart[storeIndex] = data;
				}
			}
		}

		// copy from current thread's part to next thread's padding region
		if (threadId < threadCount - 1) {
			double *lowerThreadsPart = plate0Parts[threadId + 1];
			if (totalIteration % 2 == 0) {
				lowerThreadsPart = plate1Parts[threadId + 1];
			}
			for (int x = 1; x < plateDimLength - 1; x++) {
				for (int p = 1; p <= padding; p++) {
					int myBoundaryY = padRowEnd - padRowStart + 1 - (p + padding);
					int receiverPaddingY = padding - p;
					int cellIndex = myBoundaryY * plateDimLength + x;
					double data = output[cellIndex];
					int storeIndex = receiverPaddingY * plateDimLength + x;
					lowerThreadsPart[storeIndex] = data;
				}
			}
		}
		

		// ensure that padding region copy is complete in all threads by waiting in barrier
		pthread_barrier_wait(&barrier);
	}

	// copy data back to the original plate after all computation is done
	double *finalPart = plate0Parts[threadId];
	if (iterations % 2 == 0) {
		finalPart = plate1Parts[threadId];
	}
	for (int i = rowStart; i <= rowEnd; i++) {
		for (int j = 0; j < plateDimLength; j++) {
			int plateIndex = i * plateDimLength + j;
			int storeIndex = (i - padRowStart) * plateDimLength + j;
			double data = finalPart[storeIndex];
			plate[plateIndex] = data;
		}
	}
	
	// exit pthread
        pthread_exit(NULL);
}

} // end of namespace

using namespace pthread_stencil;

// ------------------------------------------------------------ Helper function to random plate generation


void initializePlateRandomly() {

        int elementCount = plateDimLength * plateDimLength;
        plate = new double[elementCount];

        srand (time(NULL));
        for (int i = 0; i < elementCount; i++) {
                plate[i] = ((rand() % 100) / 75.00f);
        }
}


// ----------------------------------------------------------------------------------------- Main function

int main(int argc, char *argv[]) {

	struct timeval start;
        gettimeofday(&start, NULL);

	 if (argc < 6) {
                cout << "Provide the following command line arguments \n";
                cout << "\t1. A dimension length of the square plate\n";
                cout << "\t2. The number of Jaquobi iterations\n";
		cout << "\t3. The padding row+column among successive threads data\n";
                cout << "\t4. The number of threads to use\n";
                cout << "\t5. Provide the machine name tozammel|brac\n";
                exit(EXIT_FAILURE);
        }

        // parse command line arguments
        plateDimLength = atoi(argv[1]);
        iterations = atoi(argv[2]);
        padding = atoi(argv[3]);
        threadCount = atoi(argv[4]);
        char *machine = argv[5];

	// generate a random initial plate
	initializePlateRandomly();


	struct timeval end;
        gettimeofday(&end, NULL);
        double readingTime = ((end.tv_sec + end.tv_usec / 1000000.0)
                        - (start.tv_sec + start.tv_usec / 1000000.0));
        cout << "Data reading/intialization time: " << readingTime << " Seconds\n";


	// --------------------------------------------------------parallel programming starts
	gettimeofday(&start, NULL);

	// create the plate part arrays for individual threads
	plate0Parts = new double*[threadCount];
	plate1Parts = new double*[threadCount];

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
                int status = pthread_create(&threads[i], &attr, computeStencil, (void*) &threadIds[i]);
                if (status != 0) {
                        cout << "Could not create some pthreads\n";
                        exit(EXIT_FAILURE);
                }
        }


	// --------------------------------------------------------- parallel programming ends

        // join threads
        for (int i = 0; i < threadCount; i++) {
                pthread_join(threads[i], NULL);
        }
	cout << "number of threads: " << threadCount << "\n";

	gettimeofday(&end, NULL);
        double computationTime = ((end.tv_sec + end.tv_usec / 1000000.0)
                        - (start.tv_sec + start.tv_usec / 1000000.0));
        cout << "computation time: " << computationTime << " Seconds\n";

/*

	// allocate another plate variable for computation
	int plateSize = plateDims[0].length * plateDims[1].length;
	double *plate1 = new double[plateSize];

	// synchronize the alternative version of the plate with the plate read from the file
	for (int i = 0; i < plateSize; i++) {
		plate1[i] = plate0[i];
	}

	// do the iterative stencil computation
	double *oldPlate, *newPlate;
	oldPlate = NULL;
	newPlate = NULL;
	for (int i = 0; i < iterations; i++) {
		if (i % 2 == 0) {
			oldPlate = plate0;
			newPlate = plate1;
		} else {
			oldPlate = plate1;
			newPlate = plate0;
		}
		for (int y = 1; y < plateDims[0].length - 1; y++) {

			int yIndex0 = y * plateDims[1].length;
			int yIndex1 = yIndex0 - plateDims[1].length;
			int yIndex2 = yIndex0 + plateDims[1].length;

			for (int x = 1; x < plateDims[1].length - 1; x++) {

					int index0 = yIndex0 + x;
					int index1 = yIndex0 + (x - 1);
					int index2 = yIndex0 + (x + 1);
					int index3 = yIndex1 + x;
					int index4 = yIndex2 + x;
					newPlate[index0] = 0.25 * (oldPlate[index1] +
							oldPlate[index2] +
							oldPlate[index3] +
							oldPlate[index4]);
			}
		}
	}

	// compare the computed final state with the state got from file
	bool valid = true;
	int mismatchCount = 0;
	for (int i = 0; i < plateSize; i++) {
		if (newPlate[i] - plateOut[i] > 0.1 || newPlate[i] - plateOut[i] < -0.1) {
			cout << "Mismatch found: computed value " << newPlate[i];
			cout << " read value " << plateOut[i] << "\n";
			valid = false;
			mismatchCount++;
		}
	}
	if (valid == true) {
		cout << "validation successful\n";
	} else {
		cout << "outputs mismatch in " << mismatchCount << " data points\n";
	}
*/


	return 0;
}
