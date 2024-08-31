#include <string>
#include <iostream>
#include <cstdlib>
#include <cctype>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <pthread.h>
#include <math.h>
#include "../cpuCores.h"

namespace pthread_monte_simple {

//------------------------------------------------------------------- Supporting structures

typedef struct {
        int top;
        int left;
        int right;
        int bottom;
} Rectangle;

//---------------------------------------------------------------------------- program data

double *sub_area_estimates;
int cell_length;
int grid_dimension;
int points_per_cell;

//----------------------------------------------------------------- Thread Interaction data

int threadCount;

//------------------------------------------------------------------------- Thread Function

void *computeSubEstimate(void *arg) {

	int threadId = *((int*) arg);
        double cell_size = cell_length * cell_length;

	// let different threads stride through the rows of different grid blocks
	for (int y = threadId; y < grid_dimension; y += threadCount) {

		int rowIndex = y * grid_dimension;

                for (int x = 0; x < grid_dimension; x++) {

                        // determine the grid cell boundaries
                        Rectangle cell;
                        cell.top = cell_length * (y + 1) - 1;
                        cell.bottom = cell_length * y;
                        cell.left = cell_length * x;
                        cell.right = cell_length * (x + 1) - 1;

			// do appropriate number of random sampling
                        int points_inside = 0;
                        for (int s = 0; s < points_per_cell; s++) {
                                int x = rand_r((unsigned int *) &threadId) % cell_length + cell.left;
                                int y = rand_r((unsigned int *) &threadId) % cell_length + cell.bottom;
                                double result = 10 * sin(pow(x, 2)) + 50 * cos(pow(y, 3));
                                if (result <= 0.0) {
                                        points_inside++;
                                }
                        }
                        	
			// calculate the subarea estimate for the grid cell
                        double estimate = cell_size * points_inside / points_per_cell;
				
                        sub_area_estimates[rowIndex + x] = estimate;
                }
        }
	
	// exit thread
	pthread_exit(NULL);
}

} // end of namespace

//--------------------------------------------------------------------------- Main Function

using namespace pthread_monte_simple;

int main(int argc, char *argv[]) {

	if (argc < 6) {
                std::cout << "provide cell length, grid dimension in number of cells, points to";
                std::cout << " be generated per cell,\n";
		std::cout << "then specify the number of threads to be used as the last argument\n";
		std::cout << "finally, specify the machine tozammel|brac\n";
                std::exit(EXIT_FAILURE);
        }
	
	struct timeval start;
	gettimeofday(&start, NULL);

	// initialize global metadata variables
        cell_length = atoi(argv[1]);
        grid_dimension = atoi(argv[2]);
        points_per_cell = atoi(argv[3]);
	threadCount = atoi(argv[4]);
	char* machine = argv[5];

	// initialize the random number generator
        srand(time(NULL));

	// allocate necessary data structures
	int grid_size = grid_dimension * grid_dimension;
	sub_area_estimates = new double[grid_size];

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

                int status = pthread_create(&threads[i], &attr, computeSubEstimate, (void*) &threadIds[i]);
                if (status != 0) {
                        std::cout << "Could not create some pthreads\n";
                        std::exit(EXIT_FAILURE);
                }
        }

	// join threads
	for (int i = 0; i < threadCount; i++) {
                pthread_join(threads[i], NULL);
        }

	// sum up the sub-area estimates to construct the final estimate of area
	double estimate = 0.0;
	for (int i = 0; i < grid_size; i++) {
		estimate += sub_area_estimates[i];
	}

	struct timeval end;
	gettimeofday(&end, NULL);
	double executionTime = ((end.tv_sec + end.tv_usec / 1000000.0)
			- (start.tv_sec + start.tv_usec / 1000000.0));
	std::cout << "thread count: " << threadCount << "\n";
	std::cout << "execution time: " << executionTime << " Seconds\n";
	std::cout << "Estimated area under the curve: " << estimate << "\n";

	return 0;
}

