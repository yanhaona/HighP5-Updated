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

namespace mpi_stencil {

//-------------------------------------------------------------- Data Dimension Information

Dimension plateDims[2];

//--------------------------------------------------------------- Data Structure References

double *plate[2];
int maxIterations;

//-------------------------------------------------------------- Partition Config Variables

int blockSize;
int padding;
int processId;
int processCount;

//-------------------------------------------------------------------- Supporting functions

void readPlateFromFile(const char *filePath) {

	int rowsPerProcess = (plateDims[0].length + processCount - 1) / processCount;
	int rowStart = processId * rowsPerProcess;
	int rowEnd = rowStart + rowsPerProcess - 1;
	if (processId > 0) {
		rowStart = rowStart - padding;
	}
	if (processId < processCount - 1) {
		rowEnd = rowEnd + padding;
	}
	if (rowEnd >= plateDims[0].length) {
		rowEnd = plateDims[0].length - 1;
	}
	int rowCount = rowEnd - rowStart + 1;
	int localEntries = rowCount * plateDims[1].length;
	plate[0] = new double[localEntries];
	plate[1] = new double[localEntries];

	TypedInputStream<double> *stream = new TypedInputStream<double>(filePath);
	int storeIndex = 0;
	stream->open();
	List<int> *indexList = new List<int>();
        for (int i = rowStart; i <= rowEnd; i++) {
		indexList->clear();
		indexList->Append(i);
		indexList->Append(0);
		plate[0][storeIndex] = stream->readElement(indexList);
		storeIndex++;
		int count = 1;
                while (count < plateDims[1].length) {
			plate[0][storeIndex] = stream->readNextElement();
			plate[1][storeIndex] = plate[0][storeIndex];
			count++;
			storeIndex++;
		}
	}
	stream->close();

	delete indexList;
	delete stream;
}

void writePlateToFile() {

	// wait for your turn
        if (processId != 0) {
                int predecessorDone = 0;
                MPI_Status status;
                MPI_Recv(&predecessorDone, 1, MPI_INT, processId - 1, 0, MPI_COMM_WORLD, &status);
        }

        List<Dimension*> *dimLengths = new List<Dimension*>;
        dimLengths->Append(&plateDims[0]);
        dimLengths->Append(&plateDims[1]);
        TypedOutputStream<double> *stream = new TypedOutputStream<double>("plate-out.bin", dimLengths, processId == 0);
	
	int rowsPerProcess = (plateDims[0].length + processCount - 1) / processCount;
	int rowStart = processId * rowsPerProcess;
	int rowEnd = rowStart + rowsPerProcess - 1;
	if (rowEnd >= plateDims[0].length) {
		rowEnd = plateDims[0].length - 1;
	}

        List<int> *indexList = new List<int>;
        int storeIndex = 0;
	if (processId > 0) {
		storeIndex += plateDims[1].length * padding;
	}
	int plateVersion = (maxIterations * padding + 1) % 2;
        stream->open();
        for (int i = rowStart; i <= rowEnd; i++) {
        	for (int j = 0; j < plateDims[1].length; j++) {
                	indexList->clear();
                        indexList->Append(i);
                        indexList->Append(j);
                        stream->writeElement(plate[plateVersion][storeIndex], indexList);
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

//--------------------------------------------------------------------- Plate Update Function

void refinePlate() {
	
	// determine the row index range of plate belonging to the current thread;
        int rowsPerProcess = (plateDims[0].length + processCount - 1) / processCount;
        int rowStart = rowsPerProcess * processId;
        int rowEnd = rowStart + rowsPerProcess - 1;
        if (rowEnd >= plateDims[0].length) {
                rowEnd = plateDims[0].length - 1;
        }

        // determine if there is upper and lower row padding to be synchronized
        int padRowStart = rowStart;
        int padRowEnd = rowEnd;
        if (processId > 0) {
                padRowStart = rowStart - padding;
        }
        if (processId < processCount - 1) {
                padRowEnd = rowEnd + padding;
        } if (padRowEnd >= plateDims[0].length) {
		padRowEnd = plateDims[0].length - 1;
	}

	// run Jacobi iterations
        int totalIteration = 0;
        for (int i = 0; i < maxIterations; i++) {

                // setup the plate part pointers for the current and next iterations
                double *input = plate[totalIteration % 2];
		double *output = plate[(totalIteration + 1) % 2];

                // perform padding number of internal iterations before a communication attempt
                for (int p = 0; p < padding; p++) {

                        // update the output cells based on Jacobi iteration logic
                        for (int y = padRowStart + 1; y <  padRowEnd - 1; y++) {

                                int yIndex0 = (y - padRowStart) * plateDims[1].length;
                                int yIndex1 = yIndex0 - plateDims[1].length;
                                int yIndex2 = yIndex0 + plateDims[1].length;

                                for (int x = 1; x < plateDims[1].length - 1; x++) {

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

		MPI_Request sendNextReq, sendPrevReq;
        	MPI_Request recvNextReq, recvPrevReq;


		if (processId < processCount - 1) {

			// do an MPI send of current process's boundary rows to be copied in the padded rows of the next process
			int sendToNextIndexStart = rowEnd - padding + 1;
			int sendBufferIndex = (sendToNextIndexStart - padRowStart) * plateDims[1].length;
			//double *bufferToSend = output + sendBufferIndex;
			double *bufferToSend = output;
			int messageLength = padding * plateDims[1].length;	
        		MPI_Isend(bufferToSend, messageLength, MPI_DOUBLE, processId + 1, 0, MPI_COMM_WORLD, &sendNextReq);
			
			// do an MPI receive of the boundary rows from the next process to current process's padding region
			int receiveFromNextIndexStart = rowEnd + 1;
			int receiveBufferIndex = (receiveFromNextIndexStart - padRowStart) * plateDims[1].length;
			double *bufferToReceive = output + receiveBufferIndex;
        		MPI_Irecv(bufferToReceive, messageLength, MPI_DOUBLE, processId + 1, 0, MPI_COMM_WORLD, &recvNextReq);
		}
		if (processId > 0) {
			int messageLength = padding * plateDims[1].length;	
			// do an MPI send of current process's boundary rows to be copied in the padded rows of the previous process
			int sendToPrevIndexStart = rowStart;
			int sendBufferIndex = (sendToPrevIndexStart - padRowStart) * plateDims[1].length;
			double *bufferToSend = output + sendBufferIndex;
        		MPI_Isend(bufferToSend, messageLength, MPI_DOUBLE, processId - 1, 0, MPI_COMM_WORLD, &sendPrevReq);
	
			// do an MPI receive of the boundary rows from the previous process to current process's padding region 
			double *bufferToReceive = output;
        		MPI_Irecv(bufferToReceive, messageLength, MPI_DOUBLE, processId - 1, 0, MPI_COMM_WORLD, &recvPrevReq);
		}
		
		// ensure asynchronous send and receives are complete
        	MPI_Status statusSNext, statusSPrev;
        	MPI_Status statusRNext, statusRPrev;
		if (processId < processCount - 1) {
        		MPI_Wait(&sendNextReq, &statusSNext);
        		//MPI_Wait(&recvNextReq, &statusRNext);
		}
		if (processId > 0) {
        		//MPI_Wait(&sendPrevReq, &statusSPrev);
        		MPI_Wait(&recvPrevReq, &statusRPrev);
		}
	}

}

} // end of namespace

using namespace mpi_stencil;


//--------------------------------------------------------------------------- Main Function

int mainMStencil(int argc, char *argv[]) {
	
	// do MPI intialization
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &processId);
	MPI_Comm_size(MPI_COMM_WORLD, &processCount);

	if (argc < 5) {
		if (processId == 0) {
			std::cout << "Pass arguments as follows\n";
			std::cout << "\t1. path of the binary input file for the 2D plate\n";
			std::cout << "\t2. the number of jacobi iterations of plate refinements\n";
			std::cout << "\t3. amount of padding rows between successive processes\n";
			std::cout << "\t4. a 0 to disable file writing or a 1 to write the final version of the plate to a file\n";
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
	maxIterations = atoi(argv[2]);
	padding = atoi(argv[3]);
	bool fileWriteMode (atoi(argv[4]) == 1);	
	const char *filePath= argv[1];
	std::ifstream file(filePath);
        if (!file.is_open()) {
                std::cout << "could not open the specified file\n";
                std::exit(EXIT_FAILURE);
        }
	readArrayDimensionInfoFromFile(file, 2, plateDims);
	file.close();

	// read input matrices
	readPlateFromFile(filePath);
	
	struct timeval memEnd;
        gettimeofday(&memEnd, NULL);

	//------------------------------------------------------------ Computation Starts 

	refinePlate();

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
		std::cout << "Plate dimension: " << plateDims[0].length << " by " << plateDims[1].length << "\n";
		std::cout << "Padding rows: " << padding << "\n";
		std::cout << "Refinement iterations: " << (maxIterations * padding) << "\n";
	}
	
	// write outputs to files
	if (fileWriteMode) {
		if (processId == 0) {
			std::cout << "Writing the result in file plate-out.bin\n";
		}
		writePlateToFile();
	}

	// do MPI teardown
	MPI_Finalize();
	return 0;
}



