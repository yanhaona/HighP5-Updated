#ifndef _H_mm
#define _H_mm

// for error reporting and diagnostics
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

// for math functions
#include <math.h>
#include <algorithm>
#include <climits>

// for tuple definitions found in the source code
#include "tuple.h"
#include <vector>

// for LPU and PPU management data structures
#include "../../src/runtime/structure.h"
#include "../../src/runtime/lpu_management.h"

// for utility routines
#include "../../src/utils/list.h"
#include "../../src/utils/hashtable.h"
#include "../../src/utils/string_utils.h"
#include "../../src/utils/common_utils.h"

// for routines related to partition functions
#include "../../src/partition-lib/index_xform.h"
#include "../../src/partition-lib/partition_mgmt.h"

// to input-output and initialization
#include "../../src/runtime/input_prompt.h"
#include "../../src/runtime/output_prompt.h"
#include "../../src/runtime/allocator.h"

// for threading
#include <pthread.h>

// for synchronization
#include "../../src/runtime/sync.h"


namespace mm {


/*-----------------------------------------------------------------------------------
processor ordering in the hardware
------------------------------------------------------------------------------------*/

const int Processor_Order[128] = {0, 64, 24, 88, 48, 112, 8, 72, 32, 96
		, 16, 80, 40, 104, 56, 120, 4, 68, 20, 84
		, 36, 100, 52, 116, 12, 76, 28, 92, 44, 108
		, 60, 124, 58, 122, 42, 106, 26, 90, 2, 66
		, 50, 114, 34, 98, 10, 74, 62, 126, 18, 82
		, 54, 118, 30, 94, 6, 70, 38, 102, 14, 78
		, 46, 110, 22, 86, 1, 65, 17, 81, 41, 105
		, 25, 89, 9, 73, 33, 97, 49, 113, 57, 121
		, 5, 69, 21, 85, 37, 101, 53, 117, 13, 77
		, 29, 93, 45, 109, 61, 125, 59, 123, 43, 107
		, 27, 91, 3, 67, 51, 115, 35, 99, 11, 75
		, 63, 127, 19, 83, 55, 119, 31, 95, 7, 71
		, 39, 103, 15, 79, 47, 111, 23, 87};

/*-----------------------------------------------------------------------------------
constants for LPSes
------------------------------------------------------------------------------------*/
const int Space_Root = 0;
const int Space_A = 1;
const int Space_A_Sub = 2;
const int Space_Count = 3;

/*-----------------------------------------------------------------------------------
constants for PPS counts
------------------------------------------------------------------------------------*/
const int Space_4_PPUs = 1;
const int Space_3_Par_4_PPUs = 2;
const int Space_2_Par_3_PPUs = 32;
const int Space_1_Par_2_PPUs = 2;

/*-----------------------------------------------------------------------------------
constants for total and par core thread counts
------------------------------------------------------------------------------------*/
const int Space_Root_Threads = 1;
const int Space_A_Threads = 128;
const int Space_A_Sub_Threads = 128;
const int Total_Threads = 128;
const int Threads_Par_Core = 1;
const int Core_Jump = 1;
/*-----------------------------------------------------------------------------------
Data structures for Array-Metadata and Environment-Links 
------------------------------------------------------------------------------------*/

class ArrayMetadata : public Metadata {
  public:
	Dimension aDims[2];
	Dimension bDims[2];
	Dimension cDims[2];

	ArrayMetadata();
	void print(std::ofstream &stream);
};
static ArrayMetadata arrayMetadata;

class EnvironmentLinks {
  public:
	double *a;
	Dimension aDims[2];
	double *b;
	Dimension bDims[2];

	void print(std::ofstream &stream);
};
static EnvironmentLinks environmentLinks;

/*-----------------------------------------------------------------------------------
Data structures for Task-Global and Thread-Local scalar variables
------------------------------------------------------------------------------------*/

class TaskGlobals {
  public:
};

class ThreadLocals {
  public:
};

/*-----------------------------------------------------------------------------------
function to initialize the content reference objects of LPSes
------------------------------------------------------------------------------------*/
void initializeRootLPSContent(EnvironmentLinks *envLinks, ArrayMetadata *metadata);
void initializeLPSesContents(ArrayMetadata *metadata);

/*-----------------------------------------------------------------------------------
functions for retrieving partition counts in different LPSes
------------------------------------------------------------------------------------*/
int *getLPUsCountOfSpaceA(int ppuCount, Dimension cDim1, int k, Dimension cDim2, int l);
int *getLPUsCountOfSpaceA_Sub(int ppuCount, Dimension aDim2, int q);

/*-----------------------------------------------------------------------------------
functions for getting data ranges along different dimensions of an LPU
-----------------------------------------------------------------------------------*/
void getaPartForSpaceALpu(PartDimension *aLpuDims, 
		PartDimension *aParentLpuDims, 
		int *lpuCount, int *lpuId, int k);
void getbPartForSpaceALpu(PartDimension *bLpuDims, 
		PartDimension *bParentLpuDims, 
		int *lpuCount, int *lpuId, int l);
void getcPartForSpaceALpu(PartDimension *cLpuDims, 
		PartDimension *cParentLpuDims, 
		int *lpuCount, int *lpuId, int k, int l);
void getaPartForSpaceA_SubLpu(PartDimension *aLpuDims, 
		PartDimension *aParentLpuDims, 
		int *lpuCount, int *lpuId, int q);
void getbPartForSpaceA_SubLpu(PartDimension *bLpuDims, 
		PartDimension *bParentLpuDims, 
		int *lpuCount, int *lpuId, int q);

/*-----------------------------------------------------------------------------------
Data structures representing LPS and LPU contents 
------------------------------------------------------------------------------------*/

class SpaceRoot_Content {
  public:
	double *a;
	double *b;
	double *c;
};
static SpaceRoot_Content spaceRootContent;

class SpaceRoot_LPU : public LPU {
  public:
	double *a;
	PartDimension aPartDims[2];
	double *b;
	PartDimension bPartDims[2];
	double *c;
	PartDimension cPartDims[2];

	void print(std::ofstream &stream, int indent);
};

class SpaceA_Content {
  public:
	double *a;
	double *b;
	double *c;
};
static SpaceA_Content spaceAContent;

class SpaceA_LPU : public LPU {
  public:
	double *a;
	PartDimension aPartDims[2];
	double *b;
	PartDimension bPartDims[2];
	double *c;
	PartDimension cPartDims[2];
	int lpuId[2];

	void print(std::ofstream &stream, int indent);
};

class SpaceA_Sub_Content {
  public:
	double *a;
	double *b;
	double *c;
};
static SpaceA_Sub_Content spaceA_SubContent;

class SpaceA_Sub_LPU : public LPU {
  public:
	double *a;
	PartDimension aPartDims[2];
	double *b;
	PartDimension bPartDims[2];
	double *c;
	PartDimension cPartDims[2];
	int lpuId[1];

	void print(std::ofstream &stream, int indent);
};


/*-----------------------------------------------------------------------------------
function to generate PPU IDs and PPU group IDs for a thread
------------------------------------------------------------------------------------*/
ThreadIds *getPpuIdsForThread(int threadNo);

/*-----------------------------------------------------------------------------------
Thread-State implementation class for the task
------------------------------------------------------------------------------------*/

class ThreadStateImpl : public ThreadState {
  public:
	ThreadStateImpl(int lpsCount, int *lpsDimensions, 
			int *partitionArgs, 
			ThreadIds *threadIds) 
		: ThreadState(lpsCount, lpsDimensions, partitionArgs, threadIds) {}
	void setLpsParentIndexMap();
        void setRootLpu(Metadata *metadata);
        void setRootLpu(LPU *rootLpu);
	void initializeLPUs();
        int *computeLpuCounts(int lpsId);
        LPU *computeNextLpu(int lpsId, int *lpuCounts, int *nextLpuId);
};

/*-----------------------------------------------------------------------------------
function for initializing environment-links object
------------------------------------------------------------------------------------*/

EnvironmentLinks initiateEnvLinks(MMEnvironment *environment);

/*-----------------------------------------------------------------------------------
function for initializing root LPU from environment
------------------------------------------------------------------------------------*/

SpaceRoot_LPU *initiateRootLpu(MMEnvironment *environment, ArrayMetadata *metadata);

/*-----------------------------------------------------------------------------------
function for executing task
------------------------------------------------------------------------------------*/

void execute(MMEnvironment *environment, 
		MMPartition partition, 
		std::ofstream &logFile);

/*-----------------------------------------------------------------------------------
function for the initialize block
------------------------------------------------------------------------------------*/
void initializeTask(ArrayMetadata *arrayMetadata, 
		EnvironmentLinks environmentLinks, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		MMPartition partition);

/*-----------------------------------------------------------------------------------
functions for compute stages 
------------------------------------------------------------------------------------*/

int block_multiply_matrices(SpaceA_Sub_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, MMPartition partition);


/*-----------------------------------------------------------------------------------
The run method for thread simulating the task flow 
------------------------------------------------------------------------------------*/

void run(ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		MMPartition partition, ThreadStateImpl *threadState);

/*-----------------------------------------------------------------------------------
Data structure and function for Pthreads
------------------------------------------------------------------------------------*/

class PThreadArg {
  public:
	const char *taskName;
	ArrayMetadata *metadata;
	TaskGlobals *taskGlobals;
	ThreadLocals *threadLocals;
	MMPartition partition;
	ThreadStateImpl *threadState;
};

void *runPThreads(void *argument);

}
#endif
