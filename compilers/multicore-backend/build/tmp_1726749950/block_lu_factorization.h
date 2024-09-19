#ifndef _H_bluf
#define _H_bluf

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


namespace bluf {


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
const int Space_B = 2;
const int Space_C = 3;
const int Space_C_Sub = 4;
const int Space_Count = 5;

/*-----------------------------------------------------------------------------------
constants for PPS counts
------------------------------------------------------------------------------------*/
const int Space_8_PPUs = 1;
const int Space_7_Par_8_PPUs = 2;
const int Space_6_Par_7_PPUs = 2;
const int Space_5_Par_6_PPUs = 2;
const int Space_4_Par_5_PPUs = 2;
const int Space_3_Par_4_PPUs = 2;
const int Space_2_Par_3_PPUs = 2;
const int Space_1_Par_2_PPUs = 2;

/*-----------------------------------------------------------------------------------
constants for total and par core thread counts
------------------------------------------------------------------------------------*/
const int Space_Root_Threads = 1;
const int Space_A_Threads = 1;
const int Space_B_Threads = 128;
const int Space_C_Threads = 128;
const int Space_C_Sub_Threads = 128;
const int Total_Threads = 128;
const int Threads_Par_Core = 1;
const int Core_Jump = 1;
/*-----------------------------------------------------------------------------------
Data structures for Array-Metadata and Environment-Links 
------------------------------------------------------------------------------------*/

class ArrayMetadata : public Metadata {
  public:
	Dimension aDims[2];
	Dimension lDims[2];
	Dimension l_blockDims[2];
	Dimension l_columnDims[1];
	Dimension l_rowDims[1];
	Dimension pDims[1];
	Dimension p_columnDims[1];
	Dimension uDims[2];
	Dimension u_blockDims[2];

	ArrayMetadata();
	void print(std::ofstream &stream);
};
static ArrayMetadata arrayMetadata;

class EnvironmentLinks {
  public:
	double *a;
	Dimension aDims[2];

	void print(std::ofstream &stream);
};
static EnvironmentLinks environmentLinks;

/*-----------------------------------------------------------------------------------
Data structures for Task-Global and Thread-Local scalar variables
------------------------------------------------------------------------------------*/

class TaskGlobals {
  public:
	int block_size;
	int pivot;
	Range row_range;
};

class ThreadLocals {
  public:
	int k;
	int r;
};

/*-----------------------------------------------------------------------------------
function to initialize the content reference objects of LPSes
------------------------------------------------------------------------------------*/
void initializeRootLPSContent(EnvironmentLinks *envLinks, ArrayMetadata *metadata);
void initializeLPSesContents(ArrayMetadata *metadata);

/*-----------------------------------------------------------------------------------
functions for retrieving partition counts in different LPSes
------------------------------------------------------------------------------------*/
int *getLPUsCountOfSpaceB(int ppuCount, Dimension aDim2, int b);
int *getLPUsCountOfSpaceC(int ppuCount, Dimension uDim1, int b, Dimension uDim2);
int *getLPUsCountOfSpaceC_Sub(int ppuCount, Dimension u_blockDim2, int b);

/*-----------------------------------------------------------------------------------
functions for getting data ranges along different dimensions of an LPU
-----------------------------------------------------------------------------------*/
void getaPartForSpaceBLpu(PartDimension *aLpuDims, 
		PartDimension *aParentLpuDims, 
		int *lpuCount, int *lpuId, int b);
void getlPartForSpaceBLpu(PartDimension *lLpuDims, 
		PartDimension *lParentLpuDims, 
		int *lpuCount, int *lpuId, int b);
void getl_columnPartForSpaceBLpu(PartDimension *l_columnLpuDims, 
		PartDimension *l_columnParentLpuDims, 
		int *lpuCount, int *lpuId, int b);
void getuPartForSpaceBLpu(PartDimension *uLpuDims, 
		PartDimension *uParentLpuDims, 
		int *lpuCount, int *lpuId, int b);
void getu_blockPartForSpaceBLpu(PartDimension *u_blockLpuDims, 
		PartDimension *u_blockParentLpuDims, 
		int *lpuCount, int *lpuId, int b);
void getl_blockPartForSpaceCLpu(PartDimension *l_blockLpuDims, 
		PartDimension *l_blockParentLpuDims, 
		int *lpuCount, int *lpuId, int b);
void getuPartForSpaceCLpu(PartDimension *uLpuDims, 
		PartDimension *uParentLpuDims, 
		int *lpuCount, int *lpuId, int b);
void getu_blockPartForSpaceCLpu(PartDimension *u_blockLpuDims, 
		PartDimension *u_blockParentLpuDims, 
		int *lpuCount, int *lpuId, int b);
void getl_blockPartForSpaceC_SubLpu(PartDimension *l_blockLpuDims, 
		PartDimension *l_blockParentLpuDims, 
		int *lpuCount, int *lpuId, int b);
void getu_blockPartForSpaceC_SubLpu(PartDimension *u_blockLpuDims, 
		PartDimension *u_blockParentLpuDims, 
		int *lpuCount, int *lpuId, int b);

/*-----------------------------------------------------------------------------------
Data structures representing LPS and LPU contents 
------------------------------------------------------------------------------------*/

class SpaceRoot_Content {
  public:
	double *a;
	double *l;
	double *l_block;
	double *l_column;
	double *l_row;
	int *p;
	double *p_column;
	double *u;
	double *u_block;
};
static SpaceRoot_Content spaceRootContent;

class SpaceRoot_LPU : public LPU {
  public:
	double *a;
	PartDimension aPartDims[2];
	double *l;
	PartDimension lPartDims[2];
	double *l_block;
	PartDimension l_blockPartDims[2];
	double *l_column;
	PartDimension l_columnPartDims[1];
	double *l_row;
	PartDimension l_rowPartDims[1];
	int *p;
	PartDimension pPartDims[1];
	double *p_column;
	PartDimension p_columnPartDims[1];
	double *u;
	PartDimension uPartDims[2];
	double *u_block;
	PartDimension u_blockPartDims[2];

	void print(std::ofstream &stream, int indent);
};

class SpaceA_Content {
  public:
	double *a;
	double *l_block;
	double *l_column;
	double *l_row;
	int *p;
	double *p_column;
	double *u_block;
};
static SpaceA_Content spaceAContent;

class SpaceA_LPU : public LPU {
  public:
	double *a;
	PartDimension aPartDims[2];
	double *l_block;
	PartDimension l_blockPartDims[2];
	double *l_column;
	PartDimension l_columnPartDims[1];
	double *l_row;
	PartDimension l_rowPartDims[1];
	int *p;
	PartDimension pPartDims[1];
	double *p_column;
	PartDimension p_columnPartDims[1];
	double *u_block;
	PartDimension u_blockPartDims[2];

	void print(std::ofstream &stream, int indent);
};

class SpaceB_Content {
  public:
	double *a;
	double *l;
	double *l_block;
	double *l_column;
	double *l_row;
	double *p_column;
	double *u;
	double *u_block;
};
static SpaceB_Content spaceBContent;

class SpaceB_LPU : public LPU {
  public:
	double *a;
	PartDimension aPartDims[2];
	double *l;
	PartDimension lPartDims[2];
	double *l_block;
	PartDimension l_blockPartDims[2];
	double *l_column;
	PartDimension l_columnPartDims[1];
	double *l_row;
	PartDimension l_rowPartDims[1];
	double *p_column;
	PartDimension p_columnPartDims[1];
	double *u;
	PartDimension uPartDims[2];
	double *u_block;
	PartDimension u_blockPartDims[2];
	int lpuId[1];

	void print(std::ofstream &stream, int indent);
};

class SpaceC_Content {
  public:
	double *l_block;
	double *u;
	double *u_block;
};
static SpaceC_Content spaceCContent;

class SpaceC_LPU : public LPU {
  public:
	double *l_block;
	PartDimension l_blockPartDims[2];
	double *u;
	PartDimension uPartDims[2];
	double *u_block;
	PartDimension u_blockPartDims[2];
	int lpuId[2];

	void print(std::ofstream &stream, int indent);
};

class SpaceC_Sub_Content {
  public:
	double *l_block;
	double *u_block;
	double *u;
};
static SpaceC_Sub_Content spaceC_SubContent;

class SpaceC_Sub_LPU : public LPU {
  public:
	double *l_block;
	PartDimension l_blockPartDims[2];
	double *u_block;
	PartDimension u_blockPartDims[2];
	double *u;
	PartDimension uPartDims[2];
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
global synchronization primitives
------------------------------------------------------------------------------------*/

static RS *row_rangeStage4No8DSyncs[Space_Root_Threads];
static Barrier *row_rangeStage4No8ReverseSyncs[Space_Root_Threads];
static RS *uStage14No4DSyncs[Space_B_Threads];
static Barrier *uStage14No4ReverseSyncs[Space_B_Threads];
static RS *uStage20No1USyncs[Space_B_Threads];
static Barrier *uStage20No1ReverseSyncs[Space_B_Threads];
static RS *row_rangeStage4No1DSyncs[Space_Root_Threads];
static Barrier *row_rangeStage4No1ReverseSyncs[Space_Root_Threads];
static RS *u_blockStage16No1DSyncs[Space_B_Threads];
static Barrier *u_blockStage16No1ReverseSyncs[Space_B_Threads];
static RS *l_blockStage17No1DSyncs[Space_A_Threads];
static Barrier *l_blockStage17No1ReverseSyncs[Space_A_Threads];
static RS *pivotStage6No1USyncs[Space_Root_Threads];
static Barrier *pivotStage6No1ReverseSyncs[Space_Root_Threads];
static RS *pivotStage6No2RSyncs[Space_Root_Threads];
static Barrier *pivotStage6No2ReverseSyncs[Space_Root_Threads];
static RS *l_rowStage9No1RSyncs[Space_A_Threads];
static Barrier *l_rowStage9No1ReverseSyncs[Space_A_Threads];
static RS *l_columnStage12No1USyncs[Space_A_Threads];
static Barrier *l_columnStage12No1ReverseSyncs[Space_A_Threads];
static RS *p_columnStage13No1DSyncs[Space_A_Threads];
static Barrier *p_columnStage13No1ReverseSyncs[Space_A_Threads];

/*-----------------------------------------------------------------------------------
Initializer function for global synchronization primitives
------------------------------------------------------------------------------------*/
void initializeSyncPrimitives();
/*-----------------------------------------------------------------------------------
data structure and function for initializing thread's sync primitives
------------------------------------------------------------------------------------*/

class ThreadSyncPrimitive {
  public:
	RS *row_rangeStage4No8DSync;
	Barrier *row_rangeStage4No8ReverseSync;
	RS *uStage14No4DSync;
	Barrier *uStage14No4ReverseSync;
	RS *uStage20No1USync;
	Barrier *uStage20No1ReverseSync;
	RS *row_rangeStage4No1DSync;
	Barrier *row_rangeStage4No1ReverseSync;
	RS *u_blockStage16No1DSync;
	Barrier *u_blockStage16No1ReverseSync;
	RS *l_blockStage17No1DSync;
	Barrier *l_blockStage17No1ReverseSync;
	RS *pivotStage6No1USync;
	Barrier *pivotStage6No1ReverseSync;
	RS *pivotStage6No2RSync;
	Barrier *pivotStage6No2ReverseSync;
	RS *l_rowStage9No1RSync;
	Barrier *l_rowStage9No1ReverseSync;
	RS *l_columnStage12No1USync;
	Barrier *l_columnStage12No1ReverseSync;
	RS *p_columnStage13No1DSync;
	Barrier *p_columnStage13No1ReverseSync;
};

ThreadSyncPrimitive *getSyncPrimitives(ThreadIds *threadIds);
/*-----------------------------------------------------------------------------------
function for initializing environment-links object
------------------------------------------------------------------------------------*/

EnvironmentLinks initiateEnvLinks(BLUFEnvironment *environment);

/*-----------------------------------------------------------------------------------
function for initializing root LPU from environment
------------------------------------------------------------------------------------*/

SpaceRoot_LPU *initiateRootLpu(BLUFEnvironment *environment, ArrayMetadata *metadata);

/*-----------------------------------------------------------------------------------
function for executing task
------------------------------------------------------------------------------------*/

void execute(BLUFEnvironment *environment, 
		BLUFPartition partition, 
		std::ofstream &logFile);

/*-----------------------------------------------------------------------------------
function for the initialize block
------------------------------------------------------------------------------------*/
void initializeTask(ArrayMetadata *arrayMetadata, 
		EnvironmentLinks environmentLinks, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		BLUFPartition partition);

/*-----------------------------------------------------------------------------------
functions for compute stages 
------------------------------------------------------------------------------------*/

int prepare_l_and_u(SpaceB_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, BLUFPartition partition);

int calculate_row_range(SpaceA_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, BLUFPartition partition);

int select_pivot(SpaceB_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, BLUFPartition partition);

int store_pivot(SpaceA_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, BLUFPartition partition);

int interchange_rows(SpaceB_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, BLUFPartition partition);

int update_l(SpaceB_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, BLUFPartition partition);

int update_u_rows_block(SpaceB_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, BLUFPartition partition);

int collect_l_column_part(SpaceB_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, BLUFPartition partition);

int generate_pivot_column(SpaceA_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, BLUFPartition partition);

int update_u_columns(SpaceB_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, BLUFPartition partition);

int copy_updated_block_of_u(SpaceB_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, BLUFPartition partition);

int copy_updated_block_of_l(SpaceB_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, BLUFPartition partition);

int block_mm_multiply_subtract(SpaceC_Sub_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, BLUFPartition partition);


/*-----------------------------------------------------------------------------------
The run method for thread simulating the task flow 
------------------------------------------------------------------------------------*/

void run(ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		BLUFPartition partition, ThreadStateImpl *threadState);

/*-----------------------------------------------------------------------------------
Data structure and function for Pthreads
------------------------------------------------------------------------------------*/

class PThreadArg {
  public:
	const char *taskName;
	ArrayMetadata *metadata;
	TaskGlobals *taskGlobals;
	ThreadLocals *threadLocals;
	BLUFPartition partition;
	ThreadStateImpl *threadState;
};

void *runPThreads(void *argument);

}
#endif
