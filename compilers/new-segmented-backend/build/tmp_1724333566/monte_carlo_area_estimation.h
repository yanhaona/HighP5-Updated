#ifndef _H_mcae
#define _H_mcae

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

// for user defined functions found in the source code
#include "function.h"

// for LPU and PPU management data structures
#include "../../../common-libs/domain-obj/structure.h"
#include "../../src/runtime/common/lpu_management.h"

// for utility routines
#include "../../../common-libs/utils/list.h"
#include "../../../common-libs/utils/utility.h"
#include "../../../common-libs/utils/hashtable.h"
#include "../../../common-libs/utils/string_utils.h"
#include "../../../common-libs/utils/common_utils.h"
#include "../../../common-libs/utils/interval.h"
#include "../../../common-libs/utils/binary_search.h"
#include "../../../common-libs/utils/id_generation.h"

// for routines related to partition functions
#include "../../src/runtime/partition-lib/partition.h"
#include "../../src/runtime/partition-lib/index_xform.h"
#include "../../src/runtime/partition-lib/partition_mgmt.h"

// for memory management
#include "../../src/runtime/memory-management/allocation.h"
#include "../../src/runtime/memory-management/part_tracking.h"
#include "../../src/runtime/memory-management/part_generation.h"
#include "../../src/runtime/memory-management/part_management.h"

// for input-output
#include "../../src/runtime/file-io/stream.h"
#include "../../src/runtime/file-io/data_handler.h"

// for communication
#include "../../src/runtime/communication/part_folding.h"
#include "../../src/runtime/communication/part_config.h"
#include "../../src/runtime/communication/part_distribution.h"
#include "../../src/runtime/communication/confinement_mgmt.h"
#include "../../src/runtime/communication/data_transfer.h"
#include "../../src/runtime/communication/comm_buffer.h"
#include "../../src/runtime/communication/comm_statistics.h"
#include "../../src/runtime/communication/communicator.h"
#include "../../src/runtime/communication/scalar_communicator.h"
#include "../../src/runtime/communication/array_communicator.h"

// for task and program environment management and interaction
#include "../../src/runtime/environment/environment.h"
#include "../../src/runtime/environment/env_instruction.h"
#include "../../src/runtime/environment/array_transfer.h"

// for threading
#include <pthread.h>

// for MPI
#include <mpi.h>

// for synchronization
#include "../../src/runtime/common/sync.h"

// for reductions
#include "../../src/runtime/reduction/reduction_barrier.h"
#include "../../src/runtime/reduction/task_global_reduction.h"
#include "../../src/runtime/reduction/non_task_global_reduction.h"
#include "../../../common-libs/domain-obj/constant.h"

// for minimum and maximum values of numeric types
#include <limits.h>
#include <float.h>


namespace mcae {


/*--------------------------------------------------------------------------------------------------------------
processor ordering in the hardware
--------------------------------------------------------------------------------------------------------------*/

const int Processor_Order[4] = {0, 1, 2, 3};

/*--------------------------------------------------------------------------------------------------------------
constants for LPSes
--------------------------------------------------------------------------------------------------------------*/

const int Space_Root = 0;
const int Space_A = 1;
const int Space_B = 2;
const int Space_C = 3;
const int Space_Count = 4;

/*--------------------------------------------------------------------------------------------------------------
constants for PPS counts
--------------------------------------------------------------------------------------------------------------*/

const int Space_4_PPUs = 1;
const int Space_3_Par_4_PPUs = 10;
const int Space_2_Par_3_PPUs = 2;
const int Space_1_Par_2_PPUs = 2;

/*--------------------------------------------------------------------------------------------------------------
thread count constants
--------------------------------------------------------------------------------------------------------------*/

//----------------------------------------------------------------------------------------------System Constants
const int Max_Space_Root_Threads = 1;
const int Max_Space_A_Threads = 40;
const int Max_Space_B_Threads = 40;
const int Max_Space_C_Threads = 40;
const int Max_Total_Threads = 40;
static int Total_Threads = 40;

//---------------------------------------------------------------------------------------------Segment Constants
const int Threads_Per_Segment = 4;
const int Max_Segments_Count = 10;
const int Space_Root_Threads_Per_Segment = 1;
const int Space_A_Threads_Per_Segment = 4;
const int Space_B_Threads_Per_Segment = 4;
const int Space_C_Threads_Per_Segment = 4;

//--------------------------------------------------------------------------------------------Hardware Constants
const int Threads_Per_Core = 1;
const int Processors_Per_Phy_Unit = 4;
const int Core_Jump = 1;

/*--------------------------------------------------------------------------------------------------------------
Data structures for Array-Metadata and Environment-Links
--------------------------------------------------------------------------------------------------------------*/

class ArrayMetadata : public Metadata {
  public:
	Dimension estimate_diffsDims[2];
	Dimension gridDims[2];
	Dimension internal_pointsDims[2];
	Dimension local_estimatesDims[2];

	ArrayMetadata();
	void print(std::ofstream &stream);
};
static ArrayMetadata arrayMetadata;

class EnvironmentLinks {
  public:

	void print(std::ofstream &stream);
};
static EnvironmentLinks environmentLinks;


/*--------------------------------------------------------------------------------------------------------------
data structures representing LPUs
--------------------------------------------------------------------------------------------------------------*/

//----------------------------------------------------------------------------------------------------Space Root
class SpaceRoot_LPU : public LPU {
  public:
	double *estimate_diffs;
	PartDimension estimate_diffsPartDims[2];
	Rectangle *grid;
	PartDimension gridPartDims[2];
	int *internal_points;
	PartDimension internal_pointsPartDims[2];
	double *local_estimates;
	PartDimension local_estimatesPartDims[2];

	void print(std::ofstream &stream, int indent);
};

//-------------------------------------------------------------------------------------------------------Space A
class SpaceA_LPU : public LPU {
  public:
	double *local_estimates;
	PartDimension local_estimatesPartDims[2];

	void print(std::ofstream &stream, int indent);
};

//-------------------------------------------------------------------------------------------------------Space B
class SpaceB_LPU : public LPU {
  public:
	double *estimate_diffs;
	PartDimension estimate_diffsPartDims[2];
	Rectangle *grid;
	PartDimension gridPartDims[2];
	int *internal_points;
	PartDimension internal_pointsPartDims[2];
	double *local_estimates;
	PartDimension local_estimatesPartDims[2];
	int lpuId[2];

	void print(std::ofstream &stream, int indent);
};

//-------------------------------------------------------------------------------------------------------Space C
class SpaceC_LPU : public LPU {
  public:
	Rectangle *grid;
	PartDimension gridPartDims[2];
	int *internal_points;
	PartDimension internal_pointsPartDims[2];
	int lpuId[2];

	void print(std::ofstream &stream, int indent);
};

/*--------------------------------------------------------------------------------------------------------------
Data structures for Task-Global and Thread-Local scalar variables
--------------------------------------------------------------------------------------------------------------*/

class TaskGlobals {
  public:
	double area;
	int cell_length;
	int max_rounds;
	int points_per_cell;
	double precision_threshold;
};

class ThreadLocals {
  public:
	int placement_result;
	int round;
};

/*--------------------------------------------------------------------------------------------------------------
functions for generating partition configuration objects for data structures
--------------------------------------------------------------------------------------------------------------*/

//-------------------------------------------------------------------------------------------------------Space A
DataPartitionConfig *getlocal_estimatesConfigForSpaceA(ArrayMetadata *metadata, 
		MCAEPartition partition, 
		int ppuCount, 
		Hashtable<DataPartitionConfig*> *configMap);

//-------------------------------------------------------------------------------------------------------Space B
DataPartitionConfig *getestimate_diffsConfigForSpaceB(ArrayMetadata *metadata, 
		MCAEPartition partition, 
		int ppuCount, 
		Hashtable<DataPartitionConfig*> *configMap);
DataPartitionConfig *getgridConfigForSpaceB(ArrayMetadata *metadata, 
		MCAEPartition partition, 
		int ppuCount, 
		Hashtable<DataPartitionConfig*> *configMap);
DataPartitionConfig *getinternal_pointsConfigForSpaceB(ArrayMetadata *metadata, 
		MCAEPartition partition, 
		int ppuCount, 
		Hashtable<DataPartitionConfig*> *configMap);
DataPartitionConfig *getlocal_estimatesConfigForSpaceB(ArrayMetadata *metadata, 
		MCAEPartition partition, 
		int ppuCount, 
		Hashtable<DataPartitionConfig*> *configMap);

//-------------------------------------------------------------------------------------------------------Space C
DataPartitionConfig *getgridConfigForSpaceC(ArrayMetadata *metadata, 
		MCAEPartition partition, 
		int ppuCount, 
		Hashtable<DataPartitionConfig*> *configMap);
DataPartitionConfig *getinternal_pointsConfigForSpaceC(ArrayMetadata *metadata, 
		MCAEPartition partition, 
		int ppuCount, 
		Hashtable<DataPartitionConfig*> *configMap);

//---------------------------------------------------------------------------Partition Configuration Accumulator
Hashtable<DataPartitionConfig*> *getDataPartitionConfigMap(ArrayMetadata *metadata, 
		MCAEPartition partition, int *ppuCounts);


/*--------------------------------------------------------------------------------------------------------------
functions for generating memory blocks for data parts of various LPUs
--------------------------------------------------------------------------------------------------------------*/

LpsContent *genSpaceBContent(List<ThreadState*> *threads, ArrayMetadata *metadata, 
		TaskEnvironment *environment, 
		MCAEPartition partition, 
		Hashtable<DataPartitionConfig*> *partConfigMap);
LpsContent *genSpaceCContent(List<ThreadState*> *threads, ArrayMetadata *metadata, 
		TaskEnvironment *environment, 
		MCAEPartition partition, 
		Hashtable<DataPartitionConfig*> *partConfigMap);

/*--------------------------------------------------------------------------------------------------------------
functions for creating containers of non-task-global reduction results
--------------------------------------------------------------------------------------------------------------*/

void prepareSpaceBReductionResultContainers(List<ThreadState*> *threads, TaskData *taskData);

//-----------------------------------------------------------------------------------------Task Data Initializer
TaskData *initializeTaskData(List<ThreadState*> *threads, ArrayMetadata *metadata, 
		TaskEnvironment *environment, 
		SegmentState *segment, 
		MCAEPartition partition, int *ppuCounts);


/*--------------------------------------------------------------------------------------------------------------
data structure spacific part reader and writer subclasses
--------------------------------------------------------------------------------------------------------------*/

//-------------------------------------------------------------------------------------------------------Space B

//-------------------------------------------------------------------------------------------------------Space C

/*--------------------------------------------------------------------------------------------------------------
file I/O for environmental data structures
--------------------------------------------------------------------------------------------------------------*/

//------------------------------------------------------------------------------------Parts-reader map generator
Hashtable<PartReader*> *generateReadersMap(TaskData *taskData, 
			SegmentState *segment, 
			Hashtable<DataPartitionConfig*> *partConfigMap);

//------------------------------------------------------------------------------------Parts-writer map generator
Hashtable<PartWriter*> *generateWritersMap(TaskData *taskData, 
			SegmentState *segment, 
			Hashtable<DataPartitionConfig*> *partConfigMap);

/*--------------------------------------------------------------------------------------------------------------
functions for retrieving partition counts in different LPSes
--------------------------------------------------------------------------------------------------------------*/

int *getLPUsCountOfSpaceB(Hashtable<DataPartitionConfig*> *partConfigMap);
int *getLPUsCountOfSpaceC(Hashtable<DataPartitionConfig*> *partConfigMap, Dimension gridDim1, Dimension gridDim2);

/*--------------------------------------------------------------------------------------------------------------
functions for generating LPUs given LPU Ids
--------------------------------------------------------------------------------------------------------------*/

void generateSpaceALpu(ThreadState *threadState, 
			Hashtable<DataPartitionConfig*> *partConfigMap, 
			TaskData *taskData);
void generateSpaceBLpu(ThreadState *threadState, 
			Hashtable<DataPartitionConfig*> *partConfigMap, 
			TaskData *taskData);
void generateSpaceCLpu(ThreadState *threadState, 
			Hashtable<DataPartitionConfig*> *partConfigMap, 
			TaskData *taskData);

/*--------------------------------------------------------------------------------------------------------------
functions to generate PPU IDs and PPU group IDs for a thread
--------------------------------------------------------------------------------------------------------------*/

ThreadIds *getPpuIdsForThread(int threadNo);

void adjustPpuCountsAndGroupSizes(ThreadIds *threadId);


/*--------------------------------------------------------------------------------------------------------------
Thread-State implementation class for the task
--------------------------------------------------------------------------------------------------------------*/

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
        LPU *computeNextLpu(int lpsId);
	void initializeReductionResultMap();
};


/*--------------------------------------------------------------------------------------------------------------
global synchronization primitives
--------------------------------------------------------------------------------------------------------------*/

static RS *estimate_diffsStage6No1DSyncs[Space_B_Threads_Per_Segment];
static Barrier *estimate_diffsStage6No1ReverseSyncs[Space_B_Threads_Per_Segment];
static RS *estimate_diffsStage12No1DSyncs[Space_B_Threads_Per_Segment];
static Barrier *estimate_diffsStage12No1ReverseSyncs[Space_B_Threads_Per_Segment];

/*--------------------------------------------------------------------------------------------------------------
Initializer function for global synchronization primitives
--------------------------------------------------------------------------------------------------------------*/

void initializeSyncPrimitives();


/*--------------------------------------------------------------------------------------------------------------
data structure and function for initializing thread's sync primitives
--------------------------------------------------------------------------------------------------------------*/

class ThreadSyncPrimitive {
  public:
	RS *estimate_diffsStage6No1DSync;
	Barrier *estimate_diffsStage6No1ReverseSync;
	RS *estimate_diffsStage12No1DSync;
	Barrier *estimate_diffsStage12No1ReverseSync;
};

ThreadSyncPrimitive *getSyncPrimitives(ThreadIds *threadIds);

/*--------------------------------------------------------------------------------------------------------------
Reduction Primitives Management
--------------------------------------------------------------------------------------------------------------*/

//---------------------------------------------------------------------Primitive for variable 'placement_result'

class ReductionPrimitive_placement_result : public NonTaskGlobalReductionPrimitive {
  public: 
	ReductionPrimitive_placement_result(int localParticipants);
	void resetPartialResult(reduction::Result *resultVar);
  protected: 
	void updateIntermediateResult(reduction::Result *localPartialResult);
};


//---------------------------------------------------------------------------------Primitive for variable 'area'

class ReductionPrimitive_area : public TaskGlobalReductionPrimitive {
  public: 
	ReductionPrimitive_area(int localParticipants);
	void resetPartialResult(reduction::Result *resultVar);
  protected: 
	void updateIntermediateResult(reduction::Result *localPartialResult);
};


//---------------------------------------------------------------------------------Reduction Primitive Instances

static NonTaskGlobalReductionPrimitive *placement_resultReducer[Space_B_Threads_Per_Segment];
static TaskGlobalReductionPrimitive *areaReducer[Space_A_Threads_Per_Segment];

//-------------------------------------------------------------------------------Reduction Primitive Initializer

void setupReductionPrimitives(std::ofstream &logFile);

//--------------------------------------------------------------------------------Reduction Primitives Retriever

Hashtable<void*> *getReductionPrimitiveMap(ThreadIds *threadIds);

/*--------------------------------------------------------------------------------------------------------------
Task Environment Management Structures and Functions
--------------------------------------------------------------------------------------------------------------*/

//-------------------------------------------------------------------------Task Environment Implementation Class

class TaskEnvironmentImpl : public TaskEnvironment {
  public:
	double area;
  public:
	TaskEnvironmentImpl();
	void prepareItemsMap();
	void setDefaultTaskCompletionInstrs();
};

//--------------------------------------------------------------------------Environmental Links Object Generator

EnvironmentLinks initiateEnvLinks(TaskEnvironment *environment);

//---------------------------------------------------------------------------------LPS Allocation Preconfigurers

void preconfigureLpsAllocationsInEnv(TaskEnvironment *environment, 
		ArrayMetadata *metadata, 
		Hashtable<DataPartitionConfig*> *partConfigMap);

//------------------------------------------------------------------------------------Non-array variables copier

void copyBackNonArrayEnvVariables(TaskEnvironment *environment, TaskGlobals *taskGlobals);

/*--------------------------------------------------------------------------------------------------------------
task executor function
--------------------------------------------------------------------------------------------------------------*/

void execute(TaskEnvironment *environment, 
		double precision_threshold, 
		int max_rounds, 
		int cell_length, 
		int grid_dim, 
		int points_per_cell, 
		MCAEPartition partition, 
		int segmentId, 
		std::ofstream &logFile);

/*--------------------------------------------------------------------------------------------------------------
function for the initialize block
--------------------------------------------------------------------------------------------------------------*/

void initializeTask(ArrayMetadata *arrayMetadata, 
		EnvironmentLinks environmentLinks, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		MCAEPartition partition, 
		double precision_threshold, 
		int max_rounds, 
		int cell_length, 
		int grid_dim, 
		int points_per_cell);


/*--------------------------------------------------------------------------------------------------------------
functions for compute stages
--------------------------------------------------------------------------------------------------------------*/

int setupgridcells_stage_4(SpaceC_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		MCAEPartition partition, 
		std::ofstream &logFile);


int initializeestimatediffs_stage_5(SpaceB_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		MCAEPartition partition, 
		std::ofstream &logFile);


int performsampling_stage_9(SpaceC_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		Hashtable<reduction::Result*> *localReductionResultMap, 
		MCAEPartition partition, 
		std::ofstream &logFile);


int estimatesubarea_stage_10(SpaceB_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		MCAEPartition partition, 
		std::ofstream &logFile);


int estimatetotalarea_stage_11(SpaceB_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		Hashtable<reduction::Result*> *localReductionResultMap, 
		MCAEPartition partition, 
		std::ofstream &logFile);


int displayresult_stage_12(SpaceA_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		MCAEPartition partition, 
		std::ofstream &logFile);


/*--------------------------------------------------------------------------------------------------------------
run method for thread simulating the task flow
--------------------------------------------------------------------------------------------------------------*/

void run(ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		MCAEPartition partition, ThreadStateImpl *threadState);

/*-----------------------------------------------------------------------------------
Data structure and function for Pthreads
------------------------------------------------------------------------------------*/

class PThreadArg {
  public:
	const char *taskName;
	ArrayMetadata *metadata;
	TaskGlobals *taskGlobals;
	ThreadLocals *threadLocals;
	MCAEPartition partition;
	ThreadStateImpl *threadState;
};

void *runPThreads(void *argument);

}
#endif
