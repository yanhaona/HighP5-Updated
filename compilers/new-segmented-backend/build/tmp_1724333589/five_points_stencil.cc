
/*--------------------------------------------------------------------------------------------------------------
header file for the task
--------------------------------------------------------------------------------------------------------------*/
#include "five_points_stencil.h"


/*--------------------------------------------------------------------------------------------------------------
header files for different purposes
--------------------------------------------------------------------------------------------------------------*/
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


using namespace fps;


/*--------------------------------------------------------------------------------------------------------------
Functions for ArrayMetadata and EnvironmentLinks
--------------------------------------------------------------------------------------------------------------*/

fps::ArrayMetadata::ArrayMetadata() : Metadata() {
	setTaskName("Five Points Stencil");
}

void fps::ArrayMetadata::print(std::ofstream &stream) {
	stream << "Array Metadata" << std::endl;
	stream << "Array: plate";
	stream << ' ';
	plateDims[0].print(stream);
	stream << ' ';
	plateDims[1].print(stream);
	stream << std::endl;
	stream.flush();
}

/*--------------------------------------------------------------------------------------------------------------
LPU print functions
--------------------------------------------------------------------------------------------------------------*/

void fps::SpaceRoot_LPU::print(std::ofstream &stream, int indentLevel) {
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: plate" << std::endl;
	platePartDims[0].print(stream, indentLevel + 1);
	platePartDims[1].print(stream, indentLevel + 1);
	stream.flush();
}

void fps::SpaceA_LPU::print(std::ofstream &stream, int indentLevel) {
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: plate" << std::endl;
	platePartDims[0].print(stream, indentLevel + 1);
	platePartDims[1].print(stream, indentLevel + 1);
	stream.flush();
}

void fps::SpaceB_LPU::print(std::ofstream &stream, int indentLevel) {
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: plate" << std::endl;
	platePartDims[0].print(stream, indentLevel + 1);
	platePartDims[1].print(stream, indentLevel + 1);
	stream.flush();
}


/*--------------------------------------------------------------------------------------------------------------
functions for generating partition configuration objects for data structures
--------------------------------------------------------------------------------------------------------------*/

//-------------------------------------------------------------------------------------------------------Space A

DataPartitionConfig *fps::getplateConfigForSpaceA(ArrayMetadata *metadata, 
		FPSPartition partition, 
		int ppuCount, 
		Hashtable<DataPartitionConfig*> *configMap) {
	List<DimPartitionConfig*> *dimensionConfigs = new List<DimPartitionConfig*>;
	Assert(dimensionConfigs != NULL);
	int *dim0Paddings = new int[2];
	dim0Paddings[0] = partition.p1;
	dim0Paddings[1] = partition.p1;
	int *dim0Arguments = new int;
	dim0Arguments[0] = partition.k;
	dimensionConfigs->Append(new BlockCountConfig(metadata->plateDims[0], 
			dim0Arguments, dim0Paddings, ppuCount, 0));
	delete[] dim0Paddings;
	int *dim1Paddings = new int[2];
	dim1Paddings[0] = partition.p1;
	dim1Paddings[1] = partition.p1;
	int *dim1Arguments = new int;
	dim1Arguments[0] = partition.l;
	dimensionConfigs->Append(new BlockCountConfig(metadata->plateDims[1], 
			dim1Arguments, dim1Paddings, ppuCount, 1));
	delete[] dim1Paddings;
	DataPartitionConfig *config =  new DataPartitionConfig(2, dimensionConfigs);
	config->configureDimensionOrder();
	config->setLpsId(Space_A);
	return config;
}

//-------------------------------------------------------------------------------------------------------Space B

DataPartitionConfig *fps::getplateConfigForSpaceB(ArrayMetadata *metadata, 
		FPSPartition partition, 
		int ppuCount, 
		Hashtable<DataPartitionConfig*> *configMap) {
	List<DimPartitionConfig*> *dimensionConfigs = new List<DimPartitionConfig*>;
	Assert(dimensionConfigs != NULL);
	int *dim0Paddings = new int[2];
	dim0Paddings[0] = partition.p2;
	dim0Paddings[1] = partition.p2;
	int *dim0Arguments = new int;
	dim0Arguments[0] = partition.m;
	dimensionConfigs->Append(new BlockCountConfig(metadata->plateDims[0], 
			dim0Arguments, dim0Paddings, ppuCount, 0));
	delete[] dim0Paddings;
	int *dim1Paddings = new int[2];
	dim1Paddings[0] = partition.p2;
	dim1Paddings[1] = partition.p2;
	int *dim1Arguments = new int;
	dim1Arguments[0] = partition.n;
	dimensionConfigs->Append(new BlockCountConfig(metadata->plateDims[1], 
			dim1Arguments, dim1Paddings, ppuCount, 1));
	delete[] dim1Paddings;
	DataPartitionConfig *config =  new DataPartitionConfig(2, dimensionConfigs);
	config->setParent(configMap->Lookup("plateSpaceAConfig"), 1);
	config->configureDimensionOrder();
	config->setLpsId(Space_B);
	return config;
}

//---------------------------------------------------------------------------Partition Configuration Accumulator

Hashtable<DataPartitionConfig*> *fps::getDataPartitionConfigMap(ArrayMetadata *metadata, 
		FPSPartition partition, int *ppuCounts) {
	Hashtable<DataPartitionConfig*> *configMap = new Hashtable<DataPartitionConfig*>;
	DataPartitionConfig *plateSpaceAConfig = getplateConfigForSpaceA(metadata, 
		partition, ppuCounts[Space_A], configMap);
	configMap->Enter("plateSpaceAConfig", plateSpaceAConfig);
	DataPartitionConfig *plateSpaceBConfig = getplateConfigForSpaceB(metadata, 
		partition, ppuCounts[Space_B], configMap);
	configMap->Enter("plateSpaceBConfig", plateSpaceBConfig);
	return configMap;
}

/*--------------------------------------------------------------------------------------------------------------
functions for generating memory blocks for data parts of various LPUs
--------------------------------------------------------------------------------------------------------------*/

LpsContent *fps::genSpaceBContent(List<ThreadState*> *threads, ArrayMetadata *metadata, 
		TaskEnvironment *environment, 
		FPSPartition partition, 
		Hashtable<DataPartitionConfig*> *partConfigMap) {

	if(threads->NumElements() == 0) return NULL;
	LpsContent *spaceBContent = new LpsContent(Space_B);

	DataItems *plate = new DataItems("plate", 2, false);
	DataPartitionConfig *plateConfig = partConfigMap->Lookup("plateSpaceBConfig");
	plate->setPartitionConfig(plateConfig);
	DataPartsList *plateParts = plateConfig->generatePartList(2);
	plate->setPartsList(plateParts);
	std::vector<DimConfig> plateDimOrder = *(plateConfig->getDimensionOrder());
	PartIdContainer *plateContainer = NULL;
	if (plateDimOrder.size() == 1) plateContainer = new PartContainer(plateDimOrder[0]);
	else plateContainer = new PartListContainer(plateDimOrder[0]);
	Assert(plateContainer != NULL);
	List<int*> *platePartId = plateConfig->generatePartIdTemplate();
	spaceBContent->addDataItems("plate", plate);

	for (int i = 0; i < threads->NumElements(); i++) {
		ThreadState *thread = threads->Nth(i);
		int lpuId = INVALID_ID;
		while((lpuId = thread->getNextLpuId(Space_B, Space_Root, lpuId)) != INVALID_ID) {
			List<int*> *lpuIdChain = thread->getLpuIdChainWithoutCopy(
					Space_B, Space_Root);
			plateConfig->generatePartId(lpuIdChain, platePartId);
			plateContainer->insertPartId(platePartId, 2, plateDimOrder);
		}
	}

	plateParts->initializePartsList(plateConfig, plateContainer, sizeof(double));
	TaskItem *plateItem = environment->getItem("plate");
	LpsAllocation *plateAlloc = plateItem->getLpsAllocation("B");
	plateAlloc->setPartContainerTree(plateContainer);
	plateAlloc->setPartsList(new PartsList(plateParts->getPartList()));

	return spaceBContent;
}

//-----------------------------------------------------------------------------------------Task Data Initializer

TaskData *fps::initializeTaskData(List<ThreadState*> *threads, ArrayMetadata *metadata, 
		TaskEnvironment *environment, 
		SegmentState *segment, 
		FPSPartition partition, int *ppuCounts) {

	Hashtable<DataPartitionConfig*> *configMap = 
			getDataPartitionConfigMap(metadata, partition, ppuCounts);
	TaskData *taskData = new TaskData();
	Assert(taskData != NULL);
	preconfigureLpsAllocationsInEnv(environment, metadata, configMap);

	// prepare LPS contents map
	LpsContent *spaceBContent = genSpaceBContent(threads, metadata, 
			environment, partition, configMap);
	taskData->addLpsContent("B", spaceBContent);

	// prepare file I/O handlers in the environment
	Hashtable<PartReader*> *readersMap = generateReadersMap(taskData, segment, configMap);
	environment->setReadersMap(readersMap);
	Hashtable<PartWriter*> *writersMap = generateWritersMap(taskData, segment, configMap);
	environment->setWritersMap(writersMap);

	// initialize parts lists of environmental variables
	environment->preprocessProgramEnvForItems();
	environment->setupItemsPartsLists();
	environment->postprocessProgramEnvForItems();
	return taskData;
}

/*--------------------------------------------------------------------------------------------------------------
file I/O for environmental data structures
--------------------------------------------------------------------------------------------------------------*/

//------------------------------------------------------------------------------------Parts-reader map generator
Hashtable<PartReader*> *fps::generateReadersMap(TaskData *taskData, 
			SegmentState *segment, 
			Hashtable<DataPartitionConfig*> *partConfigMap) {

	Hashtable<PartReader*> *readersMap = new Hashtable<PartReader*>;

	// readers for Space B
	DataItems *plateInSpaceB = taskData->getDataItemsOfLps("B", "plate");
	if (plateInSpaceB != NULL && !plateInSpaceB->isEmpty()) {
		DataPartitionConfig *config = partConfigMap->Lookup("plateSpaceBConfig");
		plateInSpaceBReader *reader = new plateInSpaceBReader(config, plateInSpaceB->getPartsList());
		Assert(reader != NULL);
		readersMap->Enter("plateInSpaceBReader", reader);
	} else {
		readersMap->Enter("plateInSpaceBReader", NULL);
	}

	return readersMap;
}

//------------------------------------------------------------------------------------Parts-writer map generator
Hashtable<PartWriter*> *fps::generateWritersMap(TaskData *taskData, 
			SegmentState *segment, 
			Hashtable<DataPartitionConfig*> *partConfigMap) {

	Hashtable<PartWriter*> *writersMap = new Hashtable<PartWriter*>;

	int writerId;
	MPI_Comm_rank(MPI_COMM_WORLD, &writerId);
	int mpiProcessCount;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiProcessCount);
	int writersCount = min(mpiProcessCount, Max_Segments_Count);

	// writers for Space B
	DataItems *plateInSpaceB = taskData->getDataItemsOfLps("B", "plate");
	if (plateInSpaceB != NULL && !plateInSpaceB->isEmpty() && segment->computeStagesInLps(Space_B)) {
		DataPartitionConfig *config = partConfigMap->Lookup("plateSpaceBConfig");
		plateInSpaceBWriter *writer = new plateInSpaceBWriter(config, 
				writerId, plateInSpaceB->getPartsList());
		Assert(writer != NULL);
		writer->setWritersCount(writersCount);
		writersMap->Enter("plateInSpaceBWriter", writer);
	} else {
		writersMap->Enter("plateInSpaceBWriter", NULL);
	}

	return writersMap;
}

/*--------------------------------------------------------------------------------------------------------------
functions for retrieving partition counts in different LPSes
--------------------------------------------------------------------------------------------------------------*/

int *fps::getLPUsCountOfSpaceA(Hashtable<DataPartitionConfig*> *partConfigMap) {
	int *count = new int[2];
	DataPartitionConfig *plateConfig = partConfigMap->Lookup("plateSpaceAConfig");
	count[0] = plateConfig->getPartsCountAlongDimension(0);
	count[1] = plateConfig->getPartsCountAlongDimension(1);
	return count;
}

int *fps::getLPUsCountOfSpaceB(Hashtable<DataPartitionConfig*> *partConfigMap, Dimension plateDim1, Dimension plateDim2) {
	int *count = new int[2];
	DataPartitionConfig *plateConfig = partConfigMap->Lookup("plateSpaceBConfig");
	count[0] = plateConfig->getPartsCountAlongDimension(0, &plateDim1);
	count[1] = plateConfig->getPartsCountAlongDimension(1, &plateDim2);
	return count;
}

/*--------------------------------------------------------------------------------------------------------------
functions for generating LPUs given LPU Ids
--------------------------------------------------------------------------------------------------------------*/

void fps::generateSpaceALpu(ThreadState *threadState, 
			Hashtable<DataPartitionConfig*> *partConfigMap, 
			TaskData *taskData) {

	int *lpuId = threadState->getCurrentLpuId(Space_A);
	int *lpuCounts = threadState->getLpuCounts(Space_A);
	SpaceA_LPU *lpu = (SpaceA_LPU*) threadState->getCurrentLpu(Space_A, true);
	List<int*> *lpuIdChain = threadState->getLpuIdChainWithoutCopy(Space_A, Space_Root);

	lpu->lpuId[0] = lpuId[0];
	lpu->lpuId[1] = lpuId[1];

	DataPartitionConfig *plateConfig = partConfigMap->Lookup("plateSpaceAConfig");
	SpaceRoot_LPU *spaceRootLpu = (SpaceRoot_LPU*) threadState->getCurrentLpu(Space_Root);
	PartDimension *plateParentPartDims = spaceRootLpu->platePartDims;
	plateConfig->updatePartDimensionInfo(lpuId, lpuCounts, lpu->platePartDims, plateParentPartDims);
}

void fps::generateSpaceBLpu(ThreadState *threadState, 
			Hashtable<DataPartitionConfig*> *partConfigMap, 
			TaskData *taskData) {

	int *lpuId = threadState->getCurrentLpuId(Space_B);
	int *lpuCounts = threadState->getLpuCounts(Space_B);
	SpaceB_LPU *lpu = (SpaceB_LPU*) threadState->getCurrentLpu(Space_B, true);
	List<int*> *lpuIdChain = threadState->getLpuIdChainWithoutCopy(Space_B, Space_Root);

	lpu->lpuId[0] = lpuId[0];
	lpu->lpuId[1] = lpuId[1];

	DataPartitionConfig *plateConfig = partConfigMap->Lookup("plateSpaceBConfig");
	SpaceA_LPU *spaceALpu = (SpaceA_LPU*) threadState->getCurrentLpu(Space_A);
	PartDimension *plateParentPartDims = spaceALpu->platePartDims;
	plateConfig->updatePartDimensionInfo(lpuId, lpuCounts, lpu->platePartDims, plateParentPartDims);

	if (taskData != NULL) {
		PartIterator *iterator = threadState->getIterator(Space_B, "plate");
		List<int*> *partId = iterator->getPartIdTemplate();
		plateConfig->generatePartId(lpuIdChain, partId);
		DataItems *plateItems = taskData->getDataItemsOfLps("B", "plate");
		DataPart *platePart = plateItems->getDataPart(partId, iterator);
		platePart->getMetadata()->updateStorageDimension(lpu->platePartDims);
		lpu->plate = (double*) platePart->getData();
		lpu->plate_lag_1 = (double*) platePart->getData(1);
	}
}

/*--------------------------------------------------------------------------------------------------------------
functions to generate PPU IDs and PPU group IDs for a thread
--------------------------------------------------------------------------------------------------------------*/

ThreadIds *fps::getPpuIdsForThread(int threadNo)  {

	ThreadIds *threadIds = new ThreadIds;
	threadIds->threadNo = threadNo;
	threadIds->lpsCount = Space_Count;
	threadIds->ppuIds = new PPU_Ids[Space_Count];
	int idsArray[Space_Count];
	idsArray[Space_Root] = threadNo;

	// for Space Root
	threadIds->ppuIds[Space_Root].lpsName = "Root";
	threadIds->ppuIds[Space_Root].groupId = 0;
	threadIds->ppuIds[Space_Root].groupSize = Total_Threads;
	threadIds->ppuIds[Space_Root].ppuCount = 1;
	threadIds->ppuIds[Space_Root].id = (threadNo == 0) ? 0 : INVALID_ID;

	int threadCount;
	int groupSize;
	int groupThreadId;

	// for Space A;
	threadIds->ppuIds[Space_A].lpsName = "A";
	threadCount = Max_Total_Threads;
	groupSize = threadCount / 10;
	groupThreadId = idsArray[Space_Root] % groupSize;
	threadIds->ppuIds[Space_A].groupId = idsArray[Space_Root] / groupSize;
	threadIds->ppuIds[Space_A].ppuCount = 10;
	threadIds->ppuIds[Space_A].groupSize = groupSize;
	if (groupThreadId == 0) threadIds->ppuIds[Space_A].id
			= threadIds->ppuIds[Space_A].groupId;
	else threadIds->ppuIds[Space_A].id = INVALID_ID;
	idsArray[Space_A] = groupThreadId;

	// for Space B;
	threadIds->ppuIds[Space_B].lpsName = "B";
	threadCount = threadIds->ppuIds[Space_A].groupSize;
	groupSize = threadCount / 4;
	groupThreadId = idsArray[Space_A] % groupSize;
	threadIds->ppuIds[Space_B].groupId = idsArray[Space_A] / groupSize;
	threadIds->ppuIds[Space_B].ppuCount = 4;
	threadIds->ppuIds[Space_B].groupSize = groupSize;
	if (groupThreadId == 0) threadIds->ppuIds[Space_B].id
			= threadIds->ppuIds[Space_B].groupId;
	else threadIds->ppuIds[Space_B].id = INVALID_ID;
	idsArray[Space_B] = groupThreadId;

	return threadIds;
}


void fps::adjustPpuCountsAndGroupSizes(ThreadIds *threadId)  {

	int groupBegin = 0;
	int groupEnd = Total_Threads - 1;

	int groupId = 0;
	int groupSize = Total_Threads;
	int ppuCount = 1;

	groupId = threadId->ppuIds[Space_A].groupId;
	groupSize = threadId->ppuIds[Space_A].groupSize;
	ppuCount = ((groupEnd - groupBegin + 1) + (groupSize - 1)) / groupSize;
	threadId->ppuIds[Space_A].ppuCount = ppuCount;
	groupBegin = groupId * groupSize;
	groupEnd = min(groupBegin + groupSize - 1, groupEnd);
	threadId->ppuIds[Space_A].groupSize = groupEnd - groupBegin + 1;

	groupId = threadId->ppuIds[Space_B].groupId;
	groupSize = threadId->ppuIds[Space_B].groupSize;
	ppuCount = ((groupEnd - groupBegin + 1) + (groupSize - 1)) / groupSize;
	threadId->ppuIds[Space_B].ppuCount = ppuCount;
	groupBegin = groupId * groupSize;
	groupEnd = min(groupBegin + groupSize - 1, groupEnd);
	threadId->ppuIds[Space_B].groupSize = groupEnd - groupBegin + 1;

}


/*--------------------------------------------------------------------------------------------------------------
Thread-State implementation class for the task
--------------------------------------------------------------------------------------------------------------*/

//-------------------------------------------------------------------------------------------LPS Hierarchy Trace

void ThreadStateImpl::setLpsParentIndexMap() {
	lpsParentIndexMap = new int[Space_Count];
	lpsParentIndexMap[Space_Root] = INVALID_ID;
	lpsParentIndexMap[Space_A] = Space_Root;
	lpsParentIndexMap[Space_B] = Space_A;
}

//-----------------------------------------------------------------------------------------Root LPU Construction

void ThreadStateImpl::setRootLpu(Metadata *metadata) {

	ArrayMetadata *arrayMetadata = (ArrayMetadata*) metadata;

	SpaceRoot_LPU *lpu = new SpaceRoot_LPU;
	lpu->plate = NULL;
	lpu->platePartDims[0] = PartDimension();
	lpu->platePartDims[0].partition = arrayMetadata->plateDims[0];
	lpu->platePartDims[0].storage = arrayMetadata->plateDims[0].getNormalizedDimension();
	lpu->platePartDims[1] = PartDimension();
	lpu->platePartDims[1].partition = arrayMetadata->plateDims[1];
	lpu->platePartDims[1].storage = arrayMetadata->plateDims[1].getNormalizedDimension();

	lpu->setValidBit(true);
	lpsStates[Space_Root]->lpu = lpu;
}

void ThreadStateImpl::setRootLpu(LPU *lpu) {
	lpu->setValidBit(true);
	lpsStates[Space_Root]->lpu = lpu;
}

//------------------------------------------------------------------------------------------State Initialization

void ThreadStateImpl::initializeLPUs() {
	lpsStates[Space_A]->lpu = new SpaceA_LPU;
	lpsStates[Space_A]->lpu->setValidBit(false);
	lpsStates[Space_B]->lpu = new SpaceB_LPU;
	lpsStates[Space_B]->lpu->setValidBit(false);
}

//--------------------------------------------------------------------------------------------LPU Count Function


int *ThreadStateImpl::computeLpuCounts(int lpsId) {
	Hashtable<DataPartitionConfig*> *configMap = getPartConfigMap();
	if (lpsId == Space_Root) {
		return NULL;
	}
	if (lpsId == Space_A) {
		return getLPUsCountOfSpaceA(configMap);
	}
	if (lpsId == Space_B) {
		SpaceA_LPU *spaceALpu
				 = (SpaceA_LPU*) lpsStates[Space_A]->lpu;
		return getLPUsCountOfSpaceB(configMap, 
				spaceALpu->platePartDims[0].partition, 
				spaceALpu->platePartDims[1].partition);
	}
	return NULL;
}

//-------------------------------------------------------------------------------------LPU Construction Function

LPU *ThreadStateImpl::computeNextLpu(int lpsId) {
	Hashtable<DataPartitionConfig*> *partConfigMap = getPartConfigMap();
	TaskData *taskData = getTaskData();
	if (lpsId == Space_A) {
		SpaceA_LPU *currentLpu = (SpaceA_LPU*) lpsStates[Space_A]->lpu;
		generateSpaceALpu(this, partConfigMap, taskData);
		currentLpu->setValidBit(true);
		return currentLpu;
	}
	if (lpsId == Space_B) {
		SpaceB_LPU *currentLpu = (SpaceB_LPU*) lpsStates[Space_B]->lpu;
		generateSpaceBLpu(this, partConfigMap, taskData);
		currentLpu->setValidBit(true);
		return currentLpu;
	}
	return NULL;
}

//-------------------------------------------------------------------------Reduction Result Map Creator Function

void ThreadStateImpl::initializeReductionResultMap() {
	localReductionResultMap = new Hashtable<reduction::Result*>;
}

/*--------------------------------------------------------------------------------------------------------------
Initializer function for global synchronization primitives
--------------------------------------------------------------------------------------------------------------*/

void fps::initializeSyncPrimitives() {

	for (int i = 0; i < Space_A_Threads_Per_Segment; i++) {
		int participants = Space_B_Threads_Per_Segment / Space_A_Threads_Per_Segment;
		plateStage9No1GSyncs[i] = new RS(participants);
		plateStage9No1ReverseSyncs[i] = new Barrier(participants);
	}
	for (int i = 0; i < Space_B_Threads_Per_Segment; i++) {
		int participants = Space_B_Threads_Per_Segment / Space_B_Threads_Per_Segment;
		plateStage8No1GSyncs[i] = new RS(participants);
		plateStage8No1ReverseSyncs[i] = new Barrier(participants);
	}
}

/*--------------------------------------------------------------------------------------------------------------
function for initializing thread's sync primitives
--------------------------------------------------------------------------------------------------------------*/

ThreadSyncPrimitive *fps::getSyncPrimitives(ThreadIds *threadIds) {

	ThreadSyncPrimitive *threadSync = new ThreadSyncPrimitive();

	int spaceAGroup = threadIds->ppuIds[Space_A].groupId % Space_A_Threads_Per_Segment;
	threadSync->plateStage9No1GSync = plateStage9No1GSyncs[spaceAGroup];
	threadSync->plateStage9No1ReverseSync = plateStage9No1ReverseSyncs[spaceAGroup];

	int spaceBGroup = threadIds->ppuIds[Space_B].groupId % Space_B_Threads_Per_Segment;
	threadSync->plateStage8No1GSync = plateStage8No1GSyncs[spaceBGroup];
	threadSync->plateStage8No1ReverseSync = plateStage8No1ReverseSyncs[spaceBGroup];

	return threadSync;
}


/*--------------------------------------------------------------------------------------------------------------
functions to generate distribution trees for communicated variables
--------------------------------------------------------------------------------------------------------------*/

//---------------------------------------------------------------------------------------------------------plate

Container *fps::generateDistributionTreeFor_plate(List<SegmentState*> *segmentList, 
		Hashtable<DataPartitionConfig*> *configMap) {

	BranchingContainer *rootContainer = new BranchingContainer(0, LpsDimConfig());
	Assert(rootContainer != NULL);
	for (int i = 0; i < segmentList->NumElements(); i++) {
		SegmentState *segment = segmentList->Nth(i);
		int segmentTag = segment->getPhysicalId();
		List<ThreadState*> *threadList = segment->getParticipantList();
		for (int j = 0; j < threadList->NumElements(); j++) {
			ThreadState *thread = threadList->Nth(j);
			int lpuId = INVALID_ID;
			DataPartitionConfig *partConfig = NULL;
			List<int*> *partId = NULL;
			DataItemConfig *dataItemConfig = NULL;
			std::vector<LpsDimConfig> *dimOrder = NULL;

			//generating parts for: B
			if (thread->isValidPpu(Space_B)) {
				partConfig = configMap->Lookup("plateSpaceBConfig");
				partId = partConfig->generatePartIdTemplate();
				dataItemConfig = partConfig->generateStateFulVersion();
				dimOrder = dataItemConfig->generateDimOrderVector();

				while((lpuId = thread->getNextLpuId(Space_B, 
						Space_Root, lpuId)) != INVALID_ID) {
					List<int*> *lpuIdChain = thread->getLpuIdChainWithoutCopy(
							Space_B, Space_Root);
					partConfig->generatePartId(lpuIdChain, partId);
					rootContainer->insertPart(*dimOrder, segmentTag, partId);
				}
			}
		}
	}
	return rootContainer;
}

//----------------------------------------------------------------------------------------------distribution map

PartDistributionMap *fps::generateDistributionMap(List<SegmentState*> *segmentList, 
		Hashtable<DataPartitionConfig*> *configMap) {

	PartDistributionMap *distributionMap = new PartDistributionMap();
	Assert(distributionMap != NULL);
	distributionMap->setDistributionForVariable("plate", 
			generateDistributionTreeFor_plate(segmentList, configMap));
	return distributionMap;
}

/*--------------------------------------------------------------------------------------------------------------
functions to generate communication confinement configurations
--------------------------------------------------------------------------------------------------------------*/

//------------------------------------------------------------------------------------------------plateStage8No1

ConfinementConstructionConfig *fps::getConfineConstrConfigFor_plateStage8No1(TaskData *taskData, 
		Hashtable<DataPartitionConfig*> *partConfigMap, 
		int localSegmentTag, 
		PartDistributionMap *distributionMap) {

	int senderLps = Space_B;
	DataItemConfig *senderDataConfig = partConfigMap->Lookup(
			"plateSpaceBConfig")->generateStateFulVersion();
	int receiverLps = Space_B;
	DataItemConfig *receiverDataConfig = partConfigMap->Lookup(
			"plateSpaceBConfig")->generateStateFulVersion();
	int confinementLps = Space_A;
	PartIdContainer *senderPartTree = NULL;
	DataItems *senderDataItems = taskData->getDataItemsOfLps("B", "plate");
	if(senderDataItems != NULL) senderPartTree = senderDataItems->getPartIdContainer();
	PartIdContainer *receiverPartTree = NULL;
	DataItems *receiverDataItems = taskData->getDataItemsOfLps("B", "plate");
	if(receiverDataItems != NULL) receiverPartTree = receiverDataItems->getPartIdContainer();
	BranchingContainer *distributionTree = 
			(BranchingContainer*) distributionMap->getDistrubutionTree("plate");

	ConfinementConstructionConfig *confinementConfig = new ConfinementConstructionConfig(localSegmentTag, 
			senderLps, senderDataConfig, 
			receiverLps, receiverDataConfig, 
			confinementLps, 
			senderPartTree, receiverPartTree, 
			distributionTree);
	Assert(confinementConfig != NULL);
	return confinementConfig;
}

//------------------------------------------------------------------------------------------------plateStage9No1

ConfinementConstructionConfig *fps::getConfineConstrConfigFor_plateStage9No1(TaskData *taskData, 
		Hashtable<DataPartitionConfig*> *partConfigMap, 
		int localSegmentTag, 
		PartDistributionMap *distributionMap) {

	int senderLps = Space_A;
	DataItemConfig *senderDataConfig = partConfigMap->Lookup(
			"plateSpaceBConfig")->generateStateFulVersion();
	int receiverLps = Space_A;
	DataItemConfig *receiverDataConfig = partConfigMap->Lookup(
			"plateSpaceBConfig")->generateStateFulVersion();
	int confinementLps = Space_Root;
	PartIdContainer *senderPartTree = NULL;
	DataItems *senderDataItems = taskData->getDataItemsOfLps("B", "plate");
	if(senderDataItems != NULL) senderPartTree = senderDataItems->getPartIdContainer();
	PartIdContainer *receiverPartTree = NULL;
	DataItems *receiverDataItems = taskData->getDataItemsOfLps("B", "plate");
	if(receiverDataItems != NULL) receiverPartTree = receiverDataItems->getPartIdContainer();
	BranchingContainer *distributionTree = 
			(BranchingContainer*) distributionMap->getDistrubutionTree("plate");

	ConfinementConstructionConfig *confinementConfig = new ConfinementConstructionConfig(localSegmentTag, 
			senderLps, senderDataConfig, 
			receiverLps, receiverDataConfig, 
			confinementLps, 
			senderPartTree, receiverPartTree, 
			distributionTree);
	Assert(confinementConfig != NULL);
	return confinementConfig;
}

/*--------------------------------------------------------------------------------------------------------------
functions to generate data exchange lists for data dependencies
--------------------------------------------------------------------------------------------------------------*/

//------------------------------------------------------------------------------------------------plateStage8No1

List<DataExchange*> *fps::getDataExchangeListFor_plateStage8No1(TaskData *taskData, 
		Hashtable<DataPartitionConfig*> *partConfigMap, 
		int localSegmentTag, 
		PartDistributionMap *distributionMap) {

	ConfinementConstructionConfig *ccConfig = getConfineConstrConfigFor_plateStage8No1(taskData, 
			partConfigMap, localSegmentTag, distributionMap);
	List<Confinement*> *confinementList = 
			Confinement::generateAllConfinements(ccConfig, Space_Root);
	if (confinementList == NULL || confinementList->NumElements() == 0) return NULL;

	List<DataExchange*> *dataExchangeList = new List<DataExchange*>;
	Assert(dataExchangeList != NULL);
	for (int i = 0; i < confinementList->NumElements(); i++) {
		Confinement *confinement = confinementList->Nth(i);
		List<DataExchange*> *confinementExchanges = confinement->getAllDataExchanges();
		if (confinementExchanges != NULL) {
			dataExchangeList->AppendAll(confinementExchanges);
			delete confinementExchanges;
		}
	}
	if (dataExchangeList->NumElements() == 0) {
		delete dataExchangeList;
		return NULL;
	}
	return dataExchangeList;
}

//------------------------------------------------------------------------------------------------plateStage9No1

List<DataExchange*> *fps::getDataExchangeListFor_plateStage9No1(TaskData *taskData, 
		Hashtable<DataPartitionConfig*> *partConfigMap, 
		int localSegmentTag, 
		PartDistributionMap *distributionMap) {

	ConfinementConstructionConfig *ccConfig = getConfineConstrConfigFor_plateStage9No1(taskData, 
			partConfigMap, localSegmentTag, distributionMap);
	List<Confinement*> *confinementList = 
			Confinement::generateAllConfinements(ccConfig, Space_Root);
	if (confinementList == NULL || confinementList->NumElements() == 0) return NULL;

	List<DataExchange*> *dataExchangeList = new List<DataExchange*>;
	Assert(dataExchangeList != NULL);
	for (int i = 0; i < confinementList->NumElements(); i++) {
		Confinement *confinement = confinementList->Nth(i);
		List<DataExchange*> *confinementExchanges = confinement->getAllDataExchanges();
		if (confinementExchanges != NULL) {
			dataExchangeList->AppendAll(confinementExchanges);
			delete confinementExchanges;
		}
	}
	if (dataExchangeList->NumElements() == 0) {
		delete dataExchangeList;
		return NULL;
	}
	return dataExchangeList;
}

/*--------------------------------------------------------------------------------------------------------------
functions to generate communicators for exchanging data to resolve dependencies
--------------------------------------------------------------------------------------------------------------*/

//------------------------------------------------------------------------------------------------plateStage8No1

Communicator *fps::getCommunicatorFor_plateStage8No1(SegmentState *localSegment, 
		TaskData *taskData, 
		Hashtable<DataPartitionConfig*> *partConfigMap, 
		CommStatistics *commStat, 
		PartDistributionMap *distributionMap) {

	commStat->enlistDependency("plateStage8No1");
	int localSegmentTag = localSegment->getPhysicalId();
	ConfinementConstructionConfig *ccConfig = getConfineConstrConfigFor_plateStage8No1(
			taskData, 
			partConfigMap, 
			localSegmentTag, 
			distributionMap);
	ccConfig->configurePaddingInPartitionConfigsForReadWrite();
	struct timeval start;
	gettimeofday(&start, NULL);
	List<DataExchange*> *dataExchangeList = getDataExchangeListFor_plateStage8No1(taskData, 
			partConfigMap, 
			localSegmentTag, 
			distributionMap);
	if (dataExchangeList == NULL || dataExchangeList->NumElements() == 0) return NULL;

	struct timeval middle;
	gettimeofday(&middle, NULL);
	commStat->addConfinementConstrTime("plateStage8No1", start, middle);
	DataPartsList *senderDataParts = NULL;
	DataItems *senderDataItems = taskData->getDataItemsOfLps("B", "plate");
	if(senderDataItems != NULL) senderDataParts = senderDataItems->getPartsList();
	DataPartsList *receiverDataParts = NULL;
	DataItems *receiverDataItems = taskData->getDataItemsOfLps("B", "plate");
	if(receiverDataItems != NULL) receiverDataParts = receiverDataItems->getPartsList();
	SyncConfig *syncConfig = new SyncConfig(ccConfig, 
			senderDataParts, receiverDataParts, sizeof(double));
	Assert(syncConfig != NULL);

	int localSenderPpus = localSegment->getPpuCountForLps(Space_B);
	int localReceiverPpus = localSegment->getPpuCountForLps(Space_B);

	List<CommBuffer*> *bufferList = new List<CommBuffer*>;
	Assert(bufferList != NULL);
	std::vector<int> *participantTags = new std::vector<int>;
	Assert(participantTags != NULL);
	for (int i = 0; i < dataExchangeList->NumElements(); i++) {

		DataExchange *exchange = dataExchangeList->Nth(i);
		std::vector<int> senderTags = exchange->getSender()->getSegmentTags();
		binsearch::addThoseNotExist(participantTags, &senderTags);
		std::vector<int> receiverTags = exchange->getReceiver()->getSegmentTags();
		binsearch::addThoseNotExist(participantTags, &receiverTags);
		if (!exchange->involvesLocalSegment(localSegmentTag)) continue;

		CommBuffer *buffer = NULL;
		if (exchange->isIntraSegmentExchange(localSegmentTag)) {
			buffer = new SwiftIndexMappedVirtualCommBuffer(exchange, syncConfig);
		} else {
			buffer = new SwiftIndexMappedPhysicalCommBuffer(exchange, syncConfig);
		}
		Assert(buffer != NULL);
		bufferList->Append(buffer);
	}
	if (bufferList->NumElements() == 0) return NULL;
	struct timeval end;
	gettimeofday(&end, NULL);
	commStat->addBufferSetupTime("plateStage8No1", middle, end);

	Communicator *communicator = NULL;
	communicator = new GhostRegionSyncCommunicator(localSegmentTag, 
			"plateStage8No1", localSenderPpus, localReceiverPpus, bufferList);
	Assert(communicator != NULL);
	communicator->setParticipants(participantTags);
	communicator->setCommStat(commStat);
	return communicator;
}

//------------------------------------------------------------------------------------------------plateStage9No1

Communicator *fps::getCommunicatorFor_plateStage9No1(SegmentState *localSegment, 
		TaskData *taskData, 
		Hashtable<DataPartitionConfig*> *partConfigMap, 
		CommStatistics *commStat, 
		PartDistributionMap *distributionMap) {

	commStat->enlistDependency("plateStage9No1");
	int localSegmentTag = localSegment->getPhysicalId();
	ConfinementConstructionConfig *ccConfig = getConfineConstrConfigFor_plateStage9No1(
			taskData, 
			partConfigMap, 
			localSegmentTag, 
			distributionMap);
	ccConfig->configurePaddingInPartitionConfigsForReadWrite();
	struct timeval start;
	gettimeofday(&start, NULL);
	List<DataExchange*> *dataExchangeList = getDataExchangeListFor_plateStage9No1(taskData, 
			partConfigMap, 
			localSegmentTag, 
			distributionMap);
	if (dataExchangeList == NULL || dataExchangeList->NumElements() == 0) return NULL;

	struct timeval middle;
	gettimeofday(&middle, NULL);
	commStat->addConfinementConstrTime("plateStage9No1", start, middle);
	DataPartsList *senderDataParts = NULL;
	DataItems *senderDataItems = taskData->getDataItemsOfLps("B", "plate");
	if(senderDataItems != NULL) senderDataParts = senderDataItems->getPartsList();
	DataPartsList *receiverDataParts = NULL;
	DataItems *receiverDataItems = taskData->getDataItemsOfLps("B", "plate");
	if(receiverDataItems != NULL) receiverDataParts = receiverDataItems->getPartsList();
	SyncConfig *syncConfig = new SyncConfig(ccConfig, 
			senderDataParts, receiverDataParts, sizeof(double));
	Assert(syncConfig != NULL);

	int localSenderPpus = localSegment->getPpuCountForLps(Space_A);
	int localReceiverPpus = localSegment->getPpuCountForLps(Space_B);

	List<CommBuffer*> *bufferList = new List<CommBuffer*>;
	Assert(bufferList != NULL);
	std::vector<int> *participantTags = new std::vector<int>;
	Assert(participantTags != NULL);
	for (int i = 0; i < dataExchangeList->NumElements(); i++) {

		DataExchange *exchange = dataExchangeList->Nth(i);
		std::vector<int> senderTags = exchange->getSender()->getSegmentTags();
		binsearch::addThoseNotExist(participantTags, &senderTags);
		std::vector<int> receiverTags = exchange->getReceiver()->getSegmentTags();
		binsearch::addThoseNotExist(participantTags, &receiverTags);
		if (!exchange->involvesLocalSegment(localSegmentTag)) continue;

		CommBuffer *buffer = NULL;
		if (exchange->isIntraSegmentExchange(localSegmentTag)) {
			buffer = new SwiftIndexMappedVirtualCommBuffer(exchange, syncConfig);
		} else {
			buffer = new SwiftIndexMappedPhysicalCommBuffer(exchange, syncConfig);
		}
		Assert(buffer != NULL);
		bufferList->Append(buffer);
	}
	if (bufferList->NumElements() == 0) return NULL;
	struct timeval end;
	gettimeofday(&end, NULL);
	commStat->addBufferSetupTime("plateStage9No1", middle, end);

	Communicator *communicator = NULL;
	communicator = new GhostRegionSyncCommunicator(localSegmentTag, 
			"plateStage9No1", localSenderPpus, localReceiverPpus, bufferList);
	Assert(communicator != NULL);
	communicator->setParticipants(participantTags);
	communicator->setCommStat(commStat);
	return communicator;
}

//----------------------------------------------------------------------------------------------communicator map

Hashtable<Communicator*> *fps::generateCommunicators(SegmentState *localSegment, 
		List<SegmentState*> *segmentList, 
		TaskData *taskData, 
		TaskGlobals *taskGlobals, 
		Hashtable<DataPartitionConfig*> *partConfigMap, 
		PartDistributionMap *distributionMap, 
		CommStatistics *commStat, 
		std::ofstream &logFile) {

	Hashtable<Communicator*> *communicatorMap = new Hashtable<Communicator*>;
	Assert(communicatorMap != NULL);
	int segmentCount = Total_Threads / Threads_Per_Segment;

	Communicator *communicator0 = getCommunicatorFor_plateStage8No1(localSegment, 
			taskData, partConfigMap, commStat, distributionMap);
	if (communicator0 != NULL) {
		communicator0->setLogFile(&logFile);
		communicator0->setupBufferTags(1, segmentCount);
		communicator0->setupCommunicator(false);
		communicatorMap->Enter("plateStage8No1", communicator0);
	}
	Communicator *communicator1 = getCommunicatorFor_plateStage9No1(localSegment, 
			taskData, partConfigMap, commStat, distributionMap);
	if (communicator1 != NULL) {
		communicator1->setLogFile(&logFile);
		communicator1->setupBufferTags(2, segmentCount);
		communicator1->setupCommunicator(false);
		communicatorMap->Enter("plateStage9No1", communicator1);
	}

	return communicatorMap;
}
//---------------------------------------------------------------------------------communication exclude routine

void fps::excludeFromAllCommunication(int segmentId, std::ofstream &logFile) {

}

/*--------------------------------------------------------------------------------------------------------------
Task Environment Management Structures and Functions
--------------------------------------------------------------------------------------------------------------*/

//---------------------------------------------------------------------------------------------------Constructor

fps::TaskEnvironmentImpl::TaskEnvironmentImpl() : TaskEnvironment() {
	prepareItemsMap();
	resetEnvInstructions();
}

//---------------------------------------------------------------------Task Environment Function Implementations

void fps::TaskEnvironmentImpl::prepareItemsMap() {

	if (envItems->Lookup("plate") == NULL) {
		EnvironmentLinkKey *key0 = new EnvironmentLinkKey("plate", 0);
		TaskItem *item0 = new TaskItem(key0, IN_OUT, 2, sizeof(double));
		item0->setEnvironment(this);
		envItems->Enter("plate", item0);
	}
}

void fps::TaskEnvironmentImpl::setDefaultTaskCompletionInstrs() {

	TaskItem *plateItem = envItems->Lookup("plate");
	ChangeNotifyInstruction *instr0 = new ChangeNotifyInstruction(plateItem);
	addEndEnvInstruction(instr0);
}

//--------------------------------------------------------------------------Environmental Links Object Generator

EnvironmentLinks fps::initiateEnvLinks(TaskEnvironment *environment) {

	EnvironmentLinks links;
	fps::TaskEnvironmentImpl *taskEnv = (fps::TaskEnvironmentImpl *) environment;

	TaskItem *plateItem = environment->getItem("plate");
	links.plateDims[0] = plateItem->getDimension(0);
	links.plateDims[1] = plateItem->getDimension(1);

	return links;
}

//---------------------------------------------------------------------------------LPS Allocation Preconfigurers

void fps::preconfigureLpsAllocationsInEnv(TaskEnvironment *environment, 
		ArrayMetadata *metadata, 
		Hashtable<DataPartitionConfig*> *partConfigMap) {


	TaskItem *plateItem = environment->getItem("plate");
	plateItem->setDimension(0, metadata->plateDims[0]);
	plateItem->setDimension(1, metadata->plateDims[1]);
	DataPartitionConfig *plateSpaceBConfig = partConfigMap->Lookup("plateSpaceBConfig");
	LpsAllocation *plateInSpaceB = plateItem->getLpsAllocation("B");
	if (plateInSpaceB == NULL) {
		plateItem->preConfigureLpsAllocation("B", plateSpaceBConfig);
	} else {
		plateInSpaceB->setPartitionConfig(plateSpaceBConfig);
	}
}

//------------------------------------------------------------------------------------Non-array variables copier

void fps::copyBackNonArrayEnvVariables(TaskEnvironment *environment, TaskGlobals *taskGlobals) {
	fps::TaskEnvironmentImpl *taskEnv = (fps::TaskEnvironmentImpl *) environment;
}

/*--------------------------------------------------------------------------------------------------------------
task executor function
--------------------------------------------------------------------------------------------------------------*/

void fps::execute(TaskEnvironment *environment, 
		int iterations, 
		FPSPartition partition, 
		int segmentId, 
		std::ofstream &logFile) {

	environment->setLogFile(&logFile);
	if (segmentId >= Max_Segments_Count) {
		logFile << "Current segment does not participate in: Five Points Stencil\n";
		excludeFromAllCommunication(segmentId, logFile);
		return;
	}

	// setting the total-number-of-threads static variable
	int mpiProcessCount;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiProcessCount);
	int activeSegments = min(mpiProcessCount, Max_Segments_Count);
	Total_Threads = activeSegments * Threads_Per_Segment;

	// initializing environment-links object
	EnvironmentLinks envLinks = initiateEnvLinks(environment);

	// declaring other task related common variables
	TaskGlobals taskGlobals;
	ThreadLocals threadLocals;
	ArrayMetadata *metadata = new ArrayMetadata;

	// declaring and initiating segment execution timer
	struct timeval start;
	gettimeofday(&start, NULL);

	// copying partitioning parameters into an array
	int *partitionArgs = NULL;
	partitionArgs = new int[6];
	partitionArgs[0] = partition.k;
	partitionArgs[1] = partition.l;
	partitionArgs[2] = partition.m;
	partitionArgs[3] = partition.n;
	partitionArgs[4] = partition.p1;
	partitionArgs[5] = partition.p2;

	// initializing sync primitives
	initializeSyncPrimitives();

	// invoking the initializer function
	initializeTask(metadata, envLinks, &taskGlobals, &threadLocals, partition, iterations);
	metadata->plateDims[0].setLength();
	metadata->plateDims[1].setLength();
	logFile << "\ttask initialization is complete\n";
	logFile.flush();

	// declaring and initializing state variables for threads 
	ThreadLocals *threadLocalsList[Total_Threads];
	for (int i = 0; i < Total_Threads; i++) {
		threadLocalsList[i] = new ThreadLocals;
		*threadLocalsList[i] = threadLocals;
	}
	int lpsDimensions[Space_Count];
	lpsDimensions[Space_Root] = 0;
	lpsDimensions[Space_A] = 2;
	lpsDimensions[Space_B] = 2;
	ThreadIds *threadIdsList[Total_Threads];
	for (int i = 0; i < Total_Threads; i++) {
		threadIdsList[i] = getPpuIdsForThread(i);
		adjustPpuCountsAndGroupSizes(threadIdsList[i]);
	}
	int *ppuCounts = threadIdsList[0]->getAllPpuCounts();
	Hashtable<DataPartitionConfig*> *configMap = 
			getDataPartitionConfigMap(metadata, partition, ppuCounts);
	ThreadStateImpl *threadStateList[Total_Threads];
	for (int i = 0; i < Total_Threads; i++) {
		threadStateList[i] = new ThreadStateImpl(Space_Count, 
				lpsDimensions, partitionArgs, threadIdsList[i]);
		threadStateList[i]->initializeLPUs();
		threadStateList[i]->setLpsParentIndexMap();
		threadStateList[i]->setPartConfigMap(configMap);
	}

	// setting up root LPU reference in each thread's state
	for (int i = 0; i < Total_Threads; i++) {
		threadStateList[i]->setRootLpu(metadata);
	}

	// grouping threads into segments
	List<SegmentState*> *segmentList = new List<SegmentState*>;
	int segmentCount = Total_Threads / Threads_Per_Segment;
	int threadIndex = 0;
	for (int s = 0; s < segmentCount; s++) {
		SegmentState *segment = new SegmentState(s, s);
		int participantCount = 0;
		while (participantCount < Threads_Per_Segment) {
			segment->addParticipant(threadStateList[threadIndex]);
			threadIndex++;
			participantCount++;
		}
		segmentList->Append(segment);
	}
	SegmentState *mySegment = segmentList->Nth(segmentId);
	int participantStart = segmentId * Threads_Per_Segment;
	int participantEnd = participantStart + Threads_Per_Segment - 1;
	logFile << "\tsegment grouping of threads is complete\n";
	logFile.flush();

	// initializing segment memory
	TaskData *taskData = initializeTaskData(mySegment->getParticipantList(), 
			metadata, environment, mySegment, partition, ppuCounts);
	for (int i = participantStart; i <= participantEnd; i++) {
		threadStateList[i]->setTaskData(taskData);
		threadStateList[i]->setPartIteratorMap(taskData->generatePartIteratorMap());
		threadStateList[i]->initiateLogFile("fps");
		threadStateList[i]->enableLogging();
	}
	delete[] ppuCounts;
	logFile << "\tmemory allocation is complete\n";
	logFile.flush();

	// calculating memory and threads preparation time
	struct timeval end;
	gettimeofday(&end, NULL);
	double allocationTime = ((end.tv_sec + end.tv_usec / 1000000.0)
			- (start.tv_sec + start.tv_usec / 1000000.0));
	logFile << "Memory preparation time: " << allocationTime << " Seconds" << std::endl;
	double timeConsumedSoFar = allocationTime;
	logFile.flush();

	// initializing communicators map
	CommStatistics *commStat = new CommStatistics();
	PartDistributionMap *distributionMap = generateDistributionMap(segmentList, configMap);
	Hashtable<Communicator*> *communicatorMap = generateCommunicators(mySegment, 
			segmentList, taskData, &taskGlobals, configMap, distributionMap, 
			commStat, logFile);
	for (int i = participantStart; i <= participantEnd; i++) {
		threadStateList[i]->setCommunicatorMap(communicatorMap);
	}
	logFile << "\tcommunicators have been created\n";
	logFile.flush();

	// calculating communicators setup time
	gettimeofday(&end, NULL);
	double communicatorTime = ((end.tv_sec + end.tv_usec / 1000000.0)
			- (start.tv_sec + start.tv_usec / 1000000.0)) - timeConsumedSoFar;
	logFile << "Communicators setup time: " << communicatorTime << " Seconds" << std::endl;
	timeConsumedSoFar += communicatorTime;
	logFile.flush();


	// starting threads
	logFile << "\tlaunching threads\n";
	logFile.flush();
	pthread_t threads[Total_Threads];
	PThreadArg *threadArgs[Total_Threads];
	for (int i = 0; i < Total_Threads; i++) {
		threadArgs[i] = new PThreadArg;
		threadArgs[i]->taskName = "Five Points Stencil";
		threadArgs[i]->metadata = metadata;
		threadArgs[i]->taskGlobals = &taskGlobals;
		threadArgs[i]->threadLocals = threadLocalsList[i];
		threadArgs[i]->partition = partition;
		threadArgs[i]->threadState = threadStateList[i];
	}
	pthread_attr_t attr;
	cpu_set_t cpus;
	pthread_attr_init(&attr);
	int state;
	for (int i = participantStart; i <= participantEnd; i++) {
		int cpuId = (i * Core_Jump / Threads_Per_Core) % Processors_Per_Phy_Unit;
		int physicalId = Processor_Order[cpuId];
		state = pthread_create(&threads[i], &attr, runPThreads, (void *) threadArgs[i]);
		if (state) {
			std::cout << "Could not start some PThread" << std::endl;
			std::exit(EXIT_FAILURE);
		} else {
			logFile << "\t\tlaunched thread #" << i << "\n";
			logFile.flush();
		}
	}
	for (int i = participantStart; i <= participantEnd; i++) {
		pthread_join(threads[i], NULL);
	}


	// calculating computation time
	gettimeofday(&end, NULL);
	double computationTime = ((end.tv_sec + end.tv_usec / 1000000.0)
			- (start.tv_sec + start.tv_usec / 1000000.0)) - timeConsumedSoFar;
	logFile << "Computation time: " << computationTime << " Seconds" << std::endl;
	timeConsumedSoFar += computationTime;
	logFile.flush();

	double compAndOverheadTime = timeConsumedSoFar - allocationTime;
	logFile << "Computation + overhead time: " << compAndOverheadTime << " Seconds" << std::endl;
	logFile.flush();

	commStat->logStatistics(2, logFile);
	logFile.flush();

	double commTime = commStat->getTotalCommunicationTime();
	logFile << "Total communication time: " << commTime << " Seconds" << std::endl;
	logFile << "Computation without communication time: " << computationTime - commTime << " Seconds" << std::endl;
	logFile.flush();

	// doing task end environmental processing and memory cleanup
	copyBackNonArrayEnvVariables(environment, &taskGlobals);
	environment->executeTaskCompletionInstructions();
	delete taskData;
}


/*--------------------------------------------------------------------------------------------------------------
function for the initialize block
--------------------------------------------------------------------------------------------------------------*/

void fps::initializeTask(ArrayMetadata *arrayMetadata, 
		EnvironmentLinks environmentLinks, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		FPSPartition partition, 
		int iterations) {

	arrayMetadata->plateDims[0] = environmentLinks.plateDims[0];
	arrayMetadata->plateDims[1] = environmentLinks.plateDims[1];
	taskGlobals->max_iterations = iterations;
}


/*--------------------------------------------------------------------------------------------------------------
functions for compute stages
--------------------------------------------------------------------------------------------------------------*/

int fps::refineestimates_stage_7(SpaceB_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		FPSPartition partition, 
		std::ofstream &logFile) {

	//-------------------- Local Copies of Metadata -----------------------------

	// create local copies of partition and storage dimension configs of all arrays
	Dimension platePartDims[2];
	platePartDims[0] = lpu->platePartDims[0].partition;
	platePartDims[1] = lpu->platePartDims[1].partition;
	Dimension plateStoreDims[2];
	plateStoreDims[0] = lpu->platePartDims[0].storage;
	plateStoreDims[1] = lpu->platePartDims[1].storage;

	// create a local part-dimension object for later use
	PartDimension partConfig;

	// create a local transformed index variable for later use
	int xformIndex;

	//------------------- Local Variable Declarations ---------------------------

	Range localCols;
	Range localRows;
	int lpuId[2];

	//----------------------- Computation Begins --------------------------------

	localRows = platePartDims[0].range;
	localCols = platePartDims[1].range;
	{// scope entrance for parallel loop on index i
	int i;
	int iterationStart = platePartDims[0].range.min;
	int iterationBound = platePartDims[0].range.max;
	int indexIncrement = 1;
	int indexMultiplier = 1;
	if (platePartDims[0].range.min > platePartDims[0].range.max) {
		iterationBound *= -1;
		indexIncrement *= -1;
		indexMultiplier = -1;
	}
	{// scope entrance for applying index loop restrictions
	int localIterationStart = iterationStart;
	int localIterationBound = iterationBound;
	if (platePartDims[0].range.min > platePartDims[0].range.max) {
		localIterationStart = localRows.max - 1;
		localIterationBound = localRows.min + 1;
	} else {
		localIterationBound = localRows.max - 1;
		localIterationStart = localRows.min + 1;
	}
	if (platePartDims[0].range.min > platePartDims[0].range.max) {
		iterationStart = min(iterationStart, localIterationStart);
		iterationBound = max(iterationBound, localIterationBound);
	} else {
		iterationStart = max(iterationStart, localIterationStart);
		iterationBound = min(iterationBound, localIterationBound);
	}
	}// scope exit for applying index loop restrictions
	for (i = iterationStart; 
			indexMultiplier * i <= iterationBound; 
			i += indexIncrement) {
		long int iplate0 = ((long) (i - plateStoreDims[0].range.min)) * ((long) (plateStoreDims[1].length));
		{// scope entrance for parallel loop on index j
		int j;
		int iterationStart = platePartDims[1].range.min;
		int iterationBound = platePartDims[1].range.max;
		int indexIncrement = 1;
		int indexMultiplier = 1;
		if (platePartDims[1].range.min > platePartDims[1].range.max) {
			iterationBound *= -1;
			indexIncrement *= -1;
			indexMultiplier = -1;
		}
		{// scope entrance for applying index loop restrictions
		int localIterationStart = iterationStart;
		int localIterationBound = iterationBound;
		if (platePartDims[1].range.min > platePartDims[1].range.max) {
			localIterationStart = localCols.max - 1;
			localIterationBound = localCols.min + 1;
		} else {
			localIterationBound = localCols.max - 1;
			localIterationStart = localCols.min + 1;
		}
		if (platePartDims[1].range.min > platePartDims[1].range.max) {
			iterationStart = min(iterationStart, localIterationStart);
			iterationBound = max(iterationBound, localIterationBound);
		} else {
			iterationStart = max(iterationStart, localIterationStart);
			iterationBound = min(iterationBound, localIterationBound);
		}
		}// scope exit for applying index loop restrictions
		for (j = iterationStart; 
				indexMultiplier * j <= iterationBound; 
				j += indexIncrement) {
			long int jplate1 = ((long) (j - plateStoreDims[1].range.min));
			(lpu->plate[iplate0 + jplate1]) = (0.25 * ((((lpu->plate_lag_1[((long) ((i - 1) - plateStoreDims[0].range.min)) * ((long) (plateStoreDims[1].length)) + jplate1] + lpu->plate_lag_1[((long) ((i + 1) - plateStoreDims[0].range.min)) * ((long) (plateStoreDims[1].length)) + jplate1]) + lpu->plate_lag_1[iplate0 + ((long) ((j - 1) - plateStoreDims[1].range.min))]) + lpu->plate_lag_1[iplate0 + ((long) ((j + 1) - plateStoreDims[1].range.min))])));
		}
		}// scope exit for parallel loop on index j
	}
	}// scope exit for parallel loop on index i

	//------------------------- Returning Flag ----------------------------------

	return SUCCESS_RUN;
}


/*--------------------------------------------------------------------------------------------------------------
run method for thread simulating the task flow
--------------------------------------------------------------------------------------------------------------*/

void fps::run(ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		FPSPartition partition, ThreadStateImpl *threadState) {

	// log thread's affinity information
	threadState->logThreadAffinity();

	// set the root LPU in the thread state so that calculation can start
	LPU *rootLpu = threadState->getCurrentLpu(Space_Root);
	if (rootLpu == NULL) {
		threadState->setRootLpu(arrayMetadata);
	}

	// initialize thread's sync primitives holder data structure
	ThreadSyncPrimitive *threadSync = getSyncPrimitives(threadState->getThreadIds());

	// create counter variables for communicators
	int commCounter0 = 0;
	int commCounter1 = 0;

	// create a local part-dimension object for later use
	PartDimension partConfig;

	// create a local transformed index variable for later use
	int xformIndex;

	{ // scope entrance for repeat loop
	int repeatIteration = 0;
	int iterationStart = Range(1, taskGlobals->max_iterations).min;
	int iterationBound = Range(1, taskGlobals->max_iterations).max;
	int indexIncrement = 1;
	int indexMultiplier = 1;
	if (Range(1, taskGlobals->max_iterations).min > Range(1, taskGlobals->max_iterations).max) {
		iterationBound *= -1;
		indexIncrement *= -1;
		indexMultiplier = -1;
	}
	for (threadLocals->counter_1 = iterationStart; 
			indexMultiplier * threadLocals->counter_1 <= iterationBound; 
			threadLocals->counter_1 += indexIncrement) {

		// declaration of synchronization counter variables
		int plateStage9No1 = 0;

		// waiting on data reception
		if (repeatIteration > 0) {
			if (threadState->isValidPpu(Space_B)) {
				Communicator *communicator = threadState->getCommunicator("plateStage9No1");
				if (communicator != NULL) {
					communicator->receive(REQUESTING_COMMUNICATION, commCounter1);
				}
			}
			commCounter1++;
		}

		{ // scope entrance for iterating LPUs of Space A
		int spaceALpuId = INVALID_ID;
		int spaceAIteration = 0;
		SpaceA_LPU *spaceALpu = NULL;
		LPU *lpu = NULL;
		while((lpu = threadState->getNextLpu(Space_A, Space_Root, spaceALpuId)) != NULL) {
			spaceALpu = (SpaceA_LPU*) lpu;

			{ // scope entrance for repeat loop
			int repeatIteration = 0;
			int iterationStart = Range(1, partition.p1).min;
			int iterationBound = Range(1, partition.p1).max;
			int indexIncrement = partition.p2;
			int indexMultiplier = 1;
			if (Range(1, partition.p1).min > Range(1, partition.p1).max) {
				iterationBound *= -1;
				indexIncrement *= -1;
				indexMultiplier = -1;
			}
			for (threadLocals->counter_2 = iterationStart; 
					indexMultiplier * threadLocals->counter_2 <= iterationBound; 
					threadLocals->counter_2 += indexIncrement) {

				// declaration of synchronization counter variables
				int plateStage8No1 = 0;

				// waiting on data reception
				if (repeatIteration > 0) {
					if (threadState->isValidPpu(Space_B)) {
						Communicator *communicator = threadState->getCommunicator("plateStage8No1");
						if (communicator != NULL) {
							communicator->receive(REQUESTING_COMMUNICATION, commCounter0);
						}
					}
					commCounter0++;
				}

				{ // scope entrance for iterating LPUs of Space B
				int spaceBLpuId = INVALID_ID;
				int spaceBIteration = 0;
				SpaceB_LPU *spaceBLpu = NULL;
				LPU *lpu = NULL;
				while((lpu = threadState->getNextLpu(Space_B, Space_A, spaceBLpuId)) != NULL) {
					spaceBLpu = (SpaceB_LPU*) lpu;

					{ // scope entrance for repeat loop
					int repeatIteration = 0;
					int iterationStart = Range(1, partition.p2).min;
					int iterationBound = Range(1, partition.p2).max;
					int indexIncrement = 1;
					int indexMultiplier = 1;
					if (Range(1, partition.p2).min > Range(1, partition.p2).max) {
						iterationBound *= -1;
						indexIncrement *= -1;
						indexMultiplier = -1;
					}
					for (threadLocals->counter_3 = iterationStart; 
							indexMultiplier * threadLocals->counter_3 <= iterationBound; 
							threadLocals->counter_3 += indexIncrement) {

						{ // scope starts for epoch boundary
						TaskData *taskData = threadState->getTaskData();
						Hashtable<DataPartitionConfig*> *partConfigMap = threadState->getPartConfigMap();
						{ // scope starts for version updates of LPU data from Space B
						List<int*> *lpuIdChain = threadState->getLpuIdChainWithoutCopy(
								Space_B, Space_Root);
						DataPartitionConfig *config = partConfigMap->Lookup("plateSpaceBConfig");
						PartIterator *iterator = threadState->getIterator(Space_B, "plate");
						List<int*> *partId = iterator->getPartIdTemplate();
						config->generatePartId(lpuIdChain, partId);
						DataItems *items = taskData->getDataItemsOfLps("B", "plate");
						DataPart *dataPart = items->getDataPart(partId, iterator);
						dataPart->advanceEpoch();
						} // scope ends for version updates of LPU data from Space B
						generateSpaceBLpu(threadState, partConfigMap, taskData);
						spaceBLpu = (SpaceB_LPU*) threadState->getCurrentLpu(Space_B);

						if (threadState->isValidPpu(Space_B)) {
							// invoking user computation
							int stage7Executed = refineestimates_stage_7(spaceBLpu, 
									arrayMetadata, 
									taskGlobals, 
									threadLocals, 
									partition, 
									threadState->threadLog);
						}
						} // scope ends for epoch boundary

						repeatIteration++;
					}
					} // scope exit for repeat loop
					spaceBLpuId = spaceBLpu->id;
					spaceBIteration++;
				}
				threadState->removeIterationBound(Space_A);
				} // scope exit for iterating LPUs of Space B

				// barriers to ensure all readers have finished reading last update
				if (threadState->isValidPpu(Space_B)) {
					threadSync->plateStage8No1ReverseSync->wait();
				}

				// communicating updates
				if (threadState->isValidPpu(Space_B)) {
					Communicator *communicator = threadState->getCommunicator("plateStage8No1");
					if (communicator != NULL) {
						communicator->send(REQUESTING_COMMUNICATION, commCounter0);
					}
				}
				repeatIteration++;
			}
			} // scope exit for repeat loop
			spaceALpuId = spaceALpu->id;
			spaceAIteration++;
		}
		} // scope exit for iterating LPUs of Space A

		// barriers to ensure all readers have finished reading last update
		if (threadState->isValidPpu(Space_B)) {
			threadSync->plateStage9No1ReverseSync->wait();
		}

		// communicating updates
		if (threadState->isValidPpu(Space_A)) {
			Communicator *communicator = threadState->getCommunicator("plateStage9No1");
			if (communicator != NULL) {
				communicator->send(REQUESTING_COMMUNICATION, commCounter1);
			}
		}
		repeatIteration++;
	}
	} // scope exit for repeat loop

	// logging iterators' efficiency
	threadState->logIteratorStatistics();

	// close thread's log file
	threadState->closeLogFile();
}

/*-----------------------------------------------------------------------------------
PThreads run function
------------------------------------------------------------------------------------*/

void *fps::runPThreads(void *argument) {
	PThreadArg *pthreadArg = (PThreadArg *) argument;
	ThreadStateImpl *threadState = pthreadArg->threadState;
	run(pthreadArg->metadata, 
			pthreadArg->taskGlobals, 
			pthreadArg->threadLocals, 
			pthreadArg->partition, 
			threadState);
	pthread_exit(NULL);
}

