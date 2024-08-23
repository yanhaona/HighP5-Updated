
/*--------------------------------------------------------------------------------------------------------------
header file for the task
--------------------------------------------------------------------------------------------------------------*/
#include "monte_carlo_area_estimation.h"


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


// header files needed to execute external code blocks
#include <iostream>

using namespace mcae;


/*--------------------------------------------------------------------------------------------------------------
Functions for ArrayMetadata and EnvironmentLinks
--------------------------------------------------------------------------------------------------------------*/

mcae::ArrayMetadata::ArrayMetadata() : Metadata() {
	setTaskName("Monte Carlo Area Estimation");
}

void mcae::ArrayMetadata::print(std::ofstream &stream) {
	stream << "Array Metadata" << std::endl;
	stream << "Array: estimate_diffs";
	stream << ' ';
	estimate_diffsDims[0].print(stream);
	stream << ' ';
	estimate_diffsDims[1].print(stream);
	stream << std::endl;
	stream << "Array: grid";
	stream << ' ';
	gridDims[0].print(stream);
	stream << ' ';
	gridDims[1].print(stream);
	stream << std::endl;
	stream << "Array: internal_points";
	stream << ' ';
	internal_pointsDims[0].print(stream);
	stream << ' ';
	internal_pointsDims[1].print(stream);
	stream << std::endl;
	stream << "Array: local_estimates";
	stream << ' ';
	local_estimatesDims[0].print(stream);
	stream << ' ';
	local_estimatesDims[1].print(stream);
	stream << std::endl;
	stream.flush();
}

/*--------------------------------------------------------------------------------------------------------------
LPU print functions
--------------------------------------------------------------------------------------------------------------*/

void mcae::SpaceRoot_LPU::print(std::ofstream &stream, int indentLevel) {
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: estimate_diffs" << std::endl;
	estimate_diffsPartDims[0].print(stream, indentLevel + 1);
	estimate_diffsPartDims[1].print(stream, indentLevel + 1);
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: grid" << std::endl;
	gridPartDims[0].print(stream, indentLevel + 1);
	gridPartDims[1].print(stream, indentLevel + 1);
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: internal_points" << std::endl;
	internal_pointsPartDims[0].print(stream, indentLevel + 1);
	internal_pointsPartDims[1].print(stream, indentLevel + 1);
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: local_estimates" << std::endl;
	local_estimatesPartDims[0].print(stream, indentLevel + 1);
	local_estimatesPartDims[1].print(stream, indentLevel + 1);
	stream.flush();
}

void mcae::SpaceA_LPU::print(std::ofstream &stream, int indentLevel) {
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: local_estimates" << std::endl;
	local_estimatesPartDims[0].print(stream, indentLevel + 1);
	local_estimatesPartDims[1].print(stream, indentLevel + 1);
	stream.flush();
}

void mcae::SpaceB_LPU::print(std::ofstream &stream, int indentLevel) {
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: estimate_diffs" << std::endl;
	estimate_diffsPartDims[0].print(stream, indentLevel + 1);
	estimate_diffsPartDims[1].print(stream, indentLevel + 1);
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: grid" << std::endl;
	gridPartDims[0].print(stream, indentLevel + 1);
	gridPartDims[1].print(stream, indentLevel + 1);
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: internal_points" << std::endl;
	internal_pointsPartDims[0].print(stream, indentLevel + 1);
	internal_pointsPartDims[1].print(stream, indentLevel + 1);
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: local_estimates" << std::endl;
	local_estimatesPartDims[0].print(stream, indentLevel + 1);
	local_estimatesPartDims[1].print(stream, indentLevel + 1);
	stream.flush();
}

void mcae::SpaceC_LPU::print(std::ofstream &stream, int indentLevel) {
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: grid" << std::endl;
	gridPartDims[0].print(stream, indentLevel + 1);
	gridPartDims[1].print(stream, indentLevel + 1);
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: internal_points" << std::endl;
	internal_pointsPartDims[0].print(stream, indentLevel + 1);
	internal_pointsPartDims[1].print(stream, indentLevel + 1);
	stream.flush();
}


/*--------------------------------------------------------------------------------------------------------------
functions for generating partition configuration objects for data structures
--------------------------------------------------------------------------------------------------------------*/

//-------------------------------------------------------------------------------------------------------Space A

DataPartitionConfig *mcae::getlocal_estimatesConfigForSpaceA(ArrayMetadata *metadata, 
		MCAEPartition partition, 
		int ppuCount, 
		Hashtable<DataPartitionConfig*> *configMap) {
	List<DimPartitionConfig*> *dimensionConfigs = new List<DimPartitionConfig*>;
	Assert(dimensionConfigs != NULL);
	dimensionConfigs->Append(new ReplicationConfig(metadata->local_estimatesDims[0]));
	dimensionConfigs->Append(new ReplicationConfig(metadata->local_estimatesDims[1]));
	DataPartitionConfig *config =  new DataPartitionConfig(2, dimensionConfigs);
	config->configureDimensionOrder();
	config->setLpsId(Space_A);
	return config;
}

//-------------------------------------------------------------------------------------------------------Space B

DataPartitionConfig *mcae::getestimate_diffsConfigForSpaceB(ArrayMetadata *metadata, 
		MCAEPartition partition, 
		int ppuCount, 
		Hashtable<DataPartitionConfig*> *configMap) {
	List<DimPartitionConfig*> *dimensionConfigs = new List<DimPartitionConfig*>;
	Assert(dimensionConfigs != NULL);
	int *dim0Paddings = new int[2];
	dim0Paddings[0] = 0;
	dim0Paddings[1] = 0;
	int *dim0Arguments = new int;
	dim0Arguments[0] = 1;
	dimensionConfigs->Append(new BlockSizeConfig(metadata->estimate_diffsDims[0], 
			dim0Arguments, dim0Paddings, ppuCount, 0));
	delete[] dim0Paddings;
	int *dim1Paddings = new int[2];
	dim1Paddings[0] = 0;
	dim1Paddings[1] = 0;
	int *dim1Arguments = new int;
	dim1Arguments[0] = 1;
	dimensionConfigs->Append(new BlockSizeConfig(metadata->estimate_diffsDims[1], 
			dim1Arguments, dim1Paddings, ppuCount, 1));
	delete[] dim1Paddings;
	DataPartitionConfig *config =  new DataPartitionConfig(2, dimensionConfigs);
	config->configureDimensionOrder();
	config->setLpsId(Space_B);
	return config;
}

DataPartitionConfig *mcae::getgridConfigForSpaceB(ArrayMetadata *metadata, 
		MCAEPartition partition, 
		int ppuCount, 
		Hashtable<DataPartitionConfig*> *configMap) {
	List<DimPartitionConfig*> *dimensionConfigs = new List<DimPartitionConfig*>;
	Assert(dimensionConfigs != NULL);
	int *dim0Paddings = new int[2];
	dim0Paddings[0] = 0;
	dim0Paddings[1] = 0;
	int *dim0Arguments = new int;
	dim0Arguments[0] = partition.b;
	dimensionConfigs->Append(new BlockSizeConfig(metadata->gridDims[0], 
			dim0Arguments, dim0Paddings, ppuCount, 0));
	delete[] dim0Paddings;
	int *dim1Paddings = new int[2];
	dim1Paddings[0] = 0;
	dim1Paddings[1] = 0;
	int *dim1Arguments = new int;
	dim1Arguments[0] = partition.b;
	dimensionConfigs->Append(new BlockSizeConfig(metadata->gridDims[1], 
			dim1Arguments, dim1Paddings, ppuCount, 1));
	delete[] dim1Paddings;
	DataPartitionConfig *config =  new DataPartitionConfig(2, dimensionConfigs);
	config->configureDimensionOrder();
	config->setLpsId(Space_B);
	return config;
}

DataPartitionConfig *mcae::getinternal_pointsConfigForSpaceB(ArrayMetadata *metadata, 
		MCAEPartition partition, 
		int ppuCount, 
		Hashtable<DataPartitionConfig*> *configMap) {
	List<DimPartitionConfig*> *dimensionConfigs = new List<DimPartitionConfig*>;
	Assert(dimensionConfigs != NULL);
	int *dim0Paddings = new int[2];
	dim0Paddings[0] = 0;
	dim0Paddings[1] = 0;
	int *dim0Arguments = new int;
	dim0Arguments[0] = partition.b;
	dimensionConfigs->Append(new BlockSizeConfig(metadata->internal_pointsDims[0], 
			dim0Arguments, dim0Paddings, ppuCount, 0));
	delete[] dim0Paddings;
	int *dim1Paddings = new int[2];
	dim1Paddings[0] = 0;
	dim1Paddings[1] = 0;
	int *dim1Arguments = new int;
	dim1Arguments[0] = partition.b;
	dimensionConfigs->Append(new BlockSizeConfig(metadata->internal_pointsDims[1], 
			dim1Arguments, dim1Paddings, ppuCount, 1));
	delete[] dim1Paddings;
	DataPartitionConfig *config =  new DataPartitionConfig(2, dimensionConfigs);
	config->configureDimensionOrder();
	config->setLpsId(Space_B);
	return config;
}

DataPartitionConfig *mcae::getlocal_estimatesConfigForSpaceB(ArrayMetadata *metadata, 
		MCAEPartition partition, 
		int ppuCount, 
		Hashtable<DataPartitionConfig*> *configMap) {
	List<DimPartitionConfig*> *dimensionConfigs = new List<DimPartitionConfig*>;
	Assert(dimensionConfigs != NULL);
	int *dim0Paddings = new int[2];
	dim0Paddings[0] = 0;
	dim0Paddings[1] = 0;
	int *dim0Arguments = new int;
	dim0Arguments[0] = 1;
	dimensionConfigs->Append(new BlockSizeConfig(metadata->local_estimatesDims[0], 
			dim0Arguments, dim0Paddings, ppuCount, 0));
	delete[] dim0Paddings;
	int *dim1Paddings = new int[2];
	dim1Paddings[0] = 0;
	dim1Paddings[1] = 0;
	int *dim1Arguments = new int;
	dim1Arguments[0] = 1;
	dimensionConfigs->Append(new BlockSizeConfig(metadata->local_estimatesDims[1], 
			dim1Arguments, dim1Paddings, ppuCount, 1));
	delete[] dim1Paddings;
	DataPartitionConfig *config =  new DataPartitionConfig(2, dimensionConfigs);
	config->setParent(configMap->Lookup("local_estimatesSpaceAConfig"), 1);
	config->configureDimensionOrder();
	config->setLpsId(Space_B);
	return config;
}

//-------------------------------------------------------------------------------------------------------Space C

DataPartitionConfig *mcae::getgridConfigForSpaceC(ArrayMetadata *metadata, 
		MCAEPartition partition, 
		int ppuCount, 
		Hashtable<DataPartitionConfig*> *configMap) {
	List<DimPartitionConfig*> *dimensionConfigs = new List<DimPartitionConfig*>;
	Assert(dimensionConfigs != NULL);
	int *dim0Paddings = new int[2];
	dim0Paddings[0] = 0;
	dim0Paddings[1] = 0;
	int *dim0Arguments = new int;
	dim0Arguments[0] = 1;
	dimensionConfigs->Append(new BlockSizeConfig(metadata->gridDims[0], 
			dim0Arguments, dim0Paddings, ppuCount, 0));
	delete[] dim0Paddings;
	int *dim1Paddings = new int[2];
	dim1Paddings[0] = 0;
	dim1Paddings[1] = 0;
	int *dim1Arguments = new int;
	dim1Arguments[0] = 1;
	dimensionConfigs->Append(new BlockSizeConfig(metadata->gridDims[1], 
			dim1Arguments, dim1Paddings, ppuCount, 1));
	delete[] dim1Paddings;
	DataPartitionConfig *config =  new DataPartitionConfig(2, dimensionConfigs);
	config->setParent(configMap->Lookup("gridSpaceBConfig"), 1);
	config->configureDimensionOrder();
	config->setLpsId(Space_C);
	return config;
}

DataPartitionConfig *mcae::getinternal_pointsConfigForSpaceC(ArrayMetadata *metadata, 
		MCAEPartition partition, 
		int ppuCount, 
		Hashtable<DataPartitionConfig*> *configMap) {
	List<DimPartitionConfig*> *dimensionConfigs = new List<DimPartitionConfig*>;
	Assert(dimensionConfigs != NULL);
	int *dim0Paddings = new int[2];
	dim0Paddings[0] = 0;
	dim0Paddings[1] = 0;
	int *dim0Arguments = new int;
	dim0Arguments[0] = 1;
	dimensionConfigs->Append(new BlockSizeConfig(metadata->internal_pointsDims[0], 
			dim0Arguments, dim0Paddings, ppuCount, 0));
	delete[] dim0Paddings;
	int *dim1Paddings = new int[2];
	dim1Paddings[0] = 0;
	dim1Paddings[1] = 0;
	int *dim1Arguments = new int;
	dim1Arguments[0] = 1;
	dimensionConfigs->Append(new BlockSizeConfig(metadata->internal_pointsDims[1], 
			dim1Arguments, dim1Paddings, ppuCount, 1));
	delete[] dim1Paddings;
	DataPartitionConfig *config =  new DataPartitionConfig(2, dimensionConfigs);
	config->setParent(configMap->Lookup("internal_pointsSpaceBConfig"), 1);
	config->configureDimensionOrder();
	config->setLpsId(Space_C);
	return config;
}

//---------------------------------------------------------------------------Partition Configuration Accumulator

Hashtable<DataPartitionConfig*> *mcae::getDataPartitionConfigMap(ArrayMetadata *metadata, 
		MCAEPartition partition, int *ppuCounts) {
	Hashtable<DataPartitionConfig*> *configMap = new Hashtable<DataPartitionConfig*>;
	DataPartitionConfig *local_estimatesSpaceAConfig = getlocal_estimatesConfigForSpaceA(metadata, 
		partition, ppuCounts[Space_A], configMap);
	configMap->Enter("local_estimatesSpaceAConfig", local_estimatesSpaceAConfig);
	DataPartitionConfig *estimate_diffsSpaceBConfig = getestimate_diffsConfigForSpaceB(metadata, 
		partition, ppuCounts[Space_B], configMap);
	configMap->Enter("estimate_diffsSpaceBConfig", estimate_diffsSpaceBConfig);
	DataPartitionConfig *gridSpaceBConfig = getgridConfigForSpaceB(metadata, 
		partition, ppuCounts[Space_B], configMap);
	configMap->Enter("gridSpaceBConfig", gridSpaceBConfig);
	DataPartitionConfig *internal_pointsSpaceBConfig = getinternal_pointsConfigForSpaceB(metadata, 
		partition, ppuCounts[Space_B], configMap);
	configMap->Enter("internal_pointsSpaceBConfig", internal_pointsSpaceBConfig);
	DataPartitionConfig *local_estimatesSpaceBConfig = getlocal_estimatesConfigForSpaceB(metadata, 
		partition, ppuCounts[Space_B], configMap);
	configMap->Enter("local_estimatesSpaceBConfig", local_estimatesSpaceBConfig);
	DataPartitionConfig *gridSpaceCConfig = getgridConfigForSpaceC(metadata, 
		partition, ppuCounts[Space_C], configMap);
	configMap->Enter("gridSpaceCConfig", gridSpaceCConfig);
	DataPartitionConfig *internal_pointsSpaceCConfig = getinternal_pointsConfigForSpaceC(metadata, 
		partition, ppuCounts[Space_C], configMap);
	configMap->Enter("internal_pointsSpaceCConfig", internal_pointsSpaceCConfig);
	return configMap;
}

/*--------------------------------------------------------------------------------------------------------------
functions for generating memory blocks for data parts of various LPUs
--------------------------------------------------------------------------------------------------------------*/

LpsContent *mcae::genSpaceBContent(List<ThreadState*> *threads, ArrayMetadata *metadata, 
		TaskEnvironment *environment, 
		MCAEPartition partition, 
		Hashtable<DataPartitionConfig*> *partConfigMap) {

	if(threads->NumElements() == 0) return NULL;
	LpsContent *spaceBContent = new LpsContent(Space_B);

	DataItems *estimate_diffs = new DataItems("estimate_diffs", 2, true);
	DataPartitionConfig *estimate_diffsConfig = partConfigMap->Lookup("estimate_diffsSpaceBConfig");
	estimate_diffs->setPartitionConfig(estimate_diffsConfig);
	DataPartsList *estimate_diffsParts = estimate_diffsConfig->generatePartList(1);
	estimate_diffs->setPartsList(estimate_diffsParts);
	std::vector<DimConfig> estimate_diffsDimOrder = *(estimate_diffsConfig->getDimensionOrder());
	PartIdContainer *estimate_diffsContainer = NULL;
	if (estimate_diffsDimOrder.size() == 1) estimate_diffsContainer = new PartContainer(estimate_diffsDimOrder[0]);
	else estimate_diffsContainer = new PartListContainer(estimate_diffsDimOrder[0]);
	Assert(estimate_diffsContainer != NULL);
	List<int*> *estimate_diffsPartId = estimate_diffsConfig->generatePartIdTemplate();
	spaceBContent->addDataItems("estimate_diffs", estimate_diffs);

	DataItems *local_estimates = new DataItems("local_estimates", 2, true);
	DataPartitionConfig *local_estimatesConfig = partConfigMap->Lookup("local_estimatesSpaceBConfig");
	local_estimates->setPartitionConfig(local_estimatesConfig);
	DataPartsList *local_estimatesParts = local_estimatesConfig->generatePartList(1);
	local_estimates->setPartsList(local_estimatesParts);
	std::vector<DimConfig> local_estimatesDimOrder = *(local_estimatesConfig->getDimensionOrder());
	PartIdContainer *local_estimatesContainer = NULL;
	if (local_estimatesDimOrder.size() == 1) local_estimatesContainer = new PartContainer(local_estimatesDimOrder[0]);
	else local_estimatesContainer = new PartListContainer(local_estimatesDimOrder[0]);
	Assert(local_estimatesContainer != NULL);
	List<int*> *local_estimatesPartId = local_estimatesConfig->generatePartIdTemplate();
	spaceBContent->addDataItems("local_estimates", local_estimates);

	for (int i = 0; i < threads->NumElements(); i++) {
		ThreadState *thread = threads->Nth(i);
		int lpuId = INVALID_ID;
		while((lpuId = thread->getNextLpuId(Space_B, Space_Root, lpuId)) != INVALID_ID) {
			List<int*> *lpuIdChain = thread->getLpuIdChainWithoutCopy(
					Space_B, Space_Root);
			estimate_diffsConfig->generatePartId(lpuIdChain, estimate_diffsPartId);
			estimate_diffsContainer->insertPartId(estimate_diffsPartId, 2, estimate_diffsDimOrder);
			local_estimatesConfig->generatePartId(lpuIdChain, local_estimatesPartId);
			local_estimatesContainer->insertPartId(local_estimatesPartId, 2, local_estimatesDimOrder);
		}
	}

	estimate_diffsParts->initializePartsList(estimate_diffsConfig, estimate_diffsContainer, sizeof(double));
	estimate_diffsParts->allocateParts();

	local_estimatesParts->initializePartsList(local_estimatesConfig, local_estimatesContainer, sizeof(double));
	local_estimatesParts->allocateParts();

	return spaceBContent;
}

LpsContent *mcae::genSpaceCContent(List<ThreadState*> *threads, ArrayMetadata *metadata, 
		TaskEnvironment *environment, 
		MCAEPartition partition, 
		Hashtable<DataPartitionConfig*> *partConfigMap) {

	if(threads->NumElements() == 0) return NULL;
	LpsContent *spaceCContent = new LpsContent(Space_C);

	DataItems *grid = new DataItems("grid", 2, true);
	DataPartitionConfig *gridConfig = partConfigMap->Lookup("gridSpaceCConfig");
	grid->setPartitionConfig(gridConfig);
	DataPartsList *gridParts = gridConfig->generatePartList(1);
	grid->setPartsList(gridParts);
	std::vector<DimConfig> gridDimOrder = *(gridConfig->getDimensionOrder());
	PartIdContainer *gridContainer = NULL;
	if (gridDimOrder.size() == 1) gridContainer = new PartContainer(gridDimOrder[0]);
	else gridContainer = new PartListContainer(gridDimOrder[0]);
	Assert(gridContainer != NULL);
	List<int*> *gridPartId = gridConfig->generatePartIdTemplate();
	spaceCContent->addDataItems("grid", grid);

	DataItems *internal_points = new DataItems("internal_points", 2, true);
	DataPartitionConfig *internal_pointsConfig = partConfigMap->Lookup("internal_pointsSpaceCConfig");
	internal_points->setPartitionConfig(internal_pointsConfig);
	DataPartsList *internal_pointsParts = internal_pointsConfig->generatePartList(1);
	internal_points->setPartsList(internal_pointsParts);
	std::vector<DimConfig> internal_pointsDimOrder = *(internal_pointsConfig->getDimensionOrder());
	PartIdContainer *internal_pointsContainer = NULL;
	if (internal_pointsDimOrder.size() == 1) internal_pointsContainer = new PartContainer(internal_pointsDimOrder[0]);
	else internal_pointsContainer = new PartListContainer(internal_pointsDimOrder[0]);
	Assert(internal_pointsContainer != NULL);
	List<int*> *internal_pointsPartId = internal_pointsConfig->generatePartIdTemplate();
	spaceCContent->addDataItems("internal_points", internal_points);

	for (int i = 0; i < threads->NumElements(); i++) {
		ThreadState *thread = threads->Nth(i);
		int lpuId = INVALID_ID;
		while((lpuId = thread->getNextLpuId(Space_C, Space_Root, lpuId)) != INVALID_ID) {
			List<int*> *lpuIdChain = thread->getLpuIdChainWithoutCopy(
					Space_C, Space_Root);
			gridConfig->generatePartId(lpuIdChain, gridPartId);
			gridContainer->insertPartId(gridPartId, 2, gridDimOrder);
			internal_pointsConfig->generatePartId(lpuIdChain, internal_pointsPartId);
			internal_pointsContainer->insertPartId(internal_pointsPartId, 2, internal_pointsDimOrder);
		}
	}

	gridParts->initializePartsList(gridConfig, gridContainer, sizeof(Rectangle));
	gridParts->allocateParts();

	internal_pointsParts->initializePartsList(internal_pointsConfig, internal_pointsContainer, sizeof(int));
	internal_pointsParts->allocateParts();

	return spaceCContent;
}

/*--------------------------------------------------------------------------------------------------------------
functions for creating containers of non-task-global reduction results
--------------------------------------------------------------------------------------------------------------*/

void mcae::prepareSpaceBReductionResultContainers(List<ThreadState*> *threads, TaskData *taskData) {

	if(threads->NumElements() == 0) return;

	std::vector<int> lpuIdDimensions;
	lpuIdDimensions.insert(lpuIdDimensions.begin(), 2);
	lpuIdDimensions.insert(lpuIdDimensions.begin(), 1);

	ReductionResultAccessContainer *placement_resultContainer = 
			new ReductionResultAccessContainer(lpuIdDimensions);
	taskData->addReductionResultContainer("placement_result", placement_resultContainer);

	for (int i = 0; i < threads->NumElements(); i++) {
		ThreadState *thread = threads->Nth(i);
		int lpuId = INVALID_ID;
		while((lpuId = thread->getNextLpuId(Space_B, Space_Root, lpuId)) != INVALID_ID) {
			List<int*> *lpuIdChain = thread->getLpuIdChainWithoutCopy(
					Space_B, Space_Root);
			placement_resultContainer->initiateResultForLpu(lpuIdChain);
		}
	}
}

//-----------------------------------------------------------------------------------------Task Data Initializer

TaskData *mcae::initializeTaskData(List<ThreadState*> *threads, ArrayMetadata *metadata, 
		TaskEnvironment *environment, 
		SegmentState *segment, 
		MCAEPartition partition, int *ppuCounts) {

	Hashtable<DataPartitionConfig*> *configMap = 
			getDataPartitionConfigMap(metadata, partition, ppuCounts);
	TaskData *taskData = new TaskData();
	Assert(taskData != NULL);
	preconfigureLpsAllocationsInEnv(environment, metadata, configMap);

	// prepare LPS contents map
	LpsContent *spaceBContent = genSpaceBContent(threads, metadata, 
			environment, partition, configMap);
	taskData->addLpsContent("B", spaceBContent);
	LpsContent *spaceCContent = genSpaceCContent(threads, metadata, 
			environment, partition, configMap);
	taskData->addLpsContent("C", spaceCContent);

	// prepare file I/O handlers in the environment
	Hashtable<PartReader*> *readersMap = generateReadersMap(taskData, segment, configMap);
	environment->setReadersMap(readersMap);
	Hashtable<PartWriter*> *writersMap = generateWritersMap(taskData, segment, configMap);
	environment->setWritersMap(writersMap);

	// initialize parts lists of environmental variables
	environment->preprocessProgramEnvForItems();
	environment->setupItemsPartsLists();
	environment->postprocessProgramEnvForItems();

	// prepare reduction results access containers
	prepareSpaceBReductionResultContainers(threads, taskData);

	return taskData;
}

/*--------------------------------------------------------------------------------------------------------------
file I/O for environmental data structures
--------------------------------------------------------------------------------------------------------------*/

//------------------------------------------------------------------------------------Parts-reader map generator
Hashtable<PartReader*> *mcae::generateReadersMap(TaskData *taskData, 
			SegmentState *segment, 
			Hashtable<DataPartitionConfig*> *partConfigMap) {

	Hashtable<PartReader*> *readersMap = new Hashtable<PartReader*>;

	return readersMap;
}

//------------------------------------------------------------------------------------Parts-writer map generator
Hashtable<PartWriter*> *mcae::generateWritersMap(TaskData *taskData, 
			SegmentState *segment, 
			Hashtable<DataPartitionConfig*> *partConfigMap) {

	Hashtable<PartWriter*> *writersMap = new Hashtable<PartWriter*>;

	int writerId;
	MPI_Comm_rank(MPI_COMM_WORLD, &writerId);
	int mpiProcessCount;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiProcessCount);
	int writersCount = min(mpiProcessCount, Max_Segments_Count);

	return writersMap;
}

/*--------------------------------------------------------------------------------------------------------------
functions for retrieving partition counts in different LPSes
--------------------------------------------------------------------------------------------------------------*/

int *mcae::getLPUsCountOfSpaceB(Hashtable<DataPartitionConfig*> *partConfigMap) {
	int *count = new int[2];
	DataPartitionConfig *gridConfig = partConfigMap->Lookup("gridSpaceBConfig");
	count[0] = gridConfig->getPartsCountAlongDimension(0);
	count[1] = gridConfig->getPartsCountAlongDimension(1);
	return count;
}

int *mcae::getLPUsCountOfSpaceC(Hashtable<DataPartitionConfig*> *partConfigMap, Dimension gridDim1, Dimension gridDim2) {
	int *count = new int[2];
	DataPartitionConfig *gridConfig = partConfigMap->Lookup("gridSpaceCConfig");
	count[0] = gridConfig->getPartsCountAlongDimension(0, &gridDim1);
	count[1] = gridConfig->getPartsCountAlongDimension(1, &gridDim2);
	return count;
}

/*--------------------------------------------------------------------------------------------------------------
functions for generating LPUs given LPU Ids
--------------------------------------------------------------------------------------------------------------*/

void mcae::generateSpaceALpu(ThreadState *threadState, 
			Hashtable<DataPartitionConfig*> *partConfigMap, 
			TaskData *taskData) {

	int *lpuId = threadState->getCurrentLpuId(Space_A);
	int *lpuCounts = threadState->getLpuCounts(Space_A);
	SpaceA_LPU *lpu = (SpaceA_LPU*) threadState->getCurrentLpu(Space_A, true);
	List<int*> *lpuIdChain = threadState->getLpuIdChainWithoutCopy(Space_A, Space_Root);

	DataPartitionConfig *local_estimatesConfig = partConfigMap->Lookup("local_estimatesSpaceAConfig");
	SpaceRoot_LPU *spaceRootLpu = (SpaceRoot_LPU*) threadState->getCurrentLpu(Space_Root);
	PartDimension *local_estimatesParentPartDims = spaceRootLpu->local_estimatesPartDims;
	local_estimatesConfig->updatePartDimensionInfo(lpuId, lpuCounts, lpu->local_estimatesPartDims, local_estimatesParentPartDims);
}

void mcae::generateSpaceBLpu(ThreadState *threadState, 
			Hashtable<DataPartitionConfig*> *partConfigMap, 
			TaskData *taskData) {

	int *lpuId = threadState->getCurrentLpuId(Space_B);
	int *lpuCounts = threadState->getLpuCounts(Space_B);
	SpaceB_LPU *lpu = (SpaceB_LPU*) threadState->getCurrentLpu(Space_B, true);
	List<int*> *lpuIdChain = threadState->getLpuIdChainWithoutCopy(Space_B, Space_Root);

	lpu->lpuId[0] = lpuId[0];
	lpu->lpuId[1] = lpuId[1];

	DataPartitionConfig *estimate_diffsConfig = partConfigMap->Lookup("estimate_diffsSpaceBConfig");
	SpaceRoot_LPU *spaceRootLpu = (SpaceRoot_LPU*) threadState->getCurrentLpu(Space_Root);
	PartDimension *estimate_diffsParentPartDims = spaceRootLpu->estimate_diffsPartDims;
	estimate_diffsConfig->updatePartDimensionInfo(lpuId, lpuCounts, lpu->estimate_diffsPartDims, estimate_diffsParentPartDims);

	if (taskData != NULL) {
		PartIterator *iterator = threadState->getIterator(Space_B, "estimate_diffs");
		List<int*> *partId = iterator->getPartIdTemplate();
		estimate_diffsConfig->generatePartId(lpuIdChain, partId);
		DataItems *estimate_diffsItems = taskData->getDataItemsOfLps("B", "estimate_diffs");
		DataPart *estimate_diffsPart = estimate_diffsItems->getDataPart(partId, iterator);
		estimate_diffsPart->getMetadata()->updateStorageDimension(lpu->estimate_diffsPartDims);
		lpu->estimate_diffs = (double*) estimate_diffsPart->getData();
	}

	DataPartitionConfig *gridConfig = partConfigMap->Lookup("gridSpaceBConfig");
	PartDimension *gridParentPartDims = spaceRootLpu->gridPartDims;
	gridConfig->updatePartDimensionInfo(lpuId, lpuCounts, lpu->gridPartDims, gridParentPartDims);

	DataPartitionConfig *internal_pointsConfig = partConfigMap->Lookup("internal_pointsSpaceBConfig");
	PartDimension *internal_pointsParentPartDims = spaceRootLpu->internal_pointsPartDims;
	internal_pointsConfig->updatePartDimensionInfo(lpuId, lpuCounts, lpu->internal_pointsPartDims, internal_pointsParentPartDims);

	DataPartitionConfig *local_estimatesConfig = partConfigMap->Lookup("local_estimatesSpaceBConfig");
	SpaceA_LPU *spaceALpu = (SpaceA_LPU*) threadState->getCurrentLpu(Space_A);
	PartDimension *local_estimatesParentPartDims = spaceALpu->local_estimatesPartDims;
	local_estimatesConfig->updatePartDimensionInfo(lpuId, lpuCounts, lpu->local_estimatesPartDims, local_estimatesParentPartDims);

	if (taskData != NULL) {
		PartIterator *iterator = threadState->getIterator(Space_B, "local_estimates");
		List<int*> *partId = iterator->getPartIdTemplate();
		local_estimatesConfig->generatePartId(lpuIdChain, partId);
		DataItems *local_estimatesItems = taskData->getDataItemsOfLps("B", "local_estimates");
		DataPart *local_estimatesPart = local_estimatesItems->getDataPart(partId, iterator);
		local_estimatesPart->getMetadata()->updateStorageDimension(lpu->local_estimatesPartDims);
		lpu->local_estimates = (double*) local_estimatesPart->getData();
	}
}

void mcae::generateSpaceCLpu(ThreadState *threadState, 
			Hashtable<DataPartitionConfig*> *partConfigMap, 
			TaskData *taskData) {

	int *lpuId = threadState->getCurrentLpuId(Space_C);
	int *lpuCounts = threadState->getLpuCounts(Space_C);
	SpaceC_LPU *lpu = (SpaceC_LPU*) threadState->getCurrentLpu(Space_C, true);
	List<int*> *lpuIdChain = threadState->getLpuIdChainWithoutCopy(Space_C, Space_Root);

	lpu->lpuId[0] = lpuId[0];
	lpu->lpuId[1] = lpuId[1];

	DataPartitionConfig *gridConfig = partConfigMap->Lookup("gridSpaceCConfig");
	SpaceB_LPU *spaceBLpu = (SpaceB_LPU*) threadState->getCurrentLpu(Space_B);
	PartDimension *gridParentPartDims = spaceBLpu->gridPartDims;
	gridConfig->updatePartDimensionInfo(lpuId, lpuCounts, lpu->gridPartDims, gridParentPartDims);

	if (taskData != NULL) {
		PartIterator *iterator = threadState->getIterator(Space_C, "grid");
		List<int*> *partId = iterator->getPartIdTemplate();
		gridConfig->generatePartId(lpuIdChain, partId);
		DataItems *gridItems = taskData->getDataItemsOfLps("C", "grid");
		DataPart *gridPart = gridItems->getDataPart(partId, iterator);
		gridPart->getMetadata()->updateStorageDimension(lpu->gridPartDims);
		lpu->grid = (Rectangle*) gridPart->getData();
	}

	DataPartitionConfig *internal_pointsConfig = partConfigMap->Lookup("internal_pointsSpaceCConfig");
	PartDimension *internal_pointsParentPartDims = spaceBLpu->internal_pointsPartDims;
	internal_pointsConfig->updatePartDimensionInfo(lpuId, lpuCounts, lpu->internal_pointsPartDims, internal_pointsParentPartDims);

	if (taskData != NULL) {
		PartIterator *iterator = threadState->getIterator(Space_C, "internal_points");
		List<int*> *partId = iterator->getPartIdTemplate();
		internal_pointsConfig->generatePartId(lpuIdChain, partId);
		DataItems *internal_pointsItems = taskData->getDataItemsOfLps("C", "internal_points");
		DataPart *internal_pointsPart = internal_pointsItems->getDataPart(partId, iterator);
		internal_pointsPart->getMetadata()->updateStorageDimension(lpu->internal_pointsPartDims);
		lpu->internal_points = (int*) internal_pointsPart->getData();
	}
}

/*--------------------------------------------------------------------------------------------------------------
functions to generate PPU IDs and PPU group IDs for a thread
--------------------------------------------------------------------------------------------------------------*/

ThreadIds *mcae::getPpuIdsForThread(int threadNo)  {

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
	groupSize = threadCount;
	groupThreadId = idsArray[Space_Root] % groupSize;
	threadIds->ppuIds[Space_A].groupId = idsArray[Space_Root] / groupSize;
	threadIds->ppuIds[Space_A].ppuCount = 40;
	threadIds->ppuIds[Space_A].groupSize = groupSize;
	if (groupThreadId == 0) threadIds->ppuIds[Space_A].id
			= threadIds->ppuIds[Space_A].groupId;
	else threadIds->ppuIds[Space_A].id = INVALID_ID;
	idsArray[Space_A] = groupThreadId;

	// for Space B;
	threadIds->ppuIds[Space_B].lpsName = "B";
	threadCount = threadIds->ppuIds[Space_A].groupSize;
	groupSize = threadCount / 1;
	groupThreadId = idsArray[Space_A] % groupSize;
	threadIds->ppuIds[Space_B].groupId = idsArray[Space_A] / groupSize;
	threadIds->ppuIds[Space_B].ppuCount = 1;
	threadIds->ppuIds[Space_B].groupSize = groupSize;
	if (groupThreadId == 0) threadIds->ppuIds[Space_B].id
			= threadIds->ppuIds[Space_B].groupId;
	else threadIds->ppuIds[Space_B].id = INVALID_ID;
	idsArray[Space_B] = groupThreadId;

	// for Space C;
	threadIds->ppuIds[Space_C].lpsName = "C";
	threadCount = threadIds->ppuIds[Space_B].groupSize;
	groupSize = threadCount / 1;
	groupThreadId = idsArray[Space_B] % groupSize;
	threadIds->ppuIds[Space_C].groupId = idsArray[Space_B] / groupSize;
	threadIds->ppuIds[Space_C].ppuCount = 1;
	threadIds->ppuIds[Space_C].groupSize = groupSize;
	if (groupThreadId == 0) threadIds->ppuIds[Space_C].id
			= threadIds->ppuIds[Space_C].groupId;
	else threadIds->ppuIds[Space_C].id = INVALID_ID;
	idsArray[Space_C] = groupThreadId;

	return threadIds;
}


void mcae::adjustPpuCountsAndGroupSizes(ThreadIds *threadId)  {

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

	groupId = threadId->ppuIds[Space_C].groupId;
	groupSize = threadId->ppuIds[Space_C].groupSize;
	ppuCount = ((groupEnd - groupBegin + 1) + (groupSize - 1)) / groupSize;
	threadId->ppuIds[Space_C].ppuCount = ppuCount;
	groupBegin = groupId * groupSize;
	groupEnd = min(groupBegin + groupSize - 1, groupEnd);
	threadId->ppuIds[Space_C].groupSize = groupEnd - groupBegin + 1;

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
	lpsParentIndexMap[Space_C] = Space_B;
}

//-----------------------------------------------------------------------------------------Root LPU Construction

void ThreadStateImpl::setRootLpu(Metadata *metadata) {

	ArrayMetadata *arrayMetadata = (ArrayMetadata*) metadata;

	SpaceRoot_LPU *lpu = new SpaceRoot_LPU;
	lpu->estimate_diffs = NULL;
	lpu->estimate_diffsPartDims[0] = PartDimension();
	lpu->estimate_diffsPartDims[0].partition = arrayMetadata->estimate_diffsDims[0];
	lpu->estimate_diffsPartDims[0].storage = arrayMetadata->estimate_diffsDims[0].getNormalizedDimension();
	lpu->estimate_diffsPartDims[1] = PartDimension();
	lpu->estimate_diffsPartDims[1].partition = arrayMetadata->estimate_diffsDims[1];
	lpu->estimate_diffsPartDims[1].storage = arrayMetadata->estimate_diffsDims[1].getNormalizedDimension();

	lpu->grid = NULL;
	lpu->gridPartDims[0] = PartDimension();
	lpu->gridPartDims[0].partition = arrayMetadata->gridDims[0];
	lpu->gridPartDims[0].storage = arrayMetadata->gridDims[0].getNormalizedDimension();
	lpu->gridPartDims[1] = PartDimension();
	lpu->gridPartDims[1].partition = arrayMetadata->gridDims[1];
	lpu->gridPartDims[1].storage = arrayMetadata->gridDims[1].getNormalizedDimension();

	lpu->internal_points = NULL;
	lpu->internal_pointsPartDims[0] = PartDimension();
	lpu->internal_pointsPartDims[0].partition = arrayMetadata->internal_pointsDims[0];
	lpu->internal_pointsPartDims[0].storage = arrayMetadata->internal_pointsDims[0].getNormalizedDimension();
	lpu->internal_pointsPartDims[1] = PartDimension();
	lpu->internal_pointsPartDims[1].partition = arrayMetadata->internal_pointsDims[1];
	lpu->internal_pointsPartDims[1].storage = arrayMetadata->internal_pointsDims[1].getNormalizedDimension();

	lpu->local_estimates = NULL;
	lpu->local_estimatesPartDims[0] = PartDimension();
	lpu->local_estimatesPartDims[0].partition = arrayMetadata->local_estimatesDims[0];
	lpu->local_estimatesPartDims[0].storage = arrayMetadata->local_estimatesDims[0].getNormalizedDimension();
	lpu->local_estimatesPartDims[1] = PartDimension();
	lpu->local_estimatesPartDims[1].partition = arrayMetadata->local_estimatesDims[1];
	lpu->local_estimatesPartDims[1].storage = arrayMetadata->local_estimatesDims[1].getNormalizedDimension();

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
	lpsStates[Space_C]->lpu = new SpaceC_LPU;
	lpsStates[Space_C]->lpu->setValidBit(false);
}

//--------------------------------------------------------------------------------------------LPU Count Function


int *ThreadStateImpl::computeLpuCounts(int lpsId) {
	Hashtable<DataPartitionConfig*> *configMap = getPartConfigMap();
	if (lpsId == Space_Root) {
		return NULL;
	}
	if (lpsId == Space_A) {
		return NULL;
	}
	if (lpsId == Space_B) {
		return getLPUsCountOfSpaceB(configMap);
	}
	if (lpsId == Space_C) {
		SpaceB_LPU *spaceBLpu
				 = (SpaceB_LPU*) lpsStates[Space_B]->lpu;
		return getLPUsCountOfSpaceC(configMap, 
				spaceBLpu->gridPartDims[0].partition, 
				spaceBLpu->gridPartDims[1].partition);
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
	if (lpsId == Space_C) {
		SpaceC_LPU *currentLpu = (SpaceC_LPU*) lpsStates[Space_C]->lpu;
		generateSpaceCLpu(this, partConfigMap, taskData);
		currentLpu->setValidBit(true);
		return currentLpu;
	}
	return NULL;
}

//-------------------------------------------------------------------------Reduction Result Map Creator Function

void ThreadStateImpl::initializeReductionResultMap() {
	localReductionResultMap = new Hashtable<reduction::Result*>;
	localReductionResultMap->Enter("placement_result", new reduction::Result());
	localReductionResultMap->Enter("area", new reduction::Result());
}

/*--------------------------------------------------------------------------------------------------------------
Initializer function for global synchronization primitives
--------------------------------------------------------------------------------------------------------------*/

void mcae::initializeSyncPrimitives() {

	for (int i = 0; i < Space_B_Threads_Per_Segment; i++) {
		int participants = Space_C_Threads_Per_Segment / Space_B_Threads_Per_Segment;
		estimate_diffsStage6No1DSyncs[i] = new RS(participants);
		estimate_diffsStage6No1ReverseSyncs[i] = new Barrier(participants);
	}
	for (int i = 0; i < Space_B_Threads_Per_Segment; i++) {
		int participants = Space_C_Threads_Per_Segment / Space_B_Threads_Per_Segment;
		estimate_diffsStage12No1DSyncs[i] = new RS(participants);
		estimate_diffsStage12No1ReverseSyncs[i] = new Barrier(participants);
	}
}

/*--------------------------------------------------------------------------------------------------------------
function for initializing thread's sync primitives
--------------------------------------------------------------------------------------------------------------*/

ThreadSyncPrimitive *mcae::getSyncPrimitives(ThreadIds *threadIds) {

	ThreadSyncPrimitive *threadSync = new ThreadSyncPrimitive();

	int spaceBGroup = threadIds->ppuIds[Space_B].groupId % Space_B_Threads_Per_Segment;
	threadSync->estimate_diffsStage6No1DSync = estimate_diffsStage6No1DSyncs[spaceBGroup];
	threadSync->estimate_diffsStage6No1ReverseSync = estimate_diffsStage6No1ReverseSyncs[spaceBGroup];
	threadSync->estimate_diffsStage12No1DSync = estimate_diffsStage12No1DSyncs[spaceBGroup];
	threadSync->estimate_diffsStage12No1ReverseSync = estimate_diffsStage12No1ReverseSyncs[spaceBGroup];

	return threadSync;
}


/*--------------------------------------------------------------------------------------------------------------
Reduction Primitives Management
--------------------------------------------------------------------------------------------------------------*/

//---------------------------------------------------------------------Primitive for variable 'placement_result'

mcae::ReductionPrimitive_placement_result::ReductionPrimitive_placement_result(int localParticipants)
		: NonTaskGlobalReductionPrimitive(sizeof(int), SUM, localParticipants) {}

void mcae::ReductionPrimitive_placement_result::resetPartialResult(reduction::Result *resultVar) {
	resultVar->data.intValue = 0;
}

void mcae::ReductionPrimitive_placement_result::updateIntermediateResult(
		reduction::Result *localPartialResult) {
	intermediateResult->data.intValue += localPartialResult->data.intValue;
}

//---------------------------------------------------------------------------------Primitive for variable 'area'

mcae::ReductionPrimitive_area::ReductionPrimitive_area(int localParticipants)
		: TaskGlobalReductionPrimitive(sizeof(double), SUM, localParticipants) {}

void mcae::ReductionPrimitive_area::resetPartialResult(reduction::Result *resultVar) {
	resultVar->data.doubleValue = 0;
}

void mcae::ReductionPrimitive_area::updateIntermediateResult(
		reduction::Result *localPartialResult) {
	intermediateResult->data.doubleValue += localPartialResult->data.doubleValue;
}

//-------------------------------------------------------------------------------Reduction Primitive Initializer

void mcae::setupReductionPrimitives(std::ofstream &logFile) {

	int segmentId , segmentCount;
	MPI_Comm_rank(MPI_COMM_WORLD, &segmentId);
	MPI_Comm_size(MPI_COMM_WORLD, &segmentCount);


	//---------------------------------------------------------------------Primitives for 'placement_result'

	{ // scope starts
	for(int i = 0; i < Space_B_Threads_Per_Segment; i++) {
		placement_resultReducer[i] = new ReductionPrimitive_placement_result(
				Space_C_Threads_Per_Segment / Space_B_Threads_Per_Segment);
		placement_resultReducer[i]->setLogFile(&logFile);
	}
	} // scope ends


	//---------------------------------------------------------------------------------Primitives for 'area'

	{ // scope starts
	for(int i = 0; i < Space_A_Threads_Per_Segment; i++) {
		areaReducer[i] = new ReductionPrimitive_area(
				Space_B_Threads_Per_Segment / Space_A_Threads_Per_Segment);
		areaReducer[i]->setLogFile(&logFile);
	}
	} // scope ends
}

//--------------------------------------------------------------------------------Reduction Primitives Retriever

Hashtable<void*> *mcae::getReductionPrimitiveMap(ThreadIds *threadIds) {

	Hashtable<void*> *rdPrimitiveMap = new Hashtable<void*>;
	if(threadIds->ppuIds[Space_C].id != INVALID_ID) {
		int spaceBGroup = threadIds->ppuIds[Space_B].groupId % Space_B_Threads_Per_Segment;
		rdPrimitiveMap->Enter("placement_result", placement_resultReducer[spaceBGroup]);
	}
	if(threadIds->ppuIds[Space_B].id != INVALID_ID) {
		int spaceAGroup = threadIds->ppuIds[Space_A].groupId % Space_A_Threads_Per_Segment;
		rdPrimitiveMap->Enter("area", areaReducer[spaceAGroup]);
	}
	return rdPrimitiveMap;
}

/*--------------------------------------------------------------------------------------------------------------
Task Environment Management Structures and Functions
--------------------------------------------------------------------------------------------------------------*/

//---------------------------------------------------------------------------------------------------Constructor

mcae::TaskEnvironmentImpl::TaskEnvironmentImpl() : TaskEnvironment() {
	area = 0;
	prepareItemsMap();
	resetEnvInstructions();
}

//---------------------------------------------------------------------Task Environment Function Implementations

void mcae::TaskEnvironmentImpl::prepareItemsMap() {
}

void mcae::TaskEnvironmentImpl::setDefaultTaskCompletionInstrs() {
}

//--------------------------------------------------------------------------Environmental Links Object Generator

EnvironmentLinks mcae::initiateEnvLinks(TaskEnvironment *environment) {

	EnvironmentLinks links;
	mcae::TaskEnvironmentImpl *taskEnv = (mcae::TaskEnvironmentImpl *) environment;

	return links;
}

//---------------------------------------------------------------------------------LPS Allocation Preconfigurers

void mcae::preconfigureLpsAllocationsInEnv(TaskEnvironment *environment, 
		ArrayMetadata *metadata, 
		Hashtable<DataPartitionConfig*> *partConfigMap) {

}

//------------------------------------------------------------------------------------Non-array variables copier

void mcae::copyBackNonArrayEnvVariables(TaskEnvironment *environment, TaskGlobals *taskGlobals) {
	mcae::TaskEnvironmentImpl *taskEnv = (mcae::TaskEnvironmentImpl *) environment;
	taskEnv->area = taskGlobals->area;
}

/*--------------------------------------------------------------------------------------------------------------
task executor function
--------------------------------------------------------------------------------------------------------------*/

void mcae::execute(TaskEnvironment *environment, 
		double precision_threshold, 
		int max_rounds, 
		int cell_length, 
		int grid_dim, 
		int points_per_cell, 
		MCAEPartition partition, 
		int segmentId, 
		std::ofstream &logFile) {

	environment->setLogFile(&logFile);
	if (segmentId >= Max_Segments_Count) {
		logFile << "Current segment does not participate in: Monte Carlo Area Estimation\n";
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
	partitionArgs = new int[1];
	partitionArgs[0] = partition.b;

	// initializing sync primitives
	initializeSyncPrimitives();

	// initializing reduction primitives
	setupReductionPrimitives(logFile);

	// invoking the initializer function
	initializeTask(metadata, envLinks, &taskGlobals, &threadLocals, partition, precision_threshold, max_rounds, cell_length, grid_dim, points_per_cell);
	metadata->estimate_diffsDims[0].setLength();
	metadata->estimate_diffsDims[1].setLength();
	metadata->gridDims[0].setLength();
	metadata->gridDims[1].setLength();
	metadata->internal_pointsDims[0].setLength();
	metadata->internal_pointsDims[1].setLength();
	metadata->local_estimatesDims[0].setLength();
	metadata->local_estimatesDims[1].setLength();
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
	lpsDimensions[Space_A] = 0;
	lpsDimensions[Space_B] = 2;
	lpsDimensions[Space_C] = 2;
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
		threadStateList[i]->initiateLogFile("mcae");
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

	// starting threads
	logFile << "\tlaunching threads\n";
	logFile.flush();
	pthread_t threads[Total_Threads];
	PThreadArg *threadArgs[Total_Threads];
	for (int i = 0; i < Total_Threads; i++) {
		threadArgs[i] = new PThreadArg;
		threadArgs[i]->taskName = "Monte Carlo Area Estimation";
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

	// doing task end environmental processing and memory cleanup
	copyBackNonArrayEnvVariables(environment, &taskGlobals);
	environment->executeTaskCompletionInstructions();
	delete taskData;
}


/*--------------------------------------------------------------------------------------------------------------
function for the initialize block
--------------------------------------------------------------------------------------------------------------*/

void mcae::initializeTask(ArrayMetadata *arrayMetadata, 
		EnvironmentLinks environmentLinks, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		MCAEPartition partition, 
		double precision_threshold, 
		int max_rounds, 
		int cell_length, 
		int grid_dim, 
		int points_per_cell) {

	taskGlobals->precision_threshold = precision_threshold;
	taskGlobals->max_rounds = max_rounds;
	taskGlobals->cell_length = cell_length;
	taskGlobals->points_per_cell = points_per_cell;
	arrayMetadata->gridDims[0].range.min = 0;
	arrayMetadata->gridDims[0].range.max = (grid_dim - 1);
	arrayMetadata->gridDims[1] = arrayMetadata->gridDims[0];
	arrayMetadata->internal_pointsDims[0] = arrayMetadata->gridDims[0];
	arrayMetadata->internal_pointsDims[1] = arrayMetadata->gridDims[1];
	arrayMetadata->local_estimatesDims[0].range.min = 0;
	arrayMetadata->local_estimatesDims[0].range.max = ((grid_dim / partition.b) - 1);
	arrayMetadata->local_estimatesDims[1] = arrayMetadata->local_estimatesDims[0];
	arrayMetadata->estimate_diffsDims[0] = arrayMetadata->local_estimatesDims[0];
	arrayMetadata->estimate_diffsDims[1] = arrayMetadata->local_estimatesDims[1];
	init_rand_0();
}


/*--------------------------------------------------------------------------------------------------------------
functions for compute stages
--------------------------------------------------------------------------------------------------------------*/

int mcae::setupgridcells_stage_4(SpaceC_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		MCAEPartition partition, 
		std::ofstream &logFile) {

	//-------------------- Local Copies of Metadata -----------------------------

	// create local copies of partition and storage dimension configs of all arrays
	Dimension gridPartDims[2];
	gridPartDims[0] = lpu->gridPartDims[0].partition;
	gridPartDims[1] = lpu->gridPartDims[1].partition;
	Dimension gridStoreDims[2];
	gridStoreDims[0] = lpu->gridPartDims[0].storage;
	gridStoreDims[1] = lpu->gridPartDims[1].storage;

	// create a local part-dimension object for later use
	PartDimension partConfig;

	// create a local transformed index variable for later use
	int xformIndex;

	//------------------- Local Variable Declarations ---------------------------

	int cell_height;
	int cell_width;
	int lpuId[2];

	//----------------------- Computation Begins --------------------------------

	{// scope entrance for parallel loop on index i
	int i = gridPartDims[0].range.min;
	long int igrid0 = ((long) (i - gridStoreDims[0].range.min)) * ((long) (gridStoreDims[1].length));
	{// scope entrance for parallel loop on index j
	int j = gridPartDims[1].range.min;
	long int jgrid1 = ((long) (j - gridStoreDims[1].range.min));
	cell_height = taskGlobals->cell_length;
	cell_width = taskGlobals->cell_length;
	lpu->grid[igrid0 + jgrid1].left = (cell_width * i);
	lpu->grid[igrid0 + jgrid1].right = ((cell_width * (i + 1)) - 1);
	lpu->grid[igrid0 + jgrid1].bottom = (cell_height * j);
	lpu->grid[igrid0 + jgrid1].top = ((cell_height * (j + 1)) - 1);
	}// scope exit for parallel loop on index j
	}// scope exit for parallel loop on index i

	//------------------------- Returning Flag ----------------------------------

	return SUCCESS_RUN;
}


int mcae::initializeestimatediffs_stage_5(SpaceB_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		MCAEPartition partition, 
		std::ofstream &logFile) {

	//-------------------- Local Copies of Metadata -----------------------------

	// create local copies of partition and storage dimension configs of all arrays
	Dimension estimate_diffsPartDims[2];
	estimate_diffsPartDims[0] = lpu->estimate_diffsPartDims[0].partition;
	estimate_diffsPartDims[1] = lpu->estimate_diffsPartDims[1].partition;
	Dimension estimate_diffsStoreDims[2];
	estimate_diffsStoreDims[0] = lpu->estimate_diffsPartDims[0].storage;
	estimate_diffsStoreDims[1] = lpu->estimate_diffsPartDims[1].storage;

	// create a local part-dimension object for later use
	PartDimension partConfig;

	// create a local transformed index variable for later use
	int xformIndex;

	//------------------- Local Variable Declarations ---------------------------

	int lpuId[2];
	double threshold;

	//----------------------- Computation Begins --------------------------------

	threshold = taskGlobals->precision_threshold;
	{// scope entrance for parallel loop on index i
	int i = estimate_diffsPartDims[0].range.min;
	long int iestimate_diffs0 = ((long) (i - estimate_diffsStoreDims[0].range.min)) * ((long) (estimate_diffsStoreDims[1].length));
	{// scope entrance for parallel loop on index j
	int j = estimate_diffsPartDims[1].range.min;
	long int jestimate_diffs1 = ((long) (j - estimate_diffsStoreDims[1].range.min));
	lpu->estimate_diffs[iestimate_diffs0 + jestimate_diffs1] = (threshold + 1);
	}// scope exit for parallel loop on index j
	}// scope exit for parallel loop on index i

	//------------------------- Returning Flag ----------------------------------

	return SUCCESS_RUN;
}


int mcae::performsampling_stage_9(SpaceC_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		Hashtable<reduction::Result*> *localReductionResultMap, 
		MCAEPartition partition, 
		std::ofstream &logFile) {

	//-------------------- Local Copies of Metadata -----------------------------

	// create local copies of partition and storage dimension configs of all arrays
	Dimension gridPartDims[2];
	gridPartDims[0] = lpu->gridPartDims[0].partition;
	gridPartDims[1] = lpu->gridPartDims[1].partition;
	Dimension gridStoreDims[2];
	gridStoreDims[0] = lpu->gridPartDims[0].storage;
	gridStoreDims[1] = lpu->gridPartDims[1].storage;
	Dimension internal_pointsPartDims[2];
	internal_pointsPartDims[0] = lpu->internal_pointsPartDims[0].partition;
	internal_pointsPartDims[1] = lpu->internal_pointsPartDims[1].partition;
	Dimension internal_pointsStoreDims[2];
	internal_pointsStoreDims[0] = lpu->internal_pointsPartDims[0].storage;
	internal_pointsStoreDims[1] = lpu->internal_pointsPartDims[1].storage;

	// create a local part-dimension object for later use
	PartDimension partConfig;

	// create a local transformed index variable for later use
	int xformIndex;

	//------------------ Partial Results of Reductions  -------------------------

	reduction::Result *placement_result = localReductionResultMap->Lookup("placement_result");

	//------------------- Local Variable Declarations ---------------------------

	Rectangle cell;
	int lpuId[2];
	int result;
	int sample_count;
	int seed;

	//----------------------- Computation Begins --------------------------------

	sample_count = taskGlobals->points_per_cell;
	{// scope entrance for parallel loop on index i
	int i = gridPartDims[0].range.min;
	long int igrid0 = ((long) (i - gridStoreDims[0].range.min)) * ((long) (gridStoreDims[1].length));
	long int iinternal_points0 = ((long) (i - internal_pointsStoreDims[0].range.min)) * ((long) (internal_pointsStoreDims[1].length));
	{// scope entrance for parallel loop on index j
	int j = gridPartDims[1].range.min;
	long int jgrid1 = ((long) (j - gridStoreDims[1].range.min));
	long int jinternal_points1 = ((long) (j - internal_pointsStoreDims[1].range.min));
	cell = lpu->grid[igrid0 + jgrid1];
	seed = lpuId[0];
	lpu->internal_points[iinternal_points0 + jinternal_points1] = perform_sampling_0(cell, seed, taskGlobals->points_per_cell);
	}// scope exit for parallel loop on index j
	}// scope exit for parallel loop on index i
	{// scope entrance for parallel loop on index i
	int i = internal_pointsPartDims[0].range.min;
	long int iinternal_points0 = ((long) (i - internal_pointsStoreDims[0].range.min)) * ((long) (internal_pointsStoreDims[1].length));
	{// scope entrance for parallel loop on index j
	int j = internal_pointsPartDims[1].range.min;
	long int jinternal_points1 = ((long) (j - internal_pointsStoreDims[1].range.min));
	placement_result->data.intValue += lpu->internal_points[iinternal_points0 + jinternal_points1];
	}// scope exit for parallel loop on index j
	}// scope exit for parallel loop on index i

	//------------------------- Returning Flag ----------------------------------

	return SUCCESS_RUN;
}


int mcae::estimatesubarea_stage_10(SpaceB_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		MCAEPartition partition, 
		std::ofstream &logFile) {

	//-------------------- Local Copies of Metadata -----------------------------

	// create local copies of partition and storage dimension configs of all arrays
	Dimension estimate_diffsPartDims[2];
	estimate_diffsPartDims[0] = lpu->estimate_diffsPartDims[0].partition;
	estimate_diffsPartDims[1] = lpu->estimate_diffsPartDims[1].partition;
	Dimension estimate_diffsStoreDims[2];
	estimate_diffsStoreDims[0] = lpu->estimate_diffsPartDims[0].storage;
	estimate_diffsStoreDims[1] = lpu->estimate_diffsPartDims[1].storage;
	Dimension internal_pointsPartDims[2];
	internal_pointsPartDims[0] = lpu->internal_pointsPartDims[0].partition;
	internal_pointsPartDims[1] = lpu->internal_pointsPartDims[1].partition;
	Dimension internal_pointsStoreDims[2];
	internal_pointsStoreDims[0] = lpu->internal_pointsPartDims[0].storage;
	internal_pointsStoreDims[1] = lpu->internal_pointsPartDims[1].storage;
	Dimension local_estimatesPartDims[2];
	local_estimatesPartDims[0] = lpu->local_estimatesPartDims[0].partition;
	local_estimatesPartDims[1] = lpu->local_estimatesPartDims[1].partition;
	Dimension local_estimatesStoreDims[2];
	local_estimatesStoreDims[0] = lpu->local_estimatesPartDims[0].storage;
	local_estimatesStoreDims[1] = lpu->local_estimatesPartDims[1].storage;

	// create a local part-dimension object for later use
	PartDimension partConfig;

	// create a local transformed index variable for later use
	int xformIndex;

	//------------------- Local Variable Declarations ---------------------------

	int cell_size;
	int col;
	int curr_estimate;
	int internal_points;
	int lpuId[2];
	double old_estimate;
	int row;
	int sample_count;
	double updated_estimate;
	float weight;

	//----------------------- Computation Begins --------------------------------

	internal_points = threadLocals->placement_result;
	sample_count = taskGlobals->points_per_cell;
	row = local_estimatesPartDims[0].range.min;
	col = local_estimatesPartDims[1].range.min;
	cell_size = (taskGlobals->cell_length * taskGlobals->cell_length);
	curr_estimate = ((cell_size * internal_points) / sample_count);
	weight = (1 / threadLocals->round);
	old_estimate = lpu->local_estimates[((long) (row - local_estimatesStoreDims[0].range.min)) * ((long) (local_estimatesStoreDims[1].length)) + ((long) (col - local_estimatesStoreDims[1].range.min))];
	updated_estimate = ((curr_estimate * weight) + (old_estimate * (1 - weight)));
	lpu->local_estimates[((long) (row - local_estimatesStoreDims[0].range.min)) * ((long) (local_estimatesStoreDims[1].length)) + ((long) (col - local_estimatesStoreDims[1].range.min))] = updated_estimate;
	if ((updated_estimate > old_estimate)) {
		lpu->estimate_diffs[((long) (row - estimate_diffsStoreDims[0].range.min)) * ((long) (estimate_diffsStoreDims[1].length)) + ((long) (col - estimate_diffsStoreDims[1].range.min))] = (updated_estimate - old_estimate);
	} else {
		lpu->estimate_diffs[((long) (row - estimate_diffsStoreDims[0].range.min)) * ((long) (estimate_diffsStoreDims[1].length)) + ((long) (col - estimate_diffsStoreDims[1].range.min))] = (old_estimate - updated_estimate);
	}

	//------------------------- Returning Flag ----------------------------------

	return SUCCESS_RUN;
}


int mcae::estimatetotalarea_stage_11(SpaceB_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		Hashtable<reduction::Result*> *localReductionResultMap, 
		MCAEPartition partition, 
		std::ofstream &logFile) {

	//-------------------- Local Copies of Metadata -----------------------------

	// create local copies of partition and storage dimension configs of all arrays
	Dimension local_estimatesPartDims[2];
	local_estimatesPartDims[0] = lpu->local_estimatesPartDims[0].partition;
	local_estimatesPartDims[1] = lpu->local_estimatesPartDims[1].partition;
	Dimension local_estimatesStoreDims[2];
	local_estimatesStoreDims[0] = lpu->local_estimatesPartDims[0].storage;
	local_estimatesStoreDims[1] = lpu->local_estimatesPartDims[1].storage;

	// create a local part-dimension object for later use
	PartDimension partConfig;

	// create a local transformed index variable for later use
	int xformIndex;

	//------------------ Partial Results of Reductions  -------------------------

	reduction::Result *area = localReductionResultMap->Lookup("area");

	//------------------- Local Variable Declarations ---------------------------

	int lpuId[2];
	double result;

	//----------------------- Computation Begins --------------------------------

	{// scope entrance for parallel loop on index i
	int i = local_estimatesPartDims[0].range.min;
	long int ilocal_estimates0 = ((long) (i - local_estimatesStoreDims[0].range.min)) * ((long) (local_estimatesStoreDims[1].length));
	{// scope entrance for parallel loop on index j
	int j = local_estimatesPartDims[1].range.min;
	long int jlocal_estimates1 = ((long) (j - local_estimatesStoreDims[1].range.min));
	area->data.doubleValue += lpu->local_estimates[ilocal_estimates0 + jlocal_estimates1];
	}// scope exit for parallel loop on index j
	}// scope exit for parallel loop on index i

	//------------------------- Returning Flag ----------------------------------

	return SUCCESS_RUN;
}


int mcae::displayresult_stage_12(SpaceA_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		MCAEPartition partition, 
		std::ofstream &logFile) {

	//-------------------- Local Copies of Metadata -----------------------------

	// create local copies of partition and storage dimension configs of all arrays

	// create a local part-dimension object for later use
	PartDimension partConfig;

	// create a local transformed index variable for later use
	int xformIndex;

	//------------------- Local Variable Declarations ---------------------------

	double result;

	//----------------------- Computation Begins --------------------------------

	result = taskGlobals->area;

	{ // starting scope for an external code block

		// generating local variables for global scalars
		double area = taskGlobals->area;
		int cell_length = taskGlobals->cell_length;
		int max_rounds = taskGlobals->max_rounds;
		int placement_result = threadLocals->placement_result;
		int points_per_cell = taskGlobals->points_per_cell;
		double precision_threshold = taskGlobals->precision_threshold;
		int round = threadLocals->round;

		// generating local variables for array dimension metadata
		Dimension estimate_diffs_dimension[2];
		estimate_diffs_dimension[0] = arrayMetadata->estimate_diffsDims[0];
		estimate_diffs_dimension[1] = arrayMetadata->estimate_diffsDims[1];
		Dimension grid_dimension[2];
		grid_dimension[0] = arrayMetadata->gridDims[0];
		grid_dimension[1] = arrayMetadata->gridDims[1];
		Dimension internal_points_dimension[2];
		internal_points_dimension[0] = arrayMetadata->internal_pointsDims[0];
		internal_points_dimension[1] = arrayMetadata->internal_pointsDims[1];
		Dimension local_estimates_dimension[2];
		local_estimates_dimension[0] = arrayMetadata->local_estimatesDims[0];
		local_estimates_dimension[1] = arrayMetadata->local_estimatesDims[1];

		{ // external code block starts

                                        std::cout << "Estimated area under the polynomial is: " << result << "\n";
                                
		} // external code block ends

	} // ending scope for the external code block


	//------------------------- Returning Flag ----------------------------------

	return SUCCESS_RUN;
}


/*--------------------------------------------------------------------------------------------------------------
run method for thread simulating the task flow
--------------------------------------------------------------------------------------------------------------*/

void mcae::run(ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		MCAEPartition partition, ThreadStateImpl *threadState) {

	// log thread's affinity information
	threadState->logThreadAffinity();

	// set the root LPU in the thread state so that calculation can start
	LPU *rootLpu = threadState->getCurrentLpu(Space_Root);
	if (rootLpu == NULL) {
		threadState->setRootLpu(arrayMetadata);
	}

	// initialize thread's sync primitives holder data structure
	ThreadSyncPrimitive *threadSync = getSyncPrimitives(threadState->getThreadIds());

	// initializing a map for holding local, partial results of reductions
	threadState->initializeReductionResultMap();
	Hashtable<reduction::Result*> *reductionResultsMap = threadState->getLocalReductionResultMap();

	// retrieving the reduction primitives relevant to the current thread
	ThreadIds *threadIds = threadState->getThreadIds();
	Hashtable<void*> *rdPrimitiveMap = getReductionPrimitiveMap(threadIds);

	// create a local part-dimension object for later use
	PartDimension partConfig;

	// create a local transformed index variable for later use
	int xformIndex;

	// declaration of synchronization counter variables
	int estimate_diffsStage6No1 = 0;

	{ // scope entrance for iterating LPUs of Space A
	int spaceALpuId = INVALID_ID;
	int spaceAIteration = 0;
	SpaceA_LPU *spaceALpu = NULL;
	LPU *lpu = NULL;
	while((lpu = threadState->getNextLpu(Space_A, Space_Root, spaceALpuId)) != NULL) {
		spaceALpu = (SpaceA_LPU*) lpu;

		{ // beginning of a reduction boundary

		// initializing thread-local reduction result variables
		if(threadState->isValidPpu(Space_B)) {
			reduction::Result *areaLocal = reductionResultsMap->Lookup("area");
			TaskGlobalReductionPrimitive *rdPrimitive = 
					(TaskGlobalReductionPrimitive *) rdPrimitiveMap->Lookup("area");
			rdPrimitive->resetPartialResult(areaLocal);
		}


		{ // scope entrance for iterating LPUs of Space B
		int spaceBLpuId = INVALID_ID;
		int spaceBIteration = 0;
		SpaceB_LPU *spaceBLpu = NULL;
		LPU *lpu = NULL;
		while((lpu = threadState->getNextLpu(Space_B, Space_A, spaceBLpuId)) != NULL) {
			spaceBLpu = (SpaceB_LPU*) lpu;

			{ // scope entrance for reduction result preparation

			TaskData *taskData = threadState->getTaskData();
			List<int*> *lpuIdChain = threadState->getLpuIdChainWithoutCopy(Space_B, Space_Root);

			// processing for 'placement_result'
			reduction::Result *placement_result = taskData->getResultVar("placement_result", lpuIdChain);
			threadLocals->placement_result = placement_result->data.intValue;

			} // scope exit for reduction result preparation

			{ // scope entrance for iterating LPUs of Space C
			int spaceCLpuId = INVALID_ID;
			int spaceCIteration = 0;
			SpaceC_LPU *spaceCLpu = NULL;
			LPU *lpu = NULL;
			while((lpu = threadState->getNextLpu(Space_C, Space_B, spaceCLpuId)) != NULL) {
				spaceCLpu = (SpaceC_LPU*) lpu;
				if (threadState->isValidPpu(Space_C)) {
					// invoking user computation
					int stage5Executed = setupgridcells_stage_4(spaceCLpu, 
							arrayMetadata, 
							taskGlobals, 
							threadLocals, 
							partition, 
							threadState->threadLog);
				}
				spaceCLpuId = spaceCLpu->id;
				spaceCIteration++;
			}
			threadState->removeIterationBound(Space_B);
			} // scope exit for iterating LPUs of Space C

			// barriers to ensure all readers have finished reading last update
			if (threadState->isValidPpu(Space_C)) {
				threadSync->estimate_diffsStage6No1ReverseSync->wait();
			}
			if (threadState->isValidPpu(Space_B)) {
				// invoking user computation
				int stage6Executed = initializeestimatediffs_stage_5(spaceBLpu, 
						arrayMetadata, 
						taskGlobals, 
						threadLocals, 
						partition, 
						threadState->threadLog);
				estimate_diffsStage6No1 += stage6Executed;
			}

			// resolving synchronization dependencies
			if (estimate_diffsStage6No1 > 0 && threadState->isValidPpu(Space_B)) {
				threadSync->estimate_diffsStage6No1DSync->signal(0);
				estimate_diffsStage6No1 = 0;
			} else if (threadState->isValidPpu(Space_C)) {
				threadSync->estimate_diffsStage6No1DSync->wait(0);
			}

			{ // scope entrance for repeat loop
			int repeatIteration = 0;
			int iterationStart = Range(1, taskGlobals->max_rounds).min;
			int iterationBound = Range(1, taskGlobals->max_rounds).max;
			int indexIncrement = 1;
			int indexMultiplier = 1;
			if (Range(1, taskGlobals->max_rounds).min > Range(1, taskGlobals->max_rounds).max) {
				iterationBound *= -1;
				indexIncrement *= -1;
				indexMultiplier = -1;
			}
			for (threadLocals->round = iterationStart; 
					indexMultiplier * threadLocals->round <= iterationBound; 
					threadLocals->round += indexIncrement) {

				// declaration of synchronization counter variables
				int estimate_diffsStage12No1 = 0;

				// barriers to ensure all readers have finished reading last update
				if (threadState->isValidPpu(Space_C)) {
					threadSync->estimate_diffsStage12No1ReverseSync->wait();
				}
				{ // scope entrance for conditional subflow
				Dimension estimate_diffsPartDims[2];
				Dimension estimate_diffsStoreDims[2];
				estimate_diffsPartDims[0] = spaceBLpu->estimate_diffsPartDims[0].partition;
				estimate_diffsStoreDims[0] = spaceBLpu->estimate_diffsPartDims[0].storage;
				estimate_diffsPartDims[1] = spaceBLpu->estimate_diffsPartDims[1].partition;
				estimate_diffsStoreDims[1] = spaceBLpu->estimate_diffsPartDims[1].storage;
				if((spaceBLpu->estimate_diffs[((long) (spaceBLpu->lpuId[0] - estimate_diffsStoreDims[0].range.min)) * ((long) (estimate_diffsStoreDims[1].length)) + ((long) (spaceBLpu->lpuId[1] - estimate_diffsStoreDims[1].range.min))] > taskGlobals->precision_threshold)) {

					{ // beginning of a reduction boundary

					// initializing thread-local reduction result variables
					if(threadState->isValidPpu(Space_C)) {
						reduction::Result *placement_resultLocal = reductionResultsMap->Lookup("placement_result");
						NonTaskGlobalReductionPrimitive *rdPrimitive = 
								(NonTaskGlobalReductionPrimitive *) rdPrimitiveMap->Lookup("placement_result");
						rdPrimitive->resetPartialResult(placement_resultLocal);
					}


					{ // scope entrance for iterating LPUs of Space C
					int spaceCLpuId = INVALID_ID;
					int spaceCIteration = 0;
					SpaceC_LPU *spaceCLpu = NULL;
					LPU *lpu = NULL;
					while((lpu = threadState->getNextLpu(Space_C, Space_B, spaceCLpuId)) != NULL) {
						spaceCLpu = (SpaceC_LPU*) lpu;
						if (threadState->isValidPpu(Space_C)) {
							// invoking user computation
							int stage11Executed = performsampling_stage_9(spaceCLpu, 
									arrayMetadata, 
									taskGlobals, 
									threadLocals, 
									reductionResultsMap, 
									partition, 
									threadState->threadLog);
						}
						spaceCLpuId = spaceCLpu->id;
						spaceCIteration++;
					}
					threadState->removeIterationBound(Space_B);
					} // scope exit for iterating LPUs of Space C

					// executing the final step of reductions
					if(threadState->isValidPpu(Space_C)) {
						reduction::Result *localResult = reductionResultsMap->Lookup("placement_result");
						void *target = &(threadLocals->placement_result);
						List<int*> *lpuIdChain = threadState->getLpuIdChainWithoutCopy(Space_B, Space_Root);
						TaskData *taskData = threadState->getTaskData();
						reduction::Result *placement_result = taskData->getResultVar("placement_result", lpuIdChain);
						NonTaskGlobalReductionPrimitive *rdPrimitive = 
								(NonTaskGlobalReductionPrimitive*) rdPrimitiveMap->Lookup("placement_result");
						rdPrimitive->reduce(localResult, target, placement_result);
					}

					} // ending of a reduction boundary
					if (threadState->isValidPpu(Space_B)) {
						// invoking user computation
						int stage12Executed = estimatesubarea_stage_10(spaceBLpu, 
								arrayMetadata, 
								taskGlobals, 
								threadLocals, 
								partition, 
								threadState->threadLog);
						estimate_diffsStage12No1 += stage12Executed;
					}
				} // end of condition checking block
				} // scope exit for conditional subflow

				// resolving synchronization dependencies
				if (estimate_diffsStage12No1 > 0 && threadState->isValidPpu(Space_B)) {
					threadSync->estimate_diffsStage12No1DSync->signal(repeatIteration);
					estimate_diffsStage12No1 = 0;
				} else if (threadState->isValidPpu(Space_C)) {
					threadSync->estimate_diffsStage12No1DSync->wait(repeatIteration);
				}
				repeatIteration++;
			}
			} // scope exit for repeat loop
			if (threadState->isValidPpu(Space_B)) {
				// invoking user computation
				int stage13Executed = estimatetotalarea_stage_11(spaceBLpu, 
						arrayMetadata, 
						taskGlobals, 
						threadLocals, 
						reductionResultsMap, 
						partition, 
						threadState->threadLog);
			}
			spaceBLpuId = spaceBLpu->id;
			spaceBIteration++;
		}
		threadState->removeIterationBound(Space_A);
		} // scope exit for iterating LPUs of Space B

		// executing the final step of reductions
		if(threadState->isValidPpu(Space_B)) {
			reduction::Result *localResult = reductionResultsMap->Lookup("area");
			void *target = &(taskGlobals->area);
			TaskGlobalReductionPrimitive *rdPrimitive = 
					(TaskGlobalReductionPrimitive*) rdPrimitiveMap->Lookup("area");
			rdPrimitive->reduce(localResult, target);
		}

		} // ending of a reduction boundary
		if (threadState->isValidPpu(Space_A)) {
			// invoking user computation
			int stage14Executed = displayresult_stage_12(spaceALpu, 
					arrayMetadata, 
					taskGlobals, 
					threadLocals, 
					partition, 
					threadState->threadLog);
		}
		spaceALpuId = spaceALpu->id;
		spaceAIteration++;
	}
	} // scope exit for iterating LPUs of Space A

	// logging iterators' efficiency
	threadState->logIteratorStatistics();

	// close thread's log file
	threadState->closeLogFile();
}

/*-----------------------------------------------------------------------------------
PThreads run function
------------------------------------------------------------------------------------*/

void *mcae::runPThreads(void *argument) {
	PThreadArg *pthreadArg = (PThreadArg *) argument;
	ThreadStateImpl *threadState = pthreadArg->threadState;
	run(pthreadArg->metadata, 
			pthreadArg->taskGlobals, 
			pthreadArg->threadLocals, 
			pthreadArg->partition, 
			threadState);
	pthread_exit(NULL);
}

