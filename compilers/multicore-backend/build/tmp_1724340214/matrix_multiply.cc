/*-----------------------------------------------------------------------------------
header file for the task
------------------------------------------------------------------------------------*/
#include "matrix_multiply.h"

/*-----------------------------------------------------------------------------------
header files included for different purposes
------------------------------------------------------------------------------------*/
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


using namespace mm;

/*-----------------------------------------------------------------------------------
Print functions for LPUs 
------------------------------------------------------------------------------------*/

void mm::SpaceRoot_LPU::print(std::ofstream &stream, int indentLevel) {
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: a" << std::endl;
	aPartDims[0].print(stream, indentLevel + 1);
	aPartDims[1].print(stream, indentLevel + 1);
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: b" << std::endl;
	bPartDims[0].print(stream, indentLevel + 1);
	bPartDims[1].print(stream, indentLevel + 1);
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: c" << std::endl;
	cPartDims[0].print(stream, indentLevel + 1);
	cPartDims[1].print(stream, indentLevel + 1);
	stream.flush();
}

void mm::SpaceA_LPU::print(std::ofstream &stream, int indentLevel) {
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: a" << std::endl;
	aPartDims[0].print(stream, indentLevel + 1);
	aPartDims[1].print(stream, indentLevel + 1);
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: b" << std::endl;
	bPartDims[0].print(stream, indentLevel + 1);
	bPartDims[1].print(stream, indentLevel + 1);
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: c" << std::endl;
	cPartDims[0].print(stream, indentLevel + 1);
	cPartDims[1].print(stream, indentLevel + 1);
	stream.flush();
}

void mm::SpaceA_Sub_LPU::print(std::ofstream &stream, int indentLevel) {
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: a" << std::endl;
	aPartDims[0].print(stream, indentLevel + 1);
	aPartDims[1].print(stream, indentLevel + 1);
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: b" << std::endl;
	bPartDims[0].print(stream, indentLevel + 1);
	bPartDims[1].print(stream, indentLevel + 1);
	for (int i = 0; i < indentLevel; i++) stream << '\t';
	stream << "Array: c" << std::endl;
	cPartDims[0].print(stream, indentLevel + 1);
	cPartDims[1].print(stream, indentLevel + 1);
	stream.flush();
}

/*-----------------------------------------------------------------------------------
Functions for ArrayMetadata and EnvironmentLinks 
------------------------------------------------------------------------------------*/

mm::ArrayMetadata::ArrayMetadata() : Metadata() {
	setTaskName("Matrix Multiply");
}

void mm::ArrayMetadata::print(std::ofstream &stream) {
	stream << "Array Metadata" << std::endl;
	stream << "Array: a";
	stream << ' ';
	aDims[0].print(stream);
	stream << ' ';
	aDims[1].print(stream);
	stream << std::endl;
	stream << "Array: b";
	stream << ' ';
	bDims[0].print(stream);
	stream << ' ';
	bDims[1].print(stream);
	stream << std::endl;
	stream << "Array: c";
	stream << ' ';
	cDims[0].print(stream);
	stream << ' ';
	cDims[1].print(stream);
	stream << std::endl;
	stream.flush();
}

/*-----------------------------------------------------------------------------------
function to initialize the content reference objects of LPSes
------------------------------------------------------------------------------------*/

void mm::initializeRootLPSContent(EnvironmentLinks *envLinks, ArrayMetadata *metadata) {
	spaceRootContent.a = envLinks->a;
	spaceRootContent.b = envLinks->b;
	spaceRootContent.c = allocate::allocateArray <double> (2, metadata->cDims);
	allocate::zeroFillArray <double> (0, spaceRootContent.c, 2, metadata->cDims);
}

void mm::initializeLPSesContents(ArrayMetadata *metadata) {
	//Processing Space A contents
	spaceAContent.c = spaceRootContent.c;
	//Processing Space A_Sub contents
	spaceA_SubContent.a = spaceRootContent.a;
	spaceA_SubContent.b = spaceRootContent.b;
	spaceA_SubContent.c = spaceRootContent.c;
}

/*-----------------------------------------------------------------------------------
functions for retrieving partition counts in different LPSes
------------------------------------------------------------------------------------*/

int *mm::getLPUsCountOfSpaceA(int ppuCount, Dimension cDim1, int k, Dimension cDim2, int l) {
	int *count = new int[2];
	count[0] = block_size_partitionCount(cDim1, ppuCount, k);
	count[1] = block_size_partitionCount(cDim2, ppuCount, l);
	return count;
}

int *mm::getLPUsCountOfSpaceA_Sub(int ppuCount, Dimension aDim2, int q) {
	int *count = new int[1];
	count[0] = block_size_partitionCount(aDim2, ppuCount, q);
	return count;
}

/*-----------------------------------------------------------------------------------
functions for getting data ranges along different dimensions of an LPU
-----------------------------------------------------------------------------------*/

void mm::getaPartForSpaceALpu(PartDimension *aLpuDims, 
		PartDimension *aParentLpuDims, 
		int *lpuCount, int *lpuId, int k) {
	aLpuDims[0].storage = aParentLpuDims[0].storage;
	aLpuDims[0].partition = block_size_getRange(aParentLpuDims[0].partition, 
			lpuCount[0], lpuId[0], false, k, 0, 0);
	aLpuDims[0].index = lpuId[0];
	aLpuDims[0].count = lpuCount[0];
	aLpuDims[0].parent = &aParentLpuDims[0];
	aLpuDims[1] = aParentLpuDims[1];
	aLpuDims[1].count = 1;
	aLpuDims[1].index = 0;
	aLpuDims[1].parent = &aParentLpuDims[1];
}

void mm::getbPartForSpaceALpu(PartDimension *bLpuDims, 
		PartDimension *bParentLpuDims, 
		int *lpuCount, int *lpuId, int l) {
	bLpuDims[0] = bParentLpuDims[0];
	bLpuDims[0].count = 1;
	bLpuDims[0].index = 0;
	bLpuDims[0].parent = &bParentLpuDims[0];
	bLpuDims[1].storage = bParentLpuDims[1].storage;
	bLpuDims[1].partition = block_size_getRange(bParentLpuDims[1].partition, 
			lpuCount[1], lpuId[1], false, l, 0, 0);
	bLpuDims[1].index = lpuId[1];
	bLpuDims[1].count = lpuCount[1];
	bLpuDims[1].parent = &bParentLpuDims[1];
}

void mm::getcPartForSpaceALpu(PartDimension *cLpuDims, 
		PartDimension *cParentLpuDims, 
		int *lpuCount, int *lpuId, int k, int l) {
	cLpuDims[0].storage = cParentLpuDims[0].storage;
	cLpuDims[0].partition = block_size_getRange(cParentLpuDims[0].partition, 
			lpuCount[0], lpuId[0], false, k, 0, 0);
	cLpuDims[0].index = lpuId[0];
	cLpuDims[0].count = lpuCount[0];
	cLpuDims[0].parent = &cParentLpuDims[0];
	cLpuDims[1].storage = cParentLpuDims[1].storage;
	cLpuDims[1].partition = block_size_getRange(cParentLpuDims[1].partition, 
			lpuCount[1], lpuId[1], false, l, 0, 0);
	cLpuDims[1].index = lpuId[1];
	cLpuDims[1].count = lpuCount[1];
	cLpuDims[1].parent = &cParentLpuDims[1];
}

void mm::getaPartForSpaceA_SubLpu(PartDimension *aLpuDims, 
		PartDimension *aParentLpuDims, 
		int *lpuCount, int *lpuId, int q) {
	aLpuDims[0] = aParentLpuDims[0];
	aLpuDims[0].count = 1;
	aLpuDims[0].index = 0;
	aLpuDims[0].parent = &aParentLpuDims[0];
	aLpuDims[1].storage = aParentLpuDims[1].storage;
	aLpuDims[1].partition = block_size_getRange(aParentLpuDims[1].partition, 
			lpuCount[0], lpuId[0], false, q, 0, 0);
	aLpuDims[1].index = lpuId[0];
	aLpuDims[1].count = lpuCount[0];
	aLpuDims[1].parent = &aParentLpuDims[1];
}

void mm::getbPartForSpaceA_SubLpu(PartDimension *bLpuDims, 
		PartDimension *bParentLpuDims, 
		int *lpuCount, int *lpuId, int q) {
	bLpuDims[0].storage = bParentLpuDims[0].storage;
	bLpuDims[0].partition = block_size_getRange(bParentLpuDims[0].partition, 
			lpuCount[0], lpuId[0], false, q, 0, 0);
	bLpuDims[0].index = lpuId[0];
	bLpuDims[0].count = lpuCount[0];
	bLpuDims[0].parent = &bParentLpuDims[0];
	bLpuDims[1] = bParentLpuDims[1];
	bLpuDims[1].count = 1;
	bLpuDims[1].index = 0;
	bLpuDims[1].parent = &bParentLpuDims[1];
}

/*-----------------------------------------------------------------------------------
function to generate PPU IDs and PPU group IDs for a thread
------------------------------------------------------------------------------------*/

ThreadIds *mm::getPpuIdsForThread(int threadNo)  {

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
	threadCount = Total_Threads;
	groupSize = threadCount / 128;
	groupThreadId = idsArray[Space_Root] % groupSize;
	threadIds->ppuIds[Space_A].groupId = idsArray[Space_Root] / groupSize;
	threadIds->ppuIds[Space_A].ppuCount = 128;
	threadIds->ppuIds[Space_A].groupSize = groupSize;
	if (groupThreadId == 0) threadIds->ppuIds[Space_A].id
			= threadIds->ppuIds[Space_A].groupId;
	else threadIds->ppuIds[Space_A].id = INVALID_ID;
	idsArray[Space_A] = groupThreadId;

	// for Space A_Sub;
	threadIds->ppuIds[Space_A_Sub].lpsName = "A_Sub";
	threadIds->ppuIds[Space_A_Sub].groupId = 0;
	threadIds->ppuIds[Space_A_Sub].ppuCount = 1;
	threadIds->ppuIds[Space_A_Sub].groupSize = threadIds->ppuIds[Space_A].groupSize;
	threadIds->ppuIds[Space_A_Sub].id = 0;
	idsArray[Space_A_Sub] = idsArray[Space_A];

	return threadIds;
}

/*-----------------------------------------------------------------------------------
Thread-State implementation class for the task
------------------------------------------------------------------------------------*/

// Construction of task specific LPS hierarchy index map
void ThreadStateImpl::setLpsParentIndexMap() {
	lpsParentIndexMap = new int[Space_Count];
	lpsParentIndexMap[Space_Root] = INVALID_ID;
	lpsParentIndexMap[Space_A] = Space_Root;
	lpsParentIndexMap[Space_A_Sub] = Space_A;
	//threadLog << "set up parent LPS index map" << std::endl;
	//threadLog.flush();
}

// Construction of task specific root LPU
void ThreadStateImpl::setRootLpu(Metadata *metadata) {

	ArrayMetadata *arrayMetadata = (ArrayMetadata*) metadata;

	SpaceRoot_LPU *lpu = new SpaceRoot_LPU;
	lpu->a = NULL;
	lpu->aPartDims[0] = PartDimension();
	lpu->aPartDims[0].partition = arrayMetadata->aDims[0];
	lpu->aPartDims[0].storage = arrayMetadata->aDims[0].getNormalizedDimension();
	lpu->aPartDims[1] = PartDimension();
	lpu->aPartDims[1].partition = arrayMetadata->aDims[1];
	lpu->aPartDims[1].storage = arrayMetadata->aDims[1].getNormalizedDimension();

	lpu->b = NULL;
	lpu->bPartDims[0] = PartDimension();
	lpu->bPartDims[0].partition = arrayMetadata->bDims[0];
	lpu->bPartDims[0].storage = arrayMetadata->bDims[0].getNormalizedDimension();
	lpu->bPartDims[1] = PartDimension();
	lpu->bPartDims[1].partition = arrayMetadata->bDims[1];
	lpu->bPartDims[1].storage = arrayMetadata->bDims[1].getNormalizedDimension();

	lpu->c = NULL;
	lpu->cPartDims[0] = PartDimension();
	lpu->cPartDims[0].partition = arrayMetadata->cDims[0];
	lpu->cPartDims[0].storage = arrayMetadata->cDims[0].getNormalizedDimension();
	lpu->cPartDims[1] = PartDimension();
	lpu->cPartDims[1].partition = arrayMetadata->cDims[1];
	lpu->cPartDims[1].storage = arrayMetadata->cDims[1].getNormalizedDimension();

	lpu->setValidBit(true);
	lpsStates[Space_Root]->lpu = lpu;
	//threadLog << "set up root LPU" << std::endl;
	//threadLog.flush();
}

// Setting up the Root LPU reference 
void ThreadStateImpl::setRootLpu(LPU *lpu) {
	lpu->setValidBit(true);
	lpsStates[Space_Root]->lpu = lpu;
	//threadLog << "set up root LPU" << std::endl;
	//threadLog.flush();
}

// Initialization of LPU pointers of different LPSes
void ThreadStateImpl::initializeLPUs() {
	lpsStates[Space_A]->lpu = new SpaceA_LPU;
	lpsStates[Space_A]->lpu->setValidBit(false);
	lpsStates[Space_A_Sub]->lpu = new SpaceA_Sub_LPU;
	lpsStates[Space_A_Sub]->lpu->setValidBit(false);
	//threadLog << "initialized LPU pointers" << std::endl;
	//threadLog.flush();
}

// Implementation of task specific compute-LPU-Count function 
int *ThreadStateImpl::computeLpuCounts(int lpsId) {
	if (lpsId == Space_Root) {
		return NULL;
	}
	if (lpsId == Space_A) {
		int ppuCount = threadIds->ppuIds[Space_A].ppuCount;
		SpaceRoot_LPU *spaceRootLpu
				 = (SpaceRoot_LPU*) lpsStates[Space_Root]->lpu;
		return getLPUsCountOfSpaceA(ppuCount, 
				spaceRootLpu->cPartDims[0].partition, 
				partitionArgs[0], 
				spaceRootLpu->cPartDims[1].partition, 
				partitionArgs[1]);
	}
	if (lpsId == Space_A_Sub) {
		int ppuCount = threadIds->ppuIds[Space_A_Sub].ppuCount;
		SpaceA_LPU *spaceALpu
				 = (SpaceA_LPU*) lpsStates[Space_A]->lpu;
		return getLPUsCountOfSpaceA_Sub(ppuCount, 
				spaceALpu->aPartDims[1].partition, 
				partitionArgs[2]);
	}
	return NULL;
}

// Implementation of task specific compute-Next-LPU function 
LPU *ThreadStateImpl::computeNextLpu(int lpsId, int *lpuCounts, int *nextLpuId) {
	if (lpsId == Space_A) {
		SpaceRoot_LPU *spaceRootLpu
				 = (SpaceRoot_LPU*) lpsStates[Space_Root]->lpu;
		SpaceA_LPU *currentLpu
				 = (SpaceA_LPU*) lpsStates[Space_A]->lpu;
		currentLpu->lpuId[0] = nextLpuId[0];
		currentLpu->lpuId[1] = nextLpuId[1];
		currentLpu->a = spaceAContent.a;
		getaPartForSpaceALpu(currentLpu->aPartDims, 
				spaceRootLpu->aPartDims, lpuCounts, nextLpuId, 
				partitionArgs[0]);
		currentLpu->b = spaceAContent.b;
		getbPartForSpaceALpu(currentLpu->bPartDims, 
				spaceRootLpu->bPartDims, lpuCounts, nextLpuId, 
				partitionArgs[1]);
		currentLpu->c = spaceAContent.c;
		getcPartForSpaceALpu(currentLpu->cPartDims, 
				spaceRootLpu->cPartDims, lpuCounts, nextLpuId, 
				partitionArgs[0], partitionArgs[1]);
		currentLpu->setValidBit(true);
		return currentLpu;
	}
	if (lpsId == Space_A_Sub) {
		SpaceA_LPU *spaceALpu
				 = (SpaceA_LPU*) lpsStates[Space_A]->lpu;
		SpaceA_Sub_LPU *currentLpu
				 = (SpaceA_Sub_LPU*) lpsStates[Space_A_Sub]->lpu;
		currentLpu->lpuId[0] = nextLpuId[0];
		currentLpu->a = spaceA_SubContent.a;
		getaPartForSpaceA_SubLpu(currentLpu->aPartDims, 
				spaceALpu->aPartDims, lpuCounts, nextLpuId, 
				partitionArgs[2]);
		currentLpu->b = spaceA_SubContent.b;
		getbPartForSpaceA_SubLpu(currentLpu->bPartDims, 
				spaceALpu->bPartDims, lpuCounts, nextLpuId, 
				partitionArgs[2]);
		currentLpu->c = spaceA_SubContent.c;
		currentLpu->cPartDims[0] = spaceALpu->cPartDims[0];
		currentLpu->cPartDims[1] = spaceALpu->cPartDims[1];
		currentLpu->setValidBit(true);
		return currentLpu;
	}
	return NULL;
}

/*-----------------------------------------------------------------------------------
function for initializing environment-links object
------------------------------------------------------------------------------------*/

EnvironmentLinks mm::initiateEnvLinks(MMEnvironment *environment) {
	EnvironmentLinks envLinks;
	envLinks.a = environment->a;
	envLinks.aDims[0] = environment->aDims[0].partition;
	envLinks.aDims[1] = environment->aDims[1].partition;
	envLinks.b = environment->b;
	envLinks.bDims[0] = environment->bDims[0].partition;
	envLinks.bDims[1] = environment->bDims[1].partition;
	return envLinks;
}

/*-----------------------------------------------------------------------------------
function for initializing root LPU from environment
------------------------------------------------------------------------------------*/

SpaceRoot_LPU *mm::initiateRootLpu(MMEnvironment *environment, ArrayMetadata *metadata) {

	SpaceRoot_LPU *rootLpu = new SpaceRoot_LPU;
	rootLpu->a = spaceRootContent.a;
	if (environment->a != NULL) {
		rootLpu->aPartDims[0] = environment->aDims[0];
		rootLpu->aPartDims[1] = environment->aDims[1];
	} else {
		rootLpu->aPartDims[0] = PartDimension();
		rootLpu->aPartDims[0].partition = metadata->aDims[0];
		rootLpu->aPartDims[0].storage = metadata->aDims[0].getNormalizedDimension();
		rootLpu->aPartDims[1] = PartDimension();
		rootLpu->aPartDims[1].partition = metadata->aDims[1];
		rootLpu->aPartDims[1].storage = metadata->aDims[1].getNormalizedDimension();
	}

	rootLpu->b = spaceRootContent.b;
	if (environment->b != NULL) {
		rootLpu->bPartDims[0] = environment->bDims[0];
		rootLpu->bPartDims[1] = environment->bDims[1];
	} else {
		rootLpu->bPartDims[0] = PartDimension();
		rootLpu->bPartDims[0].partition = metadata->bDims[0];
		rootLpu->bPartDims[0].storage = metadata->bDims[0].getNormalizedDimension();
		rootLpu->bPartDims[1] = PartDimension();
		rootLpu->bPartDims[1].partition = metadata->bDims[1];
		rootLpu->bPartDims[1].storage = metadata->bDims[1].getNormalizedDimension();
	}

	rootLpu->c = spaceRootContent.c;
	if (environment->c != NULL) {
		rootLpu->cPartDims[0] = environment->cDims[0];
		rootLpu->cPartDims[1] = environment->cDims[1];
	} else {
		rootLpu->cPartDims[0] = PartDimension();
		rootLpu->cPartDims[0].partition = metadata->cDims[0];
		rootLpu->cPartDims[0].storage = metadata->cDims[0].getNormalizedDimension();
		rootLpu->cPartDims[1] = PartDimension();
		rootLpu->cPartDims[1].partition = metadata->cDims[1];
		rootLpu->cPartDims[1].storage = metadata->cDims[1].getNormalizedDimension();
	}
	return rootLpu;
}

/*-----------------------------------------------------------------------------------
function for executing task
------------------------------------------------------------------------------------*/

void mm::execute(MMEnvironment *environment, 
		MMPartition partition, 
		std::ofstream &logFile) {

	// initializing environment-links object
	EnvironmentLinks envLinks = initiateEnvLinks(environment);

	// declaring other task related common variables
	TaskGlobals taskGlobals;
	ThreadLocals threadLocals;
	ArrayMetadata *metadata = new ArrayMetadata;

	// copying partitioning parameters into an array
	int *partitionArgs = NULL;
	partitionArgs = new int[3];
	partitionArgs[0] = partition.k;
	partitionArgs[1] = partition.l;
	partitionArgs[2] = partition.q;

	// invoking the initializer function
	//std::cout << "invoking task initializer function\n";
	initializeTask(metadata, envLinks, &taskGlobals, &threadLocals, partition);
	metadata->aDims[0].setLength();
	metadata->aDims[1].setLength();
	metadata->bDims[0].setLength();
	metadata->bDims[1].setLength();
	metadata->cDims[0].setLength();
	metadata->cDims[1].setLength();

	// allocating memories for data structures
	initializeRootLPSContent(&envLinks, metadata);
	initializeLPSesContents(metadata);

	// initializing the root LPU reference
	SpaceRoot_LPU *rootLpu = initiateRootLpu(environment, metadata);

	// declaring and initializing state variables for threads 
	ThreadLocals *threadLocalsList[Total_Threads];
	for (int i = 0; i < Total_Threads; i++) {
		threadLocalsList[i] = new ThreadLocals;
		*threadLocalsList[i] = threadLocals;
	}
	int lpsDimensions[Space_Count];
	lpsDimensions[Space_Root] = 0;
	lpsDimensions[Space_A] = 2;
	lpsDimensions[Space_A_Sub] = 1;
	//std::cout << "generating PPU Ids for threads\n";
	ThreadIds *threadIdsList[Total_Threads];
	for (int i = 0; i < Total_Threads; i++) {
		threadIdsList[i] = getPpuIdsForThread(i);
		threadIdsList[i]->print(logFile);
	}
	//std::cout << "initiating thread-states\n";
	ThreadStateImpl *threadStateList[Total_Threads];
	for (int i = 0; i < Total_Threads; i++) {
		threadStateList[i] = new ThreadStateImpl(Space_Count, 
				lpsDimensions, partitionArgs, threadIdsList[i]);
		threadStateList[i]->initiateLogFile("mm");
		threadStateList[i]->initializeLPUs();
		threadStateList[i]->setLpsParentIndexMap();
	}

	// setting up root LPU reference in each thread's state
	for (int i = 0; i < Total_Threads; i++) {
		threadStateList[i]->setRootLpu(rootLpu);
	}

	// starting threads
	//std::cout << "starting threads\n";
	pthread_t threads[Total_Threads];
	PThreadArg *threadArgs[Total_Threads];
	for (int i = 0; i < Total_Threads; i++) {
		threadArgs[i] = new PThreadArg;
		threadArgs[i]->taskName = "Matrix Multiply";
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
	for (int i = 0; i < Total_Threads; i++) {
		int cpuId = i * Core_Jump / Threads_Par_Core;
		int physicalId = Processor_Order[cpuId];
		state = pthread_create(&threads[i], &attr, runPThreads, (void *) threadArgs[i]);
		if (state) {
			std::cout << "Could not start some PThread" << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}
	for (int i = 0; i < Total_Threads; i++) {
		pthread_join(threads[i], NULL);
	}

	// copying results of task execution into environment
	environment->a = rootLpu->a;
	environment->aDims[0] = rootLpu->aPartDims[0];
	environment->aDims[1] = rootLpu->aPartDims[1];
	environment->b = rootLpu->b;
	environment->bDims[0] = rootLpu->bPartDims[0];
	environment->bDims[1] = rootLpu->bPartDims[1];
	environment->c = rootLpu->c;
	environment->cDims[0] = rootLpu->cPartDims[0];
	environment->cDims[1] = rootLpu->cPartDims[1];
}

/*-----------------------------------------------------------------------------------
function for the initialize block
------------------------------------------------------------------------------------*/

void mm::initializeTask(ArrayMetadata *arrayMetadata, 
		EnvironmentLinks environmentLinks, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		MMPartition partition) {

	arrayMetadata->aDims[0] = environmentLinks.aDims[0];
	arrayMetadata->aDims[1] = environmentLinks.aDims[1];
	arrayMetadata->bDims[0] = environmentLinks.bDims[0];
	arrayMetadata->bDims[1] = environmentLinks.bDims[1];
	arrayMetadata->cDims[0] = arrayMetadata->aDims[0];
	arrayMetadata->cDims[1] = arrayMetadata->bDims[1];
}

/*-----------------------------------------------------------------------------------
functions for compute stages 
------------------------------------------------------------------------------------*/

int mm::block_multiply_matrices(SpaceA_Sub_LPU *lpu, 
		ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, MMPartition partition) {

	//-------------------- Local Copies of Metadata -----------------------------

	// create local copies of partition and storage dimension configs of all arrays
	Dimension aPartDims[2];
	aPartDims[0] = lpu->aPartDims[0].partition;
	aPartDims[1] = lpu->aPartDims[1].partition;
	Dimension aStoreDims[2];
	aStoreDims[0] = lpu->aPartDims[0].storage;
	aStoreDims[1] = lpu->aPartDims[1].storage;
	Dimension bPartDims[2];
	bPartDims[0] = lpu->bPartDims[0].partition;
	bPartDims[1] = lpu->bPartDims[1].partition;
	Dimension bStoreDims[2];
	bStoreDims[0] = lpu->bPartDims[0].storage;
	bStoreDims[1] = lpu->bPartDims[1].storage;
	Dimension cPartDims[2];
	cPartDims[0] = lpu->cPartDims[0].partition;
	cPartDims[1] = lpu->cPartDims[1].partition;
	Dimension cStoreDims[2];
	cStoreDims[0] = lpu->cPartDims[0].storage;
	cStoreDims[1] = lpu->cPartDims[1].storage;

	// create a local part-dimension object for later use
	PartDimension partConfig;

	// create a local transformed index variable for later use
	int xformIndex;

	//----------------------- Computation Begins --------------------------------

	{// scope entrance for parallel loop on index i
	int i;
	int iterationStart = cPartDims[0].range.min;
	int iterationBound = cPartDims[0].range.max;
	int indexIncrement = 1;
	int indexMultiplier = 1;
	if (cPartDims[0].range.min > cPartDims[0].range.max) {
		iterationBound *= -1;
		indexIncrement *= -1;
		indexMultiplier = -1;
	}
	for (i = iterationStart; 
			indexMultiplier * i <= iterationBound; 
			i += indexIncrement) {
		long int ic0 = ((long) (i)) * ((long) (cStoreDims[1].length));
		long int ia0 = ((long) (i)) * ((long) (aStoreDims[1].length));
		{// scope entrance for parallel loop on index j
		int j;
		int iterationStart = cPartDims[1].range.min;
		int iterationBound = cPartDims[1].range.max;
		int indexIncrement = 1;
		int indexMultiplier = 1;
		if (cPartDims[1].range.min > cPartDims[1].range.max) {
			iterationBound *= -1;
			indexIncrement *= -1;
			indexMultiplier = -1;
		}
		for (j = iterationStart; 
				indexMultiplier * j <= iterationBound; 
				j += indexIncrement) {
			long int jc1 = ((long) (j));
			long int jb1 = ((long) (j));
			{// scope entrance for parallel loop on index k
			int k;
			int iterationStart = aPartDims[1].range.min;
			int iterationBound = aPartDims[1].range.max;
			int indexIncrement = 1;
			int indexMultiplier = 1;
			if (aPartDims[1].range.min > aPartDims[1].range.max) {
				iterationBound *= -1;
				indexIncrement *= -1;
				indexMultiplier = -1;
			}
			for (k = iterationStart; 
					indexMultiplier * k <= iterationBound; 
					k += indexIncrement) {
				long int ka1 = ((long) (k));
				long int kb0 = ((long) (k)) * ((long) (bStoreDims[1].length));
				lpu->c[ic0 + jc1] = (lpu->c[ic0 + jc1] + (lpu->a[ia0 + ka1] * lpu->b[kb0 + jb1]));
			}
			}// scope exit for parallel loop on index k
		}
		}// scope exit for parallel loop on index j
	}
	}// scope exit for parallel loop on index i

	//------------------------- Returning Flag ----------------------------------

	return SUCCESS_RUN;
}

/*-----------------------------------------------------------------------------------
The run method for thread simulating the task flow 
------------------------------------------------------------------------------------*/

void mm::run(ArrayMetadata *arrayMetadata, 
		TaskGlobals *taskGlobals, 
		ThreadLocals *threadLocals, 
		MMPartition partition, ThreadStateImpl *threadState) {

	// log thread's affinity information
	threadState->logThreadAffinity();

	// set the root LPU in the thread state so that calculation can start
	LPU *rootLpu = threadState->getCurrentLpu(Space_Root);
	if (rootLpu == NULL) {
		threadState->setRootLpu(arrayMetadata);
	}

	// create a local part-dimension object for later use
	PartDimension partConfig;

	// create a local transformed index variable for later use
	int xformIndex;

	{ // scope entrance for iterating LPUs of Space A_Sub
	int spaceA_SubLpuId = INVALID_ID;
	int spaceA_SubIteration = 0;
	SpaceA_Sub_LPU *spaceA_SubLpu = NULL;
	LPU *lpu = NULL;
	while((lpu = threadState->getNextLpu(Space_A_Sub, Space_Root, spaceA_SubLpuId)) != NULL) {
		spaceA_SubLpu = (SpaceA_Sub_LPU*) lpu;
		if (threadState->isValidPpu(Space_A_Sub)) {
			// invoking user computation
			int stage3Executed = block_multiply_matrices(spaceA_SubLpu, 
					arrayMetadata,
					taskGlobals,
					threadLocals, partition);
		}
		spaceA_SubLpuId = spaceA_SubLpu->id;
		spaceA_SubIteration++;
	}
	} // scope exit for iterating LPUs of Space A_Sub

	// close thread's log file
	threadState->closeLogFile();
}

/*-----------------------------------------------------------------------------------
PThreads run function
------------------------------------------------------------------------------------*/

void *mm::runPThreads(void *argument) {
	PThreadArg *pthreadArg = (PThreadArg *) argument;
	ThreadStateImpl *threadState = pthreadArg->threadState;
	//std::cout << "Thread " << threadState->getThreadNo() << " has started";
	//std::cout << " executing task: " << pthreadArg->taskName << std::endl;
	run(pthreadArg->metadata, 
			pthreadArg->taskGlobals, 
			pthreadArg->threadLocals, 
			pthreadArg->partition, 
			threadState);
	//std::cout << "Thread " << threadState->getThreadNo() << " has ended" << std::endl;
	pthread_exit(NULL);
}

/*-----------------------------------------------------------------------------------
main function
------------------------------------------------------------------------------------*/

int main() {

	std::cout << "Starting Matrix Multiply Task\n";

	// declaring common task related variables
	TaskGlobals taskGlobals;
	ThreadLocals threadLocals;
	EnvironmentLinks envLinks;
	ArrayMetadata *metadata = new ArrayMetadata;
	MMEnvironment environment;
	MMPartition partition;

	// creating a program log file
	std::cout << "Creating diagnostic log: it-program.log\n";
	std::ofstream logFile;
	logFile.open("it-program.log");

	// initializing variables that are environmental links 
	std::cout << "initializing environmental links\n";
	if (outprompt::getYesNoAnswer("Want to read array \"a\" from a file?")) {
		envLinks.a = inprompt::readArrayFromFile <double> ("a", 
				2, envLinks.aDims);
	} else {
		inprompt::readArrayDimensionInfo("a", 2, envLinks.aDims);
		envLinks.a = allocate::allocateArray <double> (2, envLinks.aDims);
		allocate::randomFillPrimitiveArray <double> (envLinks.a, 
				2, envLinks.aDims);
	}
	if (outprompt::getYesNoAnswer("Want to read array \"b\" from a file?")) {
		envLinks.b = inprompt::readArrayFromFile <double> ("b", 
				2, envLinks.bDims);
	} else {
		inprompt::readArrayDimensionInfo("b", 2, envLinks.bDims);
		envLinks.b = allocate::allocateArray <double> (2, envLinks.bDims);
		allocate::randomFillPrimitiveArray <double> (envLinks.b, 
				2, envLinks.bDims);
	}

	// determining values of partition parameters
	std::cout << "determining partition parameters\n";
	int *partitionArgs = NULL;
	partitionArgs = new int[3];
	partition.k = inprompt::readPrimitive <int> ("k");
	partitionArgs[0] = partition.k;
	partition.l = inprompt::readPrimitive <int> ("l");
	partitionArgs[1] = partition.l;
	partition.q = inprompt::readPrimitive <int> ("q");
	partitionArgs[2] = partition.q;

	// determining values of initialization parameters
	std::cout << "determining initialization parameters\n";

	// invoking the initializer function
	//std::cout << "invoking task initializer function\n";
	initializeTask(metadata, envLinks, &taskGlobals, &threadLocals, partition);
	metadata->aDims[0].setLength();
	metadata->aDims[1].setLength();
	metadata->bDims[0].setLength();
	metadata->bDims[1].setLength();
	metadata->cDims[0].setLength();
	metadata->cDims[1].setLength();

	// setting the global metadata variable
	arrayMetadata = *metadata;
	metadata->print(logFile);

	// allocating memories for data structures
	std::cout << "Allocating memories\n";
	mm::initializeRootLPSContent(&envLinks, metadata);
	mm::initializeLPSesContents(metadata);

	// declaring and initializing state variables for threads 
	ThreadLocals *threadLocalsList[Total_Threads];
	for (int i = 0; i < Total_Threads; i++) {
		threadLocalsList[i] = new ThreadLocals;
		*threadLocalsList[i] = threadLocals;
	}
	int lpsDimensions[Space_Count];
	lpsDimensions[Space_Root] = 0;
	lpsDimensions[Space_A] = 2;
	lpsDimensions[Space_A_Sub] = 1;
	//std::cout << "generating PPU Ids for threads\n";
	ThreadIds *threadIdsList[Total_Threads];
	for (int i = 0; i < Total_Threads; i++) {
		threadIdsList[i] = getPpuIdsForThread(i);
		threadIdsList[i]->print(logFile);
	}
	//std::cout << "initiating thread-states\n";
	ThreadStateImpl *threadStateList[Total_Threads];
	for (int i = 0; i < Total_Threads; i++) {
		threadStateList[i] = new ThreadStateImpl(Space_Count, 
				lpsDimensions, partitionArgs, threadIdsList[i]);
		threadStateList[i]->initiateLogFile("mm");
		threadStateList[i]->initializeLPUs();
		threadStateList[i]->setLpsParentIndexMap();
	}

	// starting execution timer clock
	struct timeval start;
	gettimeofday(&start, NULL);

	// starting threads
	//std::cout << "starting threads\n";
	pthread_t threads[Total_Threads];
	PThreadArg *threadArgs[Total_Threads];
	for (int i = 0; i < Total_Threads; i++) {
		threadArgs[i] = new PThreadArg;
		threadArgs[i]->taskName = "Matrix Multiply";
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
	for (int i = 0; i < Total_Threads; i++) {
		int cpuId = i * Core_Jump / Threads_Par_Core;
		int physicalId = Processor_Order[cpuId];
		state = pthread_create(&threads[i], &attr, runPThreads, (void *) threadArgs[i]);
		if (state) {
			std::cout << "Could not start some PThread" << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}
	for (int i = 0; i < Total_Threads; i++) {
		pthread_join(threads[i], NULL);
	}

	// calculating task running time
	struct timeval end;
	gettimeofday(&end, NULL);
	double runningTime = ((end.tv_sec + end.tv_usec / 1000000.0)
			- (start.tv_sec + start.tv_usec / 1000000.0));
	logFile << "Execution Time: " << runningTime << " Seconds" << std::endl;

	// writing environment variables to files after task completion
	std::cout << "writing results to output files\n";
	if (outprompt::getYesNoAnswer("Want to save array \"a\" in a file?")) {
		outprompt::writeArrayToFile <double> ("a", 
				spaceRootContent.a, 
				2, metadata->aDims);
	}
	if (outprompt::getYesNoAnswer("Want to save array \"b\" in a file?")) {
		outprompt::writeArrayToFile <double> ("b", 
				spaceRootContent.b, 
				2, metadata->bDims);
	}
	if (outprompt::getYesNoAnswer("Want to save array \"c\" in a file?")) {
		outprompt::writeArrayToFile <double> ("c", 
				spaceRootContent.c, 
				2, metadata->cDims);
	}

	logFile.close();
	std::cout << "Parallel Execution Time: " << runningTime << " Seconds" << std::endl;
	return 0;
}
