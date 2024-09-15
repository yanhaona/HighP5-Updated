#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include "data_transfer.h"
#include "part_config.h"
#include "../memory-management/allocation.h"

#include "../../../../common-libs/utils/utility.h"
#include "../../../../common-libs/domain-obj/structure.h"

using namespace std;

//----------------------------------------------------- Data Part Index List ------------------------------------------------------/

void DataPartIndexList::clone(DataPartIndexList *source) {
	this->partIndexList->clear();
	this->partIndexList->AppendAll(source->partIndexList);
}

void DataPartIndexList::clonePartIndexList(List<DataPartIndex> *sourcePartIndexList) {
	this->partIndexList->clear();
	this->partIndexList->AppendAll(sourcePartIndexList);
}

int DataPartIndexList::read(char *destBuffer, int elementSize) {
	
	// Notice that data is read from only the first matched location in the operating memory data part storage. 
	// This is because, all locations matching a single entry in the communication buffer should have identical 
	// content.
	char *readLocation = partIndexList->Nth(0).getLocation();
        memcpy(destBuffer, readLocation, elementSize);
	return 1;
}
        
int DataPartIndexList::write(char *sourceBuffer, int elementSize) {

	// Unlike in the case for read, writing should access every single matched location in the data parts as we 
	// need to ensure that all locations matching a single entry in the communication buffer are synchronized 
	// with the same update.       
	for (int j = 0; j < partIndexList->NumElements(); j++) {
		char *writeLocation = partIndexList->Nth(j).getLocation();
		memcpy(writeLocation, sourceBuffer, elementSize);
	}

	// The return value is still once as only one element has been read from the source buffer
	return 1;
}

//-------------------------------------------------- Data Part Swift Index List ---------------------------------------------------/

DataPartSwiftIndexList::DataPartSwiftIndexList(DataPart *dataPart) : DataPartIndexList() {
	this->dataPart = dataPart;
	partIndexes = new List<long int>;
	optimizationAttempted = false;
	indexRanges = NULL;
} 

DataPartSwiftIndexList::~DataPartSwiftIndexList() {
	delete partIndexes;
	delete[] indexArray;
	delete[] indexRanges;
}

void DataPartSwiftIndexList::setupIndexArray() {

	int indexCount = partIndexes->NumElements();
       	indexArray = new long int[indexCount];
       	for (int i = 0; i < indexCount; i++) {
               indexArray[i] = partIndexes->Nth(i);
	}
	sequenceLength = indexCount;
}


void DataPartSwiftIndexList::optimizeIndexArray(int elementSize) {

	assert(indexRanges == NULL);

	totalIndices = sequenceLength;
	sequenceLength = 0;
	long int lastIndex = -30; // a negative value to account for the first index being 0
	for (int i = 0; i < totalIndices; i++) {
		long int currIndex = indexArray[i];
		if (currIndex != lastIndex + elementSize) {
			sequenceLength++; // if the current index is not the immediate next of the last
					  // then we have to resume memcopy from this new index again.
					  // So the length of the indexArray and indexRanges array should
					  // increase by 1.
		}
		lastIndex = currIndex;
	}

	long int *newIndexArray = new long int[sequenceLength];
	indexRanges = new int[sequenceLength];
	int indexesInCurrRange = 0;
	lastIndex = -30; // a negative value to account for the first index being 0
	int count = 0;
	for (int i = 0; i < totalIndices; i++) {
		long int currIndex = indexArray[i];
		if (currIndex != lastIndex + elementSize) {
			newIndexArray[count] = currIndex; // new index is a jump starting point
			if (count > 0) {
				// set the range of consecutive indices from the previous jump start
				indexRanges[count - 1] = indexesInCurrRange;
			}
			indexesInCurrRange = 1; // reset the consecutive index count
			count++; // advance the jump start point counters
		} else {
			indexesInCurrRange++; // increment the consecutive index range by on
		}
		lastIndex = currIndex;
	}
	// put the last index range value in the index range array
	indexRanges[count - 1] = indexesInCurrRange;

	// delete the old index array and attach the new one;
	delete[] indexArray;
	indexArray = newIndexArray;
}

int DataPartSwiftIndexList::read(char *destBuffer, int elementSize) {

	if (optimizationAttempted == false) {
		optimizeIndexArray(elementSize);
		optimizationAttempted = true;	
	}
	
	void *data = dataPart->getData();
        char *charData = reinterpret_cast<char*>(data);
	char *currBufferIndex = destBuffer;
	for (int i = 0; i < sequenceLength; i++) {
		char *readLocation = charData + indexArray[i];
		int consecutiveIndices = indexRanges[i];
		int copyVolume = elementSize * consecutiveIndices;
		memcpy(currBufferIndex, readLocation, copyVolume);
		currBufferIndex += copyVolume;
	}
	return totalIndices;
}
        
int DataPartSwiftIndexList::write(char *sourceBuffer, int elementSize) {
	
	if (optimizationAttempted == false) {
		optimizeIndexArray(elementSize);
		optimizationAttempted = true;	
	}

	void *data = dataPart->getData();
        char *charData = reinterpret_cast<char*>(data);
	char *currBufferIndex = sourceBuffer;
	for (int i = 0; i < sequenceLength; i++) {
		char *writeLocation = charData + indexArray[i];
		int consecutiveIndices = indexRanges[i];
		int copyVolume = elementSize * consecutiveIndices;
		memcpy(writeLocation, currBufferIndex, copyVolume);
		currBufferIndex += copyVolume;
	}
	
	return totalIndices;
}

//---------------------------------------------------- Transfer Specification -----------------------------------------------------/

TransferSpec::TransferSpec(TransferDirection direction, int elementSize) {
	this->direction = direction;
	this->elementSize = elementSize;
	this->bufferEntry = NULL;
	this->dataIndex = NULL;
	this->confinementContainerId = NULL;
}

void TransferSpec::setBufferEntry(char *bufferEntry, vector<int> *dataIndex) {
	this->bufferEntry = bufferEntry;
	this->dataIndex = dataIndex;
}

void TransferSpec::performTransfer(DataPartIndex dataPartIndex) {
	char *dataPartLocation = dataPartIndex.getLocation();
	if (direction == COMM_BUFFER_TO_DATA_PART) {
		memcpy(dataPartLocation, bufferEntry, elementSize);
	} else {
		memcpy(bufferEntry, dataPartLocation, elementSize);
	}
}

bool TransferSpec::isIncludedInTransfer(int partNo, int idDimension, int partNoIdLevel, int indexInLevel) {
	if (confinementContainerId == NULL) return true;
	int vectorSize = confinementContainerId->size();
	if (partNoIdLevel >= vectorSize) return true;
	int *idAtLevel = confinementContainerId->at(partNoIdLevel);
	return idAtLevel[indexInLevel] == partNo;
}

//--------------------------------------------------- Data Part Specification -----------------------------------------------------/

DataPartSpec::DataPartSpec(List<DataPart*> *partList, DataItemConfig *dataConfig) {
	this->partList = partList;
	this->dataConfig = dataConfig;
	dimensionality = dataConfig->getDimensionality();
	dataDimensions = new Dimension[dimensionality];
	for (int i = 0; i < dimensionality; i++) {
		dataDimensions[i] = dataConfig->getDimension(i);
	}
}

void DataPartSpec::initPartTraversalReference(vector<int> *dataIndex, vector<XformedIndexInfo*> *transformVector) {
	for (int i = 0; i < dimensionality; i++) {
		XformedIndexInfo *dimIndex = transformVector->at(i);
		dimIndex->index = dataIndex->at(i);
		dimIndex->partNo = 0;
		dimIndex->partDimension = dataDimensions[i];
	}
}

char *DataPartSpec::getUpdateLocation(PartLocator *partLocator, vector<int> *partIndex, int dataItemSize) {

	DataPartIndex dataPartIndex = getDataPartUpdateIndex(partLocator, partIndex, dataItemSize);
	return dataPartIndex.getLocation();
}

DataPartIndex DataPartSpec::getDataPartUpdateIndex(PartLocator *partLocator, 
		vector<int> *partIndex, int dataItemSize) {

	int partNo = partLocator->getPartListIndex();
        DataPart *dataPart = partList->Nth(partNo);
        PartMetadata *metadata = dataPart->getMetadata();
        Dimension *partDimensions = metadata->getBoundary();

        long int dataPointNo = 0;
        long int multiplier = 1;
        for (int i = partIndex->size() - 1; i >= 0; i--) {

                int firstIndex = partDimensions[i].range.min;
                int lastIndex = partDimensions[i].range.max;
                int dimensionIndex = partIndex->at(i);

                Assert(firstIndex <= dimensionIndex && dimensionIndex <= lastIndex);

                dataPointNo += (dimensionIndex - firstIndex) * multiplier;
                multiplier *= partDimensions[i].length;
        }

        Assert(dataPointNo < metadata->getSize());
	DataPartIndex dataPartIndex = DataPartIndex(dataPart, dataItemSize * dataPointNo);
	return dataPartIndex;
}
