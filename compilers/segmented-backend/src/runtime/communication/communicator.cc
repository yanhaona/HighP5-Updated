#include "comm_buffer.h"
#include "comm_statistics.h"
#include "communicator.h"
#include "comm_barrier.h"

#include "../../../../common-libs/utils/list.h"

#include <iostream>
#include <cstdlib>
#include <sstream>
#include <time.h>
#include <sys/time.h>

//-------------------------------------------------------------- Send Barrier ------------------------------------------------------------/

SendBarrier::SendBarrier(int participantCount, Communicator *communicator) 
		: ParallelCommBarrier(participantCount) {
	this->communicator = communicator;
}

bool SendBarrier::shouldWait(SignalType signal, int callerIterationNo) {
	return communicator->shouldWaitOnSend(signal, callerIterationNo);
}

bool SendBarrier::shouldPerformTransfer(int activeSignalsCount, int callerIterationNo) {
	return communicator->shouldSend(activeSignalsCount);
}

void SendBarrier::configureCache() {
	communicator->cacheSendBuffers();
}

void SendBarrier::beforeTransfer(int order, int participants) {
	communicator->performSendPreprocessing(order, participants);
}

void SendBarrier::transferFunction() {
	communicator->sendData();
	communicator->afterSend();
}
        
void SendBarrier::afterTransfer(int order, int participants) {
	communicator->performSendPostprocessing(order, participants);
}

bool SendBarrier::supportSingleStepTransfer() {
	return communicator->directCommunicationPossible();
}

void SendBarrier::doSingleStepTransfer(int order, int participants) {

	// let all participants join in the parallel send operation
	communicator->performDirectSend(order, participants);
	if (order == 0) {
		// let the PPU with lowest ID to do any after send reconfiguration of the communicator
		communicator->afterSend();
	}
}


void SendBarrier::recordTimingLog(TimingLogType logType, struct timeval &start, struct timeval &end) {
	CommStatistics *commStat = communicator->getCommStat();
	if (logType == BEFORE_TRANSFER_TIMING) {
		commStat->addBufferReadTime(communicator->getName(), start, end);
	} else if (logType == TRANSFER_TIMING) {
		commStat->addCommunicationTime(communicator->getName(), start, end);
	} else if (logType == AFTER_TRANSFER_TIMING) {
		commStat->addBufferWriteTime(communicator->getName(), start, end);
	}
}

void SendBarrier::executeSend() {
	struct timeval start;
        gettimeofday(&start, NULL);
	communicator->prepareBuffersForSend();
	communicator->sendData();
	communicator->afterSend();
	struct timeval end;
        gettimeofday(&end, NULL);
	CommStatistics *commStat = communicator->getCommStat();
	commStat->addCommunicationTime(communicator->getName(), start, end);
}

//------------------------------------------------------------- Receive Barrier ----------------------------------------------------------/

ReceiveBarrier::ReceiveBarrier(int participantCount, Communicator *communicator) 
		: ParallelCommBarrier(participantCount) {
	this->communicator = communicator;
}

bool ReceiveBarrier::shouldWait(SignalType signal, int callerIterationNo) {
	return communicator->shouldWaitOnReceive(signal, callerIterationNo);
}

bool ReceiveBarrier::shouldPerformTransfer(int activeSignalsCount, int callerIterationNo) {
	return communicator->shouldReceive(activeSignalsCount, callerIterationNo);
}

void ReceiveBarrier::configureCache() {
	communicator->cacheRecvBuffers();
}

void ReceiveBarrier::beforeTransfer(int order, int participants) {
	communicator->perfromRecvPreprocessing(order, participants);
}

void ReceiveBarrier::transferFunction() {
	communicator->receiveData();
	communicator->afterReceive();
}

void ReceiveBarrier::afterTransfer(int order, int participants) {
	communicator->perfromRecvPostprocessing(order, participants);
}

bool ReceiveBarrier::supportSingleStepTransfer() {
	return communicator->directCommunicationPossible();
}

void ReceiveBarrier::doSingleStepTransfer(int order, int participants) {
	
	// let all participants join in the parallel receive operation
	communicator->performDirectReceive(order, participants);
	if (order == 0) {
		// let the PPU with smallest ID to do any after receive processing of the communicator
		communicator->afterReceive();
	}
}

void ReceiveBarrier::recordTimingLog(TimingLogType logType, struct timeval &start, struct timeval &end) {
	CommStatistics *commStat = communicator->getCommStat();
	if (logType == BEFORE_TRANSFER_TIMING) {
		commStat->addBufferReadTime(communicator->getName(), start, end);
	} else if (logType == TRANSFER_TIMING) {
		commStat->addCommunicationTime(communicator->getName(), start, end);
	} else if (logType == AFTER_TRANSFER_TIMING) {
		commStat->addBufferWriteTime(communicator->getName(), start, end);
	}
}

void ReceiveBarrier::executeReceive() {
	struct timeval start;
        gettimeofday(&start, NULL);
	communicator->receiveData();
	communicator->processBuffersAfterReceive();
	communicator->afterReceive();
	struct timeval end;
        gettimeofday(&end, NULL);
	CommStatistics *commStat = communicator->getCommStat();
	commStat->addCommunicationTime(communicator->getName(), start, end);
}

//-------------------------------------------------------------- Communicator ------------------------------------------------------------/

Communicator::Communicator(int localSegmentTag, 
		const char *dependencyName, 
		int localSenderPpus, int localReceiverPpus) : CommBufferManager(dependencyName) {
	
	this->localSegmentTag = localSegmentTag;

	if (localSenderPpus > 0) {
		sendBarrier = new SendBarrier(localSenderPpus, this);
	} else {
		sendBarrier = NULL;
	}
	if (localReceiverPpus > 0) {
		receiveBarrier = new ReceiveBarrier(localReceiverPpus, this);
	} else {
		receiveBarrier = NULL;
	}
	iterationNo = 0;
	communicatorId = 0;
	commStat = NULL;
	this->localSenderPpus = localSenderPpus;
	this->localReceiverPpus = localReceiverPpus;
	cachedSendBuffers = NULL;
	cachedRecvBuffers = NULL;
}

void Communicator::describe(int indentation) {
	std::ostringstream indent;
	for (int i = 0; i < indentation; i++) indent << '\t';
	*logFile << indent.str() << "Communicator for " << dependencyName << ":\n";
	*logFile << indent.str() << '\t' << "Total Communication Buffers: ";
	*logFile << commBufferList->NumElements() << "\n";
	for (int i = 0; i < commBufferList->NumElements(); i++) {
		*logFile << indent.str() << "\tBuffer #" << i + 1 << ":\n";
		CommBuffer *buffer = commBufferList->Nth(i);
		buffer->describe(*logFile, indentation + 1);
	}
}

void Communicator::cacheSendBuffers() {
	cachedSendBuffers = getSortedList(false);
}

void Communicator::prepareBuffersForSend() {
        if (cachedSendBuffers == NULL) return;
        for (int i = 0; i < cachedSendBuffers->NumElements(); i++) {
                cachedSendBuffers->Nth(i)->readData(false, *logFile);
        }
}

void Communicator::cacheRecvBuffers() {
	cachedRecvBuffers = getSortedList(true);
}

void Communicator::processBuffersAfterReceive() {
        if (cachedRecvBuffers == NULL) return;
        for (int i = 0; i < cachedRecvBuffers->NumElements(); i++) {
                cachedRecvBuffers->Nth(i)->writeData(false, *logFile);
        }
}

void Communicator::prepareBuffersForSend(int currentPpuOrder, int participantsCount) {
        if (cachedSendBuffers == NULL) return;
        for (int i = currentPpuOrder; i < cachedSendBuffers->NumElements(); i += participantsCount) {
                cachedSendBuffers->Nth(i)->readData(false, *logFile);
        }
}
        
void Communicator::processBuffersAfterReceive(int currentPpuOrder, int participantsCount) {
        if (cachedRecvBuffers == NULL) return;
        for (int i = currentPpuOrder; i < cachedRecvBuffers->NumElements(); i += participantsCount) {
                cachedRecvBuffers->Nth(i)->writeData(false, *logFile);
        }
}

void Communicator::setupBufferTags(int communicatorId, int totalSegmentsInMachine) {
	this->communicatorId = communicatorId;
	std::ostringstream digitStr;
	digitStr << totalSegmentsInMachine;
	int digitsForSegmentId = digitStr.str().length();
	for (int i = 0; i < commBufferList->NumElements(); i++) {
		CommBuffer *buffer = commBufferList->Nth(i);
		buffer->setBufferTag(communicatorId, digitsForSegmentId);
	}
}

void Communicator::setupCommunicator(bool includeNonInteractingSegments) {
	
	*logFile << "\tSetting up communicator for " << dependencyName << "\n";
	logFile->flush();
	
	struct timeval start;
        gettimeofday(&start, NULL);
	if (includeNonInteractingSegments) {
        	segmentGroup = new SegmentGroup(*participantSegments);
	} else {
		std::vector<int> *interactingParticipants = getParticipantsTags();
        	segmentGroup = new SegmentGroup(*interactingParticipants);
		delete interactingParticipants;
	}
        segmentGroup->setupCommunicator(*logFile);
	struct timeval end;
        gettimeofday(&end, NULL);
	commStat->addCommResourcesSetupTime(dependencyName, start, end);

	*logFile << "\tSetup done for communicator for " << dependencyName << "\n";
	logFile->flush();
}

void Communicator::excludeOwnselfFromCommunication(const char *dependencyName, 
		int localSegmentTag, std::ofstream &logFile) {
	logFile << "\tExcluding myself from dependency " << dependencyName << "\n";
	logFile.flush();
	SegmentGroup::excludeSegmentFromGroupSetup(localSegmentTag, logFile);
	logFile << "\tExcluded myself from dependency " << dependencyName << "\n";
	logFile.flush();
}

void Communicator::performDirectSend(int currentPpuOrder, int participantsCount) {
	std::cout << "Direct send is not supported\n";
	std::exit(EXIT_FAILURE);
}

void Communicator::performDirectReceive(int currentPpuOrder, int participantsCount) {
	std::cout << "Direct receive is not supported\n";
	std::exit(EXIT_FAILURE);
}
