#include "../../../utils/code_constant.h"
#include "../../../../../../common-libs/utils/list.h"
#include "../../../../../../frontend/src/semantics/computation_flow.h"
#include "../../../../../../frontend/src/static-analysis/sync_stat.h"
#include "../../../../../../frontend/src/static-analysis/data_dependency.h"

#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>

List<const char*> *CompositeStage::getAllOutgoingDependencyNamesAtNestingLevel(int nestingLevel) {

        List<const char*> *arcNameList = new List<const char*>;
        List<const char*> *ownList = FlowStage::getAllOutgoingDependencyNamesAtNestingLevel(nestingLevel);
        if (ownList != NULL) {
                arcNameList->AppendAll(ownList);
        }

        for (int i = 0; i < stageList->NumElements(); i++) {
                FlowStage *stage = stageList->Nth(i);
                List<const char*> *nestedList = stage->getAllOutgoingDependencyNamesAtNestingLevel(nestingLevel);
                if (nestedList != NULL) {
                        for (int j = 0; j < nestedList->NumElements(); j++) {
                                const char *currentArcName = nestedList->Nth(j);
                                bool found = false;
                                for (int k = 0; k < arcNameList->NumElements(); k++) {
                                        if (strcmp(arcNameList->Nth(k), currentArcName) == 0) {
                                                found = true;
                                                break;
                                        }
                                }
                                if (!found) {
                                        arcNameList->Append(currentArcName);
                                }
                        }
                }
        }
        return arcNameList;
}

void CompositeStage::declareSynchronizationCounters(std::ofstream &stream, int indentation, int nestingIndex) {

        std::ostringstream indentStream;
        for (int i = 0; i < indentation; i++) indentStream << indent;
	std::string indentStr = indentStream.str();

        List<const char*> *counterNameList = getAllOutgoingDependencyNamesAtNestingLevel(nestingIndex);
        if (counterNameList != NULL && counterNameList->NumElements() > 0) {
                stream << std::endl << indentStr << "// declaration of synchronization counter variables\n";
                for (int i = 0; i < counterNameList->NumElements(); i++) {
                        const char *counter = counterNameList->Nth(i);
                        stream << indentStr << "int " << counter << " = 0";
                        stream << stmtSeparator;
                }
        }
}

void CompositeStage::generateDataReceivesForGroup(std::ofstream &stream, int indentation,
                        List<SyncRequirement*> *commDependencies) {
	
	// Note that all synchronization we are dealing here is strictly within the boundary of composite stage under
	// concern. Now if there is a synchronization dependency to a nested stage for the execution of another stage
	// that comes after it, it means that the dependency is between a latter iteration on the dependent stage on
	// an earlier iteration on the source stage. This is because, otherwise the dependent stage will execute before
	// the source stage and there should not be any dependency. Now, such dependency is only valid for iterations 
	// except the first one. So there should be a checking on iteration number before we apply waiting on such
	// dependencies. On the other hand, if the source stage executes earlier than the dependent stage then it applies
	// always, i.e., it is independent of any loop iteration, if exists, that sorrounds the composite stage.
	// Therefore, we use following two lists to partition the list of communication requirements into forward and
	// backward dependencies. 
	List<SyncRequirement*> *forwardDependencies = new List<SyncRequirement*>;
	List<SyncRequirement*> *backwardDependencies = new List<SyncRequirement*>;

	for (int i = 0; i < commDependencies->NumElements(); i++) {
		SyncRequirement *comm = commDependencies->Nth(i);
		if (!comm->isActive()) continue;

		// if the current sync requirement has been expanded to be replaced by a subsequent sync that goes even
		// deeper in the LPS hierarchy then do a waiting for the expanded sync requirement and deactivate the
		// expanded sync to avoid a subsequent receive on the same communication
		//
		// NOTE: this logic should be covered by the comm-replacement sync but it is not working. We should 
		// check if signal expansion and communication replacement can be combined.
		if (comm->isSynchingExpanded()) {
			comm = comm->getExpandedSync();
			comm->getDependencyArc()->deactivate();
		}

		DependencyArc *arc = comm->getDependencyArc();
		FlowStage *source = arc->getSource();
		FlowStage *destination = arc->getDestination();
		if (source->getIndex() < destination->getIndex()) {
			forwardDependencies->Append(comm);
		} else backwardDependencies->Append(comm);
	}

	std::ostringstream indentStr;
	for (int i = 0; i < indentation; i++) indentStr << indent;
	
	if (commDependencies->NumElements() > 0) {
		stream << std::endl << indentStr.str() << "// waiting on data reception\n";
	}
	
	// Write the code for forwards communication requirements.
	for (int i = 0; i < forwardDependencies->NumElements(); i++) {
		
		SyncRequirement *comm = forwardDependencies->Nth(i);

		// if the data send for this communicator has been replaced with some earlier communicator then the data 
		// receive is done using that earlier communicator too 
		SyncRequirement *commReplacement = comm->getReplacementSync();
		bool signalReplaced = false;
		if (commReplacement != NULL) {
			comm = commReplacement;
			signalReplaced = true;
		}	

		int commIndex = comm->getIndex();
		Space *dependentLps = comm->getDependentLps();
		stream << indentStr.str() << "if (threadState->isValidPpu(Space_" << dependentLps->getName();
		stream << ")) {\n";
		stream << indentStr.str() << indent << "Communicator *communicator = threadState->getCommunicator(\"";
		stream << comm->getDependencyArc()->getArcName() << "\")" << stmtSeparator;
		stream << indentStr.str() << indent << "if (communicator != NULL) {\n";
		stream << indentStr.str() << doubleIndent << "communicator->receive(REQUESTING_COMMUNICATION";
		stream << paramSeparator << "commCounter" << commIndex << ")" << stmtSeparator; 
		stream << indentStr.str() << indent << "}\n";
		stream << indentStr.str() << "}\n";

		// The counter should be advanced regardless of this PPU's participation in communication to keep the 
		// counter value uniform across all PPUs and segments. The exeception is when the data send signal has 
		// been replaced by some other communicator. This is because, we will have two receive calls for the 
		// replacement communicator then and we would want to make the second call to bypass any data processing. 
		if (!signalReplaced) {
			stream << indentStr.str() << "commCounter" << commIndex << "++" << stmtSeparator;
		}
	}

	// write the code for backword sync requirements within an if block
	if (backwardDependencies->NumElements() > 0) {
		stream << indentStr.str() << "if (repeatIteration > 0) {\n";
		for (int i = 0; i < backwardDependencies->NumElements(); i++) {
			
			SyncRequirement *comm = backwardDependencies->Nth(i);
			
			SyncRequirement *commReplacement = comm->getReplacementSync();
			bool signalReplaced = false;
			if (commReplacement != NULL) {
				comm = commReplacement;
				signalReplaced = true;
			}	

			int commIndex = comm->getIndex();
			Space *dependentLps = comm->getDependentLps();
			stream << indentStr.str() << indent; 
			stream << "if (threadState->isValidPpu(Space_" << dependentLps->getName();
			stream << ")) {\n";
			stream << indentStr.str() << doubleIndent;
			stream << "Communicator *communicator = threadState->getCommunicator(\"";
			stream << comm->getDependencyArc()->getArcName() << "\")" << stmtSeparator;
			stream << indentStr.str() << doubleIndent << "if (communicator != NULL) {\n";
			stream << indentStr.str() << tripleIndent;
			stream << "communicator->receive(REQUESTING_COMMUNICATION";
			stream << paramSeparator << "commCounter" << commIndex << ")" << stmtSeparator; 
			stream << indentStr.str() << doubleIndent << "}\n";
			stream << indentStr.str() << indent << "}\n";

			// the counter is advanced outside the if condition to keep it in sync with all other PPUs
			if (!signalReplaced) {
				stream << indentStr.str() << indent << "commCounter";
				stream << commIndex << "++" << stmtSeparator;
			}
		}
		stream << indentStr.str() << "}\n";
	}

	// finally deactive all sync dependencies as they are already been taken care of here	
	for (int i = 0; i < commDependencies->NumElements(); i++) {
		SyncRequirement *comm = commDependencies->Nth(i);
		comm->deactivate();
	}
}

void CompositeStage::genSimplifiedWaitingForReactivationCode(std::ofstream &stream, int indentation,
		List<SyncRequirement*> *syncRequirements, Space *lastWaitingLps) {

	std::ostringstream indentStream;
	for (int i = 0; i < indentation; i++) indentStream << indent;
	std::string indentStr = indentStream.str();

	// Find the lowest LPS in the sync requirements that should wait on a barrier. If we make all the PPUs
        // handling code from that LPS waiting once then it will, by the sake of the hierarchy, synchronize all
        // parent PPUs also.
        Space *lowestWaitingLps = NULL;
        for (int i = 0; i < syncRequirements->NumElements(); i++) {
                SyncRequirement *sync = syncRequirements->Nth(i);
                Space *syncSpanLps = sync->getSyncSpan();
                if (lowestWaitingLps == NULL) {
                        if (lastWaitingLps == NULL
                                        || strcmp(lastWaitingLps->getName(), syncSpanLps->getName()) != 0) {
                                lowestWaitingLps = syncSpanLps;
                        }
                } else if (syncSpanLps->isParentSpace(lowestWaitingLps)) {
                        lowestWaitingLps = syncSpanLps;
                }
        }

        if (lowestWaitingLps == NULL) return;

	if (syncRequirements->NumElements() > 0) {
		stream << std::endl << indentStr;
		stream << "// barriers to ensure all readers have finished reading last update\n";
	}

	for (int i = 0; i < syncRequirements->NumElements(); i++) {
		SyncRequirement *sync = syncRequirements->Nth(i);
		Space *syncSpanLps = sync->getSyncSpan();
		stream << indentStr << "if (threadState->isValidPpu(Space_";
		stream << syncSpanLps->getName();
		stream << ")) {\n";	
		stream << indentStr << indent;
		stream << "threadSync->" << sync->getReverseSyncName() << "->wait()";
		stream << stmtSeparator;
		stream << indentStr << "}\n";
	}
}

Space *CompositeStage::genSimplifiedSignalsForGroupTransitionsCode(std::ostream &stream, int indentation,
		List<SyncRequirement*> *syncRequirements) {

	std::ostringstream indentStream;
	for (int i = 0; i < indentation; i++) indentStream << indent;
	std::string indentStr = indentStream.str();

	if (syncRequirements->NumElements() > 0) {
		stream << std::endl << indentStr << "// resolving synchronization dependencies\n";
	}

	// iterate over all the synchronization signals and then issue signals and waits in a lock-step fasion
	Space *lastSignalingLps = NULL;
        Space *lastWaitingLps = NULL;
	for (int i = 0; i < syncRequirements->NumElements(); i++) {
		
		// Within the shared memory environment, keeping PPUs at a space waiting for a certain data is
                // sufficient to synchronize them about other data also as long as the PPU doing signaling operates
                // at the same or a descendent level. Therefore, this checking is added to skip some redundant
                // synchronizations.
                SyncRequirement *currentSync = syncRequirements->Nth(i);
                FlowStage *sourceStage = currentSync->getDependencyArc()->getSource();
                Space *syncSpanLps = currentSync->getSyncSpan();
                Space *signalingLps = sourceStage->getSpace();
                if (lastWaitingLps != NULL
                                && strcmp(syncSpanLps->getName(), lastWaitingLps->getName()) == 0
                                && (
                                        (strcmp(lastSignalingLps->getName(), signalingLps->getName()) == 0)
                                         || lastWaitingLps->isParentSpace(signalingLps)
                                         || lastSignalingLps->isParentSpace(signalingLps)
                                )) {
                        continue;
                }
                lastWaitingLps = syncSpanLps;
                lastSignalingLps = signalingLps;

		const char *counterVarName = currentSync->getDependencyArc()->getArcName();
	
		// Check if the current synchronization is conditional, i.e., it only gets signaled by threads that
		// executed a certain execution flow-stage. In that case, there will be a counter variable set to a
		// non-zero value. For synchronization signal issued by compiler injected sync-stages, there will be
		// no such counter variable. 
		bool needCounter = currentSync->getCounterRequirement();
		
		stream << indentStr << "if (";
		// check if the concerned update did take place
		if (needCounter) {
			stream << counterVarName << " > 0 && ";
		}
		// also check if the current PPU is a valid candidate for signaling update
		stream << "threadState->isValidPpu(Space_" << signalingLps->getName();
		stream << ")) {\n";
		// then signal synchronization
		stream << indentStr << indent;
		stream << "threadSync->" << currentSync->getSyncName() << "->signal(";
		FlowStage *signalSource = currentSync->getDependencyArc()->getSignalSrc();
		if (signalSource->getRepeatIndex() > 0) stream << "repeatIteration";
		else stream << "0";
		stream << ")" << stmtSeparator;
		// then reset the counter
		if (needCounter) {	 
			stream << indentStr << indent;
			stream << counterVarName << " = 0" << stmtSeparator;
		}

		// the waiting is in an else-block coupled with the signaling if-block as the current implementation
		// of synchronization primitives does not support the PPU (or PPUs) in the signaling block to also be
		// among the list of waiting PPUs.  
		stream << indentStr << "} else if (";
		FlowStage *waitingStage = currentSync->getWaitingComputation();
		stream << "threadState->isValidPpu(Space_" << syncSpanLps->getName();
		stream << ")) {\n";
		stream << indentStr << indent;
		stream << "threadSync->" << currentSync->getSyncName() << "->wait(";
		FlowStage *signalSink = currentSync->getDependencyArc()->getSignalSink();
		if (signalSink->getRepeatIndex() > 0) stream << "repeatIteration";
		else stream << "0";
		stream << ")" << stmtSeparator;
		stream << indentStr << "}\n";
	}

	return lastWaitingLps;
}

void CompositeStage::generateDataSendsForGroup(std::ofstream &stream, int indentation, 
		List<SyncRequirement*> *commRequirements) {
	
	std::ostringstream indentStr;
	for (int i = 0; i < indentation; i++) indentStr << indent;

	if (commRequirements->NumElements() > 0) {
		stream << std::endl << indentStr.str() << "// communicating updates\n";
	}

	// iterate over all the update signals
	for (int i = 0; i < commRequirements->NumElements(); i++) {
		
		SyncRequirement *currentComm = commRequirements->Nth(i);
		if (!currentComm->isActive()) continue;

		const char *counterVarName = currentComm->getDependencyArc()->getArcName();
		
		// Check if there are other communication requests to the same data receiver LPS for the same data item.
		// If there are such requests are redundant and deactivate them.
		const char *varName = currentComm->getDependencyArc()->getVarName();
		Space *syncSpanLps = currentComm->getSyncSpan();
		for (int j = i + 1; j < commRequirements->NumElements(); j++) {
			SyncRequirement *nextComm = commRequirements->Nth(j);
			if (nextComm->isActive() && strcmp(varName, nextComm->getDependencyArc()->getVarName()) == 0) {
                		Space *nextSyncSpanLps = nextComm->getSyncSpan();
				if (strcmp(syncSpanLps->getName(), nextSyncSpanLps->getName()) == 0) {
					const char *nextCounterVarName = currentComm->getDependencyArc()->getArcName();
					nextComm->getDependencyArc()->deactivate();
				} else if (nextSyncSpanLps->isParentSpace(syncSpanLps)) {

					// on the other hand, if the subsequent  communication request on the same variable
					// has a deeper receiving LPS than the current one then we should rather deactive
					// the current communication signal and setup information to move the reception
					// on the next communication request forward.
					currentComm->expandSynchingTo(nextComm);	
				}
			}
		}
		if (currentComm->isSynchingExpanded()) continue;

		int commIndex = currentComm->getIndex();
		
		// check if the current PPU is a valid candidate for signaling update
		Space *signalingLps = currentComm->getDependencyArc()->getSource()->getSpace();
		stream << indentStr.str() << "if (threadState->isValidPpu(Space_" << signalingLps->getName();
		stream << ")) {\n";
		
		// retrieve the communicator for this dependency
		stream << indentStr.str() << indent;
		stream << "Communicator *communicator = threadState->getCommunicator(\"";
		stream << currentComm->getDependencyArc()->getArcName() << "\")" << stmtSeparator;
		stream << indentStr.str() << indent << "if (communicator != NULL) {\n";
		
		// Check if the current communication is conditional, i.e., it only gets signaled by threads that
		// executed a certain execution flow-stage. In that case, there will be a counter variable set to a
		// non-zero value. For communication signals issued by compiler injected sync-stages, there will be
		// no such counter variable. 
		bool needCounter = currentComm->getCounterRequirement();
		
		if (needCounter) {
			// If the counter variable for the sync has been updated then the current PPU controller has 
			// data to send so it should indicate that fact in its call to the communicator. Otherwise
			// it should only report that it has reached this particular execution point.
			stream << indentStr.str() << doubleIndent << "if (" << counterVarName << " > 0) {\n";
			stream << indentStr.str() << tripleIndent << "communicator->send(";
			stream << "REQUESTING_COMMUNICATION" << paramSeparator;
			stream << "commCounter" << commIndex <<  ")" << stmtSeparator;
			stream << indentStr.str() << doubleIndent << "} else communicator->send(";
			stream << "PASSIVE_REPORTING" << paramSeparator;
			stream << "commCounter" << commIndex << ")" << stmtSeparator;
		} else {
			stream << indentStr.str() << doubleIndent << "communicator->send(";
			stream << "REQUESTING_COMMUNICATION" << paramSeparator;
			stream << "commCounter" << commIndex <<  ")" << stmtSeparator;
		}

		stream << indentStr.str() << indent << "}\n";
		
		// then reset the counter
		if (needCounter) {	 
			stream << indentStr.str() << indent << counterVarName << " = 0" << stmtSeparator;
		}
		stream << indentStr.str() << "}\n";
	}
}

void CompositeStage::generateInvocationCode(std::ofstream &stream, int indentation, Space *containerSpace) {

        // if the index is 0 then it is the first composite stage representing the entire compution. We declare any 
	// synchronization counter that is applicable outside all repeat-control-block boundaries
        if (this->index == 0) {
                declareSynchronizationCounters(stream, indentation, this->repeatIndex + 1);
        }

	// Iterate over groups of flow stages where each group executes within a single LPS. This scheme has the
        // consequence of generating LPU only one time for all stages of a group then execute all of them before
        // proceed to the next LPU 
        List<List<FlowStage*>*> *stageGroups = getConsecutiveNonLPSCrossingStages();
	Space *lastWaitingLps = NULL;
        for (int i = 0; i < stageGroups->NumElements(); i++) {
		
		List<FlowStage*> *currentGroup = stageGroups->Nth(i);

		//------------------------------------------------------------------------ Dependency Handling Starts
		
                // retrieve all data dependencies, and sort them to ensure waitings for updates happen in order
                List<SyncRequirement*> *dataDependencies = getDataDependeciesOfGroup(currentGroup);
                dataDependencies = SyncRequirement::sortList(dataDependencies);

		// we should separate dependencies into communication and synchronization dependencies then take 
		// appropriate actions based on type 
		// (the current simplified implementation does nothing with the synchronization dependencies)
                int segmentedPPS = space->getSegmentedPPS();
                List<SyncRequirement*> *commDependencies = new List<SyncRequirement*>;
                List<SyncRequirement*> *syncDependencies = new List<SyncRequirement*>;
                SyncRequirement::separateCommunicationFromSynchronizations(segmentedPPS,
                		dataDependencies, commDependencies, syncDependencies);
                generateDataReceivesForGroup(stream, indentation, commDependencies);

		//-------------------------------------------------------------------------- Dependency Handling Ends
		
		//---------------------------------------------------------------- Write After Read Activation Starts

		// retrieve all shared data update signals that need to be activated if stages in the group execute
                List<SyncRequirement*> *updateSignals = getUpdateSignalsOfGroup(currentGroup);
                // mark these signals as signaled so that they are not reactivated within the nested code
                for (int j = 0; j < updateSignals->NumElements(); j++) {
                        updateSignals->Nth(j)->signal();
                }
                // sort the update signals to ensure waiting for signal clearance (equivalent to get signals from the
                // readers that all of them have finished reading the last update) happens in proper order
                updateSignals = SyncRequirement::sortList(updateSignals);

                // divide the signals between those issuing communications and those that do not  
                List<SyncRequirement*> *commSignals = new List<SyncRequirement*>;
                List<SyncRequirement*> *syncSignals = new List<SyncRequirement*>;
                SyncRequirement::separateCommunicationFromSynchronizations(segmentedPPS,
                                updateSignals, commSignals, syncSignals);

		// If there is any reactivating condition that need be checked before we let the flow of control 
                // enter the nested stages then we wait for those condition clearance.
		// 
		// If the last group of stages of the current composite stage issue a sync signal for an LPS
                // that the first group of stages operate in then we can ignore a backward reactivation barrier
                // before the first group of stages as a special condition.
                if (i == 0) {
                        std::ostringstream tempStream;
                        List<FlowStage*> *lastGroup = stageGroups->Nth(stageGroups->NumElements() - 1);
                        List<SyncRequirement*> *lastUpdateSignals = getUpdateSignalsOfGroup(lastGroup);
                        lastUpdateSignals = SyncRequirement::sortList(lastUpdateSignals);
                	List<SyncRequirement*> *lastCommSignals = new List<SyncRequirement*>;
                	List<SyncRequirement*> *lastSyncSignals = new List<SyncRequirement*>;
                	SyncRequirement::separateCommunicationFromSynchronizations(segmentedPPS,
                                	lastUpdateSignals, lastCommSignals, lastSyncSignals);
                        lastWaitingLps = genSimplifiedSignalsForGroupTransitionsCode(tempStream, 0, lastSyncSignals);
                }

                genSimplifiedWaitingForReactivationCode(stream, indentation, syncSignals, lastWaitingLps);

		//------------------------------------------------------------------ Write After Read Activation Ends

		//----------------------------------------------- Signal Management for Sync Stage Only Groups Starts

		// Sync stages that dictate additional data movement needs are not needed during code generation. so 
		// we filter them out.     
                currentGroup = filterOutSyncStages(currentGroup);
                if (currentGroup->NumElements() == 0) {
                        // before rulling sync stages out, we need to ensure that whatever signals they were supposed
                        // issue are by-default issued and whatever data they were supposed to send are sent
                        genSimplifiedSignalsForGroupTransitionsCode(stream, indentation, syncSignals);
                        generateDataSendsForGroup(stream, indentation, commSignals);
                        continue;
                }		

		//------------------------------------------------- Signal Management for Sync Stage Only Groups Ends

		//------------------------------------------------------ Code Generation for the Current Group Starts

		for (int j = 0; j < currentGroup->NumElements(); j++) {
			FlowStage *stage = currentGroup->Nth(j);
			stage->generateInvocationCode(stream, indentation, space);
		}		

		//-------------------------------------------------------- Code Generation for the Current Group ends

		//--------------------------------------------------------------------- Update Signal Handling Starts

		Space *waitingLps = genSimplifiedSignalsForGroupTransitionsCode(stream, indentation, syncSignals);
		if (lastWaitingLps == NULL) {
                        lastWaitingLps = waitingLps;
                } else if (waitingLps != NULL) {
                        if (waitingLps->isParentSpace(lastWaitingLps)) {
                                lastWaitingLps = waitingLps;
                        }
                }
                generateDataSendsForGroup(stream, indentation, commSignals);	

		//----------------------------------------------------------------------- Update Signal Handling Ends
	}	
}	


