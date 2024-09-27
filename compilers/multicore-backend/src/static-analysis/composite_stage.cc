#include "data_flow.h"
#include "data_access.h"
#include "sync_stat.h"
#include "../syntax/ast.h"
#include "../syntax/ast_expr.h"
#include "../syntax/ast_stmt.h"
#include "../syntax/ast_task.h"
#include "../syntax/ast_type.h"
#include "../utils/list.h"
#include "../utils/hashtable.h"
#include "../semantics/task_space.h"
#include "../semantics/scope.h"
#include "../semantics/symbol.h"
#include "../codegen/name_transformer.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <string.h>

//------------------------------------------------ Composite Stage -------------------------------------------------------/

CompositeStage::CompositeStage(int index, Space *space, Expr *executeCond) : FlowStage(index, space, executeCond) {
	stageList = new List<FlowStage*>;	
}

void CompositeStage::addStageAtBeginning(FlowStage *stage) {
	stageList->InsertAt(stage, 0);
	stage->setParent(this);
}

void CompositeStage::addStageAtEnd(FlowStage *stage) {
	stageList->Append(stage);
	SyncStage *syncStage = dynamic_cast<SyncStage*>(stage);
	if (syncStage != NULL) {
		syncStage->setIndex(stageList->NumElements());
	}
	stage->setParent(this);
}

void CompositeStage::insertStageAt(int index, FlowStage *stage) {
	stageList->InsertAt(stage, index);
	stage->setIndex(index);
	stage->setParent(this);
}

void CompositeStage::removeStageAt(int stageIndex) { stageList->RemoveAt(stageIndex); }

void CompositeStage::setStageList(List<FlowStage*> *stageList) { 
        this->stageList = stageList;
        for (int i = 0; i < stageList->NumElements(); i++) {
		stageList->Nth(i)->setParent(this);
	} 
}

bool CompositeStage::isStageListEmpty() {
	return getLastNonSyncStage() == NULL;
}

Space *CompositeStage::getLastNonSyncStagesSpace() {
	Space *space = this->space;
	FlowStage *lastNonSyncStage = getLastNonSyncStage();
	if (lastNonSyncStage != NULL) return lastNonSyncStage->getSpace();
	return space;
}

FlowStage *CompositeStage::getLastNonSyncStage() {
	int i = stageList->NumElements() - 1; 
	for (; i >= 0; i--) {
		FlowStage *stage = stageList->Nth(i);
		SyncStage *syncStage = dynamic_cast<SyncStage*>(stage);
		if (syncStage == NULL) return stage;
	}
	return NULL;

}

void CompositeStage::addSyncStagesBeforeExecution(FlowStage *nextStage, List<FlowStage*> *stageList) {
	
	Space *previousSpace = getLastNonSyncStagesSpace();
	Space *nextSpace = nextStage->getSpace();
	
	List<Space*> *spaceTransitionChain = Space::getConnetingSpaceSequenceForSpacePair(previousSpace, nextSpace);
	if (spaceTransitionChain == NULL || spaceTransitionChain->NumElements() == 0) {
		if (!isStageListEmpty()) {
			spaceTransitionChain = new List<Space*>;
			spaceTransitionChain->Append(previousSpace);
			spaceTransitionChain->Append(nextSpace);
		}
	} 

	int nextStageIndex = nextStage->getIndex();
	if (spaceTransitionChain != NULL) {
		for (int i = 1; i < spaceTransitionChain->NumElements(); i++) {
			Space *oldSpace = spaceTransitionChain->Nth(i - 1);
			Space *newSpace = spaceTransitionChain->Nth(i);
			if (oldSpace->isParentSpace(newSpace)) {
				// new space is higher in the space hierarchy; so an exit from the old space 
				// should be recorded along with an entry to the new space
				SpaceEntryCheckpoint *oldCheckpoint = SpaceEntryCheckpoint::getCheckpoint(oldSpace);
				SyncStage *oldEntrySyncStage = oldCheckpoint->getEntrySyncStage();
				Hashtable<VariableAccess*> *accessLogs = getAccessLogsForSpaceInIndexLimit(oldSpace,
                        			stageList, oldCheckpoint->getStageIndex(), 
						nextStageIndex - 1, true);
				// If there is an entry sync stage for the old space then we need to populate
				// its accessMap currectly.
				if (oldEntrySyncStage != NULL) {
					SyncStageGenerator::populateAccessMapOfEntrySyncStage(oldEntrySyncStage, accessLogs);
				}
				// generate and add to the list all possible sync stages that are required due
				// to the exit from the old space
				SpaceEntryCheckpoint::removeACheckpoint(oldSpace);
				List<SyncStage*> *exitSyncs = SyncStageGenerator::generateExitSyncStages(oldSpace, accessLogs);
				for (int i = 0; i < exitSyncs->NumElements(); i++) {
					SyncStage *exitSyncStage = exitSyncs->Nth(i);
					addStageAtEnd(exitSyncStage);
				}
				// generate and add any potential return sync stage to the new space
				accessLogs = getAccessLogsForReturnToSpace(newSpace, stageList, nextStageIndex - 1);
				SyncStage *returnSync = SyncStageGenerator::generateReturnSyncStage(newSpace, accessLogs);
				if (returnSync != NULL) {
					addStageAtEnd(returnSync);
				}
			} else if (newSpace->isParentSpace(oldSpace)) {
				// Old space is higher in the space hierarchy; so an entry to the new space
				// should be recorded. The entry sync stage here, if present, is just a place holder.
				// Later on during exit its accessLog is filled with appropriate data
				SyncStage *entrySyncStage = SyncStageGenerator::generateEntrySyncStage(newSpace);
				SpaceEntryCheckpoint *checkpoint = 
						SpaceEntryCheckpoint::addACheckpointIfApplicable(newSpace, nextStageIndex);
				checkpoint->setEntrySyncStage(entrySyncStage);
				if (entrySyncStage != NULL) {
					addStageAtEnd(entrySyncStage);
				}
			} else if (oldSpace == newSpace) {
				// transition is made from some computation of the same space to another; we
				// may need to synchronize overlapping boundary regions of data structures
				FlowStage *lastStageOnSpace = getLastNonSyncStage();
				Hashtable<VariableAccess*>  *accessLogs = lastStageOnSpace->getAccessMap();
				SyncStage *reappearanceSync = 
						SyncStageGenerator::generateReappearanceSyncStage(newSpace, accessLogs);
				if (reappearanceSync != NULL) {
					addStageAtEnd(reappearanceSync);
				}
			} else {
				printf("We should never have a disjoint space transition chain\n");
				exit(1);
			}
		}
	}
}

void CompositeStage::addSyncStagesOnReturn(List<FlowStage*> *stageList) {
        
	Space *previousSpace = getLastNonSyncStagesSpace();
	Space *currentSpace = getSpace();
	List<Space*> *spaceTransitionChain = Space::getConnetingSpaceSequenceForSpacePair(previousSpace, currentSpace);
	if (spaceTransitionChain == NULL || spaceTransitionChain->NumElements() == 0) return;
	
	int lastStageIndex = getLastNonSyncStage()->getIndex();
	for (int i = 1; i < spaceTransitionChain->NumElements(); i++) {
        	Space *oldSpace = spaceTransitionChain->Nth(i - 1);
        	Space *newSpace = spaceTransitionChain->Nth(i);
		SpaceEntryCheckpoint *oldCheckpoint = SpaceEntryCheckpoint::getCheckpoint(oldSpace);
                SyncStage *oldEntrySyncStage = oldCheckpoint->getEntrySyncStage();
                Hashtable<VariableAccess*> *accessLogs = getAccessLogsForSpaceInIndexLimit(oldSpace,
                		stageList, oldCheckpoint->getStageIndex(), lastStageIndex, true);
                // If there is an entry sync stage for the old space then we need to populate
                // its accessMap currectly.
                if (oldEntrySyncStage != NULL) {
                	SyncStageGenerator::populateAccessMapOfEntrySyncStage(oldEntrySyncStage, accessLogs);
                }
                // generate and add to the list all possible sync stages that are required due
                // to the exit from the old space
                SpaceEntryCheckpoint::removeACheckpoint(oldSpace);
                List<SyncStage*> *exitSyncs = SyncStageGenerator::generateExitSyncStages(oldSpace, accessLogs);
                for (int i = 0; i < exitSyncs->NumElements(); i++) {
                	SyncStage *exitSyncStage = exitSyncs->Nth(i);
                	addStageAtEnd(exitSyncStage);
                }
                // generate and add any potential return sync stage to the new space
                accessLogs = getAccessLogsForReturnToSpace(newSpace, stageList, lastStageIndex);
                SyncStage *returnSync = SyncStageGenerator::generateReturnSyncStage(newSpace, accessLogs);
                if (returnSync != NULL) {
                	addStageAtEnd(returnSync);
                }	
	}	
}

void CompositeStage::print(int indent) {
	for (int i = 0; i < indent; i++) printf("\t");
	printf("Flow Stage: %s (Space: %s)", name, space->getName());
	if (executeCond != NULL) printf(" Conditionally Executed");
	printf("\n");
	for (int i = 0; i < stageList->NumElements(); i++) {
		stageList->Nth(i)->print(indent + 1);
	}
	dataDependencies->print(indent + 1);
}

void CompositeStage::performDependencyAnalysis(PartitionHierarchy *hierarchy) {
	for (int i = 0; i < stageList->NumElements(); i++) {
		stageList->Nth(i)->performDependencyAnalysis(hierarchy);
	}
}

void CompositeStage::reorganizeDynamicStages() {
	
	int currentStageIndex = 0;
	while (currentStageIndex < stageList->NumElements()) {
		
		FlowStage *stage = stageList->Nth(currentStageIndex);
		
		// do the reorganization for nested stages wherever necessary
		stage->reorganizeDynamicStages();

		stage->setIndex(currentStageIndex);
		currentStageIndex++;

		if (!stage->getSpace()->isDynamic()) continue;
		SyncStage *syncStage = dynamic_cast<SyncStage*>(stage);
		if (syncStage == NULL || !syncStage->isLoaderSync()) continue;
		
		// gather the stages that should be included in the generated composite stage
		List<Expr*> *activationConditions = new List<Expr*>;
		List<FlowStage*> *toBeNestedStages = new List<FlowStage*>;
		int stageListStart = currentStageIndex;
		Space *dynamicSpace = stage->getSpace();
		while (stageList->NumElements() > stageListStart) {
			FlowStage *newStage = stageList->Nth(stageListStart);
			Space *newSpace = newStage->getSpace();
			if (newSpace != dynamicSpace && !newSpace->isParentSpace(dynamicSpace)) break;
			if (newSpace == dynamicSpace && newStage->getExecuteCondition() != NULL) {
				activationConditions->Append(newStage->getExecuteCondition());
			}
			toBeNestedStages->Append(newStage);
			// do the reorganization for nested stages wherever necessary
			newStage->reorganizeDynamicStages();
			// remove the stage from current list
			removeStageAt(stageListStart);
		}
		
		// generate an updated execute condition for the new to be added composite stage
		Expr *overallCondition = NULL;
		if (activationConditions->NumElements() == 1) {
			for (int i = 0; i < toBeNestedStages->NumElements(); i++) {
				FlowStage *newStage = toBeNestedStages->Nth(i);
				if (newStage->getSpace() == dynamicSpace) {
					// TODO probably the correct solution is to wrap up the computation
					// inside within an conditional block that test the validity of 
					// entering into this stage. That way we will avoid the problem of
					// reevaluating active LPUs for this stage and at the same time ensure
					// that we are not going to execute any unintended code.	
					newStage->setExecuteCondition(NULL);
				}
			}
			overallCondition = activationConditions->Nth(0);
		} else if (activationConditions->NumElements() > 1) {
			overallCondition = activationConditions->Nth(0);
			for (int i = 1; i < activationConditions->NumElements(); i++) {
				Expr *nextCondition = activationConditions->Nth(i);
				overallCondition = new LogicalExpr(overallCondition, OR, 
						nextCondition, *overallCondition->GetLocation());
			}	
		}

		// generate the composite stage
		CompositeStage *generatedStage = new CompositeStage(currentStageIndex - 1, 
				dynamicSpace, overallCondition);
		generatedStage->setName("\"Dynamic Computation\"");
		generatedStage->setStageList(toBeNestedStages);
		generatedStage->addStageAtBeginning(syncStage);
		
		// remove the sync stage from current list and add the composite stage in its place
		removeStageAt(currentStageIndex - 1);
		insertStageAt(currentStageIndex - 1, generatedStage);
	}
}

List<List<FlowStage*>*> *CompositeStage::getConsecutiveNonLPSCrossingStages() {
	
	List<List<FlowStage*>*> *stageGroups = new List<List<FlowStage*>*>;
	Space *currentSpace = stageList->Nth(0)->getSpace();
	List<FlowStage*> *currentGroup = new List<FlowStage*>;
	currentGroup->Append(stageList->Nth(0));

	for (int i = 1; i < stageList->NumElements(); i++) {
		FlowStage *stage = stageList->Nth(i);
		RepeatCycle *repeatCycle = dynamic_cast<RepeatCycle*>(stage);
		// repeat cycles are never put into any compiler generated group to keep the semantics of the
		// program intact
		if (repeatCycle != NULL) {
			// add the so far group into the stage groups list if it is not empty
			if (currentGroup->NumElements() > 0) stageGroups->Append(currentGroup);
			// create a new, isolated group for repeat and add that in the list too 
			List<FlowStage*> *repeatGroup = new List<FlowStage*>;
			repeatGroup->Append(repeatCycle);
			stageGroups->Append(repeatGroup);
			// then reset the current group
			currentGroup = new List<FlowStage*>;
		// if the stage is executing in a different LPS then create a new group
		}	
		else if (stage->getSpace() != currentSpace) {
			currentSpace = stage->getSpace();
			if (currentGroup->NumElements() > 0) stageGroups->Append(currentGroup);
			currentGroup = new List<FlowStage*>;
			currentGroup->Append(stage);
		// otherwise add the stage in the current group if there is no synchronization dependency from
		// current stage to any stage already within the group. If there is a dependency then create a
		// new group. Note that this policy helps us to drag down LPU-LPU synchronization dependencies 
		// as transition has been made between computation stages into PPU-PPU dependencies. 
		} else {
			bool dependencyExists = false;
			for (int j = 0; j < currentGroup->NumElements(); j++) {
				FlowStage *earlierStage = currentGroup->Nth(j);
				if (earlierStage->isDependentStage(stage) 
						|| stage->isDependentStage(earlierStage)) {
					dependencyExists = true;
					break;
				}
			}
			if (dependencyExists) {
				stageGroups->Append(currentGroup);
				currentGroup = new List<FlowStage*>;
				currentGroup->Append(stage);
			} else 	currentGroup->Append(stage);
		}
	}
	if (currentGroup->NumElements() > 0) stageGroups->Append(currentGroup);
	return stageGroups;
}

void CompositeStage::generateInvocationCode(std::ofstream &stream, int indentation, Space *containerSpace) {
	
	std::string stmtSeparator = ";\n";
	std::string paramSeparator = ", ";
	std::ostringstream indent;
	for (int i = 0; i < indentation; i++) indent << '\t';
	int nextIndentation = indentation;
	std::ostringstream nextIndent;
	nextIndent << indent.str();

	// if the index is 0 then it is the first composite stage representing the entire compution. We declare
	// any synchronization counter that is applicable outside all repeat-cycle boundaries
	if (this->index == 0) {
		declareSynchronizationCounters(stream, indentation, this->repeatIndex);
	}

	// if their is an LPS transition due to entering this stage then create a while loop traversing LPUs
	// of newly entered LPS
	if (this->space != containerSpace) {
		const char *spaceName = space->getName();
		nextIndentation++;
		nextIndent << '\t';
		// create a new local scope for traversing LPUs of this new scope
		stream << std::endl;
		stream << indent.str() << "{ // scope entrance for iterating LPUs of Space ";
		stream << spaceName << "\n";
		// declare a new variable for tracking the last LPU id
		stream << indent.str() << "int space" << spaceName << "LpuId = INVALID_ID" << stmtSeparator;
		// declare another variable to track the iteration number of the while loop
		stream << indent.str() << "int space" << spaceName << "Iteration = 0" << stmtSeparator;
		// declare a new variable to hold on to current LPU of this LPS
		stream << indent.str() << "Space" << spaceName << "_LPU *space" << spaceName << "Lpu = NULL";
		stream << stmtSeparator;
		// declare another variable to assign the value of get-Next-LPU call
		stream << indent.str() << "LPU *lpu = NULL" << stmtSeparator;
		stream << indent.str() << "while((lpu = threadState->getNextLpu(";
		stream << "Space_" << spaceName << paramSeparator << "Space_" << containerSpace->getName();
		stream << paramSeparator << "space" << spaceName << "LpuId)) != NULL) {\n";
		// cast the common LPU variable to LPS specific LPU		
		stream << nextIndent.str() << "space" << spaceName << "Lpu = (Space" << spaceName;
		stream  << "_LPU*) lpu" << stmtSeparator;
	}

	// Check if there is any activating condition attached with this container stage. If there is such a
	// condition then the stages inside should be executed only if the active condition evaluates to true.
	if (executeCond != NULL) {
		
		// first get a hold of the LPU reference
		std::ostringstream lpuVarName;
		lpuVarName << "space" << space->getName() << "Lpu->";
		
		// then generate local variables for any array been accessed within the condition
		List<FieldAccess*> *fields = executeCond->getTerminalFieldAccesses();
		for (int i = 0; i < fields->NumElements(); i++) {
			DataStructure *structure = space->getLocalStructure(fields->Nth(i)->getField()->getName());
			if (structure == NULL) continue;
			ArrayDataStructure *array = dynamic_cast<ArrayDataStructure*>(structure);
			if (array == NULL) continue;
			const char *arrayName = array->getName();
        		int dimensions = array->getDimensionality();
        		stream << nextIndent.str() << "Dimension ";
                	stream  << arrayName << "PartDims[" << dimensions << "];\n";
        		stream << nextIndent.str() << "Dimension ";
                	stream  << arrayName << "StoreDims[" << dimensions << "];\n";
                	for (int j = 0; j < dimensions; j++) {
                		stream << indent.str();
               			stream << arrayName << "PartDims[" << j << "] = " << lpuVarName.str();
                        	stream << arrayName << "PartDims[" << j << "].partition;\n";
                		stream << indent.str();
               			stream << arrayName << "StoreDims[" << j << "] = " << lpuVarName.str();
                        	stream << arrayName << "PartDims[" << j << "].storage;\n";
        		}
		}
		
		// now set the name transformer's LPU prefix properly so that if the activate condition involves
		// accessing elements of the LPU, it works correctly
		ntransform::NameTransformer::transformer->setLpuPrefix(lpuVarName.str().c_str());
	
		// then generate an if condition for condition checking
		stream << nextIndent.str() << "if(!(";
		std::ostringstream conditionStream;
		executeCond->translate(conditionStream, nextIndentation, 0, space);
		stream << conditionStream.str();
		stream << ")) {\n";
		// we skip the current LPU if the condition evaluates to false	
		stream << nextIndent.str() << "\tcontinue" << stmtSeparator;
		stream << nextIndent.str() << "}\n";	
	}
	
	// Iterate over groups of flow stages where each group executes within a single LPS. This scheme has the
	// consequence of generating LPU only one time for all stages of a group then execute all of them before
	// proceed to the next LPU 
	List<List<FlowStage*>*> *stageGroups = getConsecutiveNonLPSCrossingStages();

	Space *lastWaitingLps = NULL;
	List<Space*> *growingWaitingList = new List<Space*>;	
	for (int i = 0; i < stageGroups->NumElements(); i++) {
		
		List<FlowStage*> *currentGroup = stageGroups->Nth(i);

		// commented out this code as the barrier based priliminary implementation does not need this
		/*-------------------------------------------------------------------------------------------------
		// first generate code for waiting on all synchronization dependencies before diving into 
		// multiplexed LPU iterations.
		List<SyncRequirement*> *syncDependencies = getSyncDependeciesOfGroup(currentGroup);
		// sort the sync dependencies to ensure that their is a consistent order of the waiting primitives
		syncDependencies = SyncRequirement::sortList(syncDependencies);
		generateSyncCodeForGroupTransitions(stream, nextIndentation, syncDependencies);
		//-----------------------------------------------------------------------------------------------*/
		
		// retrieve all synchronization signals that need to be activated if the stages in the group execute
		List<SyncRequirement*> *syncSignals = getSyncSignalsOfGroup(currentGroup);
		// mark these signals as signaled so that they are not reactivated within the nested code
		for (int j = 0; j < syncSignals->NumElements(); j++) {
			syncSignals->Nth(j)->signal();
		}
		// sort the sync signals to ensure waiting for signal clearance (equivalent to get signals from the
		// readers that all of them have finished reading the last update) happens in proper order
		syncSignals = SyncRequirement::sortList(syncSignals);

		// simplified implementation	
		//-------------------------------------------------------------------------------------------------
		// If there is any reactivating condition that need be checked before we let the flow of control 
		// enter the nested stages then we wait for those condition clearance.
		//
		// if the last group of stages of the current composite stage issue a sync signal for an LPS
		// that the first group of stages operate in then we can ignore a backward reactivation barrier
		// before the first group of stages as a special condition.
		if (i == 0) {
			std::ostringstream tempStream;
			List<FlowStage*> *lastGroup = stageGroups->Nth(stageGroups->NumElements() - 1);
			List<SyncRequirement*> *lastSyncSignals = getSyncSignalsOfGroup(lastGroup);
			lastSyncSignals = SyncRequirement::sortList(lastSyncSignals);
			lastWaitingLps = genSimplifiedSignalsForGroupTransitionsCode(tempStream, 0, lastSyncSignals);
		}
		List<Space*> *followingWaitList = getWaitingLpsInFollowingStages(i, stageGroups);
		followingWaitList->AppendAll(growingWaitingList);
		List<Space*> *lpsesWaited = genSimplifiedWaitingForReactivationCode(stream, nextIndentation, 
				syncSignals, lastWaitingLps, followingWaitList);
		growingWaitingList->AppendAll(lpsesWaited);
		//-------------------------------------------------------------------------------------------------

		// Sync stages -- not synchronization dependencies -- that dictate additional data movement 
		// operations are not needed (according to our current (date: Jan-30-2015) understanding); so we
		// filter them out.	
		currentGroup = filterOutSyncStages(currentGroup);
		if (currentGroup->NumElements() == 0) {
			// before rulling sync stages out, we need to ensure that whatever signals they were supposed
			// issue are by-default issued
			generateSignalCodeForGroupTransitions(stream, nextIndentation, syncSignals); 
			continue;
		}

		// there should be a special checking for repeat loops. Repeat loops will have the while loop that
		// iterates over LPUs of the LPS under concern inside its body if it does not use LPS dependent 
		// variables in repeat evaluation process. Otherwise, it should follow the normal code generation
		// procedure followed for other stages 
		bool groupIsLpsIndRepeat = false;
		if (currentGroup->NumElements() == 1) {
			FlowStage *stage = currentGroup->Nth(0);
			RepeatCycle *repeatCycle = dynamic_cast<RepeatCycle*>(stage);
			if (repeatCycle != NULL && !repeatCycle->isLpsDependent()) {
				Space *repeatSpace = repeatCycle->getSpace();
				// temporarily the repeat cycle is assigned to be executed in current LPS to ensure
				// that LPU iterations for stages inside it proceed accurately
				repeatCycle->changeSpace(this->space);
				repeatCycle->generateInvocationCode(stream, nextIndentation, space);
				// reset the repeat cycle's LPS to the previous one once code generation is done
				repeatCycle->changeSpace(repeatSpace);
				groupIsLpsIndRepeat = true;
			}	
		}

		if (!groupIsLpsIndRepeat) {
			Space *groupSpace = currentGroup->Nth(0)->getSpace();
			// if the LPS of the group is not the same of this one then we need to consider LPS entry, 
			// i.e., include a while loop for the group to traverse over the LPUs
			if (groupSpace != space) {
				// create a local composite stage to apply the logic of this function recursively
				CompositeStage *tempStage = new CompositeStage(-1, groupSpace, NULL);
				tempStage->setStageList(currentGroup);
				tempStage->generateInvocationCode(stream, nextIndentation, space);
				delete tempStage;
			// otherwise the current LPU of this composite stage will suffice and we execute all the 
			// nested stages one after one.	 
			} else {
				for (int j = 0; j < currentGroup->NumElements(); j++) {
					FlowStage *stage = currentGroup->Nth(j);
					stage->generateInvocationCode(stream, nextIndentation, space);
				}
			}
		}

		// simplified implementation	
		//-------------------------------------------------------------------------------------------------
		// generate code for signaling updates and waiting on those updates for any shared variable change
		// made by stages within current group  
		Space *waitingLps = genSimplifiedSignalsForGroupTransitionsCode(stream, nextIndentation, syncSignals);
		if (lastWaitingLps == NULL) {
			lastWaitingLps = waitingLps;
		} else if (waitingLps != NULL) {
			if (waitingLps->isParentSpace(lastWaitingLps)) {
				lastWaitingLps = waitingLps;
			}
		}

		//-------------------------------------------------------------------------------------------------

		// commented out this code as the barrier based priliminary implementation does not need this
		/*-------------------------------------------------------------------------------------------------
		// finally if some updater stages can be executed again because current group has finished execution
		// then we need to enable the reactivation signals
		generateCodeForReactivatingDataModifiers(stream, nextIndentation, syncDependencies);
		//-----------------------------------------------------------------------------------------------*/
	}

	// close the while loop if applicable
	if (space != containerSpace) {
		const char *spaceName = space->getName();
		// update the iteration number and next LPU id
		stream << nextIndent.str() << "space" << spaceName << "LpuId = space" << spaceName;
		stream << "Lpu->id" << stmtSeparator;	
		stream << nextIndent.str() << "space" << spaceName << "Iteration++" << stmtSeparator;
		stream << indent.str() << "}\n";
		// at the end remove checkpoint if the container LPS is not the root LPS
		if (!containerSpace->isRoot()) {
			stream << indent.str() << "threadState->removeIterationBound(Space_";
			stream << containerSpace->getName() << ')' << stmtSeparator;
		}
		// exit from the scope
		stream << indent.str() << "} // scope exit for iterating LPUs of Space ";
		stream << space->getName() << "\n";
	}	
}

// A composite stage is a group entry if it has flow stages of multiple LPSes inside or any stage inside it
// is a group entry.
bool CompositeStage::isGroupEntry() {
	for (int i = 0; i < stageList->NumElements(); i++) {
		FlowStage *stage = stageList->Nth(i);
		if (stage->isGroupEntry()) return true;
	}
	List<List<FlowStage*>*> *stageGroups = getConsecutiveNonLPSCrossingStages();
	return (stageGroups->NumElements() > 1 || space != stageGroups->Nth(0)->Nth(0)->getSpace());
}

void CompositeStage::calculateLPSUsageStatistics() {
	for (int i = 0; i < stageList->NumElements(); i++) {
		FlowStage *stage = stageList->Nth(i);
		stage->calculateLPSUsageStatistics();
	}	
}

void CompositeStage::analyzeSynchronizationNeeds() {	
	for (int i = 0; i < stageList->NumElements(); i++) {
		FlowStage *stage = stageList->Nth(i);
		stage->analyzeSynchronizationNeeds();
	}
}

void CompositeStage::retriveExternCodeBlocksConfigs(IncludesAndLinksMap *externConfigMap) {
        for (int i = 0; i < stageList->NumElements(); i++) {
                FlowStage *stage = stageList->Nth(i);
                stage->retriveExternCodeBlocksConfigs(externConfigMap);
        }
}

int CompositeStage::getHighestNestedStageIndex() {
        int stageCount = stageList->NumElements();
        FlowStage *lastStage = stageList->Nth(stageCount - 1);
        CompositeStage *nestedCompositeStage = dynamic_cast<CompositeStage*>(lastStage);
        if (nestedCompositeStage == NULL) {
                return lastStage->getIndex();
        } else {
                return nestedCompositeStage->getHighestNestedStageIndex();
        }
}

void CompositeStage::analyzeSynchronizationNeedsForComposites() {	
	
	// evaluate synchronization requirements and dependencies within its nested composite stages
	for (int i = 0; i < stageList->NumElements(); i++) {
		FlowStage *stage = stageList->Nth(i);
		CompositeStage *compositeStage = dynamic_cast<CompositeStage*>(stage);
		if (compositeStage != NULL) {
			compositeStage->analyzeSynchronizationNeedsForComposites();
		}
	}
	
	// get the first and last flow index within the group to know the boundary of this composite stage
	int beginIndex = stageList->Nth(0)->getIndex();
	int endIndex = getHighestNestedStageIndex();
	
	// now start extracting boundary crossing synchronizations out
	synchronizationReqs = new StageSyncReqs(this);

	// now iterate over the nested stages and extract any boundary crossing synchronization founds within
	// and assign that to current composite stage
	for (int i = 0; i < stageList->NumElements(); i++) {
		FlowStage *stage = stageList->Nth(i);
		StageSyncReqs *nestedSyncs = stage->getAllSyncRequirements();
		List<VariableSyncReqs*> *varSyncList = nestedSyncs->getVarSyncList();
		for (int j = 0; j < varSyncList->NumElements(); j++) {
			VariableSyncReqs *varSyncs = varSyncList->Nth(j);
			List<SyncRequirement*> *syncReqList = varSyncs->getSyncList();
			for (int k =  0; k < syncReqList->NumElements(); k++) {
				SyncRequirement *syncReq = syncReqList->Nth(k);
				//FlowStage *waitingStage = syncReq->getWaitingComputation();
				FlowStage *waitingStage = syncReq->getDependencyArc()->getSignalSink();
				int sinkIndex = waitingStage->getIndex();
				if (sinkIndex < beginIndex || sinkIndex > endIndex) {
					synchronizationReqs->addVariableSyncReq(varSyncs->getVarName(), 
							syncReq, false);
					// update the nesting index of the arc so that it can be discovered later
					// by its nesting index
					syncReq->getDependencyArc()->setNestingIndex(this->repeatIndex);
					// then update the signal source for the arc to indicate inside that it 
					// has been lifted up
					syncReq->getDependencyArc()->setSignalSrc(this);
				}
			}
		}
	}
}

void CompositeStage::deriveSynchronizationDependencies() {
	
	// derive dependencies within its nested composite stages
	for (int i = 0; i < stageList->NumElements(); i++) {
		FlowStage *stage = stageList->Nth(i);
		CompositeStage *compositeStage = dynamic_cast<CompositeStage*>(stage);
		if (compositeStage != NULL) {
			compositeStage->deriveSynchronizationDependencies();
		}
	}

	// get the first and last flow index within the group to know the boundary of this composite stage
	int beginIndex = stageList->Nth(0)->getIndex();
	int endIndex = getHighestNestedStageIndex();
	
	// then check those stage and draw out any synchronization dependencies they have due to changes made in
	// stages outside this composite stage as the dependencies of the current composite stage itself.
	for (int i = 0; i < stageList->NumElements(); i++) {
		FlowStage *stage = stageList->Nth(i);
		StageSyncDependencies *nestedDependencies = stage->getAllSyncDependencies();
		List<SyncRequirement*> *dependencyList = nestedDependencies->getDependencyList();
		for (int j = 0; j < dependencyList->NumElements(); j++) {
			SyncRequirement *sync = dependencyList->Nth(j);
			FlowStage *syncSource = sync->getDependencyArc()->getSignalSrc();
			//FlowStage *syncSource = sync->getDependencyArc()->getSource();
			// if the source's index is within the nesting boundary then the dependency is internal
			// otherwise, it is external and we have to pull it out.
			int syncIndex = syncSource->getIndex();
			if (syncIndex < beginIndex || syncIndex > endIndex) {
				this->syncDependencies->addDependency(sync);
				// a pulled out dependency's arc should be updated to reflect that the current
				// composite stage is its sink
				sync->getDependencyArc()->setSignalSink(this);
			}
		}
	}
}

void CompositeStage::printSyncRequirements() {
	for (int i = 0; i < stageList->NumElements(); i++) {
		FlowStage *stage = stageList->Nth(i);
		stage->printSyncRequirements();
	}	
}

int CompositeStage::assignIndexAndGroupNo(int currentIndex, int currentGroupNo, int currentRepeatCycle) {
	int nextIndex = FlowStage::assignIndexAndGroupNo(currentIndex, currentGroupNo, currentRepeatCycle);
	for (int i = 0; i < stageList->NumElements(); i++) {
		FlowStage *stage = stageList->Nth(i);
		nextIndex = stage->assignIndexAndGroupNo(nextIndex, this->index, currentRepeatCycle);
	}
	return nextIndex;	
}

List<FlowStage*> *CompositeStage::filterOutSyncStages(List<FlowStage*> *originalList) {
	List<FlowStage*> *filteredList = new List<FlowStage*>;
	for (int i = 0; i < originalList->NumElements(); i++) {
		FlowStage *stage = originalList->Nth(i);
		SyncStage *sync = dynamic_cast<SyncStage*>(stage);
		if (sync == NULL) filteredList->Append(stage);
	}
	return filteredList;
}

List<SyncRequirement*> *CompositeStage::getSyncDependeciesOfGroup(List<FlowStage*> *group) {
	List<SyncRequirement*> *syncList = new List<SyncRequirement*>;
	for (int i = 0; i < group->NumElements(); i++) {
		FlowStage *stage = group->Nth(i);
		StageSyncDependencies *sync = stage->getAllSyncDependencies();
		List<SyncRequirement*> *activeDependencies = sync->getActiveDependencies();
		syncList->AppendAll(activeDependencies);
	}
	return syncList;
}

void CompositeStage::generateSyncCodeForGroupTransitions(std::ofstream &stream, int indentation,
                        List<SyncRequirement*> *syncDependencies) {

	// Note that all synchronization we are dealing here is strictly within the boundary of composite stage under
	// concern. Now if there is a synchronization dependency to a nested stage for the execution of another stage
	// that comes after it, it means that the dependency is between a latter iteration on the dependent stage on
	// an earlier iteration on the source stage. This is because, otherwise the dependent stage will execute before
	// the source stage and there should not be any dependency. Now, such dependency is only valid for iterations 
	// except the first one. So there should be a checking on iteration number before we apply waiting on such
	// dependencies. On the other hand, if the source stage executes earlier than the dependent stage then it applies
	// always, i.e., it is independent of any loop iteration, if exists, that sorrounds the composite stage.
	// Therefore, we use following two lists to partition the list of synchronization requirements into forward and
	// backword dependencies. 
	List<SyncRequirement*> *forwardSyncs = new List<SyncRequirement*>;
	List<SyncRequirement*> *backwordSyncs = new List<SyncRequirement*>;

	for (int i = 0; i < syncDependencies->NumElements(); i++) {
		SyncRequirement *sync = syncDependencies->Nth(i);
		DependencyArc *arc = sync->getDependencyArc();
		FlowStage *source = arc->getSource();
		FlowStage *destination = arc->getDestination();
		if (source->getIndex() < destination->getIndex()) {
			forwardSyncs->Append(sync);
		} else backwordSyncs->Append(sync);
	}

	std::ostringstream indentStr;
	for (int i = 0; i < indentation; i++) indentStr << '\t';
	
	if (syncDependencies->NumElements() > 0) {
		stream << std::endl << indentStr.str() << "// waiting for synchronization signals\n";
	}
	
	// Write the code for forwards sync requirements. Now we are just printing the requirements that should be replaced
	// with actual waiting logic later. Note the if condition that ensures that the thread is supposed to wait for the 
	// signal by checking that if it has a valid ppu-id corresponding to the LPS under concern.
	for (int i = 0; i < forwardSyncs->NumElements(); i++) {
		SyncRequirement *sync = forwardSyncs->Nth(i);
		Space *dependentLps = sync->getDependentLps();
		stream << indentStr.str() << "if (threadState->isValidPpu(Space_" << dependentLps->getName();
		stream << ")) {\n";
		stream << indentStr.str() << '\t';
		sync->writeDescriptiveComment(stream, true);
		stream << indentStr.str() << "}\n";
	}

	// write the code for backword sync requirements within an if block
	if (backwordSyncs->NumElements() > 0) {
		stream << indentStr.str() << "if (repeatIteration > 0) {\n";
		for (int i = 0; i < backwordSyncs->NumElements(); i++) {
			SyncRequirement *sync = backwordSyncs->Nth(i);
			Space *dependentLps = sync->getDependentLps();
			stream << indentStr.str() << "\tif (threadState->isValidPpu(Space_" << dependentLps->getName();
			stream << ")) {\n";
			stream << indentStr.str() << "\t\t";
			sync->writeDescriptiveComment(stream, true);
			stream << indentStr.str() << "\t}\n";
		}
		stream << indentStr.str() << "}\n";
	}

	// finally deactive all sync dependencies as they are already been taken care of here	
	for (int i = 0; i < syncDependencies->NumElements(); i++) {
		SyncRequirement *sync = syncDependencies->Nth(i);
		sync->deactivate();
	}
}

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

	std::string stmtSeparator = ";\n";
	std::ostringstream indent;
	for (int i = 0; i < indentation; i++) indent << '\t';

	List<const char*> *counterNameList = getAllOutgoingDependencyNamesAtNestingLevel(nestingIndex);
	if (counterNameList != NULL && counterNameList->NumElements() > 0) {
		stream << std::endl << indent.str() << "// declaration of synchronization counter variables\n";
		for (int i = 0; i < counterNameList->NumElements(); i++) {
			const char *counter = counterNameList->Nth(i);
			stream << indent.str() << "int " << counter << " = 0";
			stream << stmtSeparator;
		}
	}
}

List<SyncRequirement*> *CompositeStage::getSyncSignalsOfGroup(List<FlowStage*> *group) {
	List<SyncRequirement*> *syncList = new List<SyncRequirement*>;
	for (int i = 0; i < group->NumElements(); i++) {
		FlowStage *stage = group->Nth(i);
		List<SyncRequirement*> *activeStageSignals 
				= stage->getAllSyncRequirements()->getAllNonSignaledSyncReqs();
		if (activeStageSignals != NULL) {
			syncList->AppendAll(activeStageSignals);
		}
	}
	return syncList;
}

void CompositeStage::generateSignalCodeForGroupTransitions(std::ofstream &stream, int indentation,
                        List<SyncRequirement*> *syncRequirements) {

	std::string stmtSeparator = ";\n";
	std::string paramSeparator = ", ";
	std::ostringstream indent;
	for (int i = 0; i < indentation; i++) indent << '\t';

	if (syncRequirements->NumElements() > 0) {
		stream << std::endl << indent.str() << "// issuing synchronization signals\n";
	}

	// iterate over all the synchronization signals
	for (int i = 0; i < syncRequirements->NumElements(); i++) {
		SyncRequirement *currentSync = syncRequirements->Nth(i);
		const char *counterVarName = currentSync->getDependencyArc()->getArcName();
		// check if the concerned update did take place
		stream << indent.str() << "if (" << counterVarName << " > 0 && ";
		// also check if the current PPU is a valid candidate for signaling update
		Space *signalingLps = currentSync->getDependencyArc()->getSource()->getSpace();
		stream << "threadState->isValidPpu(Space_" << signalingLps->getName();
		stream << ")) {\n";
		// then signal synchronization (which is now just writing a comment)
		stream << indent.str() << '\t';
		currentSync->writeDescriptiveComment(stream, false);
		// then reset the counter	 
		stream << indent.str() << '\t';
		stream << counterVarName << " = 0" << stmtSeparator;
		stream << indent.str() << "}\n";
	}
}

void CompositeStage::setReactivatorFlagsForSyncReqs() {
	
	// set reactivator flags for all sync primitives operating on this composite stage level
	FlowStage::setReactivatorFlagsForSyncReqs();

	// then set reactivator flags for nested computations 
	for (int i = 0; i < stageList->NumElements(); i++) {
		FlowStage *stage = stageList->Nth(i);
		stage->setReactivatorFlagsForSyncReqs(); 
	} 
}

void CompositeStage::generateCodeForWaitingForReactivation(std::ofstream &stream, int indentation, 
		List<SyncRequirement*> *syncRequirements) {
	
	std::ostringstream indent;
	for (int i = 0; i < indentation; i++) indent << '\t';

	List<DependencyArc*> *reactivatorSyncs = new List<DependencyArc*>;
	for (int i = 0; i < syncRequirements->NumElements(); i++) {
		SyncRequirement *sync = syncRequirements->Nth(i);
		DependencyArc *arc = sync->getDependencyArc();
		if (arc->isReactivator()) {
			reactivatorSyncs->Append(arc);	
		}
	}

	if (reactivatorSyncs->NumElements() > 0) {
		stream << std::endl << indent.str() << "// waiting for signal clearance from readers\n";
	}

	for (int i = 0; i < reactivatorSyncs->NumElements(); i++) {
		DependencyArc *arc = reactivatorSyncs->Nth(i);
		Space *signalingLps = arc->getSignalSrc()->getSpace();
		stream << indent.str() << "if (threadState->isValidPpu(Space_";
		stream << signalingLps->getName();
		stream << ")) {\n";
		stream << indent.str() << '\t';
		stream << "// waiting the signal for " << arc->getArcName() << " to be cleared\n";
		stream << indent.str() << "}\n";
	}			
}

void CompositeStage::generateCodeForReactivatingDataModifiers(std::ofstream &stream, int indentation,
		List<SyncRequirement*> *syncDependencies) {
	
	std::ostringstream indent;
	for (int i = 0; i < indentation; i++) indent << '\t';

	List<DependencyArc*> *reactivatorSyncs = new List<DependencyArc*>;
	for (int i = 0; i < syncDependencies->NumElements(); i++) {
		SyncRequirement *sync = syncDependencies->Nth(i);
		DependencyArc *arc = sync->getDependencyArc();
		if (arc->isReactivator()) {
			reactivatorSyncs->Append(arc);	
		}
	}

	if (reactivatorSyncs->NumElements() > 0) {
		stream << std::endl << indent.str() << "// sending clearance signals to writers\n";
	}

	for (int i = 0; i < reactivatorSyncs->NumElements(); i++) {
		DependencyArc *arc = reactivatorSyncs->Nth(i);
		Space *sinkLps = arc->getSignalSink()->getSpace();
		stream << indent.str() << "if (threadState->isValidPpu(Space_";
		stream << sinkLps->getName();
		stream << ")) {\n";
		stream << indent.str() << '\t';
		stream << "// sending clearance for " << arc->getArcName() << " signal\n";
		stream << indent.str() << "}\n";
	}			
}

List<Space*> *CompositeStage::getWaitingLpsInFollowingStages(int currentIndex, List<List<FlowStage*>*> *stageGroups) {
	
	List<Space*> *waitingLpsList = new List<Space*>;
	int totalGroups = stageGroups->NumElements();
	int groupIndex = (currentIndex + 1) % totalGroups;
	while (groupIndex != currentIndex) {
		List<FlowStage*> *currentGroup = stageGroups->Nth(groupIndex);
		List<SyncRequirement*> *syncSignals = getSyncSignalsOfGroup(currentGroup);
		for (int i = 0; i < syncSignals->NumElements(); i++) {
			SyncRequirement *sync = syncSignals->Nth(i);
			Space *syncSpanLps = sync->getSyncSpan();
			bool lpsIncluded = false;
			for (int j = 0; j < waitingLpsList->NumElements(); j++) {
				Space *waitingLps = waitingLpsList->Nth(j);
				if (strcmp(syncSpanLps->getName(), waitingLps->getName()) == 0) {
					lpsIncluded = true;
					break;
				}
			}
			if (lpsIncluded == false) {
				waitingLpsList->Append(syncSpanLps);
			}
		}
		groupIndex = (groupIndex + 1) % totalGroups;
	}
	return waitingLpsList;
}

List<Space*> *CompositeStage::genSimplifiedWaitingForReactivationCode(std::ostream &stream, int indentation,
		List<SyncRequirement*> *syncRequirements, 
		Space *lastWaitingLps, List<Space*> *allWaitingLpses) {

	std::string stmtSeparator = ";\n";	
	std::ostringstream indent;
	for (int i = 0; i < indentation; i++) indent << '\t';

	List<Space*> *newWaitingLpses = new List<Space*>;


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

	if (lowestWaitingLps == NULL) return newWaitingLpses;

	if (syncRequirements->NumElements() > 0) {
		stream << std::endl << indent.str();
		stream << "// barriers to ensure all readers have finished reading last update\n";
	}

	// even after finding a lowest waiting LPS, we may need to put more than one sync barrier as there may
	// be cross LPS syncing requirements. That is, PPUs of two unrelated LPSes need to wait.
	bool syncBeingPut = false;
	for (int i = 0; i < syncRequirements->NumElements(); i++) {

		SyncRequirement *sync = syncRequirements->Nth(i);
		Space *syncSpanLps = sync->getSyncSpan();
		if (lowestWaitingLps->isParentSpace(syncSpanLps)) continue;
		bool redundentWait = false;
		for (int j = 0; j < allWaitingLpses->NumElements(); j++) {
			Space *waitingLps = allWaitingLpses->Nth(j);
			if (strcmp(syncSpanLps->getName(), waitingLps->getName()) == 0 
					||  syncSpanLps->getMappedPpsId() == waitingLps->getMappedPpsId()) {
				redundentWait = true;
				break;
			}
		}
		if (redundentWait == true) continue;
		else if (strcmp(lowestWaitingLps->getName(), syncSpanLps->getName()) == 0) {
			if (syncBeingPut == false) {
				syncBeingPut = true;
			} else {
				continue;
			}
		}
		
		FlowStage *sourceStage = sync->getDependencyArc()->getSource();
		Space *signalingLps = sourceStage->getSpace();
		if (syncSpanLps->isParentSpace(signalingLps) 
				&& syncSpanLps->getMappedPpsId() == signalingLps->getMappedPpsId()) {
			continue;
		}

		newWaitingLpses->Append(syncSpanLps);

		stream << indent.str() << "if (threadState->isValidPpu(Space_";
		stream << syncSpanLps->getName();
		stream << ")) {\n";	
		stream << indent.str() << '\t';
		stream << "threadSync->" << sync->getReverseSyncName() << "->wait()";
		stream << stmtSeparator;
		stream << indent.str() << "}\n";
	}

	return newWaitingLpses;
}

Space *CompositeStage::genSimplifiedSignalsForGroupTransitionsCode(std::ostream &stream, int indentation,
		List<SyncRequirement*> *syncRequirements) {

	std::string stmtSeparator = ";\n";
	std::string paramSeparator = ", ";
	std::ostringstream indent;
	for (int i = 0; i < indentation; i++) indent << '\t';

	if (syncRequirements->NumElements() > 0) {
		stream << std::endl << indent.str() << "// resolving synchronization dependencies\n";
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
		
		Assert(syncSpanLps->getMappedPpsId() > -1 && signalingLps->getMappedPpsId() > -1);
		if (syncSpanLps->isParentSpace(signalingLps) 
				&& syncSpanLps->getMappedPpsId() == signalingLps->getMappedPpsId()) {
			currentSync->getDependencyArc()->deactivate();
			continue;
		}

		lastWaitingLps = syncSpanLps;
		lastSignalingLps = signalingLps;

		const char *counterVarName = currentSync->getDependencyArc()->getArcName();
		
		// check if the concerned update did take place
		stream << indent.str() << "if (" << counterVarName << " > 0 && ";
		// also check if the current PPU is a valid candidate for signaling update
		stream << "threadState->isValidPpu(Space_" << signalingLps->getName();
		stream << ")) {\n";
		// then signal synchronization
		stream << indent.str() << '\t';
		stream << "threadSync->" << currentSync->getSyncName() << "->signal(";
		FlowStage *signalSource = currentSync->getDependencyArc()->getSignalSrc();
		if (signalSource->getRepeatIndex() > 0) stream << "repeatIteration";
		else stream << "0";
		stream << ")" << stmtSeparator;
		// then reset the counter	 
		stream << indent.str() << '\t';
		stream << counterVarName << " = 0" << stmtSeparator;

		// the waiting is in an else block coupled with the signaling if block as the current implementation
		// of synchronization primitives does not support the PPU (or PPUs) in the signaling block to also be
		// among the list of waiting PPUs.  
		stream << indent.str() << "} else if (";
		FlowStage *waitingStage = currentSync->getWaitingComputation();
		stream << "threadState->isValidPpu(Space_" << syncSpanLps->getName();
		stream << ")) {\n";
		stream << indent.str() << '\t';
		stream << "threadSync->" << currentSync->getSyncName() << "->wait(";
		FlowStage *signalSink = currentSync->getDependencyArc()->getSignalSink();
		if (signalSink->getRepeatIndex() > 0) stream << "repeatIteration";
		else stream << "0";
		stream << ")" << stmtSeparator;
		stream << indent.str() << "}\n";
	}

	return lastWaitingLps;
}

//------------------------------------------------- Repeat Cycle ------------------------------------------------------/

RepeatCycle::RepeatCycle(int index, Space *space, RepeatCycleType type, Expr *executeCond) 
		: CompositeStage(index, space, NULL) {
	this->type = type;
	this->repeatCond = executeCond;
}

void RepeatCycle::addSyncStagesOnReturn(List<FlowStage*> *stageList) {
	CompositeStage::addSyncStagesOnReturn(stageList);
	// For repeat cycles iterating on a space that has data structures with overlapped/ghost regions, extra
	// reappearance sync may be needed. Therefore, repeat cycle overrides this method and add an additional
	// clause.
        SyncStage *reappearanceSync = SyncStageGenerator::generateReappearanceSyncStage(space, accessMap);
        if (reappearanceSync != NULL) {
        	addStageAtEnd(reappearanceSync);
        }
}

// Here we do the dependency analysis twice to take care of any new dependencies that may arise due to the
// return from the last stage of repeat cycle to the first one. Note that the consequence of this double iteration
// may seems to be addition of superfluous dependencies when we have nesting of repeat cycles. However that should
// not happen as there is a redundancy checking mechanism in place where we add dependency arcs to flow stages. So
// unwanted arcs will be dropped off. 
void RepeatCycle::performDependencyAnalysis(PartitionHierarchy *hierarchy) {
	CompositeStage::performDependencyAnalysis(hierarchy);	
	FlowStage::performDependencyAnalysis(repeatConditionAccessMap, hierarchy);
	CompositeStage::performDependencyAnalysis(hierarchy);	
}

void RepeatCycle::calculateLPSUsageStatistics() {
	Iterator<VariableAccess*> iterator = repeatConditionAccessMap->GetIterator();
	VariableAccess *accessLog;
	while ((accessLog = iterator.GetNextValue()) != NULL) {
		if (!accessLog->isContentAccessed()) continue;
		const char *varName = accessLog->getName();
		DataStructure *structure = space->getLocalStructure(varName);
		AccessFlags *accessFlags = accessLog->getContentAccessFlags();
		LPSVarUsageStat *usageStat = structure->getUsageStat();
		// As repeat condition is supposed to be evaluated multiple times, any data structure been used
		// here should be marked as multiple access. Thus, we add access information twice.
		if (accessFlags->isRead() || accessFlags->isWritten()) {
			usageStat->addAccess();
			usageStat->addAccess();
		}
	}
	// Two calls have been made for evaluating the usage statistics in nested computation stages assuming that
	// a repeat cycle will in most cases executes at least twice.
	CompositeStage::calculateLPSUsageStatistics();
	CompositeStage::calculateLPSUsageStatistics();
}

int RepeatCycle::assignIndexAndGroupNo(int currentIndex, int currentGroupNo, int currentRepeatCycle) {
	
	this->index = currentIndex;
	this->groupNo = currentGroupNo;
	this->repeatIndex = currentRepeatCycle;

	int nextIndex = currentIndex + 1;
	int nextRepeatIndex = currentRepeatCycle + 1;

	for (int i = 0; i < stageList->NumElements(); i++) {
		FlowStage *stage = stageList->Nth(i);
		nextIndex = stage->assignIndexAndGroupNo(nextIndex, this->index, nextRepeatIndex);
	}
	return nextIndex;
	
}

void RepeatCycle::generateInvocationCode(std::ofstream &stream, int indentation, Space *containerSpace) {
	
	std::string stmtSeparator = ";\n";

	// if it is not a sub-partition repeat block then we have to add a while or for loop depending on the
	// condition
	if (type == Conditional_Repeat) {
		
		std::ostringstream indent;
		for (int i = 0; i < indentation; i++) indent << '\t';
		
		// create a scope for repeat loop
		stream << std::endl << indent.str() << "{ // scope entrance for repeat loop\n";

		// declare a repeat iteration number tracking variable
		stream << indent.str() << "int repeatIteration = 0" << stmtSeparator;

		// get the name of the lpu for the execution LPS
		std::ostringstream lpuName;
		lpuName << "space" << space->getName() << "Lpu->";

		// If the repeat condition involves accessing metadata of some task global array then we need 
		// to create local copies of its metadata so that name transformer can work properly 
        	if (isLpsDependent()) {
			List<const char*> *localArrays = filterInArraysFromAccessMap(repeatConditionAccessMap);
			for (int i = 0; i < localArrays->NumElements(); i++) {
        			const char *arrayName = localArrays->Nth(i);
        			ArrayDataStructure *array = (ArrayDataStructure*) space->getStructure(arrayName);
        			int dimensions = array->getDimensionality();
        			stream << indent.str() << "Dimension ";
                		stream  << arrayName << "PartDims[" << dimensions << "];\n";
        			stream << indent.str() << "Dimension ";
                		stream  << arrayName << "StoreDims[" << dimensions << "];\n";
                		for (int j = 0; j < dimensions; j++) {
                			stream << indent.str();
               				stream << arrayName << "PartDims[" << j << "] = " << lpuName.str();
                        		stream << arrayName << "PartDims[" << j << "].partition;\n";
                			stream << indent.str();
               				stream << arrayName << "StoreDims[" << j << "] = " << lpuName.str();
                        		stream << arrayName << "PartDims[" << j << "].storage;\n";
        			}
        		}
		}

		// update the name transformer for probable array access within repeat condition
		ntransform::NameTransformer::transformer->setLpuPrefix(lpuName.str().c_str());
		
		// if the repeat condition is a logical expression then it is a while loop in the source code
		if (dynamic_cast<LogicalExpr*>(repeatCond) != NULL) {
			std::ostringstream condition;
			repeatCond->translate(condition, indentation, 0, space);
			stream << indent.str() << "while (" << condition.str() << ") {\n";
			// declare all synchronization counter variables here that will be updated inside 
			// repeat loop 
			declareSynchronizationCounters(stream, indentation + 1, this->repeatIndex + 1);
			// invoke the repeat body
			CompositeStage::generateInvocationCode(stream, indentation + 1, containerSpace);
			// increase the loop iteration counter
			stream << indent.str() << "\trepeatIteration++" << stmtSeparator;
			stream << indent.str() << "}\n";
		// otherwise it is a for loop based on a range condition that we need to translate
		} else {
			// translate the range expression into a for loop
			RangeExpr *rangeExpr = dynamic_cast<RangeExpr*>(repeatCond);
			std::ostringstream rangeLoop;
			rangeExpr->generateLoopForRangeExpr(rangeLoop, indentation, space);
			stream << rangeLoop.str();
			// declare all synchronization counter variables here that will be updated inside 
			// repeat loop 
			declareSynchronizationCounters(stream, indentation + 1, this->repeatIndex + 1);
			// translate the repeat body
			CompositeStage::generateInvocationCode(stream, indentation + 1, containerSpace);
			// increase the loop iteration counter
			stream << indent.str() << "\trepeatIteration++" << stmtSeparator;
			// close the range loop
			stream << indent.str() << "}\n";
		}
		// exit the scope created for the repeat loop 
		stream << indent.str() << "} // scope exit for repeat loop\n";
	} else {
		// declare any synchronization counter variable needed inside the repeat loop
		declareSynchronizationCounters(stream, indentation, this->repeatIndex + 1);
		
		// TODO probably we need to maintain a repeat iteration counter in this case two. Then we should
		// change this straightforward superclass's code execution strategy. We should investigate this
		// in the future 
		CompositeStage::generateInvocationCode(stream, indentation, containerSpace);
	}
}

bool RepeatCycle::isLpsDependent() {
	VariableAccess *accessLog;
	Iterator<VariableAccess*> iterator = repeatConditionAccessMap->GetIterator();
	while ((accessLog = iterator.GetNextValue()) != NULL) {
		if (!(accessLog->isContentAccessed() 
			|| (accessLog->isMetadataAccessed() 
				&& accessLog->isLocalAccess()))) continue;
		const char *varName = accessLog->getName();
		DataStructure *structure = space->getStructure(varName);
		ArrayDataStructure *array = dynamic_cast<ArrayDataStructure*>(structure);
		if (array != NULL) return true;
	}
	return false;
}

