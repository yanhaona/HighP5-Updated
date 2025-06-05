#include "../../../utils/code_constant.h"
#include "../../../utils/name_transformer.h"
#include "../../../../../../common-libs/domain-obj/constant.h"
#include "../../../../../../frontend/src/syntax/ast_expr.h"
#include "../../../../../../frontend/src/semantics/task_space.h"
#include "../../../../../../frontend/src/semantics/computation_flow.h"
#include "../../../../../../frontend/src/static-analysis/sync_stat.h"
#include "../../../../../../frontend/src/static-analysis/data_dependency.h"

#include <fstream>
#include <sstream>

void RepeatControlBlock::generateInvocationCode(std::ofstream &stream, int indentation, Space *containerSpace) {

	std::ostringstream indentStream;
	for (int i = 0; i < indentation; i++) indentStream << indent;
	std::string indentStr = indentStream.str();

	// create a scope for repeat loop
	stream << std::endl << indentStr << "{ // scope entrance for repeat loop\n";

	// declare a repeat iteration number tracking variable
	stream << indentStr << "int repeatIteration = 0" << stmtSeparator;

	// get the name of the lpu for the execution LPS
	std::ostringstream lpuName;
	lpuName << "space" << space->getName() << "Lpu->";

	// If the repeat condition involves accessing metadata of some task global array then we need to create 
	// local copies of its metadata so that name transformer can work properly 
	if (isLpsDependent()) {
		List<const char*> *localArrays = filterInArraysFromAccessMap();
		for (int i = 0; i < localArrays->NumElements(); i++) {
			const char *arrayName = localArrays->Nth(i);
			ArrayDataStructure *array = (ArrayDataStructure*) space->getStructure(arrayName);
			int dimensions = array->getDimensionality();
			stream << indentStr << "Dimension ";
			stream  << arrayName << "PartDims[" << dimensions << "]" << stmtSeparator;
			stream << indentStr << "Dimension ";
			stream  << arrayName << "StoreDims[" << dimensions << "]" << stmtSeparator;
			for (int j = 0; j < dimensions; j++) {
				stream << indentStr;
				stream << arrayName << "PartDims[" << j << "] = " << lpuName.str();
				stream << arrayName << "PartDims[" << j << "].partition" << stmtSeparator;
				stream << indentStr;
				stream << arrayName << "StoreDims[" << j << "] = " << lpuName.str();
				stream << arrayName << "PartDims[" << j << "].storage" << stmtSeparator;
			}
		}
	}

	// update the name transformer for probable array access within repeat condition
	ntransform::NameTransformer::transformer->setLpuPrefix(lpuName.str().c_str());

	if (this->type == Condition_Repeat) {
		std::ostringstream condStream;
		condition->translate(condStream, indentation, 0, space);
		stream << indentStr << "while (" << condStream.str() << ") {\n";
	} else {
		// translate the range expression into a for loop
		RangeExpr *rangeExpr = dynamic_cast<RangeExpr*>(condition);
		std::ostringstream rangeLoop;
		rangeExpr->generateLoopForRangeExpr(rangeLoop, indentation, space);
		stream << rangeLoop.str();
	}
	// declare all synchronization counter variables here that will be updated inside repeat loop 
	declareSynchronizationCounters(stream, indentation + 1, this->repeatIndex + 1);

	// accumulate the backward dependencies from the repeat loop
        List<SyncRequirement*> *backwardDependencies = new List<SyncRequirement*>;
	List<List<FlowStage*>*> *stageGroups = getConsecutiveNonLPSCrossingStages();
        for (int i = 0; i < stageGroups->NumElements(); i++) {

                List<FlowStage*> *currentGroup = stageGroups->Nth(i);

                // retrieve all data dependencies, and sort them to ensure waitings for updates happen in order
                List<SyncRequirement*> *dataDependencies = getDataDependeciesOfGroup(currentGroup);
                dataDependencies = SyncRequirement::sortList(dataDependencies);

                int segmentedPPS = space->getSegmentedPPS();
                List<SyncRequirement*> *commDependencies = new List<SyncRequirement*>;
                List<SyncRequirement*> *syncDependencies = new List<SyncRequirement*>;
                SyncRequirement::separateCommunicationFromSynchronizations(segmentedPPS,
                                dataDependencies, commDependencies, syncDependencies);

		// filter out the backward dependencies from a later stage of the repeat to an earlier stage
        	for (int i = 0; i < commDependencies->NumElements(); i++) {
                	SyncRequirement *comm = commDependencies->Nth(i);
                	if (!comm->isActive()) continue;
                	DependencyArc *arc = comm->getDependencyArc();
                	FlowStage *source = arc->getSource();
                	FlowStage *destination = arc->getDestination();
                	if (source->getIndex() >= destination->getIndex()) {
                		backwardDependencies->Append(comm);
			}
        	}
	}

	// translate the repeat body
	CompositeStage::generateInvocationCode(stream, indentation + 1, containerSpace);
	// increase the loop iteration counter
	stream << indentStr << "\trepeatIteration++" << stmtSeparator;
	// close the range loop
	stream << indentStr << "}\n";

	// exit the scope created for the repeat loop 
	stream << indentStr << "} // scope exit for repeat loop\n";

	// generate data receive code for all inverse data dependencies at the end of the repeat loop
        if (backwardDependencies->NumElements() > 0) {
                stream << std::endl << indentStr << "// waiting on backward dependencies at repeat loop end\n";
		for (int i = 0; i < backwardDependencies->NumElements(); i++) {

                        SyncRequirement *comm = backwardDependencies->Nth(i);

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
                        stream << indentStr;
                        stream << "if (threadState->isValidPpu(Space_" << dependentLps->getName();
                        stream << ")) {\n";
                        stream << indentStr << indent;
                        stream << "Communicator *communicator = threadState->getCommunicator(\"";
                        stream << comm->getDependencyArc()->getArcName() << "\")" << stmtSeparator;
                        stream << indentStr << indent << "if (communicator != NULL) {\n";
                        stream << indentStr << doubleIndent;
                        stream << "communicator->receive(REQUESTING_COMMUNICATION";
                        stream << paramSeparator << "commCounter" << commIndex << ")" << stmtSeparator;
                        stream << indentStr << indent << "}\n";
                        stream << indentStr << "}\n";

			// The counter should be advanced regardless of this PPU's participation in communication to keep the
                	// counter value uniform across all PPUs and segments. The exeception is when the data send signal has
                	// been replaced by some other communicator. This is because, we will have two receive calls for the
                	// replacement communicator then and we would want to make the second call to bypass any data processing.
                	if (!signalReplaced) {
                        	stream << indentStr << "commCounter" << commIndex << "++" << stmtSeparator;
                	}
        	}
	}
}
