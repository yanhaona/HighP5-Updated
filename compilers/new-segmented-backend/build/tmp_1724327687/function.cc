
/*--------------------------------------------------------------------------------------------------------------
common header files for different purposes
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



/*--------------------------------------------------------------------------------------------------------------
header files needed to execute external code blocks
--------------------------------------------------------------------------------------------------------------*/

#include <time.h>
#include <cstdlib>
#include <math.h>


/*--------------------------------------------------------------------------------------------------------------
init_rand instances
--------------------------------------------------------------------------------------------------------------*/

void init_rand_0() {

	//---------------------------------------------------------------------------Local Variable Declarations


	//-----------------------------------------------------------------------------------------Function Body

	{ // starting scope for an external code block

		{ // external code block starts
 srand(time(NULL)); 
		} // external code block ends

	} // ending scope for the external code block

}


/*--------------------------------------------------------------------------------------------------------------
perform_sampling instances
--------------------------------------------------------------------------------------------------------------*/

int perform_sampling_0(Rectangle cell, int seed, int trial_count) {

	//---------------------------------------------------------------------------Local Variable Declarations
	int cell_height;
	int cell_width;
	int internal_points;
	int trial;


	//-----------------------------------------------------------------------------------------Function Body
	cell_height = ((cell.top - cell.bottom) + 1);
	cell_width = ((cell.right - cell.left) + 1);
	internal_points = 0;
	trial = 0;
	do {

		{ // starting scope for an external code block

			{ // external code block starts

				int x = rand_r((unsigned int *) &seed) % cell_width + cell.left;
				int y = rand_r((unsigned int *) &seed) % cell_height + cell.bottom;
				
				// tested polynomial is 10 sin x^2 + 50 cos y^3
				double result = 10 * sin(pow(x, 2)) + 50 * cos(pow(y, 3));
				if (result <= 0.0) {
					internal_points++; 
				}
			
			} // external code block ends

		} // ending scope for the external code block

		trial = (trial + 1);
	} while((trial < trial_count));
	return internal_points;
}

