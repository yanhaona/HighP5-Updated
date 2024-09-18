// A modified PCubeS Description of a 128 Core Machine in Tozammel's University
// The modification is done so that we can test 4, 8, and 16 way parallel versions
// of HighP5 sample codes in that machine. Note that this modification does not
// indicate any weakness in PCubeS type architecture. Rather, it compensates for
// a compiler limitation at this moment that does not allow us to select a subset
// of PPS from a level to restrict parallelism.

// processor model name: Intel(R) Xeon(R) Gold 6438Y+ 
// Total Shared Memory 251 GB
// Data width 64 bit

//--------------------------------------------------------------------------------------
//Space #Number : 	$Space-Name	(#PPU-Count)	// Comment
//--------------------------------------------------------------------------------------
Space 	8<unit><segment>:  	CPU 		(1)		// 251 GB RAM Total			
Space 	7: 		        NUMA-Node 	(2) 		// 60 MB L-3 Cache		
Space 	6: 		        Added-PPS4 	(2)		// 
Space 	5: 		        Added-PPS3 	(2)		//
Space 	4: 		        Added-PPS2 	(2)		//
Space 	3: 		        Added-PPS1 	(2)		//
Space 	2: 		        Core-Pair 	(2)		// 2 MB L-2 Cache (1 floating point unit per core-pair ????????)
Space 	1<core>:	        Core		(2)		// 48 KB L-1 Cache 
