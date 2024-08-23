// PCubeS Description of a 128 Core Machine in Tozammel's University
// processor model name: Intel(R) Xeon(R) Gold 6438Y+ 
// Total Shared Memory 251 GB
// Data width 64 bit

//--------------------------------------------------------------------------------------
//Space #Number : 	$Space-Name	(#PPU-Count)	// Comment
//--------------------------------------------------------------------------------------
Space 	4:  		CPU 		(1)		// 251 GB RAM Total			
Space 	3: 		NUMA-Node 	(2) 		// 60 MB L-3 Cache		
Space 	2: 		Core-Pair 	(32)		// 2 MB L-2 Cache (1 floating point unit per core-pair ????????)
Space 	1<core>:	Core		(2)		// 48 KB L-1 Cache 
