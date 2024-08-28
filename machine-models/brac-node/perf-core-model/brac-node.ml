// 2nd PCubeS Description of a 14-Core Intel Machine with 6 2-hyperthreaded cores
// This model does not consider the 8 efficient single-threaded cores
// processor model name    : 13th Gen Intel(R) Core(TM) i5-13600KF
// total memory: 15 GB

//--------------------------------------------------------------------------------------
//Space #Number : 	$Space-Name	(#PPU-Count)	// Comment
//--------------------------------------------------------------------------------------
Space 	4:  		CPU 		(1)		// 15 GB RAM in CPU			
Space 	3: 		L3-Cache 	(1) 		// 24 MB L3 cache shared by all cores		
Space   2: 	        Core-Pair     	(6)		// 2 MB Cache shared between core pair
Space   1<core>:        Core            (2)             // 24 KB Cache per core
