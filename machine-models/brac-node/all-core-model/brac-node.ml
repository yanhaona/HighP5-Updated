// PCubeS Description of a 14-Core Intel Machine with 6 2-hyperthreaded cores
// processor model name    : 13th Gen Intel(R) Core(TM) i5-13600KF
// total memory: 15 GB

//--------------------------------------------------------------------------------------
//Space #Number : 	$Space-Name	(#PPU-Count)	// Comment
//--------------------------------------------------------------------------------------
Space 	3:  		CPU 		(1)		// 15 GB RAM in CPU			
Space 	2: 		L3-Cache 	(1) 		// 24 MB L3 cache shared by all cores		
Space 	1<core>: 	Core     	(20)		// 1 MB Cache per core 
