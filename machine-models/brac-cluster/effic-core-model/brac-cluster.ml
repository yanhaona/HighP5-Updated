// 3nd PCubeS Description of a 4 node cluster of  14-Core Intel Machine with 8 single-threaded 
// cores. This model does not consider the 6 dual-threaded performance cores
// processor model name    : 13th Gen Intel(R) Core(TM) i5-13600KF
// total memory: 15 GB

//--------------------------------------------------------------------------------------
//Space #Number : 	        $Space-Name	(#PPU-Count)	// Comment
//--------------------------------------------------------------------------------------
Space 	5<unit><segment>:  	Cluster		(1)		// 			
Space 	4:  		        Node 		(1)		// 15 GB RAM per CPU			
Space 	3: 		        L3-Cache 	(1) 		// 24 MB L3 cache shared by all cores		
Space   2: 	                Core-Group     	(2)		// 4 MB Cache shared among 4 cores
Space   1<core>:                Core            (4)             // 24 KB Cache per core
