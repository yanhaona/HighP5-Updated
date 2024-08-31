#include <stdlib.h>
#include <map>
#include <string.h>

using namespace std;

int Brac_All[20] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9
                , 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};

int Brac_Perf[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

int Brac_Effic[8] = {12, 13, 14, 15, 16, 17, 18, 19};

int Tozammel[128] = {0, 64, 24, 88, 48, 112, 8, 72, 32, 96
                , 16, 80, 40, 104, 56, 120, 4, 68, 20, 84
                , 36, 100, 52, 116, 12, 76, 28, 92, 44, 108
                , 60, 124, 58, 122, 42, 106, 26, 90, 2, 66
                , 50, 114, 34, 98, 10, 74, 62, 126, 18, 82
                , 54, 118, 30, 94, 6, 70, 38, 102, 14, 78
                , 46, 110, 22, 86, 1, 65, 17, 81, 41, 105
                , 25, 89, 9, 73, 33, 97, 49, 113, 57, 121
                , 5, 69, 21, 85, 37, 101, 53, 117, 13, 77
                , 29, 93, 45, 109, 61, 125, 59, 123, 43, 107
                , 27, 91, 3, 67, 51, 115, 35, 99, 11, 75
                , 63, 127, 19, 83, 55, 119, 31, 95, 7, 71
                , 39, 103, 15, 79, 47, 111, 23, 87};



int *getCoreNumberingArray(const char *machine, int threadCount) {

		
	map<int, int*> bracCpuMap = {
		{1, Brac_All},
		{20, Brac_All},
		{2, Brac_Effic},
		{8, Brac_Effic},
		{6, Brac_Perf},
		{12, Brac_Perf}
	};

	if (strcmp(machine, "tozammel") == 0) {
		return Tozammel;
	} else if (strcmp(machine, "brac") == 0) {
		return bracCpuMap.at(threadCount);
	} else {
		return NULL;
	}
}

int getCoreCount(const char *machine, int threadCount) {
	
	map<int, int> bracCpuMap = {
		{1, 20},
		{20, 20},
		{2, 8},
		{8, 8},
		{6, 12},
		{12, 12}
	};

	if (strcmp(machine, "tozammel") == 0) {
		return 128;
	} else if (strcmp(machine, "brac") == 0) {
		return bracCpuMap.at(threadCount);
	} else {
		return 0;
	}
}
