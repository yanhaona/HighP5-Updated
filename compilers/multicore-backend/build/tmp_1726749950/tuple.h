#ifndef _H_tuple
#define _H_tuple

#include <iostream>
#include <vector>

#include "../../src/runtime/structure.h"

class BLUFEnvironment;
class BLUFPartition;

class BLUFEnvironment {
  public:
	double* a;
	PartDimension aDims[2];
	double* u;
	PartDimension uDims[2];
	double* l;
	PartDimension lDims[2];
	int* p;
	PartDimension pDims[1];
	const char* name;
	BLUFEnvironment() {
		a = NULL;
		u = NULL;
		l = NULL;
		p = NULL;
		name = NULL;
	}
};

class BLUFPartition {
  public:
	int b;
	BLUFPartition() {
		b = 0;
	}
};

#endif
