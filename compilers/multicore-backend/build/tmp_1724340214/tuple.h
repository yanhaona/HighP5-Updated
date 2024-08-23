#ifndef _H_tuple
#define _H_tuple

#include <iostream>
#include <vector>

#include "../../src/runtime/structure.h"

class MMEnvironment;
class MMPartition;

class MMEnvironment {
  public:
	double* a;
	PartDimension aDims[2];
	double* b;
	PartDimension bDims[2];
	double* c;
	PartDimension cDims[2];
	const char* name;
	MMEnvironment() {
		a = NULL;
		b = NULL;
		c = NULL;
		name = NULL;
	}
};

class MMPartition {
  public:
	int k;
	int l;
	int q;
	MMPartition() {
		k = 0;
		l = 0;
		q = 0;
	}
};

#endif
