#ifndef _H_tuple
#define _H_tuple

#include <iostream>
#include <vector>

#include "../../../common-libs/domain-obj/structure.h"

class ProgramArgs;
class FPSPartition;

class ProgramArgs {
  public:
	const char* input_file;
	int iterations;
	int k;
	int l;
	int m;
	int n;
	const char* output_file;
	int p1;
	int p2;
	ProgramArgs() {
		input_file = NULL;
		iterations = 0;
		k = 0;
		l = 0;
		m = 0;
		n = 0;
		output_file = NULL;
		p1 = 0;
		p2 = 0;
	}
};

class FPSPartition {
  public:
	int k;
	int l;
	int m;
	int n;
	int p1;
	int p2;
	FPSPartition() {
		k = 0;
		l = 0;
		m = 0;
		n = 0;
		p1 = 0;
		p2 = 0;
	}
};

#endif
