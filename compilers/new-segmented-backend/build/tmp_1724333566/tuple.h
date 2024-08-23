#ifndef _H_tuple
#define _H_tuple

#include <iostream>
#include <vector>

#include "../../../common-libs/domain-obj/structure.h"

class ProgramArgs;
class Rectangle;
class MCAEPartition;

class ProgramArgs {
  public:
	int b;
	int cell_length;
	int grid_dim;
	int max_rounds;
	int points_per_cell;
	double precision;
	ProgramArgs() {
		b = 0;
		cell_length = 0;
		grid_dim = 0;
		max_rounds = 0;
		points_per_cell = 0;
		precision = 0;
	}
};

class Rectangle {
  public:
	int top;
	int bottom;
	int left;
	int right;
	Rectangle() {
		top = 0;
		bottom = 0;
		left = 0;
		right = 0;
	}
};

class MCAEPartition {
  public:
	int b;
	MCAEPartition() {
		b = 0;
	}
};

#endif
