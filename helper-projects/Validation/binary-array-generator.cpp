#include <cstdlib>
#include <cctype>
#include <stdlib.h>
#include <string.h>
#include <deque>
#include <sys/time.h>
#include <time.h>
#include <string>
#include <random>

#include "sparse_matrix_gen.h"
#include "fileUtility.h"
#include "structures.h"
#include "stream.h"
#include "list.h"


using namespace std;

namespace bin_array_gen {

int getElementCount(List<Dimension*> *dimList) {
	
	int length = 1;
	for (int i = 0; i < dimList->NumElements(); i++) {
		length = length * dimList->Nth(i)->length;
	}
	return length;
}

void saveIntArrayInFile(const char* fileName, List<Dimension*> *dimList, int *data) {
	
	TypedOutputStream<int> stream = TypedOutputStream<int>(fileName, dimList, true);
	stream.open();
	int elementCount = getElementCount(dimList);
	for (int i = 0; i < elementCount; i++) {
		stream.writeNextElement(data[i]);
	}
	stream.close();
}

void saveFloatArrayInFile(const char* fileName, List<Dimension*> *dimList, float *data) {
	
	TypedOutputStream<float> stream = TypedOutputStream<float>(fileName, dimList, true);
	stream.open();
	int elementCount = getElementCount(dimList);
	for (int i = 0; i < elementCount; i++) {
		stream.writeNextElement(data[i]);
	}
	stream.close();
}

void saveDoubleArrayInFile(const char* fileName, List<Dimension*> *dimList, double *data) {
	
	TypedOutputStream<double> stream = TypedOutputStream<double>(fileName, dimList, true);
	stream.open();
	int elementCount = getElementCount(dimList);
	for (int i = 0; i < elementCount; i++) {
		stream.writeNextElement(data[i]);
	}
	stream.close();
}


}// end of namespace

using namespace bin_array_gen;

void generateArray(int dataType, List<Dimension*> *dimList, const char* fileName) {

        srand(time(NULL));
	int elementCount = getElementCount(dimList);
	Dimension total = Dimension();
	total.length = elementCount;

	if (dataType == 1) {
		int *data = allocate<int>(1, &total);
		for (int i = 0; i < elementCount; i++) {
			data[i] = rand() * 11 % 1000;	
		}
		saveIntArrayInFile(fileName, dimList, data);
	} else if (dataType == 2) {
		float *data = allocate<float>(1, &total);
		for (int i = 0; i < elementCount; i++) {
			data[i] = (rand() % 100 * 1.0) / (rand() % 100 + 1.0);	
		}
		saveFloatArrayInFile(fileName, dimList, data);
	} else if (dataType == 3) {
		double *data = allocate<double>(1, &total);
		for (int i = 0; i < elementCount; i++) {
			data[i] = (rand() % 100 * 1.0) / (rand() % 100 + 1.0);	
		}
		saveDoubleArrayInFile(fileName, dimList, data);
	}
}

int mainBAGen(int argc, const char* argv[]) {

	if (argc < 4) {
                std::cout << "provide the following information\n";
		std::cout << "1. data type of the array\n";
		std::cout << "\t1 = int\n";
		std::cout << "\t2 = float\n";
		std::cout << "\t3 = double\n";
		std::cout << "2. number of dimensions in the array\n";
		std::cout << "3. name of the file to store the array\n";
		std::cout << "4. then provide the length of individual dimensions one after another\n";
                std::exit(EXIT_FAILURE);
        }

	int dataType = atoi(argv[1]);
	int dimCount = atoi(argv[2]);
	const char* fileName = argv[3];
	if (argc < dimCount + 4) {
		std::cout << "3. you haven't provided the length of all the dimensions\n";
                std::exit(EXIT_FAILURE);
	}
	List<Dimension*> *dimList = new List<Dimension*>;
	for (int d = 0; d < dimCount; d++) {
		Dimension *dim = new Dimension();
		dim->length = atoi(argv[4 + d]);
		dim->range.min = 0;
		dim->range.max = dim->length - 1;
		dimList->Append(dim);
	}

	generateArray(dataType, dimList, fileName);
	std::cout << "Saved generated array in file: " << fileName << "\n";
	return 0;
}

