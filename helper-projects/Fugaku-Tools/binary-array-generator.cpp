#include <cstdlib>
#include <cctype>
#include <stdlib.h>
#include <string.h>
#include <deque>
#include <sys/time.h>
#include <time.h>
#include <string>
#include <random>
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

int main(int argc, const char* argv[]) {

	if (argc < 3) {
                std::cout << "provide the following information\n";
		std::cout << "1. number of dimensions in the double type array\n";
		std::cout << "3. name of the file to store the array\n";
		std::cout << "4. then provide the length of individual dimensions one after another\n";
                std::exit(EXIT_FAILURE);
        }

	int dimCount = atoi(argv[1]);
	const char* fileName = argv[2];
	if (argc < dimCount + 3) {
		std::cout << "3. you haven't provided the length of all the dimensions\n";
                std::exit(EXIT_FAILURE);
	}
	int *dimLengths = new int[dimCount];
	long int elementCount = 1;
	for (int d = 0; d < dimCount; d++) {
		dimLengths[d] = atoi(argv[3 + d]);
		elementCount = elementCount * dimLengths[d];
	}
        
	srand(time(NULL));
	ofstream stream(fileName, ios_base::binary);

        for (int i = 0; i < dimCount; i++) {
        	ostringstream str;
                if (i > 0) str << "*";
                str << dimLengths[i];
                int strLength = str.str().length();
                stream.write(str.str().c_str(), sizeof(char) * strLength);
                stream.flush();
        }
        char lineEnd = '\n';
        stream.write(&lineEnd, sizeof(char));
        stream.flush();


	for (int i = 0; i < elementCount; i++) {
		double item = (rand() % 100 * 1.0) / (rand() % 100 + 1.0);	
		stream.write(reinterpret_cast<char*>(&item), sizeof(item));
	}
	stream.close();

	std::cout << "Saved generated array in file: " << fileName << "\n";
	return 0;
}

