void readAFromFile(const char *filePath) {

	int rowsPerProcess = (aDims[0].length + processCount - 1) / processCount;
	int rowStart = processId * rowsPerProcess;
	int rowEnd = rowStart + rowsPerProcess - 1;
	if (rowEnd >= aDims[0].length) {
		rowEnd = aDims[0].length - 1;
	}
	int rowCount = rowEnd - rowStart + 1;
	int localEntries = rowCount * aDims[1].length;
	a = new double[localEntries];

	TypedInputStream<double> *stream = new TypedInputStream<double>(filePath);
	int storeIndex = 0;
	stream->open();
	List<int> *indexList = new List<int>();
        for (int i = rowStart; i <= rowEnd; i++) {
		indexList->clear();
		indexList->Append(i);
		indexList->Append(0);
		a[storeIndex] = stream->readElement(indexList);
		storeIndex++;
		int count = 1;
                while (count < aDims[1].length) {
			a[storeIndex] = stream->readNextElement();
			storeIndex++;
			count++;
		}
	}
	stream->close();

	delete indexList;
	delete stream;
}

void readBFromFile(const char *filePath) {

	int bSize = bDims[0].length * bDims[1].length;
	b = new double[bSize];
	TypedInputStream<double> *stream = new TypedInputStream<double>(filePath);
	int storeIndex = 0;
	stream->open();
	for (int i = 0; i < bSize; i++) {
		b[storeIndex] = stream->readNextElement();
	}
	stream->close();
	delete stream;
}


	// parse command line arguments
	blockSize = atoi(argv[1]);
	bool fileWriteMode (atoi(argv[4]) == 1);	
	const char *filePathA = argv[2];
	std::ifstream fileA(filePathA);
        if (!fileA.is_open()) {
                std::cout << "could not open the specified file\n";
                std::exit(EXIT_FAILURE);
        }
	readArrayDimensionInfoFromFile(fileA, 2, aDims);
	fileA.close();
	const char *filePathB = argv[3];
	std::ifstream fileB(filePathB);
        if (!fileB.is_open()) {
                std::cout << "could not open the specified file\n";
                std::exit(EXIT_FAILURE);
        }
	readArrayDimensionInfoFromFile(fileB, 2, bDims);
	fileB.close();

	// read input matrices
	readAFromFile(filePathA);
	readBFromFile(filePathB);

	struct timeval memEnd;
        gettimeofday(&memEnd, NULL);




