all: test_IsoThresholdGenerator test_IsoOrderedGenerator

test_IsoThresholdGenerator:
	clang++ -std=c++17 -x c++ test_IsoThresholdGenerator.c -o isoThreshold -lpthread

test_IsoOrderedGenerator:
	clang++ -std=c++17 -x c++ test_IsoOrderedGenerator.c -o isoOrdered -lpthread
clean:
	rm -f isoOrdered isoThreshold
