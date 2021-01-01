CC = g++

CFLAGS = -O0 -g -Wall -pedantic -Wno-write-strings # -fsanitize=address -fno-omit-frame-pointer

TARGET = PhyloCSFpp

all: $(TARGET) automatic_tests

$(TARGET): prototype.cpp src/alignment_reader.hpp src/additional_scores.hpp src/maf_parser.hpp src/fit.hpp src/ecm.hpp src/models.hpp src/fixed_lik.hpp src/omega.hpp src/instance.hpp src/newick.hpp src/translation.hpp
	$(CC) $(CFLAGS) -o $(TARGET) -lgsl -lgslcblas prototype.cpp

automatic_tests: automatic_tests.cpp src/alignment_reader.hpp src/additional_scores.hpp src/maf_parser.hpp src/fit.hpp src/ecm.hpp src/models.hpp src/fixed_lik.hpp src/omega.hpp src/instance.hpp src/newick.hpp src/translation.hpp
	$(CC) $(CFLAGS) -o automatic_tests -lgsl -lgslcblas automatic_tests.cpp

clean:
	$(RM) $(TARGET) automatic_tests
