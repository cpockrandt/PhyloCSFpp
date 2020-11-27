CC = g++

CFLAGS = -O0 -g -fsanitize=address -fno-omit-frame-pointer -Wall -pedantic

TARGET = PhyloCSFpp

all: $(TARGET)

$(TARGET): prototype.cpp src/alignment_reader.hpp src/maf_parser.hpp src/fit.hpp src/ecm.hpp src/models.hpp src/fixed_lik.hpp src/omega.hpp src/instance.hpp src/newick.hpp src/translation.hpp
	$(CC) $(CFLAGS) -o $(TARGET) -lgsl -lgslcblas prototype.cpp

clean:
	$(RM) $(TARGET)
