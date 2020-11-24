CC = g++

CFLAGS = -O0 -g -fsanitize=address -fno-omit-frame-pointer -Wall -pedantic # -static-libasan

TARGET = PhyloCSF++

# SOURCES := $(shell find ./src -name '*.hpp')

all: $(TARGET)

$(TARGET): prototype.cpp src/alignment_reader.hpp src/ecm.hpp src/fixed_lik.hpp src/instance.hpp src/newick.hpp src/translation.hpp
	$(CC) $(CFLAGS) -o $(TARGET) -lgsl -lgslcblas prototype.cpp

clean:
	$(RM) $(TARGET)
