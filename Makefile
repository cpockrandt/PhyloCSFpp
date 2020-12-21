CFLAGS = -O0 -g -std=c++11 -fsanitize=address -fno-omit-frame-pointer -Wall -pedantic -Wno-write-strings
LDLIBS = `pkg-config --libs gsl` `pkg-config --cflags gsl`
TARGET = PhyloCSFpp

all: $(TARGET)

$(TARGET): prototype.cpp src/estimate_hmm_parameter.hpp src/create_tracks.hpp src/alignment_reader.hpp src/maf_parser.hpp src/fit.hpp src/ecm.hpp src/models.hpp src/fixed_lik.hpp src/omega.hpp src/instance.hpp src/newick.hpp src/translation.hpp
	$(CXX) $(CFLAGS) -o $(TARGET) $(LDLIBS) prototype.cpp

clean:
	$(RM) $(TARGET)
