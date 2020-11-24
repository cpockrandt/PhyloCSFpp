CC = g++

CFLAGS = -O0 -g -fsanitize=address -fno-omit-frame-pointer -Wall -pedantic # -static-libasan

TARGET = PhyloCSF++

all: $(TARGET)

$(TARGET): prototype.cpp
	$(CC) $(CFLAGS) -o $(TARGET) -lgsl -lgslcblas prototype.cpp

clean:
	$(RM) $(TARGET)
