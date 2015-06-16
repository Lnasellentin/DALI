# You may need to adjust the compiler, compiler flags or library paths. -lgsl and -lgslcblas refer to the Gnu Scientific Library, which must be installed.

SHELL = /bin/sh

S = ../source/

# Settings for C++ compiler.
CPP = g++
CPPFLAGS = -g -O0 
CPPINCLUDES =

DALI_SOURCE = $(S)DaliBase.h $(S)DaliBase.cpp $(S)DaliGen.cpp $(S)Utils.cpp $(S)DaliPaint.cpp $(S)DaliPaint.h
DALI_OBJ = $(S)DaliBase.o $(S)DaliGen.o $(S)Utils.o $(S)DaliPaint.o


# Linker settings
LINKER_CMD = g++
AR_CMD = ar
LDFLAGS = $(CPPFLAGS)
LDLIBS = -lgsl -lgslcblas 

GEN_DEP = TestGen.cpp $(DALI_SOURCE)
GEN_OBJ = TestGen.o $(DALI_OBJ)

TUT_DEP = Tutorial.cpp $(DALI_SOURCE)
TUT_OBJ = Tutorial.o $(DALI_OBJ)

LOAD_DEP = LoadGenFromFile.cpp $(DALI_SOURCE)
LOAD_OBJ = LoadGenFromFile.o $(DALI_OBJ)



.PHONY: all clean test

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c $*.cpp -o $*.o

gen: $(GEN_OBJ)
	$(LINKER_CMD) $(LDFLAGS) -o gen $(GEN_OBJ) $(LDLIBS)

tutorial: $(TUT_OBJ)
	  $(LINKER_CMD) $(LDFLAGS) -o tutorial $(TUT_OBJ) $(LDLIBS)

load: $(LOAD_OBJ)
	  $(LINKER_CMD) $(LDFLAGS) -o load $(LOAD_OBJ) $(LDLIBS)

clean:
	rm -f $(S)*.o; rm -f *.o 
