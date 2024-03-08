CXX=g++
CXXFLAGS=-std=c++11 -Wall
# Include directories
INCDIR=-I/opt/homebrew/Cellar/openblas/0.3.26/include -I/opt/homebrew/Cellar/boost/1.84.0_1/include
# Library directories
LIBDIR=-L/opt/homebrew/Cellar/openblas/0.3.26/lib -L/opt/homebrew/Cellar/boost/1.84.0_1/lib
# Libraries to link against
LIBS=-lblas -llapack -lcblas -lboost_program_options

default: solver

LidDrivenCavitySolver.o: LidDrivenCavitySolver.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) -o LidDrivenCavitySolver.o -c LidDrivenCavitySolver.cpp 

LidDrivenCavity.o: LidDrivenCavity.cpp LidDrivenCavity.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -o LidDrivenCavity.o -c LidDrivenCavity.cpp

SolverCG.o: SolverCG.cpp SolverCG.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -o SolverCG.o -c SolverCG.cpp

solver: LidDrivenCavitySolver.o LidDrivenCavity.o SolverCG.o
	$(CXX) -o solver LidDrivenCavitySolver.o LidDrivenCavity.o SolverCG.o $(LIBDIR) $(LIBS)

doc:
	doxygen Doxyfile

lidDrivenCavity-test: lidDrivenCavity-test.cpp LidDrivenCavity.o SolverCG.o
	$(CXX) $(CXXFLAGS) $(INCDIR) -o lidDrivenCavity-test lidDrivenCavity-test.cpp LidDrivenCavity.o SolverCG.o $(LIBDIR) $(LIBS) -lboost_unit_test_framework

SolverCG-test: SolverCG-test.cpp SolverCG.o
	$(CXX) $(CXXFLAGS) $(INCDIR) -o SolverCG-test SolverCG-test.cpp SolverCG.o $(LIBDIR) $(LIBS) -lboost_unit_test_framework

run-tests: lidDrivenCavity-test SolverCG-test
	./lidDrivenCavity-test
	./SolverCG-test
