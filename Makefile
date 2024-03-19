CXX= OMPI_CXX=g++-10 mpicxx
CXXFLAGS= -std=c++11 -Wall -O3 -g
# CXX= mpicxx
# CXXFLAGS= -std=c++11 -Wall -O3

# Include directories
INCDIR=-I/opt/homebrew/Cellar/openblas/0.3.26/include -I/opt/homebrew/Cellar/boost/1.84.0_1/include
# Library directories
LIBDIR=-L/opt/homebrew/Cellar/openblas/0.3.26/lib -L/opt/homebrew/Cellar/boost/1.84.0_1/lib
# Libraries to link against
LIBS=-lblas -llapack  -lboost_program_options

default: solver

.PHONY: clean, run, runOMP, doc, test

LidDrivenCavitySolver.o: LidDrivenCavitySolver.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) -fopenmp -o LidDrivenCavitySolver.o -c LidDrivenCavitySolver.cpp 

LidDrivenCavity.o: LidDrivenCavity.cpp LidDrivenCavity.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -fopenmp -o LidDrivenCavity.o -c LidDrivenCavity.cpp

SolverCG.o: SolverCG.cpp SolverCG.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -fopenmp -o SolverCG.o -c SolverCG.cpp

solver: LidDrivenCavitySolver.o LidDrivenCavity.o SolverCG.o
	$(CXX) -fopenmp -o solver LidDrivenCavitySolver.o LidDrivenCavity.o SolverCG.o $(LIBDIR) $(LIBS)

run: solver
	OMP_NUM_THREADS=3 mpiexec -n 1 ./solver --Lx 1 --Ly 1 --Nx 201 --Ny 201 --Re 1000 --dt 0.005 --T 0.5

doc:
	doxygen Doxyfile

lidDrivenCavity-test: lidDrivenCavity-test.cpp LidDrivenCavity.o SolverCG.o
	$(CXX) $(CXXFLAGS) $(INCDIR) -fopenmp -o lidDrivenCavity-test lidDrivenCavity-test.cpp LidDrivenCavity.o SolverCG.o $(LIBDIR) $(LIBS) -lboost_unit_test_framework

SolverCG-test: SolverCG-test.cpp SolverCG.o
	$(CXX) $(CXXFLAGS) $(INCDIR) -fopenmp -o SolverCG-test SolverCG-test.cpp SolverCG.o $(LIBDIR) $(LIBS) -lboost_unit_test_framework

test: lidDrivenCavity-test SolverCG-test
	OMP_NUM_THREADS=1 mpiexec -n 4 ./lidDrivenCavity-test
	OMP_NUM_THREADS=1 mpiexec -n 4 ./SolverCG-test

clean:
	-rm -f *.o solver lidDrivenCavity-test SolverCG-test