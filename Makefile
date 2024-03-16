CXX=mpicxx
CXXFLAGS= -std=c++11 -Wall -O3
# Include directories
INCDIR=-I/opt/homebrew/Cellar/openblas/0.3.26/include -I/opt/homebrew/Cellar/boost/1.84.0_1/include
# Library directories
LIBDIR=-L/opt/homebrew/Cellar/openblas/0.3.26/lib -L/opt/homebrew/Cellar/boost/1.84.0_1/lib
# Libraries to link against
LIBS=-lblas -llapack -lcblas -lboost_program_options

default: solver

.PHONY: clean, run, runOMP, doc, run-tests, run-testSolver

LidDrivenCavitySolver.o: LidDrivenCavitySolver.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) -o LidDrivenCavitySolver.o -c LidDrivenCavitySolver.cpp 

LidDrivenCavity.o: LidDrivenCavity.cpp LidDrivenCavity.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -o LidDrivenCavity.o -c LidDrivenCavity.cpp

SolverCG.o: SolverCG.cpp SolverCG.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -o SolverCG.o -c SolverCG.cpp

solver: LidDrivenCavitySolver.o LidDrivenCavity.o SolverCG.o
	$(CXX) -o solver LidDrivenCavitySolver.o LidDrivenCavity.o SolverCG.o $(LIBDIR) $(LIBS)

run: solver
	mpiexec -n 9 ./solver --Lx 1 --Ly 1 --Nx 201 --Ny 201 --Re 1000 --dt 0.005 --T 5

LidDrivenCavitySolverOMP.o: LidDrivenCavitySolver.cpp
	$(CXX) $(CXXFLAGS) -fopenmp $(INCDIR) -o LidDrivenCavitySolverOMP.o -c LidDrivenCavitySolver.cpp 

LidDrivenCavityOMP.o: LidDrivenCavity.cpp LidDrivenCavity.h
	$(CXX) $(CXXFLAGS) -fopenmp $(INCDIR) -o LidDrivenCavityOMP.o -c LidDrivenCavity.cpp

SolverCGOMP.o: SolverCG.cpp SolverCG.h
	$(CXX) $(CXXFLAGS) -fopenmp $(INCDIR) -o SolverCGOMP.o -c SolverCG.cpp

solverOMP: LidDrivenCavitySolverOMP.o LidDrivenCavityOMP.o SolverCGOMP.o
	$(CXX) -fopenmp -o solverOMP LidDrivenCavitySolverOMP.o LidDrivenCavityOMP.o SolverCGOMP.o $(LIBDIR) $(LIBS)

runOMP: solverOMP
	export OMP_NUM_THREADS=4
	mpiexec ./solverOMP --Lx 1 --Ly 1 --Nx 201 --Ny 201 --Re 1000 --dt 0.005 --T 0.5

doc:
	doxygen Doxyfile

lidDrivenCavity-test: lidDrivenCavity-test.cpp LidDrivenCavity.o SolverCG.o
	$(CXX) $(CXXFLAGS) $(INCDIR) -o lidDrivenCavity-test lidDrivenCavity-test.cpp LidDrivenCavity.o SolverCG.o $(LIBDIR) $(LIBS) -lboost_unit_test_framework

SolverCG-test: SolverCG-test.cpp SolverCG.o
	$(CXX) $(CXXFLAGS) $(INCDIR) -o SolverCG-test SolverCG-test.cpp SolverCG.o $(LIBDIR) $(LIBS) -lboost_unit_test_framework

run-tests: lidDrivenCavity-test SolverCG-test
	./lidDrivenCavity-test
	./SolverCG-test

solverCG-test-parallel: SolverCG-test.cpp SolverCG.cpp SolverCG.h
	mpicxx $(INCDIR) -o solverCG-test-parallel SolverCG-test.cpp SolverCG.cpp $(LIBDIR) $(LIBS) -lboost_unit_test_framework -std=c++11

run-testSolver: solverCG-test-parallel
	mpiexec -n 4 ./solverCG-test-parallel

clean:
	-rm -f *.o solver solverOMP lidDrivenCavity-test SolverCG-test solverCG-test-parallel LidDrivenCavitySolverOMP.o LidDrivenCavityOMP.o SolverCGOMP.o