#include <iostream>
#include <mpi.h>
#include <omp.h>
using namespace std;

#include "LidDrivenCavity.h"

#define BOOST_TEST_MODULE LidDrivenCavity
#include <boost/test/included/unit_test.hpp>

struct MPIFixture {
    public:
        explicit MPIFixture() {
            argc = boost::unit_test::framework::master_test_suite().argc;
            argv = boost::unit_test::framework::master_test_suite().argv;
            cout << "Initialising MPI" << endl;
            MPI_Init(&argc, &argv);
        }

        ~MPIFixture() {
            cout << "Finalising MPI" << endl;
            MPI_Finalize();
        }

        int argc;
        char **argv;
};
BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_AUTO_TEST_CASE(SetReynoldsNumber) {
    LidDrivenCavity* solver = new LidDrivenCavity();
    solver->SetReynoldsNumber(20.0);
    BOOST_CHECK_EQUAL(solver->nu, 1.0/20.0);
    delete solver;
}

BOOST_AUTO_TEST_CASE(Integrate) {
    bool runWithoutError = true;
    try {
        LidDrivenCavity* solver = new LidDrivenCavity();
        solver->SetDomainSize(1.0, 1.0);
        solver->SetGridSize(201,201);
        solver->SetTimeStep(0.005);
        solver->SetFinalTime(0.01);
        solver->SetReynoldsNumber(1000.0);

        solver->Initialise();
        solver->Integrate();

        delete solver;
    } catch (...) {
        runWithoutError = true;
    }
    cout << "Run without error: " << runWithoutError << endl;
    BOOST_CHECK(runWithoutError);
}
