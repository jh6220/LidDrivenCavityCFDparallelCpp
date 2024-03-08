#include "SolverCG.h"

#define BOOST_TEST_MODULE SolverCG
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE(Solve) {
    int Nx = 4;
    int Ny = 4;
    double dx = 0.1;
    double dy = 0.1;
    SolverCG* cg  = new SolverCG(Nx, Ny, dx, dy);
    double w[] = {0,0,0,0,
                0,1,2,0,
                0,3,4,0,
                0,0,0,0};
    double* s = new double[Nx*Ny]();
    cg->Solve(w, s);
    double s_true[] = {0,0,0,0,
                    0,0.00875, 0.01125,0,
                    0,0.01375, 0.01625,0,
                    0,0,0,0};

    double tol = 0.0001;
    
    for (int i = 0; i < Nx*Ny; i++) 
    {
        BOOST_CHECK_CLOSE(s[i], s_true[i], tol);
    }
}