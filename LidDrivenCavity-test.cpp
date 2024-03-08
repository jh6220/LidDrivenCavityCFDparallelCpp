#include "LidDrivenCavity.h"

#define BOOST_TEST_MODULE LidDrivenCavity
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE(SetReynoldsNumber) {
    LidDrivenCavity* solver = new LidDrivenCavity();
    solver->SetReynoldsNumber(20.0);
    BOOST_CHECK_EQUAL(solver->nu, 1.0/20.0);
    }