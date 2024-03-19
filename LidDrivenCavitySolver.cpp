#include <iostream>
#include <mpi.h>
#include <omp.h>
using namespace std;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "LidDrivenCavity.h"

/*! \mainpage Lid driven cavity flow solver documentation
 *
 * \section makefile commands
 *
 * you can compile and run the code with the following command:
 * \code make run \endcode
 * 
 * you can run the tests with the following command:
 * \code make test \endcode
 * 
 * you can clean the execution and linked files with:
 * \code make clean \endcode
 * 
 * you can generate the documentation with:
 * \code make docs \endcode
 *
 * \section Code structure
 * 
 * The code is structured in 3 files:
 *
 * \subsection LidDrivenCavitySolver 
 * 
 * This file contains the main function. It reads the input parameters and
 * creates an instance of the LidDrivenCavity class. It then calls the
 * Initialise and Integrate methods of the LidDrivenCavity class.
 * It also writes the initial and final solutions to files.
 *  
 * \subsection LidDrivenCavity
 * 
 * This file contains the implementation of the LidDrivenCavity class. This
 * objects stores the simulation parameters and the solution and has the methos
 * initialise, integrate and write the solution to a file.
 * 
 * \subsection SolverCG
 * 
 * This file contains the implementation of the SolverCG class. This class
 * implements the Conjugate Gradient method to solve a linear system of equations.
 * 
 */

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    /// Parse input
    po::options_description opts(
        "Solver for the 2D lid-driven cavity incompressible flow problem");
    opts.add_options()
        ("Lx",  po::value<double>()->default_value(1.0),
                 "Length of the domain in the x-direction.")
        ("Ly",  po::value<double>()->default_value(1.0),
                 "Length of the domain in the y-direction.")
        ("Nx",  po::value<int>()->default_value(9),
                 "Number of grid points in x-direction.")
        ("Ny",  po::value<int>()->default_value(9),
                 "Number of grid points in y-direction.")
        ("dt",  po::value<double>()->default_value(0.01),
                 "Time step size.")
        ("T",   po::value<double>()->default_value(1.0),
                 "Final time.")
        ("Re",  po::value<double>()->default_value(10),
                 "Reynolds number.")
        ("verbose",    "Be more verbose.")
        ("help",       "Print help message.");


    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << opts << endl;
        return 0;
    }

    /// Initialise the solver and set the simulation parameters
    LidDrivenCavity* solver = new LidDrivenCavity();
    solver->SetDomainSize(vm["Lx"].as<double>(), vm["Ly"].as<double>());
    solver->SetGridSize(vm["Nx"].as<int>(),vm["Ny"].as<int>());
    solver->SetTimeStep(vm["dt"].as<double>());
    solver->SetFinalTime(vm["T"].as<double>());
    solver->SetReynoldsNumber(vm["Re"].as<double>());

    solver->PrintConfiguration();

    solver->Initialise();

    solver->WriteSolution("ic.txt");

    /// Run the simulation
    solver->Integrate();

    solver->WriteSolution("final.txt");

    MPI_Finalize();
	return 0;
}
