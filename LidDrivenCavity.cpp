#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>
using namespace std;

#include <cblas.h>
#include <mpi.h>

#define IDX(I,J) ((J)*Nx + (I))

#include "LidDrivenCavity.h"
#include "SolverCG.h"

LidDrivenCavity::LidDrivenCavity()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &world_size);
    MPI_Comm_size(MPI_COMM_WORLD, &world_rank);
    
    // Check if the world size is p^2
    world_size_root = sqrt(world_size);
    if (world_size != world_size_root*world_size_root)
    {
        if world_rank == 0 {
            cout << "The number of processes must be a perfect square" << endl;
        }
        MPI_Finalize();
        exit(-1);
    }
}

LidDrivenCavity::~LidDrivenCavity()
{
    CleanUp();
}

void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
    this->Lx = xlen;
    this->Ly = ylen;
    UpdateDxDy();
}

void LidDrivenCavity::SetGridSize(int nx, int ny)
{
    this->Nx = nx;
    this->Ny = ny;
    UpdateDxDy();

    // Define 2D Cartesian topology 
    // Needs to be commented 
    int dims[2] = {world_size_root, world_size_root};
    int periods[2] = {0, 0};
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &mygrid);

    int keep[2];

    MPI_Comm_rank(mygrid, &mygrid_rank);
    MPI_Cart_coords(mygrid, mygrid_rank, 2, coords);
    keep[0] = 0;
    keep[1] = 1; // keep rows in subgrid
    MPI_Cart_sub(mygrid, keep, &yCoordComm);
    keep[0] = 1; // keep columns in subgrid
    keep[1] = 0;
    MPI_Cart_sub(mygrid, keep, &xCoordComm);

    // Define the grid of the local MPI process
    kx = (Nx-2)%world_size_root;
    Nx_local = (Nx-2-kx)/world_size_root + 2;
    if (coords[0] < kx) {
            Nx_local = Nx_local + 1;
    }
    ky = Ny%world_size_root;
    Ny_local = (Ny-2-ky)/world_size_root + 2;
    if (coords[1] < ky) {
        Ny_local = Ny_local + 1;
    }
}

void LidDrivenCavity::SetTimeStep(double deltat)
{
    this->dt = deltat;
}

void LidDrivenCavity::SetFinalTime(double finalt)
{
    this->T = finalt;
}

void LidDrivenCavity::SetReynoldsNumber(double re)
{
    this->Re = re;
    this->nu = 1.0/re;
}

void LidDrivenCavity::Initialise()
{
    CleanUp();

    v   = new double[Npts_local]();
    s   = new double[Npts_local]();
    tmp = new double[Npts_local]();
    cg  = new SolverCG(Nx_local, Ny_local, dx, dy);

    dataB_bottom_sent = new double[Nx_local];
    dataB_top_sent = new double[Nx_local];
    dataB_left_sent = new double[Ny_local];
    dataB_right_sent = new double[Ny_local];

    dataB_bottom_recv = new double[Nx_local];
    dataB_top_recv = new double[Nx_local];
    dataB_left_recv = new double[Ny_local];
    dataB_right_recv = new double[Ny_local];
}

void LidDrivenCavity::Integrate()
{
    int NSteps = ceil(T/dt);
    for (int t = 0; t < NSteps; ++t)
    {
        std::cout << "Step: " << setw(8) << t
                  << "  Time: " << setw(8) << t*dt
                  << std::endl;
        Advance(t);
    }
}

/**
 * @brief Writes the solution into a file
 * @param file  name of the file to which the solution is written
*/
void LidDrivenCavity::WriteSolution(std::string file)
{
    double* u0 = new double[Nx*Ny]();
    double* u1 = new double[Nx*Ny]();
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            u0[IDX(i,j)] =  (s[IDX(i,j+1)] - s[IDX(i,j)]) / dy;
            u1[IDX(i,j)] = -(s[IDX(i+1,j)] - s[IDX(i,j)]) / dx;
        }
    }
    for (int i = 0; i < Nx; ++i) {
        u0[IDX(i,Ny-1)] = U;
    }

    std::ofstream f(file.c_str());
    std::cout << "Writing file " << file << std::endl;
    int k = 0;
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            k = IDX(i, j);
            f << i * dx << " " << j * dy << " " << v[k] <<  " " << s[k] 
              << " " << u0[k] << " " << u1[k] << std::endl;
        }
        f << std::endl;
    }
    f.close();

    delete[] u0;
    delete[] u1;
}


void LidDrivenCavity::PrintConfiguration()
{
    cout << "Grid size: " << Nx << " x " << Ny << endl;
    cout << "Spacing:   " << dx << " x " << dy << endl;
    cout << "Length:    " << Lx << " x " << Ly << endl;
    cout << "Grid pts:  " << Npts << endl;
    cout << "Timestep:  " << dt << endl;
    cout << "Steps:     " << ceil(T/dt) << endl;
    cout << "Reynolds number: " << Re << endl;
    cout << "Linear solver: preconditioned conjugate gradient" << endl;
    cout << endl;
    if (nu * dt / dx / dy > 0.25) {
        cout << "ERROR: Time-step restriction not satisfied!" << endl;
        cout << "Maximum time-step is " << 0.25 * dx * dy / nu << endl;
        exit(-1);
    }
}


void LidDrivenCavity::CleanUp()
{
    if (v) {
        delete[] v;
        delete[] s;
        delete[] tmp;
        delete cg;
    }
}


void LidDrivenCavity::UpdateDxDy()
{
    dx = Lx / (Nx-1);
    dy = Ly / (Ny-1);
    Npts = Nx * Ny;
    Npts_local = Nx_local * Ny_local;
}

void LidDrivenCavity::UpdateDataWithParallelProcesses(double* data, int tag) {

    // Collect boundary data to buffers and sent to neighbours
    if (coords[0] != 0) {
        // left
        for (int i = 0; i < Ny_local; ++i) {
            dataB_left_sent[i] = data[IDX(1, i)];
        }
        MPI_Isent(dataB_left_sent, Ny_local, MPI_DOUBLE, coords[0]-1, tag, xCoordComm, MPI_STATUS_IGNORE);
    }
    if (coords[0] != world_size_root-1) {
        // right
        for (int i = 0; i < Ny_local; ++i) {
            dataB_right_sent[i] = data[IDX(Nx_local-2, j)];
        }
        MPI_Isent(dataB_right_sent, Ny_local, MPI_DOUBLE, coords[0]+1, tag, xCoordComm, MPI_STATUS_IGNORE);
    }
    if (coords[1] != 0) {
        // top
        for (int j = 0; j < Nx_local; ++j) {
            dataB_top_sent[j] = data[IDX(j,1)];
        }
        MPI_Isent(dataB_top_sent, Nx_local, MPI_DOUBLE, coords[1]-1, tag, yCoordComm, MPI_STATUS_IGNORE);
    }
    if (coords[1] != world_size_root-1) {
        // bottom
        for (int j = 0; j < Nx_local; ++j) {
            dataB_bottom_sent[j] = data[IDX(j, Ny_local-2)];
        }
        MPI_Isent(dataB_bottom_sent, Nx_local, MPI_DOUBLE, coords[1]+1, tag, yCoordComm, MPI_STATUS_IGNORE);
    }

    // Receive boundary data from neighbours
    if (coords[0] != 0) {
        MPI_Recv(dataB_left_recv, Ny_local, MPI_DOUBLE, coords[0]-1, tag, xCoordComm, MPI_STATUS_IGNORE);
        for (int i = 0; i < Ny_local; ++i) {
            data[IDX(0,i)] = dataB_left_recv[i];
        }
    }
    if (coords[0] != world_size_root-1) {
        MPI_Recv(dataB_right_recv, Ny_local, MPI_DOUBLE, coords[0]+1, tag, xCoordComm, MPI_STATUS_IGNORE);
        for (int i = 0; i < Ny_local; ++i) {
            data[IDX(Nx_local-1,i)] = dataB_right_recv[i];
        }
    }
    if (coords[1] != 0) {
        MPI_Recv(dataB_top_recv, Nx_local, MPI_DOUBLE, coords[1]-1, tag, yCoordComm, MPI_STATUS_IGNORE);
        for (int j = 0; j < Nx_local; ++j) {
            data[IDX(j,0)] = dataB_top_recv[j];
        }
    }
    if (coords[1] != world_size_root-1) {
        MPI_Recv(dataB_bottom_recv, Nx_local, MPI_DOUBLE, coords[1]+1, tag, yCoordComm, MPI_STATUS_IGNORE);
        for (int j = 0; j < Nx_local; ++j) {
            data[IDX(j,Ny_local-1)] = dataB_bottom_recv[j];
        }
    }
}


void LidDrivenCavity::Advance(int idxT)
{
    double dxi  = 1.0/dx;
    double dyi  = 1.0/dy;
    double dx2i = 1.0/dx/dx;
    double dy2i = 1.0/dy/dy;

    int n_tags = 4;
    
    // Vorticity boundary conditions
    if (coords[0] == 0) {
        for (int i = 1; i < Ny_local-1; ++i) {
            // left
            v[IDX(0,j)]    = 2.0 * dx2i * (s[IDX(0,j)]    - s[IDX(1,j)]);
        }
    }

    if (coords[0] == world_size_root-1) {
        for (int i = 1; i < Ny_local-1; ++i) {
            // right
            v[IDX(Nx_local-1,j)] = 2.0 * dx2i * (s[IDX(Nx_local-1,j)] - s[IDX(Nx_local-2,j)]);
        }
    }

    if (coords[1] == 0) {
        for (int j = 1; j < Nx_local-1; ++j) {
            // top
            v[IDX(i,0)]    = 2.0 * dy2i * (s[IDX(i,0)]    - s[IDX(i,1)]);
        }
    }

    if (coords[1] == world_size_root-1) {
        for (int j = 1; j < Nx_local-1; ++j) {
            // bottom
            v[IDX(i,Ny_local-1)] = 2.0 * dy2i * (s[IDX(i,Ny_local-1)] - s[IDX(i,Ny_local-2)])
                           - 2.0 * dyi*U;
        }
    }

    // Compute interior vorticity
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            v[IDX(i,j)] = dx2i*(
                    2.0 * s[IDX(i,j)] - s[IDX(i+1,j)] - s[IDX(i-1,j)])
                        + 1.0/dy/dy*(
                    2.0 * s[IDX(i,j)] - s[IDX(i,j+1)] - s[IDX(i,j-1)]);
        }
    }

    //Exchange vorticity data with parallel processes
    UpdateDataWithParallelProcesses(v, n_tags*idxT+0);

    // Time advance vorticity
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            v[IDX(i,j)] = v[IDX(i,j)] + dt*(
                ( (s[IDX(i+1,j)] - s[IDX(i-1,j)]) * 0.5 * dxi
                 *(v[IDX(i,j+1)] - v[IDX(i,j-1)]) * 0.5 * dyi)
              - ( (s[IDX(i,j+1)] - s[IDX(i,j-1)]) * 0.5 * dyi
                 *(v[IDX(i+1,j)] - v[IDX(i-1,j)]) * 0.5 * dxi)
              + nu * (v[IDX(i+1,j)] - 2.0 * v[IDX(i,j)] + v[IDX(i-1,j)])*dx2i
              + nu * (v[IDX(i,j+1)] - 2.0 * v[IDX(i,j)] + v[IDX(i,j-1)])*dy2i);
        }
    }
    

    // Sinusoidal test case with analytical solution, which can be used to test
    // the Poisson solver
    /*
    const int k = 3;
    const int l = 3;
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            v[IDX(i,j)] = -M_PI * M_PI * (k * k + l * l)
                                       * sin(M_PI * k * i * dx)
                                       * sin(M_PI * l * j * dy);
        }
    }
    */

    // Solve Poisson problem
    cg->Solve(v, s);
}
