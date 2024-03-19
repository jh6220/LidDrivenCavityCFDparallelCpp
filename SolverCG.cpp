#include <iostream>
#include <algorithm>
#include <cstring>
using namespace std;

#include <cblas.h>
#include <mpi.h>
#include <omp.h>
#include <cmath>

#include "SolverCG.h"

#define IDX(I,J) ((J)*Nx + (I))

/**
 * @brief Construct a new SolverCG::SolverCG object
 * @param pNx   Number of grid points in x-direction
 * @param pNy   Number of grid points in y-direction
 * @param pdx   Grid spacing in x-direction
 * @param pdy   Grid spacing in y-direction
*/
SolverCG::SolverCG(int pNx, int pNy, double pdx, double pdy)
{
    dx = pdx;
    dy = pdy;
    Nx = pNx;
    Ny = pNy;
    // Initialise solution arrays
    n = Nx*Ny;
    r = new double[n];
    p = new double[n];
    z = new double[n];
    t = new double[n]; //temp
}

/**
 * @brief Destroy the SolverCG::SolverCG object
*/
SolverCG::~SolverCG()
{
    delete[] r;
    delete[] p;
    delete[] z;
    delete[] t;
}

/**
 * @brief Set the parameters regarding MPI parallelism
 * @param pcoords           Coordinates of the process in the 2D grid (int[2])
 * @param pworld_size_root  Number of processes in each direction
 * @param pxCoordComm       Communicator for the x-direction of 2D Cartesian topology
 * @param pyCoordComm       Communicator for the y-direction of 2D Cartesian topology
 * @param pworld_rank       Rank of the process in the global communicator
*/
void SolverCG::SetParallelParams(int* pcoords, int pworld_size_root, MPI_Comm pxCoordComm, MPI_Comm pyCoordComm, int pworld_rank)
{
    coords[0] = pcoords[0];
    coords[1] = pcoords[1];
    world_size_root = pworld_size_root;
    world_size = world_size_root*world_size_root;
    xCoordComm = pxCoordComm;
    yCoordComm = pyCoordComm;
    world_rank = pworld_rank;

    // Allocate memory for the buffers array which are used to update the boundary data with parallel processes
    dataB_left_sent = new double[Ny];
    dataB_right_sent = new double[Ny];
    dataB_top_sent = new double[Nx];
    dataB_bottom_sent = new double[Nx];

    dataB_left_recv = new double[Ny];
    dataB_right_recv = new double[Ny];
    dataB_top_recv = new double[Nx];
    dataB_bottom_recv = new double[Nx];
}

/**
 * @brief Impose boundary conditions on the local stram function array, applied only to the global boundaries of the local domain
 * @param inout     Input array
*/
void SolverCG::ImposeBC(double* inout) {
    int i,j;
    if (coords[0] == 0) {
        #pragma omp parallel for default(shared) private(j) schedule(static)
        for (j = 0; j < Ny; ++j) {
            inout[IDX(0,j)]    = 0.0;
        }
    }

    if (coords[0] == world_size_root-1) {
        #pragma omp parallel for default(shared) private(j) schedule(static)
        for (j = 0; j < Ny; ++j) {
            inout[IDX(Nx-1,j)] = 0.0;        }
    }

    if (coords[1] == 0) {
        #pragma omp parallel for default(shared) private(i) schedule(static)
        for (i = 0; i < Nx; ++i) {
            inout[IDX(i,0)]    = 0.0;
        }
    }

    if (coords[1] == world_size_root-1) {
        #pragma omp parallel for default(shared) private(i) schedule(static)
        for (i = 0; i < Nx; ++i) {
            inout[IDX(i,Ny-1)] = 0.0;
        }
    }

}

/**
 * @brief Update the boundary data of the local stream function array with the data from the neighbouring processes
 * @param data     Input array
 * @param tag      Tag for the MPI communication
*/
void SolverCG::UpdateDataWithParallelProcesses(double* data, int tag) {
    MPI_Request request_left, request_right, request_top, request_bottom;
    int i,j;
    /// Collect boundary data to buffers and sent to neighbours
    if (coords[0] != 0) {
        // #pragma omp parallel for default(shared) private(j) schedule(static)
        for (j = 0; j < Ny; ++j) {
            dataB_left_sent[j] = data[IDX(1, j)];
        }
        MPI_Isend(dataB_left_sent, Ny, MPI_DOUBLE, coords[0]-1, tag, xCoordComm, &request_left);
    }
    if (coords[0] != world_size_root-1) {
        #pragma omp parallel for default(shared) private(j) schedule(static)
        for (j = 0; j < Ny; ++j) {
            dataB_right_sent[j] = data[IDX(Nx-2, j)];
        }
        MPI_Isend(dataB_right_sent, Ny, MPI_DOUBLE, coords[0]+1, tag, xCoordComm, &request_right);
    }
    if (coords[1] != 0) {
        // #pragma omp parallel for default(shared) private(i) schedule(static)
        for (i = 0; i < Nx; ++i) {
            dataB_top_sent[i] = data[IDX(i,1)];
        }
        MPI_Isend(dataB_top_sent, Nx, MPI_DOUBLE, coords[1]-1, tag, yCoordComm, &request_top);
    }
    if (coords[1] != world_size_root-1) {
        // #pragma omp parallel for default(shared) private(i) schedule(static)
        for (i = 0; i < Nx; ++i) {
            dataB_bottom_sent[i] = data[IDX(i, Ny-2)];
        }
        MPI_Isend(dataB_bottom_sent, Nx, MPI_DOUBLE, coords[1]+1, tag, yCoordComm, &request_bottom);
    }

    /// Receive boundary data from neighbours
    if (coords[0] != 0) {
        MPI_Recv(dataB_left_recv, Ny, MPI_DOUBLE, coords[0]-1, tag, xCoordComm, MPI_STATUS_IGNORE);
        // #pragma omp parallel for default(shared) private(j) schedule(static)
        for (j = 0; j < Ny; ++j) {
            data[IDX(0,j)] = dataB_left_recv[j];
        }
    }
    if (coords[0] != world_size_root-1) {
        MPI_Recv(dataB_right_recv, Ny, MPI_DOUBLE, coords[0]+1, tag, xCoordComm, MPI_STATUS_IGNORE);
        // #pragma omp parallel for default(shared) private(j) schedule(static)
        for (j = 0; j < Ny; ++j) {
            data[IDX(Nx-1,j)] = dataB_right_recv[j];
        }
    }
    if (coords[1] != 0) {
        MPI_Recv(dataB_top_recv, Nx, MPI_DOUBLE, coords[1]-1, tag, yCoordComm, MPI_STATUS_IGNORE);
        // #pragma omp parallel for default(shared) private(i) schedule(static)
        for (i = 0; i < Nx; ++i) {
            data[IDX(i,0)] = dataB_top_recv[i];
        }
    }
    if (coords[1] != world_size_root-1) {
        MPI_Recv(dataB_bottom_recv, Nx, MPI_DOUBLE, coords[1]+1, tag, yCoordComm, MPI_STATUS_IGNORE);
        // #pragma omp parallel for default(shared) private(i) schedule(static)
        for (i = 0; i < Nx; ++i) {
            data[IDX(i,Ny-1)] = dataB_bottom_recv[i];
        }
    }
}

/**
 * @brief Calculate the global residual of the data from all processes
 * @param r     Input array
 * @return double 
*/
double SolverCG::CalculateEpsGlobalParallel(double* r) {
    eps = 0;
    int j,Nj,kk;

    /// Calculate the local contribution to the residual
    #pragma omp parallel for private(j,Nj,kk) reduction(+:eps)
    for (j=1; j<Ny-1; j++) {
        Nj = j*Nx;
        for (kk = Nj+1; kk < Nj+Nx-1; ++kk) {
            eps = eps + r[kk]*r[kk];
        }
    }

    /// Calculate the global residual by summing over all processes
    MPI_Allreduce(&eps, &eps_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    eps_global = sqrt(eps_global);
    return eps_global;
}

/**
 * @brief Solve the linear system Ax = b using the Conjugate Gradient method with MPI and OMP parallelism
 * @param b     Right-hand side of the linear system
 * @param x     Solution of the linear system
*/
void SolverCG::Solve(double* b, double* x) {

    int k,i;
    int j,Nj,kk;
    double alpha_num, alpha_den, alpha_num_global, alpha_den_global, alpha_global;
    double beta, beta_global;
    double tol = 0.001;
    double dx2i = 1.0/dx/dx;
    double dy2i = 1.0/dy/dy;
    /// Calculate the inverse factor for the preconditioner
    double factor = 1.0/(2.0*(1.0/(dx*dx) + 1.0/(dy*dy)));

    /// Calculate the global residual and terminate b=0
    eps_global = CalculateEpsGlobalParallel(b);
    if (eps_global < tol) {
        std::fill(x, x+n, 0.0);
        cout << "Norm is " << eps_global << endl;
        return;
    }

    /// Applting the operator to the x_0
    UpdateDataWithParallelProcesses(p, 0);
    ApplyOperator(x, t);
    ImposeBC(r);

    
    #pragma omp parallel for private(i)
    for (i=0;i<Nx*Ny;i++) {
        r[i] = b[i]-t[i]; // r_0 = b - A x_0
        z[i] = r[i]*factor; // preconditioning z_0
        p[i] = z[i];
    }

    k = 0;
    int k_max = 5000;
    /// Main loop interates until the maximum number of iterations is reached or the residual is below the tolerance
    do {
        k++;

        /// Compute the cocal contribution to alpha_k numerator and denominator
        UpdateDataWithParallelProcesses(p, k+1);

        alpha_num = 0;
        alpha_den = 0;
        #pragma omp parallel for private(j,Nj,k) reduction(+:alpha_num, alpha_den)
        for (j = 1; j < Ny - 1; ++j) {
            Nj = j*Nx;
            for (kk = Nj+1; kk < Nj+Nx-1; ++kk) {
                t[kk] = ( -     p[kk-1]
                                + 2.0*p[kk]
                                -     p[kk+1])*dx2i
                            + ( -     p[kk-Nx]
                                + 2.0*p[kk]
                                -     p[kk+Nx])*dy2i;
                alpha_num = alpha_num + r[kk]*z[kk];
                alpha_den = alpha_den + t[kk]*p[kk];
            }
        }

        /// Compute the global alpha_k by summing the denominators and numerators values over all processes
        MPI_Allreduce(&alpha_num, &alpha_num_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&alpha_den, &alpha_den_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        alpha_global = alpha_num_global/alpha_den_global;

        /// Update the x_{k+1}, r_{k+1} and compute the local contributions to the beta_k and residual eps
        eps = 0;
        beta = 0;
        #pragma omp parallel for private(j,Nj) reduction(+:eps, beta)
        for (j = 1; j < Ny - 1; ++j) {
            Nj = j*Nx;
            for (kk = Nj+1; kk < Nj+Nx-1; ++kk) {
                x[kk] = x[kk] + alpha_global*p[kk];
                r[kk] = r[kk] - alpha_global*t[kk];
                z[kk] = r[kk]*factor;
                beta = beta + r[kk]*z[kk];
                eps = eps + r[kk]*r[kk];
            }
        }

        /// Compute the global residual by summing the local contributions over all processes and check for convergence
        MPI_Allreduce(&eps, &eps_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        eps_global = sqrt(eps_global);

        if (eps_global < tol*tol) {
            if (world_rank == 0) {
                cout << "Converged in " << k << " iterations. eps = " << eps_global << endl;
            }
            break;
        }

        /// Compute the global beta_k by summing the local contributions over all processes
        MPI_Allreduce(&beta, &beta_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        beta_global = beta_global/alpha_num_global;

        /// Update the p_{k+1}
        #pragma omp parallel for private(i)
        for (i=0; i<Nx*Ny; i++) {
            p[i]=z[i] + beta_global*p[i];
        }
        

    } while (k < k_max); // Set a maximum number of iterations

    if (k == k_max) {
        if (world_rank == 0) {
        cout << "FAILED TO CONVERGE" << endl;
        }
        exit(-1);
    }
}

/**
 * @brief Apply the linear operator to the input array
 * @param in    Input array
 * @param out   Output array
*/
void SolverCG::ApplyOperator(double* in, double* out) {
    double dx2i = 1.0/dx/dx;
    double dy2i = 1.0/dy/dy;
    int j,Nj,k;
    #pragma omp parallel for private(j,Nj,k)
    for (j = 1; j < Ny - 1; ++j) {
        Nj = j*Nx;
        for (k = Nj+1; k < Nj+Nx-1; ++k) {
            out[k] = ( -     in[k-1]
                              + 2.0*in[k]
                              -     in[k+1])*dx2i
                          + ( -     in[k-Nx]
                              + 2.0*in[k]
                              -     in[k+Nx])*dy2i;
        }
    }
}