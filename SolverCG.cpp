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

SolverCG::SolverCG(int pNx, int pNy, double pdx, double pdy)
{
    dx = pdx;
    dy = pdy;
    Nx = pNx;
    Ny = pNy;
    n = Nx*Ny;
    r = new double[n];
    p = new double[n];
    z = new double[n];
    t = new double[n]; //temp
}


SolverCG::~SolverCG()
{
    delete[] r;
    delete[] p;
    delete[] z;
    delete[] t;
}

void SolverCG::SetParallelParams(int* pcoords, int pworld_size_root, MPI_Comm pxCoordComm, MPI_Comm pyCoordComm, int pworld_rank)
{
    coords[0] = pcoords[0];
    coords[1] = pcoords[1];
    world_size_root = pworld_size_root;
    world_size = world_size_root*world_size_root;
    xCoordComm = pxCoordComm;
    yCoordComm = pyCoordComm;
    world_rank = pworld_rank;

    dataB_left_sent = new double[Ny];
    dataB_right_sent = new double[Ny];
    dataB_top_sent = new double[Nx];
    dataB_bottom_sent = new double[Nx];

    dataB_left_recv = new double[Ny];
    dataB_right_recv = new double[Ny];
    dataB_top_recv = new double[Nx];
    dataB_bottom_recv = new double[Nx];
}

void SolverCG::ImposeBCParallel(double* inout) {
    int i,j;
    if (coords[0] == 0) {
        #pragma omp parallel for default(shared) private(j) schedule(static)
        for (j = 0; j < Ny; ++j) {
            // left
            inout[IDX(0,j)]    = 0.0;
        }
    }

    if (coords[0] == world_size_root-1) {
        #pragma omp parallel for default(shared) private(j) schedule(static)
        for (j = 0; j < Ny; ++j) {
            // right
            inout[IDX(Nx-1,j)] = 0.0;        }
    }

    if (coords[1] == 0) {
        #pragma omp parallel for default(shared) private(i) schedule(static)
        for (i = 0; i < Nx; ++i) {
            // top
            inout[IDX(i,0)]    = 0.0;
        }
    }

    if (coords[1] == world_size_root-1) {
        #pragma omp parallel for default(shared) private(i) schedule(static)
        for (i = 0; i < Nx; ++i) {
            // bottom
            inout[IDX(i,Ny-1)] = 0.0;
        }
    }

}

// void SolverCG::ImposeBCParallel(double* inout) {

//     #pragma omp parallel
//     {
//         #pragma omp single
//         {    
//             if (coords[0] == 0) {
//                 #pragma omp task
//                 for (int j = 0; j < Ny; ++j) {
//                     // left
//                     inout[IDX(0,j)]    = 0.0;
//                 }
//             }

//             if (coords[0] == world_size_root-1) {
//                 #pragma omp task
//                 for (int j = 0; j < Ny; ++j) {
//                     // right
//                     inout[IDX(Nx-1,j)] = 0.0;        }
//             }

//             if (coords[1] == 0) {
//                 #pragma omp task
//                 for (int i = 0; i < Nx; ++i) {
//                     // top
//                     inout[IDX(i,0)]    = 0.0;
//                 }
//             }

//             if (coords[1] == world_size_root-1) {
//                 #pragma omp task
//                 for (int i = 0; i < Nx; ++i) {
//                     // bottom
//                     inout[IDX(i,Ny-1)] = 0.0;
//                 }
//             }
//         }
//         #pragma omp taskwait // Ensure all tasks are completed
//     }

// }

void SolverCG::UpdateDataWithParallelProcesses(double* data, int tag) {
    MPI_Request request_left, request_right, request_top, request_bottom;
    // Collect boundary data to buffers and sent to neighbours
    if (coords[0] != 0) {
        // left
        for (int j = 0; j < Ny; ++j) {
            dataB_left_sent[j] = data[IDX(1, j)];
        }
        MPI_Isend(dataB_left_sent, Ny, MPI_DOUBLE, coords[0]-1, tag, xCoordComm, &request_left);
    }
    if (coords[0] != world_size_root-1) {
        // right
        for (int j = 0; j < Ny; ++j) {
            dataB_right_sent[j] = data[IDX(Nx-2, j)];
        }
        MPI_Isend(dataB_right_sent, Ny, MPI_DOUBLE, coords[0]+1, tag, xCoordComm, &request_right);
    }
    if (coords[1] != 0) {
        // top
        for (int i = 0; i < Nx; ++i) {
            dataB_top_sent[i] = data[IDX(i,1)];
        }
        MPI_Isend(dataB_top_sent, Nx, MPI_DOUBLE, coords[1]-1, tag, yCoordComm, &request_top);
    }
    if (coords[1] != world_size_root-1) {
        // bottom
        for (int i = 0; i < Nx; ++i) {
            dataB_bottom_sent[i] = data[IDX(i, Ny-2)];
        }
        MPI_Isend(dataB_bottom_sent, Nx, MPI_DOUBLE, coords[1]+1, tag, yCoordComm, &request_bottom);
    }

    // Receive boundary data from neighbours
    if (coords[0] != 0) {
        MPI_Recv(dataB_left_recv, Ny, MPI_DOUBLE, coords[0]-1, tag, xCoordComm, MPI_STATUS_IGNORE);
        for (int j = 0; j < Ny; ++j) {
            data[IDX(0,j)] = dataB_left_recv[j];
        }
    }
    if (coords[0] != world_size_root-1) {
        MPI_Recv(dataB_right_recv, Ny, MPI_DOUBLE, coords[0]+1, tag, xCoordComm, MPI_STATUS_IGNORE);
        for (int j = 0; j < Ny; ++j) {
            data[IDX(Nx-1,j)] = dataB_right_recv[j];
        }
    }
    if (coords[1] != 0) {
        MPI_Recv(dataB_top_recv, Nx, MPI_DOUBLE, coords[1]-1, tag, yCoordComm, MPI_STATUS_IGNORE);
        for (int i = 0; i < Nx; ++i) {
            data[IDX(i,0)] = dataB_top_recv[i];
        }
    }
    if (coords[1] != world_size_root-1) {
        MPI_Recv(dataB_bottom_recv, Nx, MPI_DOUBLE, coords[1]+1, tag, yCoordComm, MPI_STATUS_IGNORE);
        for (int i = 0; i < Nx; ++i) {
            data[IDX(i,Ny-1)] = dataB_bottom_recv[i];
        }
    }
}

void SolverCG::PrintMatrix(int Ny, int Nx, double* M) {
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            cout << M[j*Nx+i] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

double SolverCG::CalculateEpsGlobalParallel(double* r) {
    eps = 0;
    int j = 0;

    #pragma omp parallel for default(shared) private(j) reduction(+:eps) schedule(static)
    for (j=1; j<Ny-1; j++) {
        eps += cblas_ddot(Nx-2, r+j*Nx+1, 1, r+j*Nx+1, 1);
    }
    MPI_Allreduce(&eps, &eps_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    eps_global = sqrt(eps_global);
    return eps_global;
}

void SolverCG::PreconditionParallel(double* in, double* out, double factor) {
    cblas_dcopy(n, in, 1, out, 1);
    cblas_dscal(n, factor, out, 1);
}


void SolverCG::SolveParallel(double* b, double* x) {
    int k,i;
    double alpha_num, alpha_den, alpha_num_global, alpha_den_global, alpha_global;
    double beta, beta_global;
    double tol = 0.001;
    double factor = 1.0/(2.0*(1.0/(dx*dx) + 1.0/(dy*dy)));

    eps_global = CalculateEpsGlobalParallel(b);
    if (eps_global < tol) {
        std::fill(x, x+n, 0.0);
        cout << "Norm is " << eps_global << endl;
        return;
    }



    UpdateDataWithParallelProcesses(p, 0);

    ApplyOperator(x, t);
    cblas_dcopy(n, b, 1, r, 1);        // r_0 = b (i.e. b)
    ImposeBCParallel(r);

    cblas_daxpy(n, -1.0, t, 1, r, 1);
    PreconditionParallel(r, z, factor);
    cblas_dcopy(n, z, 1, p, 1);        // p_0 = r_0

    k = 0;
    int k_max = 5000;
    do {
        k++;
        UpdateDataWithParallelProcesses(p, k+1);
        ApplyOperator(p, t);

        alpha_num = 0;
        alpha_den = 0;
        #pragma omp parallel for default(shared) private(i) reduction(+:alpha_num, alpha_den) schedule(static)
        for (int i=1; i<Ny-1; i++) {
            alpha_num += cblas_ddot(Nx-2, p+i*Nx+1, 1, t+i*Nx+1, 1);
            alpha_den += cblas_ddot(Nx-2, t+i*Nx+1, 1, p+i*Nx+1, 1);
        }

        MPI_Allreduce(&alpha_num, &alpha_num_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&alpha_den, &alpha_den_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        alpha_global = alpha_num_global/alpha_den_global;

        cblas_daxpy(n,  alpha_global, p, 1, x, 1);  // x_{k+1} = x_k + alpha_k p_k
        cblas_daxpy(n, -alpha_global, t, 1, r, 1); // r_{k+1} = r_k - alpha_k A p_k

        eps_global = CalculateEpsGlobalParallel(r);
        // if (world_rank==0) {
        //     cout << "eps_global = " << eps_global << endl;}

        if (eps_global < tol*tol) {
            if (world_rank == 0) {
                cout << "Converged in " << k << " iterations. eps = " << eps_global << endl;
            }
            break;
        }

        PreconditionParallel(r, z, factor);
        beta = 0;
        #pragma omp parallel for default(shared) private(i) reduction(+:beta) schedule(static)
        for (i=1; i<Ny-1; i++) {
            beta += cblas_ddot(Nx-2, r+i*Nx+1, 1, z+i*Nx+1, 1);
        }

        MPI_Allreduce(&beta, &beta_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        beta_global = beta_global/alpha_den_global;

        cblas_dcopy(n, z, 1, t, 1);
        cblas_daxpy(n, beta_global, p, 1, t, 1);
        cblas_dcopy(n, t, 1, p, 1);

    } while (k < k_max); // Set a maximum number of iterations

    if (k == k_max) {
        if (world_rank == 0) {
        cout << "FAILED TO CONVERGE" << endl;
        }
        exit(-1);
    }
}

void SolverCG::Solve(double* b, double* x) {
    unsigned int n = Nx*Ny;
    int k;
    double alpha;
    double beta;
    double eps;
    double tol = 0.001;

    eps = cblas_dnrm2(n, b, 1);
    if (eps < tol*tol) {
        std::fill(x, x+n, 0.0);
        cout << "Norm is " << eps << endl;
        return;
    }

    ApplyOperator(x, t);
    cblas_dcopy(n, b, 1, r, 1);        // r_0 = b (i.e. b)
    ImposeBC(r);

    cblas_daxpy(n, -1.0, t, 1, r, 1);
    Precondition(r, z);
    cblas_dcopy(n, z, 1, p, 1);        // p_0 = r_0

    k = 0;
    do {
        k++;
        // Perform action of Nabla^2 * p
        ApplyOperator(p, t);

        alpha = cblas_ddot(n, t, 1, p, 1);  // alpha = p_k^T A p_k
        alpha = cblas_ddot(n, r, 1, z, 1) / alpha; // compute alpha_k
        beta  = cblas_ddot(n, r, 1, z, 1);  // z_k^T r_k

        cblas_daxpy(n,  alpha, p, 1, x, 1);  // x_{k+1} = x_k + alpha_k p_k
        cblas_daxpy(n, -alpha, t, 1, r, 1); // r_{k+1} = r_k - alpha_k A p_k

        eps = cblas_dnrm2(n, r, 1);

        if (eps < tol*tol) {
            break;
        }
        Precondition(r, z);
        beta = cblas_ddot(n, r, 1, z, 1) / beta;

        cblas_dcopy(n, z, 1, t, 1);
        cblas_daxpy(n, beta, p, 1, t, 1);
        cblas_dcopy(n, t, 1, p, 1);

    } while (k < 5000); // Set a maximum number of iterations

    if (k == 5000) {
        cout << "FAILED TO CONVERGE" << endl;
        exit(-1);
    }

    cout << "Converged in " << k << " iterations. eps = " << eps << endl;
}

// void SolverCG::ApplyOperator(double* in, double* out) {
//     // Assume ordered with y-direction fastest (column-by-column)
//     double dx2i = 1.0/dx/dx;
//     double dy2i = 1.0/dy/dy;
//     int jm1 = 0, jp1 = 2;

//     for (int j = 1; j < Ny - 1; ++j) {
//         for (int i = 1; i < Nx - 1; ++i) {
//             out[IDX(i,j)] = ( -     in[IDX(i-1, j)]
//                               + 2.0*in[IDX(i,   j)]
//                               -     in[IDX(i+1, j)])*dx2i
//                           + ( -     in[IDX(i, jm1)]
//                               + 2.0*in[IDX(i,   j)]
//                               -     in[IDX(i, jp1)])*dy2i;
//         }
//         jm1++;
//         jp1++;
//     }

void SolverCG::ApplyOperator(double* in, double* out) {
    // Assume ordered with y-direction fastest (column-by-column)
    double dx2i = 1.0/dx/dx;
    double dy2i = 1.0/dy/dy;
    int jm1 = 0, jp1 = 2;
    int i,j;
    #pragma omp parallel for default(shared) private(i,j,jm1,jp1) schedule(static)
    for (j = 1; j < Ny - 1; ++j) {
        jm1 = j-1;
        jp1 = j+1;
        for (i = 1; i < Nx - 1; ++i) {
            out[IDX(i,j)] = ( -     in[IDX(i-1, j)]
                              + 2.0*in[IDX(i,   j)]
                              -     in[IDX(i+1, j)])*dx2i
                          + ( -     in[IDX(i, jm1)]
                              + 2.0*in[IDX(i,   j)]
                              -     in[IDX(i, jp1)])*dy2i;
        }
    }
}


void SolverCG::Precondition(double* in, double* out) {
    int i, j;
    double dx2i = 1.0/dx/dx;
    double dy2i = 1.0/dy/dy;
    double factor = 2.0*(dx2i + dy2i);
    for (i = 1; i < Nx - 1; ++i) {
        for (j = 1; j < Ny - 1; ++j) {
            out[IDX(i,j)] = in[IDX(i,j)]/factor; // OPTIMIZE by multiplying by 1/factor
        }
    }
    // Boundaries
    for (i = 0; i < Nx; ++i) {
        out[IDX(i, 0)] = in[IDX(i,0)];
        out[IDX(i, Ny-1)] = in[IDX(i, Ny-1)];
    }

    for (j = 0; j < Ny; ++j) {
        out[IDX(0, j)] = in[IDX(0, j)];
        out[IDX(Nx - 1, j)] = in[IDX(Nx - 1, j)];
    }
}

void SolverCG::ImposeBC(double* inout) {
        // Boundaries
    for (int i = 0; i < Nx; ++i) {
        inout[IDX(i, 0)] = 0.0;
        inout[IDX(i, Ny-1)] = 0.0;
    }

    for (int j = 0; j < Ny; ++j) {
        inout[IDX(0, j)] = 0.0;
        inout[IDX(Nx - 1, j)] = 0.0;
    }

}
