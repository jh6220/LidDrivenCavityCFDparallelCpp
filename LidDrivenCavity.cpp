#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>
using namespace std;

#include <cblas.h>
#include <mpi.h>
#include <omp.h>

#define IDX(I,J) ((J)*Nx + (I))
#define IDX_local(I,J) ((J)*Nx_local + (I))
#define IDX_write(I,J) ((J)*Nx_write + (I))

#include "LidDrivenCavity.h"
#include "SolverCG.h"

LidDrivenCavity::LidDrivenCavity()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    // Check if the world size is p^2
    world_size_root = sqrt(world_size);
    if (world_size != world_size_root*world_size_root)
    {
        if (world_rank == 0) {
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

    // Define 2D Cartesian topology for the MPI processes
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
    i_start_global = (Nx_local-2)*coords[0] + min(coords[0], kx); //where the index in the x direction starts of the domain of this process (including the boundary)
    if (coords[0] < kx) {
        Nx_local = Nx_local + 1;
    }
    ky = (Ny-2)%world_size_root;
    Ny_local = (Ny-2-ky)/world_size_root + 2;
    j_start_global = (Ny_local-2)*coords[1] + min(coords[1], ky); //where the index in the y direction starts of the domain of this process (including the boundary)
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
    Npts_local = Nx_local * Ny_local;
    v   = new double[Npts_local]();
    vnew= new double[Npts_local]();
    s   = new double[Npts_local]();
    tmp = new double[Npts_local]();
    cg  = new SolverCG(Nx_local, Ny_local, dx, dy);
    cg->SetParallelParams(coords, world_size_root, xCoordComm, yCoordComm, world_rank);

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
        if (world_rank == 0) {
            std::cout << "Step: " << setw(8) << t
                    << "  Time: " << setw(8) << t*dt
                    << std::endl;
        }
        Advance(t);
    }
}

// /**
//  * @brief Writes the solution into a file
//  * @param file  name of the file to which the solution is written
// */
// void LidDrivenCavity::WriteSolution(std::string file)
// {
//     int Nx_write = Nx_local-2;
//     int Ny_write = Ny_local-2;
//     int i_start_write = 1;
//     int i_end_write = Nx_local-2;
//     int j_start_write = 1;
//     int j_end_write = Ny_local-2;

//     if (coords[0] == 0) {
//         Nx_write = Nx_write + 1;
//         i_start_write = 0;
//     }
//     if (coords[0] == world_size_root-1) {
//         Nx_write = Nx_write + 1;
//         i_end_write = Nx_local-1;
//     }
//     if (coords[1] == 0) {
//         Ny_write = Ny_write + 1;
//         j_start_write = 0;
//     }
//     if (coords[1] == world_size_root-1) {
//         Ny_write = Ny_write + 1;
//         j_end_write = Ny_local-1;
//     }

//     double* u0_write = new double[Nx_write*Ny_write]();
//     double* u1_write = new double[Nx_write*Ny_write]();
//     double* s_write = new double[Nx_write*Ny_write]();
//     double* v_write = new double[Nx_write*Ny_write]();
//     double* x_write = new double[Nx_write*Ny_write]();
//     double* y_write = new double[Nx_write*Ny_write]();
//     int i2,j2;
//     // double* x_local = new double[]
//     for (int i = 1; i < Nx_local-1; ++i) {
//         i2 = i-i_start_write;
//         for (int j = 1; j < Nx_local-1; ++j) {
//             j2 = j-j_start_write;
//             u0_write[IDX_write(i2,j2)] =  (s[IDX_local(i,j+1)] - s[IDX_local(i,j)]) / dy;
//             u1_write[IDX_write(i2,j2)] = -(s[IDX_local(i+1,j)] - s[IDX_local(i,j)]) / dx;
//         }
//     }

//     if (coords[1] == world_size_root-1) {   
//         int JN = (Ny_write-1)*Nx_write;
//         for (int i = 0; i < Nx_write; ++i) {
//             u0[JN+i] = U;
//         }
//     }

//     for (int i = 0; i < Nx_write; ++i) {
//         i2 = i+i_start_write;
//         ig = i2+i_start_global;
//         for (int j = 0; j < Ny_write; ++j) {
//             j2 = j+j_start_write;
//             s_write[IDX_write(i,j)] = s[IDX_local(i2,j2)];
//             v_write[IDX_write(i,j)] = v[IDX_local(i2,j2)];
//             x_write[IDX_write(i,j)] = (ig) * dx;
//             y_write[IDX_write(i,j)] = (j2+j_start_write) * dy;
//         }
//     }

//     double* u0_global, *u1_global, *s_global, *v_global, *x_global, *y_global;
//     if (world_rank == 0) {
//         u0_global = new double[Nx*Ny]();
//         u1_global = new double[Nx*Ny]();
//         s_global = new double[Nx*Ny]();
//         v_global = new double[Nx*Ny]();
//         x_global = new double[Nx*Ny]();
//         y_global = new double[Nx*Ny]();
//     }
//     MPI_Gather()

//     std::ofstream f(file.c_str());
//     std::cout << "Writing file " << file << std::endl;
//     int k = 0;
//     for (int i = 0; i < Nx_local; ++i)
//     {
//         for (int j = 0; j < Ny_local; ++j)
//         {
//             k = IDX_local(i, j);
//             f << i * dx << " " << j * dy << " " << v[k] <<  " " << s[k] 
//               << " " << u0[k] << " " << u1[k] << std::endl;
//         }
//         f << std::endl;
//     }
//     f.close();

//     delete[] u0;
//     delete[] u1;
// }

/**
 * @brief Writes the solution into a file
 * @param file  name of the file to which the solution is written
*/
void LidDrivenCavity::WriteSolutionParallel(std::string file)
{
    int Nx_write = Nx_local-2;
    int Ny_write = Ny_local-2;
    int i_start_write = 1;
    int i_end_write = Nx_local-2;
    int j_start_write = 1;
    int j_end_write = Ny_local-2;

    if (coords[0] == 0) {
        Nx_write = Nx_write + 1;
        i_start_write = 0;
    }
    if (coords[0] == world_size_root-1) {
        Nx_write = Nx_write + 1;
        i_end_write = Nx_local-1;
    }
    if (coords[1] == 0) {
        Ny_write = Ny_write + 1;
        j_start_write = 0;
    }
    if (coords[1] == world_size_root-1) {
        Ny_write = Ny_write + 1;
        j_end_write = Ny_local-1;
    }

    double* u0 = new double[Npts_local]();
    double* u1 = new double[Npts_local]();
    for (int i = 1; i < Nx_local - 1; ++i) {
        for (int j = 1; j < Ny_local - 1; ++j) {
            u0[IDX_local(i,j)] =  (s[IDX_local(i,j+1)] - s[IDX_local(i,j)]) / dy;
            u1[IDX_local(i,j)] = -(s[IDX_local(i+1,j)] - s[IDX_local(i,j)]) / dx;
        }
    }

    if (coords[1] == world_size_root-1) {    
        for (int i = 0; i < Nx_local; ++i) {
            u0[IDX_local(i,Ny_local-1)] = U;
        }
    }

    if (world_rank==0) {
        std::cout << "Writing file " << file << std::endl;
    }
    
    stringstream ss;
    int k = 0;
    for (int i = i_start_write ; i < i_end_write+1; ++i)
    {
        for (int j = j_start_write; j < j_end_write+1; ++j)
        {
            k = IDX_local(i, j);
            ss << (i + i_start_global) * dx << " " << (j + j_start_global) * dy 
                << " " << v[k] <<  " " << s[k] 
                << " " << u0[k] << " " << u1[k] << endl;
        }
        // ss << std::endl;
    }

    string data = ss.str();
    int size = data.size();
    int all_sizes[world_size];
    MPI_Allgather(&size, 1, MPI_INT, all_sizes, 1, MPI_INT, MPI_COMM_WORLD);

    MPI_Offset offset = 0;
    for (int i = 0; i < world_rank; ++i) {
        offset += all_sizes[i];
    }

    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, file.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    // Write the data
    MPI_File_write_at_all(fh, offset, data.c_str(), size, MPI_CHAR, MPI_STATUS_IGNORE);

    // Close the file
    MPI_File_close(&fh);

    delete[] u0;
    delete[] u1;
}


void LidDrivenCavity::PrintConfiguration()
{
    if (world_rank == 0) {        
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
}


void LidDrivenCavity::CleanUp()
{
    if (v) {
        delete[] v;
        delete[] vnew;
        delete[] s;
        delete[] tmp;
        delete cg;

        delete[] dataB_bottom_sent;
        delete[] dataB_top_sent;
        delete[] dataB_left_sent;
        delete[] dataB_right_sent;
        delete[] dataB_bottom_recv;
        delete[] dataB_top_recv;
        delete[] dataB_left_recv;
        delete[] dataB_right_recv;
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

    MPI_Request request_left, request_right, request_top, request_bottom;
    // Collect boundary data to buffers and sent to neighbours
    if (coords[0] != 0) {
        // left
        for (int j = 0; j < Ny_local; ++j) {
            dataB_left_sent[j] = data[IDX_local(1, j)];
        }
        MPI_Isend(dataB_left_sent, Ny_local, MPI_DOUBLE, coords[0]-1, tag, xCoordComm, &request_left);
    }
    if (coords[0] != world_size_root-1) {
        // right
        for (int j = 0; j < Ny_local; ++j) {
            dataB_right_sent[j] = data[IDX_local(Nx_local-2, j)];
        }
        MPI_Isend(dataB_right_sent, Ny_local, MPI_DOUBLE, coords[0]+1, tag, xCoordComm, &request_right);
    }
    if (coords[1] != 0) {
        // top
        for (int i = 0; i < Nx_local; ++i) {
            dataB_top_sent[i] = data[IDX_local(i,1)];
        }
        MPI_Isend(dataB_top_sent, Nx_local, MPI_DOUBLE, coords[1]-1, tag, yCoordComm, &request_top);
    }
    if (coords[1] != world_size_root-1) {
        // bottom
        for (int i = 0; i < Nx_local; ++i) {
            dataB_bottom_sent[i] = data[IDX_local(i, Ny_local-2)];
        }
        MPI_Isend(dataB_bottom_sent, Nx_local, MPI_DOUBLE, coords[1]+1, tag, yCoordComm, &request_bottom);
    }

    // Receive boundary data from neighbours
    if (coords[0] != 0) {
        MPI_Recv(dataB_left_recv, Ny_local, MPI_DOUBLE, coords[0]-1, tag, xCoordComm, MPI_STATUS_IGNORE);
        for (int j = 0; j < Ny_local; ++j) {
            data[IDX_local(0,j)] = dataB_left_recv[j];
        }
    }
    if (coords[0] != world_size_root-1) {
        MPI_Recv(dataB_right_recv, Ny_local, MPI_DOUBLE, coords[0]+1, tag, xCoordComm, MPI_STATUS_IGNORE);
        for (int j = 0; j < Ny_local; ++j) {
            data[IDX_local(Nx_local-1,j)] = dataB_right_recv[j];
        }
    }
    if (coords[1] != 0) {
        MPI_Recv(dataB_top_recv, Nx_local, MPI_DOUBLE, coords[1]-1, tag, yCoordComm, MPI_STATUS_IGNORE);
        for (int i = 0; i < Nx_local; ++i) {
            data[IDX_local(i,0)] = dataB_top_recv[i];
        }
    }
    if (coords[1] != world_size_root-1) {
        MPI_Recv(dataB_bottom_recv, Nx_local, MPI_DOUBLE, coords[1]+1, tag, yCoordComm, MPI_STATUS_IGNORE);
        for (int i = 0; i < Nx_local; ++i) {
            data[IDX_local(i,Ny_local-1)] = dataB_bottom_recv[i];
        }
    }
}


void LidDrivenCavity::Advance(int idxT)
{
    double dxi  = 1.0/dx;
    double dyi  = 1.0/dy;
    double dx2i = 1.0/dx/dx;
    double dy2i = 1.0/dy/dy;
    int i, j;

    int n_tags = 3;
    
    // Vorticity boundary conditions
    #pragma omp parallel
    {
        #pragma omp single
        {
            if (coords[0] == 0) {
                #pragma omp task
                for (j = 1; j < Ny_local-1; ++j) {
                    // left
                    v[IDX_local(0,j)]    = 2.0 * dx2i * (s[IDX_local(0,j)]    - s[IDX_local(1,j)]);
                }
            }

            if (coords[0] == world_size_root-1) {
                #pragma omp task
                for (j = 1; j < Ny_local-1; ++j) {
                    // right
                    v[IDX_local(Nx_local-1,j)] = 2.0 * dx2i * (s[IDX_local(Nx_local-1,j)] - s[IDX_local(Nx_local-2,j)]);
                }
            }

            if (coords[1] == 0) {
                #pragma omp task
                for (i = 1; i < Nx_local-1; ++i) {
                    // top
                    v[IDX_local(i,0)]    = 2.0 * dy2i * (s[IDX_local(i,0)]    - s[IDX_local(i,1)]);
                }
            }

            if (coords[1] == world_size_root-1) {
                #pragma omp task
                for (i = 1; i < Nx_local-1; ++i) {
                    // bottom
                    v[IDX_local(i,Ny_local-1)] = 2.0 * dy2i * (s[IDX_local(i,Ny_local-1)] - s[IDX_local(i,Ny_local-2)])
                                - 2.0 * dyi*U;
                }
            }
        }
        #pragma omp taskwait // Ensure all tasks are completed
    }

    // if (coords[0] == 0) {
    //     #pragma omp parallel for default(shared) private(j) schedule(static)
    //     for (j = 1; j < Ny_local-1; ++j) {
    //         // left
    //         v[IDX_local(0,j)]    = 2.0 * dx2i * (s[IDX_local(0,j)]    - s[IDX_local(1,j)]);
    //     }
    // }

    // if (coords[0] == world_size_root-1) {
    //     #pragma omp parallel for default(shared) private(j) schedule(static)
    //     for (j = 1; j < Ny_local-1; ++j) {
    //         // right
    //         v[IDX_local(Nx_local-1,j)] = 2.0 * dx2i * (s[IDX_local(Nx_local-1,j)] - s[IDX_local(Nx_local-2,j)]);
    //     }
    // }

    // if (coords[1] == 0) {
    //     #pragma omp parallel for default(shared) private(i) schedule(static)
    //     for (i = 1; i < Nx_local-1; ++i) {
    //         // top
    //         v[IDX_local(i,0)]    = 2.0 * dy2i * (s[IDX_local(i,0)]    - s[IDX_local(i,1)]);
    //     }
    // }

    // if (coords[1] == world_size_root-1) {
    //     #pragma omp parallel for default(shared) private(i) schedule(static)
    //     for (i = 1; i < Nx_local-1; ++i) {
    //         // bottom
    //         v[IDX_local(i,Ny_local-1)] = 2.0 * dy2i * (s[IDX_local(i,Ny_local-1)] - s[IDX_local(i,Ny_local-2)])
    //                        - 2.0 * dyi*U;
    //     }
    // }

    //Exchange vorticity data with parallel processes
    UpdateDataWithParallelProcesses(v, n_tags*idxT+0);

    // int threadid;
    // #pragma omp parallel private(threadid)
    // {
    //     threadid = omp_get_thread_num();
    //     #pragma omp critical
    //     cout << "Thread: " << threadid <<endl;
    // }

    // Compute interior vorticity
    #pragma omp parallel for collapse(2) default(shared) private(i,j) schedule (static)
    for (j = 1; j < Ny_local - 1; ++j) {
        for (i = 1; i < Nx_local - 1; ++i) {
            v[IDX_local(i,j)] = dx2i*(
                    2.0 * s[IDX_local(i,j)] - s[IDX_local(i+1,j)] - s[IDX_local(i-1,j)])
                        + 1.0/dy/dy*(
                    2.0 * s[IDX_local(i,j)] - s[IDX_local(i,j+1)] - s[IDX_local(i,j-1)]);
        }
    }

    //Exchange vorticity data with parallel processes
    UpdateDataWithParallelProcesses(v, n_tags*idxT+1);

    // Time advance vorticity
    #pragma omp parallel for collapse(2) default(shared) private(i,j) schedule (static)
    for (j = 1; j < Ny_local - 1; ++j) {
        for (i = 1; i < Nx_local - 1; ++i) {
            vnew[IDX_local(i,j)] = v[IDX_local(i,j)] + dt*(
                ( (s[IDX_local(i+1,j)] - s[IDX_local(i-1,j)]) * 0.5 * dxi
                 *(v[IDX_local(i,j+1)] - v[IDX_local(i,j-1)]) * 0.5 * dyi)
              - ( (s[IDX_local(i,j+1)] - s[IDX_local(i,j-1)]) * 0.5 * dyi
                 *(v[IDX_local(i+1,j)] - v[IDX_local(i-1,j)]) * 0.5 * dxi)
              + nu * (v[IDX_local(i+1,j)] - 2.0 * v[IDX_local(i,j)] + v[IDX_local(i-1,j)])*dx2i
              + nu * (v[IDX_local(i,j+1)] - 2.0 * v[IDX_local(i,j)] + v[IDX_local(i,j-1)])*dy2i);
        }
    }

    // Solve Poisson problem
    cg->SolveParallel(vnew, s);

    UpdateDataWithParallelProcesses(s, n_tags*idxT+2);
}
