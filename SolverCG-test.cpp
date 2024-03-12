#include <iostream>
#include "SolverCG.h"

using namespace std;

#define BOOST_TEST_MODULE SolverCG
#include <boost/test/included/unit_test.hpp>

#define IDX(I,J) ((J)*Nx + (I))

void PrintMatrix(int Ny, int Nx, double* M) {
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            cout << M[j*Nx+i] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

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

BOOST_AUTO_TEST_CASE(Solve2) {

    int Nx = 101;
    int Ny = 101;
    double dx = 1.0/(Nx-1);
    double dy = 1.0/(Nx-1);
    double* v = new double[Nx*Ny]();
    double* s = new double[Nx*Ny]();
    double* s_true = new double[Nx*Ny]();
    SolverCG* cg  = new SolverCG(Nx, Ny, dx, dy);

    double tol = 0.001;

    const int k = 3;
    const int l = 3;
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            v[IDX(i,j)] = -M_PI * M_PI * (k * k + l * l)
                                       * sin(M_PI * k * i * dx)
                                       * sin(M_PI * l * j * dy);
            s_true[IDX(i,j)] = -sin(M_PI * k * i * dx) * sin(M_PI * l * j * dy);
        }
    }

    cg->Solve(v, s);

    for (int i = 0; i < Nx*Ny; i++) 
    {
        BOOST_CHECK_SMALL(s[i]-s_true[i], tol);
    }   
}


BOOST_AUTO_TEST_CASE(SolveParallel) {
    MPI_Init(NULL, NULL);

    int world_size, world_rank, world_size_root;
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
    int Nx = 101;
    int Ny = 101;
    double dx = 1.0/(Nx-1);
    double dy = 1.0/(Nx-1);

    // Define 2D Cartesian topology for the MPI processes
    int dims[2] = {world_size_root, world_size_root};
    int periods[2] = {0, 0};
    int coords[2];
    int reorder = 1;
    int mygrid_rank, Nx_local, Ny_local, kx, ky, i_start_global, j_start_global;
    MPI_Comm mygrid, xCoordComm, yCoordComm;

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

    // cout << "Rank: " << world_rank << " xCoord: " << coords[0] << " yCoord: " << coords[1] << " Nx_local: " << Nx_local << " Ny_local: " << Ny_local << " i_start_global: " << i_start_global << " j_start_global " << j_start_global << endl;

    // MPI_Barrier(MPI_COMM_WORLD);

    int Npts_local = Nx_local*Ny_local;
    double* v = new double[Npts_local]();
    double* s = new double[Npts_local]();
    double* s_true = new double[Npts_local]();

    double tol = 0.001;

    const int k = 1;
    const int l = 1;
    for (int i = 0; i < Nx_local; ++i) {
        for (int j = 0; j < Ny_local; ++j) {
            v[j*Nx_local+i] = -M_PI * M_PI * (k * k + l * l)
                                       * sin(M_PI * k * (i+i_start_global) * dx)
                                       * sin(M_PI * l * (j+j_start_global) * dy);
            s_true[j*Nx_local+i] = -sin(M_PI * k * (i+i_start_global) * dx) * sin(M_PI * l * (j+j_start_global) * dy);
        }
    }

    // if (world_rank != 0) {
    //     // If not the first process, receive a token from the previous process
    //     int token;
    //     MPI_Recv(&token, 1, MPI_INT, world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // }

    // // Perform the print operation
    // cout << "Process " << world_rank << " is printing its matrix." << endl;
    // cout << "Nx_local: " << Nx_local << " Ny_local: " << Ny_local << " i_start_global: " << i_start_global << " j_start_global " << j_start_global << endl;
    // PrintMatrix(Ny_local, Nx_local, s_true);

    // if (world_rank != world_size - 1) {
    //     // If not the last process, send a token to the next process
    //     int token = 0; // The value of the token doesn't matter in this case
    //     MPI_Send(&token, 1, MPI_INT, world_rank + 1, 0, MPI_COMM_WORLD);
    // }


    SolverCG* cg  = new SolverCG(Nx_local, Ny_local, dx, dy);
    cg->SetParallelParams(coords, world_size_root, xCoordComm, yCoordComm, world_rank);
    cg->SolveParallel(v, s);

    for (int i = 1; i < Nx_local-1; i++) 
    {
        for (int j = 1; j < Ny_local-1; j++) 
        {
            BOOST_CHECK_SMALL(s[j*Nx_local+i]-s_true[j*Nx_local+i], tol);
        }
    }

    MPI_Finalize();
    
}