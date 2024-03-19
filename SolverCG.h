#pragma once

#include <mpi.h>

class SolverCG
{
public:
    SolverCG(int pNx, int pNy, double pdx, double pdy);
    ~SolverCG();

    void Solve(double* b, double* x);
    void SetParallelParams(int* pcoords, int pworld_size_root, MPI_Comm pxCoordComm, MPI_Comm pyCoordComm, int pworld_rank);

private:
    double dx;
    double dy;
    double eps, eps_global;
    int Nx;
    int Ny;
    unsigned int n;
    double* r;
    double* p;
    double* z;
    double* t;
    double* dataB_left_sent, *dataB_right_sent, *dataB_top_sent, *dataB_bottom_sent;
    double* dataB_left_recv, *dataB_right_recv, *dataB_top_recv, *dataB_bottom_recv;
    int coords[2];
    int world_size_root, world_rank, world_size;
    MPI_Comm xCoordComm, yCoordComm;

    void ImposeBC(double* p);
    void ApplyOperator(double* p, double* t);
    void UpdateDataWithParallelProcesses(double* data, int tag);
    double CalculateEpsGlobalParallel(double* r);
};

