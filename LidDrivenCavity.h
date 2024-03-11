#pragma once

#include <string>
#include <mpi.h>
using namespace std;

class SolverCG;

class LidDrivenCavity
{
public:
    LidDrivenCavity();
    ~LidDrivenCavity();

    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);

    void Initialise();
    void Integrate();
    void WriteSolution(std::string file);
    void PrintConfiguration();
    double nu   = 0.1;

private:
    double* v   = nullptr;
    double* s   = nullptr;
    double* tmp = nullptr;
    double* dataB_left_sent = nullptr;
    double* dataB_right_sent = nullptr;
    double* dataB_top_sent = nullptr;
    double* dataB_bottom_sent = nullptr;
    double* dataB_left_recv = nullptr;
    double* dataB_right_recv = nullptr;
    double* dataB_top_recv = nullptr;
    double* dataB_bottom_recv = nullptr;

    double dt   = 0.01;
    double T    = 1.0;
    double dx;
    double dy;
    int    Nx   = 9;
    int    Ny   = 9;
    int    Npts = 81;
    int    world_rank;
    int    world_size;
    int    world_size_root, mygrid_rank;
    MPI_Comm mygrid, xCoordsComm, yCoordComm;
    int  kx, ky, Nx_local, Ny_local, Npts_local;
    int coords[2];
    double Lx   = 1.0;
    double Ly   = 1.0;
    double Re   = 10;
    double U    = 1.0;
    // double nu   = 0.1;

    SolverCG* cg = nullptr;

    void CleanUp();
    void UpdateDxDy();
    void Advance();
};

