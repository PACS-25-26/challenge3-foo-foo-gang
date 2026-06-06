#include "LaplaceProblem.hpp"

#include <fstream>
#include <mpi.h>

// used properly only by rank 0
Real Domain::compute_L2error(const Function &u_ex) const
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank != 0) return 0.;

    Real err_L2 = 0.;

    for(size_t i = 1; i < n-1; ++i){
        for(size_t j = 1; j < n-1; ++j){
            Real x = getCoord(i);
            Real y = getCoord(j);
            Real exact_val = u_ex(x,y);
            err_L2 += (U(i,j) - exact_val) * (U(i,j) - exact_val);
        }
    }
    return sqrt(h * err_L2);
}

// used properly only by rank 0
void Domain::exportVTK(const std::string &filename) const
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank != 0) return;

    std::ofstream file(filename);
    if(!file.is_open()){
        std::cerr << "Impossible to create the solution file" << std::endl;
        return;
    }

    // Header
    file << "# vtk DataFile Version 2.0\n";
    file << "Laplace Problem solution, obtained with a Jacobi Solver on a structured squared mesh\n";
    file << "ASCII\n";

    // Dataset type
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS " << n << " " << n << " 1\n";
    file << "ORIGIN " << problem.lb << " " << problem.lb << " 0.0\n";
    file << "SPACING " << h << " " << h << " 1\n";

    // Solution Points
    file << "POINT_DATA " << n*n << "\n";
    file << "SCALARS u(xh,yh) double 1\n";
    file << "LOOKUP_TABLE default\n";

    for(int j = 0; j < n; ++j)
        for(int i = 0; i < n; ++i)
            file << U(i,j) << "\n";

    file.close();
}
