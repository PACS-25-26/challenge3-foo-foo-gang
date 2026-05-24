#include "LaplaceProblem.hpp"
#include <fstream>
#include <mpi.h>

void Domain::exportVTK(const std::string &filename) const
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank != 0) return;

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
    file << "POINT DATA " << n*n << "\n";
    file << "SCALARS u(xh,yh) double 1\n";
    file << "LOOKUP_TABLE default\n";

    for(int j = 0; j < n; ++j)
        for(int i = 0; i < n; ++i)
            file << U(i,j) << "\n";

    file.close();
}
