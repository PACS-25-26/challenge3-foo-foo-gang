#include "JacobiSolver.hpp"
#include "muParserInterface.hpp"

#include <mpi.h>

int main (int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 3){
        if(rank == 0)
            std::cerr << "You must pass through command line in order: "
                      << "Dimension of the grid n and Forcing term f!" << std::endl;        
        MPI_Finalize();
        return -1;
    }
    
    unsigned int n = std::stoi(argv[1]);
    
    std::string fun_expr = argv[2];
    muParserWrapper fun_parser(fun_expr);
    Function f = [&](Real x, Real y) -> Real {return fun_parser(x,y);};

    LaplaceProblem problem;
    problem.force = f;
    Domain domain(n, problem);
    JacobiSolver solver(1000, 1e-5);

    solver.solve(domain);

    domain.exportVTK("Laplace_Solution.vtk");

    MPI_Finalize();
}