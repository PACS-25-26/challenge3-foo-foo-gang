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
    problem.dirichlet_bc = [](Real x, Real y) -> Real {return 0.;};

    Domain domain(n, problem);
    JacobiSolver solver(100000, 1e-8);

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    size_t iter_count = solver.solve(domain);
    
    Function u_ex = [](Real x, Real y) -> Real {
        return sin(2*M_PI*x) * sin(2*M_PI*y);
    };
    Real err_L2 = domain.compute_L2error(u_ex);
    
    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();
    double local_elapsed = end_time - start_time;
    double max_elapsed = 0.0;

    MPI_Reduce(&local_elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(rank == 0) {
        std::cout << "N = " << n << " | Procs = " << size 
                  << " | Time = " << max_elapsed << "s | L2 Error = " << err_L2 
                  << " | Iters = " << iter_count << std::endl;
    }

    domain.exportVTK("Laplace_Solution.vtk");

    MPI_Finalize();
}