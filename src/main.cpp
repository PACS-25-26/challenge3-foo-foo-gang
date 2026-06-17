#include "JacobiSolver.hpp"
#include "muParserInterface.hpp"

#include <mpi.h>

/**
 * @file main.cpp
 * @brief Initialize MPI communication and solve the desired Laplace Problem.
 * Beware that the number of iterations is set to 100,000 in order to allow that Jacobi Algorithm
 * converges for the grid dimension up to n = 256, while the tolerance is set to 1e-8. Reducing it
 * allows to reach more accurate results, although being more computationally demanding.
 * 
 * @param argc number of arguments passed by command line (must be 3)
 * @param argv array containing arguments, in order: 
 *             - executable name
 *             - domain-grid dimensions
 *             - forcing term f(x,y)
 */
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(argc != 3){
        if(rank == 0)
            std::cerr << "You must pass through command line in order: "
                      << "Dimension of the grid n and Forcing term f!" << std::endl;        
        MPI_Finalize();
        return -1;
    }
    
    unsigned int n = std::stoi(argv[1]); // string to integer
    
    // muParser utility usage
    std::string fun_expr = argv[2];
    muParserWrapper fun_parser(fun_expr);
    Function f = [&](Real x, Real y) -> Real {return fun_parser(x,y);};

    // definition of the desired problem
    LaplaceProblem problem;
    problem.force = f;
    problem.dirichlet_bc = [](Real x, Real y) -> Real {return 0.;}; // change here the BCs

    // construction of the computational grid
    Domain domain(n, problem);
    // construction of the Jacobi solver
    JacobiSolver solver(100000, 1e-8);

    // Start of the clocking algorithm using MPI clocking functions:
    // - MPI_Barrier() synchronyzes all the processes before starting chronometer
    // - MPI_Wtime() computes the wall-clock time used for convergence
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    // solve the problem
    size_t iter_count = solver.solve(domain);
    
    // exact solution definition: u_ex(x,y) = sin(2*pi*x)*sin(2*pi*y)
    Function u_ex = [](Real x, Real y) -> Real {
        return sin(2*M_PI*x) * sin(2*M_PI*y);
    };
    // L2 error computation
    Real err_L2 = domain.compute_L2error(u_ex);
    
    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();
    double local_elapsed = end_time - start_time;
    double max_elapsed = 0.0;
    // End of the clocking algorithm

    // take only the maximum time elapsed among all the ranks
    MPI_Reduce(&local_elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(rank == 0) {
        std::cout << "N = " << n << " | Procs = " << size 
                  << " | Time = " << max_elapsed << "s | L2 Error = " << err_L2 
                  << " | Iters = " << iter_count << std::endl;
    }

    // export the solution in the corresponfing VTK file
    domain.exportVTK("Laplace_Solution.vtk");

    MPI_Finalize();
}