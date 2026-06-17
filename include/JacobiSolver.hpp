#pragma once

#include "LaplaceProblem.hpp"

/**
 * @class JacobiSolver
 * @brief Implement a parallel solver for the Laplace Equation, specifically using
 * the Jacobi iterative method, exploiting hybrid parallel implementation - MPI + OpenMP.
 */
class JacobiSolver
{
    private:
    size_t max_iter;    // maximum number of iterations
    Real tolerance;     // tolerance for the convergence criterion

    public:
    /**
     * @brief Solver constructor
     * @param max_it maximum iterations number
     * @param tol tolerance
     */
    JacobiSolver(size_t max_it, Real tol) : max_iter(max_it), tolerance(tol) {}
    
    /**
     * @brief Core function of the class: solve the problem on the assigned domain. Execute
     * Jacobi Algorithm sharing rows workloads over different MPI processes, with a further
     * parallelization through an OpenMP directive over 'for' loops.
     * 
     * @param d reference computational domain, which includes the grid and the forcing term 
     * @return number of iterations used to reach convergence
     */
    size_t solve(Domain &d) const;
};
