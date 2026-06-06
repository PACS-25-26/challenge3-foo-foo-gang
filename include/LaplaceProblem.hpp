#pragma once

#include "TypeTraits.hpp"

using namespace Traits;

/**
 * @struct LaplaceProblem
 * @brief Contains mathematical parameters of Laplace Problem for a squared domain.
 */
struct LaplaceProblem
{
    Real lb = 0.;           // lower bound for the grid
    Real ub = 1.;           // upper bound for the grid
    Function force;         // forcing term, f(x,y) in the equation
    Function dirichlet_bc;  // Dirichlet Boundary Conditions
};

/**
 * @class Domain
 * @brief Represents the computational grid and save the global solution U as an Eigen/Dense matrix.
 */
class Domain
{
    public:
    unsigned int n;                // number of nodes alongside each direction (grid dimensions)
    Real h;                        // grid spacing between two nodes
    const LaplaceProblem &problem; // Laplace problem's reference
    Matrix U;                      // global solution U of the numerical problem

    /**
     * @brief Domain constructor
     * @param dim dimension of the computational grid
     * @param prob reference to the problem to solve
     */
    Domain(unsigned int dim, const LaplaceProblem &prob) : n(dim), problem(prob) 
    {
        h = (problem.ub - problem.lb) / (n - 1);
        U = Matrix::Zero(n, n);
    }

    /**
     * @brief Mapping from indexes to physical coordinates of the domain
    */
    inline Real getCoord(size_t i) const
    {return problem.lb + i * h;}

    /**
     * @brief Compute the L2-norm Error wrt the exact solution u_ex 
     */
    Real compute_L2error(const Function &u_ex) const;

    /**
     * @brief Print the solution in VTK format for the visualization
     */
    void exportVTK(const std::string &filename) const;
};

