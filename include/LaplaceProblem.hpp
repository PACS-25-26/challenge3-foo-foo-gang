#pragma once

#include "TypeTraits.hpp"

using namespace Traits;

struct LaplaceProblem
{
    Real lb = 0.;
    Real ub = 1.;
    Function force;
    Function dirichlet_bc;
};

class Domain
{
    public:
    unsigned int n;
    Real h;
    const LaplaceProblem &problem;
    Matrix U;

    Domain(unsigned int dim, const LaplaceProblem &prob) : n(dim), problem(prob) 
    {
        h = (problem.ub - problem.lb) / (n - 1);
        U = Matrix::Zero(n, n);
    }

    inline Real getCoord(size_t i) const
    {return problem.lb + i * h;}

    Real compute_L2error(const Function &u_ex) const;

    void exportVTK(const std::string &filename) const;
};

