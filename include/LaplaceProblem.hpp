#pragma once

#include "TypeTraits.hpp"

using namespace Traits;

struct LaplaceProblem
{
    Real lb;
    Real ub;
    Function force;
    // Function bound_conds;
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
};

