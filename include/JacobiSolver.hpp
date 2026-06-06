#pragma once

#include "LaplaceProblem.hpp"


class JacobiSolver
{
    private:
    size_t max_iter;
    Real tolerance;

    public:
    JacobiSolver(size_t max_it, Real tol) : max_iter(max_it), tolerance(tol) {}
    size_t solve(Domain &d) const;
};
