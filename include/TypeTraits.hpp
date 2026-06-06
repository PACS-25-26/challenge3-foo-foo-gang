#pragma once

#include <functional>
#include <Eigen/Dense>

/**
 * @file TypeTraits.hpp
 * @brief Definition of type aliases used in the project in order to ease
 * readbility and code modifications.
 */

namespace Traits
{
    using Real = double;
    using ReturnType = double;
    using Function = std::function <ReturnType(Real, Real)>;                              // Function: R^2 -> R
    using Matrix = Eigen::Matrix <Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>; // dynamic Eigen matrix with Row-Major default ordering
};
