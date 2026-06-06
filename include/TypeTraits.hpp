#pragma once

#include <functional>
#include <Eigen/Dense>

namespace Traits
{
    using Real = double;
    using ReturnType = double;
    using Function = std::function <ReturnType(Real, Real)>;
    using Matrix = Eigen::Matrix <Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
};
