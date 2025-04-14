#pragma once

#include <Eigen/Dense>

namespace wmtk::components::internal::utils {

/**
 * @brief Unsigned mean ratio metric.
 */
inline double mean_ratio_metric(
    const Eigen::Ref<const Eigen::Vector3d> p0,
    const Eigen::Ref<const Eigen::Vector3d> p1,
    const Eigen::Ref<const Eigen::Vector3d> p2)
{
    const Eigen::Vector3d a = (p1 - p0);
    const Eigen::Vector3d b = (p2 - p1);
    const Eigen::Vector3d c = (p0 - p2);

    const double sq_length_sum = a.squaredNorm() + b.squaredNorm() + c.squaredNorm();
    const double area = a.cross(b).norm();

    const double prefactor = 2 * std::sqrt(3);

    return prefactor * area / sq_length_sum;
}

} // namespace wmtk::components::internal::utils