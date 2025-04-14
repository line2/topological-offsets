#pragma once

#include <Eigen/Dense>

namespace wmtk::components::internal::utils {

/**
 * @brief Computes the angle between two lines represented by normalized vectors in degree.
 *
 * Max angle is 90 degrees.
 */
inline double normal_angle(const Eigen::Vector3d& a, const Eigen::Vector3d& b)
{
    const double ab_norm = std::clamp(a.cross(b).norm(), -1., 1.);
    const double angle_rad = std::asin(ab_norm);
    return (180. / M_PI) * angle_rad;
}
inline double normal_angle(const Eigen::Vector2d& a, const Eigen::Vector2d& b)
{
    const double dot = std::abs(std::clamp(a.dot(b), -1., 1.));
    const double angle_rad = std::acos(dot);
    return (180. / M_PI) * angle_rad;
}

/**
 * @brief Computes the angle between two normalized vectors in degree.
 *
 * Max angle is 180 degrees.
 */
inline double normal_angle_180(
    const Eigen::Ref<const Eigen::Vector3d> a,
    const Eigen::Ref<const Eigen::Vector3d> b)
{
    const double ab_dot = std::clamp(a.dot(b), -1., 1.);
    const double angle_rad = std::acos(ab_dot);
    return (180. / M_PI) * angle_rad;
}
} // namespace wmtk::components::internal::utils