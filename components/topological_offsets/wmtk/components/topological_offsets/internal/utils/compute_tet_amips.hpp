#pragma once

#include <Eigen/Dense>

namespace wmtk::components::internal::utils {

Eigen::VectorXd compute_tet_amips(const Eigen::MatrixXd& P);

Eigen::VectorXd compute_tri_amips(const Eigen::MatrixXd& P);

} // namespace wmtk::components::internal::utils