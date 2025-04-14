#pragma once

#include <Eigen/Dense>

namespace wmtk::components::internal::utils {

Eigen::VectorXd compute_edge_length(const Eigen::MatrixXd& P);

} // namespace wmtk::components::internal::utils