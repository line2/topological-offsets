#include "compute_edge_length.hpp"

namespace wmtk::components::internal::utils {

Eigen::VectorXd compute_edge_length(const Eigen::MatrixXd& P)
{
    assert(P.cols() == 2);
    assert(P.rows() == 2 || P.rows() == 3);
    return Eigen::VectorXd::Constant(1, (P.col(0) - P.col(1)).norm());
};

} // namespace wmtk::components::internal::utils