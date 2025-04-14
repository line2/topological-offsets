#include "compute_tet_amips.hpp"

#include <wmtk/function/utils/amips.hpp>

namespace wmtk::components::internal::utils {

Eigen::VectorXd compute_tet_amips(const Eigen::MatrixXd& P)
{
    assert(P.rows() == 3); // rows --> attribute dimension
    assert(P.cols() == 4);
    // tet
    std::array<double, 12> pts;
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            pts[3 * i + j] = P(j, i);
        }
    }
    const double a = function::utils::Tet_AMIPS_energy(pts);
    return Eigen::VectorXd::Constant(1, a);
};

Eigen::VectorXd compute_tri_amips(const Eigen::MatrixXd& P)
{
    assert(P.rows() == 2); // rows --> attribute dimension
    assert(P.cols() == 3);
    // tet
    std::array<double, 6> pts;
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            pts[2 * i + j] = P(j, i);
        }
    }
    const double a = function::utils::Tri_AMIPS_energy(pts);
    return Eigen::VectorXd::Constant(1, a);
};

} // namespace wmtk::components::internal::utils