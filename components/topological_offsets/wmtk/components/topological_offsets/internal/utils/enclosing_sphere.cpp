#include "enclosing_sphere.hpp"

#include <wmtk/simplex/faces_single_dimension.hpp>

namespace wmtk::components::internal::utils {

std::tuple<Eigen::VectorXd, double> enclosing_sphere(
    const attribute::MeshAttributeHandle& pos_handle,
    const simplex::Simplex& s)
{
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor, 3> VecType;

    const Mesh& m = pos_handle.mesh();

    const auto p_acc = m.create_const_accessor(pos_handle.as<double>());

    if (m.top_simplex_type() == PrimitiveType::Vertex) {
        Eigen::VectorXd p = p_acc.const_vector_attribute(s);
        return std::make_tuple(p, 0);
    }

    std::vector<VecType> pts;
    VecType center;
    center.setZero(p_acc.dimension(), 1);
    for (const simplex::Simplex& v : simplex::faces_single_dimension(m, s, PrimitiveType::Vertex)) {
        pts.emplace_back(p_acc.const_vector_attribute(v));
        center += pts.back();
    }
    center /= pts.size();

    assert(pts.size() > 1);

    // find the longest edge
    double l_squared = -1;
    for (const VecType& p : pts) {
        const double ll = (p - center).squaredNorm();
        l_squared = std::max(l_squared, ll);
    }

    return std::make_tuple(center, std::sqrt(l_squared));
}

} // namespace wmtk::components::internal::utils
