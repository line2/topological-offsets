#include "bbox_from_mesh.hpp"

#include <wmtk/simplex/faces_single_dimension.hpp>

namespace wmtk::utils {


Eigen::MatrixXd bbox_from_mesh(const attribute::MeshAttributeHandle& position_handle)
{
    const Mesh& m = position_handle.mesh();

    const int64_t dim = position_handle.dimension();

    Eigen::MatrixXd bbox;
    bbox.resize(2, dim);
    for (int64_t d = 0; d < dim; ++d) {
        bbox(0, d) = std::numeric_limits<double>::max();
        bbox(1, d) = std::numeric_limits<double>::lowest();
    }

    const auto p_acc = m.create_const_accessor<double>(position_handle);

    for (const Tuple& t : m.get_all(PrimitiveType::Vertex)) {
        const auto p = p_acc.const_vector_attribute(t);
        for (int64_t d = 0; d < dim; ++d) {
            bbox(0, d) = std::min(bbox(0, d), p[d]);
            bbox(1, d) = std::max(bbox(1, d), p[d]);
        }
    }

    return bbox;
}

Eigen::MatrixXd bbox_from_mesh(
    const attribute::MeshAttributeHandle& position_handle,
    const attribute::MeshAttributeHandle& tag_handle,
    const int64_t tag_value)
{
    const Mesh& m = position_handle.mesh();

    const int64_t dim = position_handle.dimension();

    Eigen::MatrixXd bbox;
    bbox.resize(2, dim);
    for (int64_t d = 0; d < dim; ++d) {
        bbox(0, d) = std::numeric_limits<double>::max();
        bbox(1, d) = std::numeric_limits<double>::lowest();
    }

    const auto p_acc = m.create_const_accessor<double>(position_handle);
    const auto t_acc = m.create_const_accessor<int64_t>(tag_handle);

    for (const Tuple& s : m.get_all(t_acc.primitive_type())) {
        if (t_acc.const_scalar_attribute(s) != tag_value) {
            continue;
        }

        const auto vs = simplex::faces_single_dimension(
            m,
            simplex::Simplex(m, t_acc.primitive_type(), s),
            PrimitiveType::Vertex);

        for (const simplex::Simplex& t : vs) {
            const Eigen::Vector3d p = p_acc.const_vector_attribute(t);
            for (int64_t d = 0; d < dim; ++d) {
                bbox(0, d) = std::min(bbox(0, d), p[d]);
                bbox(1, d) = std::max(bbox(1, d), p[d]);
            }
        }
    }

    return bbox;
}

double bbox_diagonal_from_mesh(const attribute::MeshAttributeHandle& position_handle)
{
    const Eigen::MatrixXd bbox = bbox_from_mesh(position_handle);

    const double d = (bbox.row(0) - bbox.row(1)).norm();

    return d;
}

double bbox_diagonal_from_mesh(
    const attribute::MeshAttributeHandle& position_handle,
    const attribute::MeshAttributeHandle& tag_handle,
    const int64_t tag_value)
{
    const Eigen::MatrixXd bbox = bbox_from_mesh(position_handle, tag_handle, tag_value);

    const double d = (bbox.row(0) - bbox.row(1)).norm();

    return d;
}

} // namespace wmtk::utils