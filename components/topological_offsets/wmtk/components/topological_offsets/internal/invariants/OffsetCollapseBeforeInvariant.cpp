#include "OffsetCollapseBeforeInvariant.hpp"

#include <wmtk/Mesh.hpp>
#include <wmtk/components/topological_offsets/internal/utils/normal_angle.hpp>
#include <wmtk/simplex/top_dimension_cofaces.hpp>
#include <wmtk/utils/Logger.hpp>

namespace wmtk::components::internal::invariants {


OffsetCollapseBeforeInvariant::OffsetCollapseBeforeInvariant(
    const attribute::MeshAttributeHandle& normal_samples_handle,
    const attribute::MeshAttributeHandle& position_handle,
    const double max_normal_deviation,
    const int64_t collapse_type)
    : Invariant(normal_samples_handle.mesh(), true, false, false)
    , m_normal_samples_handle(normal_samples_handle)
    , m_position_handle(position_handle)
    , m_max_normal_deviation(max_normal_deviation)
    , m_collapse_type(collapse_type)
{
    assert(
        normal_samples_handle.primitive_type() == normal_samples_handle.mesh().top_simplex_type());
    assert(normal_samples_handle.mesh() == position_handle.mesh());
}

bool OffsetCollapseBeforeInvariant::before(const simplex::Simplex& s) const
{
    const Tuple v0 = s.tuple();
    const Tuple v1 = mesh().switch_tuple(v0, PrimitiveType::Vertex);

    const auto pos_acc = mesh().create_const_accessor(m_position_handle.as<double>());

    const auto normal_samples_acc =
        mesh().create_const_accessor(m_normal_samples_handle.as<double>());

    const Eigen::Vector3d p0 = pos_acc.const_vector_attribute(v0);
    const Eigen::Vector3d p1 = pos_acc.const_vector_attribute(v1);

    const Eigen::Vector3d e_dir = (p1 - p0).normalized();

    std::vector<Tuple> tris;
    switch (m_collapse_type) {
    case 1: {
        // copy other --> v1
        tris = simplex::top_dimension_cofaces_tuples(mesh(), simplex::Simplex::vertex(mesh(), v1));
        break;
    }
    case 2: {
        // copy tuple --> v0
        tris = simplex::top_dimension_cofaces_tuples(mesh(), simplex::Simplex::vertex(mesh(), v0));
        break;
    }
    default: return false;
    }

    double min_angle = std::numeric_limits<double>::max();
    double max_angle = std::numeric_limits<double>::lowest();

    // compute angles between normals and edge
    for (const Tuple& tri : tris) {
        const auto normal_samples_vec = normal_samples_acc.const_vector_attribute(tri);
        std::array<Eigen::Vector3d, 4> normal_samples;
        normal_samples[0] = normal_samples_vec.block(3, 0, 3, 1);
        normal_samples[1] = normal_samples_vec.block(6, 0, 3, 1);
        normal_samples[2] = normal_samples_vec.block(9, 0, 3, 1);
        normal_samples[3] = normal_samples_vec.block(12, 0, 3, 1);

        for (const auto& n : normal_samples) {
            const double angle = utils::normal_angle_180(e_dir, n);
            min_angle = std::min(min_angle, angle);
            max_angle = std::max(max_angle, angle);
        }
    }

    const double nd = max_angle - min_angle;

    return nd < m_max_normal_deviation;
}

} // namespace wmtk::components::internal::invariants