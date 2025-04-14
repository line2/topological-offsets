#include "LongestEdgeInvariant.hpp"

#include <wmtk/Mesh.hpp>
#include <wmtk/utils/Logger.hpp>

namespace wmtk::components::internal::invariants {

LongestEdgeInvariant::LongestEdgeInvariant(const attribute::MeshAttributeHandle& position_handle)
    : Invariant(position_handle.mesh(), true, false, false)
    , m_position_handle(position_handle)
{}

bool LongestEdgeInvariant::before(const simplex::Simplex& s) const
{
    if (s.primitive_type() != PrimitiveType::Edge) {
        return true;
    }

    assert(!mesh().is_boundary(s));
    assert(mesh().top_simplex_type() == PrimitiveType::Triangle);

    const Tuple v0 = s.tuple();
    const Tuple v1 = mesh().switch_tuple(v0, PrimitiveType::Vertex);
    const Tuple v2 = mesh().switch_tuples(v0, {PrimitiveType::Edge, PrimitiveType::Vertex});
    const Tuple v3 = mesh().switch_tuples(
        v0,
        {PrimitiveType::Triangle, PrimitiveType::Edge, PrimitiveType::Vertex});

    const auto pos_acc = mesh().create_const_accessor(m_position_handle.as<double>());

    const Eigen::Vector3d p0 = pos_acc.const_vector_attribute(v0);
    const Eigen::Vector3d p1 = pos_acc.const_vector_attribute(v1);
    const Eigen::Vector3d p2 = pos_acc.const_vector_attribute(v2);
    const Eigen::Vector3d p3 = pos_acc.const_vector_attribute(v3);

    const double e01 = (p1 - p0).squaredNorm();

    const double e02 = (p2 - p0).squaredNorm();
    const double e03 = (p3 - p0).squaredNorm();

    const double e12 = (p2 - p1).squaredNorm();
    const double e13 = (p2 - p1).squaredNorm();

    if (e01 > e02 && e01 > e12) {
        return true;
    }

    if (e01 > e03 && e01 > e13) {
        return true;
    }

    return false;
}

} // namespace wmtk::components::internal::invariants