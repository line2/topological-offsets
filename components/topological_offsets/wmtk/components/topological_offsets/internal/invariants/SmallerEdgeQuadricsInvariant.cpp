#include "SmallerEdgeQuadricsInvariant.hpp"

#include <wmtk/components/layered_offsets/utils/Quadrics.hpp>

namespace wmtk::components::internal::invariants {

SmallerEdgeQuadricsInvariant::SmallerEdgeQuadricsInvariant(
    const attribute::MeshAttributeHandle& quadrics_handle,
    const attribute::MeshAttributeHandle& position_handle,
    const double comp_dist)
    : Invariant(quadrics_handle.mesh(), true, false, false)
    , m_quadrics_handle(quadrics_handle)
    , m_position_handle(position_handle)
    , m_comp_dist(comp_dist)
{
    assert(quadrics_handle.primitive_type() == quadrics_handle.mesh().top_simplex_type());
    assert(quadrics_handle.mesh() == position_handle.mesh());
}

bool SmallerEdgeQuadricsInvariant::before(const simplex::Simplex& s) const
{
    if (s.primitive_type() != PrimitiveType::Edge) {
        return true;
    }

    using utils::Quadrics;

    const auto q_acc = mesh().create_const_accessor(m_quadrics_handle.as<double>());

    const double sq_dist = Quadrics::squared_distance_of_projected_simplex_to_optimum(
        m_position_handle,
        m_quadrics_handle,
        s);

    return sq_dist < m_comp_dist * m_comp_dist;
}

} // namespace wmtk::components::internal::invariants