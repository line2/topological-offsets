#pragma once

#include <wmtk/attribute/MeshAttributeHandle.hpp>
#include <wmtk/invariants/Invariant.hpp>

namespace wmtk::components::internal::invariants {

/**
 * Compute the distance from the center of the quadric-projected edge (both vertices moved to the
 * quadric minimum) to the edge's quadric minimum. Compare this distance with the given value
 */
class SmallerEdgeQuadricsInvariant : public Invariant
{
public:
    SmallerEdgeQuadricsInvariant(
        const attribute::MeshAttributeHandle& quadrics_handle,
        const attribute::MeshAttributeHandle& position_handle,
        const double comp_dist);

    bool before(const simplex::Simplex& s) const override;

private:
    const attribute::MeshAttributeHandle m_quadrics_handle;
    const attribute::MeshAttributeHandle m_position_handle;
    const double m_comp_dist;
};

} // namespace wmtk::components::internal::invariants