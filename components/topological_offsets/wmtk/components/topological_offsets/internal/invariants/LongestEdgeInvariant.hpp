#pragma once

#include <wmtk/attribute/MeshAttributeHandle.hpp>
#include <wmtk/invariants/Invariant.hpp>

namespace wmtk::components::internal::invariants {

/**
 * Compute the distance from the center of the quadric-projected edge (both vertices moved to the
 * quadric minimum) to the edge's quadric minimum. Compare this distance with the given value
 */
class LongestEdgeInvariant : public Invariant
{
public:
    LongestEdgeInvariant(const attribute::MeshAttributeHandle& position_handle);

    bool before(const simplex::Simplex& s) const override;

private:
    const attribute::MeshAttributeHandle m_position_handle;
};

} // namespace wmtk::components::internal::invariants