#pragma once

#include <wmtk/attribute/MeshAttributeHandle.hpp>
#include <wmtk/invariants/Invariant.hpp>

namespace wmtk::components::internal::invariants {

/**
 * Compute the distance from the center of the quadric-projected edge (both vertices moved to the
 * quadric minimum) to the edge's quadric minimum. Compare this distance with the given value
 */
class OffsetSwapInvariant : public Invariant
{
public:
    OffsetSwapInvariant(
        const attribute::MeshAttributeHandle& normal_samples_handle,
        const attribute::MeshAttributeHandle& position_handle,
        const double max_normal_deviation);

    bool before(const simplex::Simplex& s) const override;

private:
    const attribute::MeshAttributeHandle m_normal_samples_handle;
    const attribute::MeshAttributeHandle m_position_handle;
    const double m_max_normal_deviation;
};

} // namespace wmtk::components::internal::invariants