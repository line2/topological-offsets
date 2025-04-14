#pragma once

#include <wmtk/attribute/MeshAttributeHandle.hpp>
#include <wmtk/invariants/Invariant.hpp>

namespace wmtk::components::internal::invariants {

/**
 * Compute the normal deviation after an operation and compare it to a user defined max.
 */
class AmipsBeforeInvariant : public Invariant
{
public:
    AmipsBeforeInvariant(
        const attribute::MeshAttributeHandle& amips_handle,
        const double min_amips);

    bool before(const simplex::Simplex& s) const override;

private:
    const attribute::MeshAttributeHandle m_amips_handle;
    const double m_min_amips;
};

} // namespace wmtk::components::internal::invariants