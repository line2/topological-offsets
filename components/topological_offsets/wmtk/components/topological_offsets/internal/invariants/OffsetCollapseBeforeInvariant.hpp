#pragma once

#include <wmtk/attribute/MeshAttributeHandle.hpp>
#include <wmtk/invariants/Invariant.hpp>

namespace wmtk::components::internal::invariants {

class OffsetCollapseBeforeInvariant : public Invariant
{
public:
    OffsetCollapseBeforeInvariant(
        const attribute::MeshAttributeHandle& normal_samples_handle,
        const attribute::MeshAttributeHandle& position_handle,
        const double max_normal_deviation,
        const int64_t collapse_type);

    bool before(const simplex::Simplex& s) const override;

private:
    const attribute::MeshAttributeHandle m_normal_samples_handle;
    const attribute::MeshAttributeHandle m_position_handle;
    const double m_max_normal_deviation;
    const int64_t m_collapse_type;
};
} // namespace wmtk::components::internal::invariants
