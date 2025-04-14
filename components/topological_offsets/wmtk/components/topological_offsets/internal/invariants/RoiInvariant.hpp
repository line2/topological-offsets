#pragma once

#include <wmtk/attribute/MeshAttributeHandle.hpp>
#include <wmtk/invariants/Invariant.hpp>

namespace wmtk::components::internal::invariants {

/**
 * Returns true if any adjacent top-simplex is tagged.
 */
class RoiInvariant : public Invariant
{
public:
    RoiInvariant(const attribute::MeshAttributeHandle& roi_handle, const int64_t tag_value = 1);

    bool before(const simplex::Simplex& s) const override;

private:
    const attribute::MeshAttributeHandle m_roi_handle;
    const int64_t m_tag_value;
};

} // namespace wmtk::components::internal::invariants