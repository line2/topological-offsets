#include "RoiInvariant.hpp"

#include <wmtk/Mesh.hpp>
#include <wmtk/simplex/top_dimension_cofaces_iterable.hpp>
#include <wmtk/utils/Logger.hpp>

namespace wmtk::components::internal::invariants {

RoiInvariant::RoiInvariant(
    const attribute::MeshAttributeHandle& roi_handle,
    const int64_t tag_value)
    : Invariant(roi_handle.mesh(), true, false, false)
    , m_roi_handle(roi_handle)
    , m_tag_value(tag_value)
{
    assert(m_roi_handle.dimension() == 1);
}

bool RoiInvariant::before(const simplex::Simplex& s) const
{
    const auto acc = mesh().create_const_accessor(m_roi_handle.as<int64_t>());

    auto tuples = simplex::top_dimension_cofaces_iterable(mesh(), s);
    for (const Tuple& t : tuples) {
        if (acc.const_scalar_attribute(t) == m_tag_value) {
            return true;
        }
    }

    return false;
}

} // namespace wmtk::components::internal::invariants