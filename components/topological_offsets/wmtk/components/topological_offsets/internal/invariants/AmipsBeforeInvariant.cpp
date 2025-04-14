#include "AmipsBeforeInvariant.hpp"

#include <wmtk/Mesh.hpp>
#include <wmtk/simplex/top_dimension_cofaces.hpp>
#include <wmtk/utils/Logger.hpp>

namespace wmtk::components::internal::invariants {

AmipsBeforeInvariant::AmipsBeforeInvariant(
    const attribute::MeshAttributeHandle& amips_handle,
    const double min_amips)
    : Invariant(amips_handle.mesh(), true, false, false)
    , m_amips_handle(amips_handle)
    , m_min_amips(min_amips)
{
    assert(m_amips_handle.dimension() == 1);
}

bool AmipsBeforeInvariant::before(const simplex::Simplex& s) const
{
    const auto tuples = simplex::top_dimension_cofaces_tuples(mesh(), s);

    const auto acc = mesh().create_const_accessor(m_amips_handle.as<double>());

    for (const Tuple& t : tuples) {
        if (acc.const_scalar_attribute(t) > m_min_amips) {
            return true;
        }
    }

    return false;
}

} // namespace wmtk::components::internal::invariants