#include "ConvergenceInvariant.hpp"

#include <wmtk/Mesh.hpp>
#include <wmtk/simplex/top_dimension_cofaces.hpp>
#include <wmtk/utils/Logger.hpp>

namespace wmtk::components::internal::invariants {

ConvergenceInvariant::ConvergenceInvariant(const attribute::MeshAttributeHandle& convergence_handle)
    : Invariant(convergence_handle.mesh(), true, false, false)
    , m_convergence_handle(convergence_handle)
{
    assert(m_convergence_handle.dimension() == 1);
    assert(
        mesh().top_simplex_type() == PrimitiveType::Triangle ||
        mesh().top_simplex_type() == PrimitiveType::Edge);
}

bool ConvergenceInvariant::before(const simplex::Simplex& s) const
{
    const auto tuples = simplex::top_dimension_cofaces_tuples(mesh(), s);

    const auto acc = mesh().create_const_accessor(m_convergence_handle.as<int64_t>());

    for (const Tuple& t : tuples) {
        if (acc.const_scalar_attribute(t) == 0) {
            return true;
        }
    }

    return false;
}

} // namespace wmtk::components::internal::invariants