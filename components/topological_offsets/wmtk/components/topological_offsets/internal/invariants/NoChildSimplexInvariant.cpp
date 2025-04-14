#include "NoChildSimplexInvariant.hpp"

#include <wmtk/Mesh.hpp>

namespace wmtk::components::internal::invariants {
NoChildSimplexInvariant::NoChildSimplexInvariant(const Mesh& parent_mesh, const Mesh& child_mesh)
    : Invariant(parent_mesh, true, false, false)
    , m_child_mesh(child_mesh)
{}
bool NoChildSimplexInvariant::before(const simplex::Simplex& s) const
{
    return !mesh().simplex_is_in_child(m_child_mesh, s);
}
} // namespace wmtk::components::internal::invariants
