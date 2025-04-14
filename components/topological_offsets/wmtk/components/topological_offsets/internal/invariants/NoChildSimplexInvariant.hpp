#pragma once

#include <wmtk/invariants/Invariant.hpp>

namespace wmtk::components::internal::invariants {
/**
 * This invariant returns true if the given simplex is not mapping to any child mesh.
 */
class NoChildSimplexInvariant : public Invariant
{
public:
    NoChildSimplexInvariant(const Mesh& parent_mesh, const Mesh& child_mesh);
    bool before(const simplex::Simplex& s) const override;

private:
    const Mesh& m_child_mesh;
};
} // namespace wmtk::components::internal::invariants
