#pragma once

#include <wmtk/attribute/MeshAttributeHandle.hpp>
#include <wmtk/invariants/Invariant.hpp>

namespace wmtk::components::internal::invariants {

/**
 * Returns false (!) if all adjacent triangles are converged.
 */
class ConvergenceInvariant : public Invariant
{
public:
    ConvergenceInvariant(const attribute::MeshAttributeHandle& convergence_handle);

    bool before(const simplex::Simplex& s) const override;

private:
    const attribute::MeshAttributeHandle m_convergence_handle;
};

} // namespace wmtk::components::internal::invariants