#pragma once

#include <wmtk/attribute/MeshAttributeHandle.hpp>
#include <wmtk/invariants/Invariant.hpp>

namespace wmtk::components::internal::invariants {

/**
 * Compute the normal deviation after an operation and compare it to a user defined max.
 */
class NormalDeviationAfterInvariant : public Invariant
{
public:
    NormalDeviationAfterInvariant(
        const attribute::MeshAttributeHandle& normal_deviation_handle,
        const double max_normal_deviation,
        const bool compare_with_before);

    bool after(
        const std::vector<Tuple>& top_dimension_tuples_before,
        const std::vector<Tuple>& top_dimension_tuples_after) const override;

private:
    const attribute::MeshAttributeHandle m_normal_deviation_handle;
    const double m_max_normal_deviation;
    const bool m_compare_with_before;
};

} // namespace wmtk::components::internal::invariants