#pragma once

#include <wmtk/attribute/MeshAttributeHandle.hpp>
#include <wmtk/invariants/Invariant.hpp>

namespace wmtk::components::internal::invariants {

/**
 * Compute the mean ratio metric after an operation and compare it to a user defined min.
 */
class MeanRatioAfterInvariant : public Invariant
{
public:
    MeanRatioAfterInvariant(
        const attribute::MeshAttributeHandle& mean_ratio_metric_handle,
        const double min_mean_ratio_metric);

    bool after(
        const std::vector<Tuple>& top_dimension_tuples_before,
        const std::vector<Tuple>& top_dimension_tuples_after) const override;

private:
    const attribute::MeshAttributeHandle m_mean_ratio_metric_handle;
    const double m_min_mean_ratio_metric;
};

} // namespace wmtk::components::internal::invariants