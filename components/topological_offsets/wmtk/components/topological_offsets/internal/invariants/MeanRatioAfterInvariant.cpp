#include "MeanRatioAfterInvariant.hpp"

#include <wmtk/Mesh.hpp>
#include <wmtk/utils/Logger.hpp>

namespace wmtk::components::internal::invariants {

MeanRatioAfterInvariant::MeanRatioAfterInvariant(
    const attribute::MeshAttributeHandle& mean_ratio_metric_handle,
    const double min_mean_ratio_metric)
    : Invariant(mean_ratio_metric_handle.mesh(), false, true, true)
    , m_mean_ratio_metric_handle(mean_ratio_metric_handle)
    , m_min_mean_ratio_metric(min_mean_ratio_metric)
{
    assert(mean_ratio_metric_handle.primitive_type() == mesh().top_simplex_type());
}

bool MeanRatioAfterInvariant::after(
    const std::vector<Tuple>& top_dimension_tuples_before,
    const std::vector<Tuple>& top_dimension_tuples_after) const
{
    assert(top_dimension_tuples_before.empty());

    const double min_mrm_before =
        mesh().parent_scope([this, &top_dimension_tuples_before]() -> double {
            const auto acc = mesh().create_const_accessor(m_mean_ratio_metric_handle.as<double>());
            double min_mrm = std::numeric_limits<double>::max();
            for (const Tuple& t : top_dimension_tuples_before) {
                min_mrm = std::min(min_mrm, acc.const_scalar_attribute(t));
            }
            return min_mrm;
        });


    const double mrm_min_min = std::min(min_mrm_before, m_min_mean_ratio_metric);

    {
        const auto acc = mesh().create_const_accessor(m_mean_ratio_metric_handle.as<double>());
        for (const Tuple& t : top_dimension_tuples_after) {
            if (acc.const_scalar_attribute(t) <= mrm_min_min) {
                return false;
            }
        }
    }

    return true;
}

} // namespace wmtk::components::internal::invariants