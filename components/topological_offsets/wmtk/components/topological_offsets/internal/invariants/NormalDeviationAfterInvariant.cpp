#include "NormalDeviationAfterInvariant.hpp"

#include <wmtk/Mesh.hpp>
#include <wmtk/utils/Logger.hpp>

namespace wmtk::components::internal::invariants {

NormalDeviationAfterInvariant::NormalDeviationAfterInvariant(
    const attribute::MeshAttributeHandle& normal_deviation_handle,
    const double max_normal_deviation,
    const bool compare_with_before)
    : Invariant(normal_deviation_handle.mesh(), false, true, true)
    , m_normal_deviation_handle(normal_deviation_handle)
    , m_max_normal_deviation(max_normal_deviation)
    , m_compare_with_before(compare_with_before)
{}

bool NormalDeviationAfterInvariant::after(
    const std::vector<Tuple>& top_dimension_tuples_before,
    const std::vector<Tuple>& top_dimension_tuples_after) const
{
    // assert(top_dimension_tuples_before.empty());

    if (m_compare_with_before) {
        const double max_nd_before =
            mesh().parent_scope([this, &top_dimension_tuples_before]() -> double {
                const auto acc =
                    mesh().create_const_accessor(m_normal_deviation_handle.as<double>());
                double max_nd = std::numeric_limits<double>::lowest();
                for (const Tuple& t : top_dimension_tuples_before) {
                    max_nd = std::max(max_nd, acc.const_scalar_attribute(t));
                }
                return max_nd;
            });

        // only prohibit moves that degrade already good quality
        if (max_nd_before > m_max_normal_deviation) {
            return true;
        }
    }

    {
        const auto acc = mesh().create_const_accessor(m_normal_deviation_handle.as<double>());

        for (const Tuple& t : top_dimension_tuples_after) {
            if (acc.const_scalar_attribute(t) >= m_max_normal_deviation) {
                return false;
            }
        }
    }

    return true;
}

} // namespace wmtk::components::internal::invariants