#include "OffsetSwapInvariant.hpp"

#include <wmtk/Mesh.hpp>
#include <wmtk/components/topological_offsets/internal/utils/mean_ratio_metric.hpp>
#include <wmtk/components/topological_offsets/internal/utils/normal_angle.hpp>
#include <wmtk/utils/Logger.hpp>

namespace wmtk::components::internal::invariants {

OffsetSwapInvariant::OffsetSwapInvariant(
    const attribute::MeshAttributeHandle& normal_samples_handle,
    const attribute::MeshAttributeHandle& position_handle,
    const double max_normal_deviation)
    : Invariant(normal_samples_handle.mesh(), true, false, false)
    , m_normal_samples_handle(normal_samples_handle)
    , m_position_handle(position_handle)
    , m_max_normal_deviation(max_normal_deviation)
{
    assert(
        normal_samples_handle.primitive_type() == normal_samples_handle.mesh().top_simplex_type());
    assert(&(normal_samples_handle.mesh()) == &(position_handle.mesh()));
}

bool OffsetSwapInvariant::before(const simplex::Simplex& s) const
{
    if (s.primitive_type() != PrimitiveType::Edge) {
        return true;
    }

    assert(!mesh().is_boundary(s));

    const Tuple v0 = s.tuple();
    const Tuple v1 = mesh().switch_tuple(v0, PrimitiveType::Vertex);
    const Tuple v2 = mesh().switch_tuples(v0, {PrimitiveType::Edge, PrimitiveType::Vertex});
    const Tuple v3 = mesh().switch_tuples(
        v0,
        {PrimitiveType::Triangle, PrimitiveType::Edge, PrimitiveType::Vertex});

    const auto pos_acc = mesh().create_const_accessor(m_position_handle.as<double>());

    const auto normal_samples_acc =
        mesh().create_const_accessor(m_normal_samples_handle.as<double>());

    const Eigen::Vector3d p0 = pos_acc.const_vector_attribute(v0);
    const Eigen::Vector3d p1 = pos_acc.const_vector_attribute(v1);
    const Eigen::Vector3d p2 = pos_acc.const_vector_attribute(v2);
    const Eigen::Vector3d p3 = pos_acc.const_vector_attribute(v3);

    const Eigen::Vector3d e_old_dir = (p1 - p0).normalized();
    const Eigen::Vector3d e_new_dir = (p2 - p3).normalized();

    const auto f0_normal_samples = normal_samples_acc.const_vector_attribute(v0);
    const auto f1_normal_samples = normal_samples_acc.const_vector_attribute(v3);

    std::array<Eigen::Vector3d, 8> normal_samples;
    normal_samples[0] = f0_normal_samples.block(3, 0, 3, 1);
    normal_samples[1] = f0_normal_samples.block(6, 0, 3, 1);
    normal_samples[2] = f0_normal_samples.block(9, 0, 3, 1);
    normal_samples[3] = f0_normal_samples.block(12, 0, 3, 1);
    normal_samples[4] = f1_normal_samples.block(3, 0, 3, 1);
    normal_samples[5] = f1_normal_samples.block(6, 0, 3, 1);
    normal_samples[6] = f1_normal_samples.block(9, 0, 3, 1);
    normal_samples[7] = f1_normal_samples.block(12, 0, 3, 1);

    // compute angles betwee normals and current or swapped edge
    auto edge_normal_deviation =
        [&normal_samples](const Eigen::Ref<const Eigen::Vector3d> e) -> double {
        double min_angle = std::numeric_limits<double>::max();
        double max_angle = std::numeric_limits<double>::lowest();

        for (const auto& n : normal_samples) {
            const double angle = utils::normal_angle_180(e, n);
            min_angle = std::min(min_angle, angle);
            max_angle = std::max(max_angle, angle);
        }

        return max_angle - min_angle;
    };

    const double nd_old = edge_normal_deviation(e_old_dir);
    const double nd_new = edge_normal_deviation(e_new_dir);

    if (nd_old < m_max_normal_deviation && nd_new >= m_max_normal_deviation) {
        return false;
    }

    //const double nd_metric_old = (nd_old / 90.0);
    //const double nd_metric_new = (nd_new / 90.0);

    const double mrm_old_0 = utils::mean_ratio_metric(p0, p1, p2);
    const double mrm_old_1 = utils::mean_ratio_metric(p0, p1, p3);
    const double mrm_old = std::min(mrm_old_0, mrm_old_1);

    const double mrm_new_0 = utils::mean_ratio_metric(p2, p3, p0);
    const double mrm_new_1 = utils::mean_ratio_metric(p2, p3, p1);
    const double mrm_new = std::min(mrm_new_0, mrm_new_1);

    // const double q_old = mrm_old_inv * nd_metric_old * nd_metric_old * nd_metric_old;
    // const double q_new = mrm_new_inv * nd_metric_new * nd_metric_new * nd_metric_new;

    // return q_new < q_old;
    return mrm_new > mrm_old;
}

} // namespace wmtk::components::internal::invariants