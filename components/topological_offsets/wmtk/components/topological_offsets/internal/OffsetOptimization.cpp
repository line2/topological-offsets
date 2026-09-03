#include "OffsetOptimization.hpp"

#include <fstream>
#include <wmtk/Scheduler.hpp>
#include <wmtk/function/LocalNeighborsSumFunction.hpp>
#include <wmtk/function/simplex/AMIPS.hpp>
#include <wmtk/function/utils/amips.hpp>
#include <wmtk/invariants/EdgeValenceInvariant.hpp>
#include <wmtk/invariants/InteriorEdgeInvariant.hpp>
#include <wmtk/invariants/MaxFunctionInvariant.hpp>
#include <wmtk/invariants/MultiMeshLinkConditionInvariant.hpp>
#include <wmtk/invariants/SeparateSubstructuresInvariant.hpp>
#include <wmtk/invariants/SimplexInversionInvariant.hpp>
#include <wmtk/invariants/Swap23EnergyBeforeInvariant.hpp>
#include <wmtk/invariants/Swap32EnergyBeforeInvariant.hpp>
#include <wmtk/invariants/Swap44EnergyBeforeInvariant.hpp>
#include <wmtk/invariants/Swap56EnergyBeforeInvariant.hpp>
#include <wmtk/invariants/TodoInvariant.hpp>
#include <wmtk/io/ParaviewWriter.hpp>
#include <wmtk/multimesh/consolidate.hpp>
#include <wmtk/operations/EdgeCollapse.hpp>
#include <wmtk/operations/EdgeSplit.hpp>
#include <wmtk/operations/MinOperationSequence.hpp>
#include <wmtk/operations/OptimizationSmoothing.hpp>
#include <wmtk/operations/OrOperationSequence.hpp>
#include <wmtk/operations/attribute_new/CollapseNewAttributeStrategy.hpp>
#include <wmtk/operations/attribute_new/NewAttributeStrategy.hpp>
#include <wmtk/operations/composite/TetEdgeSwap.hpp>
#include <wmtk/operations/composite/TetFaceSwap.hpp>
#include <wmtk/operations/composite/TriEdgeSwap.hpp>
#include <wmtk/simplex/cofaces_single_dimension_iterable.hpp>
#include <wmtk/simplex/faces_single_dimension.hpp>
#include <wmtk/simplex/link_single_dimension_iterable.hpp>
#include <wmtk/simplex/top_dimension_cofaces.hpp>
#include <wmtk/simplex/top_dimension_cofaces_iterable.hpp>
#include <wmtk/utils/Logger.hpp>
#include <wmtk/utils/TupleInspector.hpp>
#include <wmtk/utils/orient.hpp>

#include <wmtk/components/topological_offsets/internal/invariants/NoChildSimplexInvariant.hpp>
#include <wmtk/components/topological_offsets/internal/utils/bvh_from_mesh.hpp>
#include <wmtk/components/topological_offsets/internal/utils/bvh_to_mesh_index_map.hpp>

#include "functions/AmipsOffsetOptimization.hpp"
#include "functions/AmipsOptimization.hpp"
#include "invariants/AmipsBeforeInvariant.hpp"
#include "invariants/ConvergenceInvariant.hpp"
#include "invariants/LongestEdgeInvariant.hpp"
#include "invariants/NormalDeviationAfterInvariant.hpp"
#include "invariants/OffsetCollapseBeforeInvariant.hpp"
#include "invariants/OffsetSwapInvariant.hpp"
#include "invariants/RoiInvariant.hpp"
#include "utils/Quadrics.hpp"
#include "utils/compute_edge_length.hpp"
#include "utils/compute_tet_amips.hpp"
#include "utils/is_tet_in_offset.hpp"
#include "utils/mean_ratio_metric.hpp"
#include "utils/normal_angle.hpp"
#include "utils/write_mesh.hpp"


namespace wmtk::components::internal {

std::shared_ptr<operations::CollapseNewAttributeStrategy<double>> keep_substructure_strategy(
    Mesh& m,
    const attribute::MeshAttributeHandle& attr)
{
    auto strat = std::make_shared<operations::CollapseNewAttributeStrategy<double>>(attr);
    strat->set_simplex_predicate([&m](const simplex::Simplex& s) -> bool {
        for (const auto& child : m.get_child_meshes()) {
            if (m.simplex_is_in_child(*child, s)) {
                return true;
            }
            // if (!m.map_to_child_tuples(*child, s).empty()) {
            //     return true;
            // }
        }
        return false;
    });
    strat->set_strategy(operations::CollapseBasicStrategy::Default);

    return strat;
}

OffsetOptimization::OffsetOptimization(Mesh& mesh, attribute::MeshAttributeHandle& pos_handle)
    : m_mesh(mesh)
    , m_pos_handle(pos_handle)
{
    m_metrics_report["metrics"] = nlohmann::json::array();
}

void OffsetOptimization::set_offset_mesh(
    const std::shared_ptr<Mesh> offset_mesh,
    const std::shared_ptr<Mesh> input_mesh,
    attribute::MeshAttributeHandle offset_pos_handle,
    attribute::MeshAttributeHandle input_pos_handle)
{
    assert(input_mesh != offset_mesh);

    bool o_is_child = false;
    bool i_is_child = false;
    for (const auto child : m_mesh.get_all_child_meshes()) {
        if (child == offset_mesh) {
            o_is_child = true;
            continue;
        }
        if (child == input_mesh) {
            i_is_child = true;
            continue;
        }
    }
    assert(o_is_child);
    assert(i_is_child);

    m_offset_mesh = offset_mesh;
    m_input_mesh = input_mesh;
    m_offset_pos_handle = offset_pos_handle;
    m_input_pos_handle = input_pos_handle;
}

void OffsetOptimization::set_input_bvh(const std::shared_ptr<utils::OffsetBvh>& input_bvh)
{
    m_input_bvh = input_bvh;
}


void OffsetOptimization::set_offset_distance(const double offset_distance)
{
    if (offset_distance < 0) {
        logger().warn(
            "Setting a negative offset distance might cause undefined behavior. Offset distance is "
            "{}.",
            offset_distance);
    }
    logger().info("offset_distance = {}", offset_distance);
    m_offset_distance = offset_distance;
}

void OffsetOptimization::set_normal_deviation(
    const double min_normal_deviation,
    const double max_normal_deviation)
{
    logger().info("min_normal_deviation = {}", min_normal_deviation);
    logger().info("max_normal_deviation = {}", max_normal_deviation);
    m_min_normal_deviation = min_normal_deviation;
    m_max_normal_deviation = max_normal_deviation;
}

void OffsetOptimization::set_edge_length_constraints(
    const double min_edge_length,
    const double max_edge_length)
{
    logger().info("min_edge_length = {}", min_edge_length);
    logger().info("max_edge_length = {}", max_edge_length);
    m_min_edge_length = min_edge_length;
    m_max_edge_length = max_edge_length;
}

void OffsetOptimization::set_tetrahedron_tag_handle(
    const attribute::MeshAttributeHandle& tet_tag_handle)
{
    m_tet_tag_handle = tet_tag_handle;
}

void OffsetOptimization::init_embedding_optimization()
{
    logger().info("Init embedding optimization");

    //////////////////////////////////
    // Region of interest (roi)
    m_embedding_roi_attribute =
        m_mesh.register_attribute<int64_t>("offset_optimization_roi", m_mesh.top_simplex_type(), 1);

    //////////////////////////////////
    // AMIPS
    m_embedding_amips_attribute = m_mesh.register_attribute<double>(
        "offset_optimization_amips",
        m_mesh.top_simplex_type(),
        1);

    if (m_mesh.top_simplex_type() == PrimitiveType::Tetrahedron) {
        m_embedding_amips_transfer =
            std::make_shared<wmtk::operations::SingleAttributeTransferStrategy<double, double>>(
                m_embedding_amips_attribute,
                m_pos_handle,
                utils::compute_tet_amips);
    } else if (m_mesh.top_simplex_type() == PrimitiveType::Triangle) {
        m_embedding_amips_transfer =
            std::make_shared<wmtk::operations::SingleAttributeTransferStrategy<double, double>>(
                m_embedding_amips_attribute,
                m_pos_handle,
                utils::compute_tri_amips);
    }
    m_embedding_amips_transfer->run_on_all();

    m_embedding_amips = std::make_shared<function::AMIPS>(m_mesh, m_pos_handle);

    m_amips_energy = std::make_shared<function::LocalNeighborsSumFunction>(
        m_mesh,
        m_pos_handle,
        *m_embedding_amips);

    //////////////////////////////////
    // Storing edge lengths
    m_embedding_edge_length_attribute = m_mesh.register_attribute<double>(
        "offset_optimization_edge_length",
        PrimitiveType::Edge,
        1);
    m_embedding_edge_length_transfer =
        std::make_shared<wmtk::operations::SingleAttributeTransferStrategy<double, double>>(
            m_embedding_edge_length_attribute,
            m_pos_handle,
            utils::compute_edge_length);
    m_embedding_edge_length_transfer->run_on_all();

    //////////////////////////////////
    // Target edge length
    m_embedding_target_edge_length_attribute = m_mesh.register_attribute<double>(
        "offset_optimization_embedding_target_edge_length",
        PrimitiveType::Edge,
        1,
        false,
        m_max_edge_length); // defaults to max edge length

    auto compute_target_edge_length =
        [this](const Eigen::MatrixXd& P, const std::vector<Tuple>& tuples) -> Eigen::VectorXd {
        assert(P.rows() == 1); // rows --> attribute dimension

        auto target_edge_length_accessor =
            m_mesh.create_accessor(m_embedding_target_edge_length_attribute.as<double>());

        const double current_target_edge_length =
            target_edge_length_accessor.const_scalar_attribute(tuples[0]);

        const auto el_acc =
            m_mesh.create_const_accessor(m_embedding_edge_length_attribute.as<double>());
        const double el_current = el_acc.const_scalar_attribute(tuples[0]);

        const double max_amips = P.maxCoeff();

        double new_target_edge_length = el_current;
        if (max_amips > 100) {
            new_target_edge_length *= 0.5;
        } else {
            new_target_edge_length *= 2;
        }

        // at maximum twice the incident edge length
        {
            const auto vertices = simplex::faces_single_dimension(
                m_mesh,
                simplex::Simplex::edge(m_mesh, tuples[0]),
                PrimitiveType::Vertex);

            for (const simplex::Simplex& v : vertices) {
                auto edges =
                    simplex::cofaces_single_dimension_iterable(m_mesh, v, PrimitiveType::Edge);

                double min_el = std::numeric_limits<double>::max();
                for (const Tuple& e : edges) {
                    min_el = std::min(min_el, el_acc.const_scalar_attribute(e));
                }
                new_target_edge_length = std::min(new_target_edge_length, 3 * min_el);
            }
        }

        new_target_edge_length = std::min(new_target_edge_length, m_max_edge_length); // upper bound
        new_target_edge_length = std::max(new_target_edge_length, m_min_edge_length); // lower bound

        return Eigen::VectorXd::Constant(1, new_target_edge_length);
    };
    m_embedding_target_edge_length_transfer =
        std::make_shared<wmtk::operations::SingleAttributeTransferStrategy<double, double>>(
            m_embedding_target_edge_length_attribute,
            m_embedding_amips_attribute,
            compute_target_edge_length);

    // set target length to current length
    {
        auto target_edge_length_acc =
            m_mesh.create_accessor(m_embedding_target_edge_length_attribute.as<double>());
        const auto edge_length_acc =
            m_mesh.create_const_accessor(m_embedding_edge_length_attribute.as<double>());

        for (const Tuple& t : m_mesh.get_all(PrimitiveType::Edge)) {
            target_edge_length_acc.scalar_attribute(t) = edge_length_acc.const_scalar_attribute(t);
        }
    }
    m_embedding_target_edge_length_transfer->run_on_all();

    //////////////////////////////////
    // Lambdas for priority
    //////////////////////////////////
    m_embedding_prio_long_edges_first = [this](const simplex::Simplex& s) -> double {
        assert(
            s.primitive_type() == PrimitiveType::Edge ||
            s.primitive_type() == PrimitiveType::Triangle);

        const auto acc = m_mesh.create_const_accessor<double>(m_embedding_edge_length_attribute);

        if (s.primitive_type() == PrimitiveType::Edge) {
            return -acc.const_scalar_attribute(s.tuple());
        } else if (s.primitive_type() == PrimitiveType::Triangle) {
            // get max edge length
            auto edges = simplex::faces_single_dimension_tuples(m_mesh, s, PrimitiveType::Edge);
            double max_el = std::numeric_limits<double>::lowest();
            for (const Tuple& e : edges) {
                max_el = std::max(max_el, acc.const_scalar_attribute(e));
            }
            return -max_el;
        }
        return 0;
    };
    m_embedding_prio_short_edges_first = [this](const simplex::Simplex& s) -> double {
        assert(s.primitive_type() == PrimitiveType::Edge);
        auto acc = m_mesh.create_accessor<double>(m_embedding_edge_length_attribute);
        return acc.scalar_attribute(s.tuple());
    };

    //////////////////////////////////
    // Invariants
    m_embedding_inversion_invariant =
        std::make_shared<SimplexInversionInvariant<double>>(m_mesh, m_pos_handle.as<double>());

    m_link_conditions = std::make_shared<wmtk::invariants::InvariantCollection>(m_mesh);
    m_link_conditions->add(std::make_shared<MultiMeshLinkConditionInvariant>(m_mesh));
    // Link conditions for children should not be necessary to check
    // for (const auto& child : m_mesh.get_child_meshes()) {
    //    m_link_conditions->add(std::make_shared<MultiMeshLinkConditionInvariant>(*child));
    //}

    m_embedding_no_child_invariants =
        std::make_shared<wmtk::invariants::InvariantCollection>(m_mesh);
    for (const auto& child : m_mesh.get_child_meshes()) {
        m_embedding_no_child_invariants->add(
            std::make_shared<internal::invariants::NoChildSimplexInvariant>(m_mesh, *child));
    }

    m_embedding_separate_substructures_invariant =
        std::make_shared<wmtk::invariants::SeparateSubstructuresInvariant>(m_mesh);

    m_embedding_todo_larger = std::make_shared<TodoLargerInvariant>(
        m_mesh,
        m_embedding_edge_length_attribute.as<double>(),
        m_embedding_target_edge_length_attribute.as<double>(),
        4.0 / 3.0);

    m_embedding_todo_smaller = std::make_shared<TodoSmallerInvariant>(
        m_mesh,
        m_embedding_edge_length_attribute.as<double>(),
        m_embedding_target_edge_length_attribute.as<double>(),
        3.0 / 4.0);

    m_embedding_roi_invariant =
        std::make_shared<internal::invariants::RoiInvariant>(m_embedding_roi_attribute);

    m_embedding_valence_3_invariant =
        std::make_shared<wmtk::invariants::EdgeValenceInvariant>(m_mesh, 3);
    m_embedding_valence_4_invariant =
        std::make_shared<wmtk::invariants::EdgeValenceInvariant>(m_mesh, 4);
    m_embedding_swap44_invariant = std::make_shared<wmtk::Swap44EnergyBeforeInvariant<double>>(
        m_mesh,
        m_pos_handle.as<double>(),
        0);
    m_embedding_swap44_2_invariant = std::make_shared<wmtk::Swap44EnergyBeforeInvariant<double>>(
        m_mesh,
        m_pos_handle.as<double>(),
        1);
    m_embedding_swap32_invariant = std::make_shared<wmtk::Swap32EnergyBeforeInvariant<double>>(
        m_mesh,
        m_pos_handle.as<double>());
    m_embedding_swap23_invariant = std::make_shared<wmtk::Swap23EnergyBeforeInvariant<double>>(
        m_mesh,
        m_pos_handle.as<double>());

    //////////////////////////////////
    // Collapse strategies
    m_embedding_pos_collapse_strat = keep_substructure_strategy(m_mesh, m_pos_handle);

    m_embedding_pos_collapse_strat_1 = keep_substructure_strategy(m_mesh, m_pos_handle);
    m_embedding_pos_collapse_strat_1->set_strategy(operations::CollapseBasicStrategy::CopyOther);
    m_embedding_pos_collapse_strat_2 = keep_substructure_strategy(m_mesh, m_pos_handle);
    m_embedding_pos_collapse_strat_2->set_strategy(operations::CollapseBasicStrategy::CopyTuple);
}

void OffsetOptimization::init_offset_optimization()
{
    logger().info("Init offset optimization");

    const PrimitiveType pt_face = get_primitive_type_from_id(m_mesh.top_cell_dimension() - 1);

    //////////////////////////////////
    // Storing edge lengths
    m_offset_edge_length_attribute = m_offset_mesh->register_attribute<double>(
        "offset_optimization_offset_edge_length",
        PrimitiveType::Edge,
        1);
    m_offset_edge_length_transfer =
        std::make_shared<wmtk::operations::SingleAttributeTransferStrategy<double, double>>(
            m_offset_edge_length_attribute,
            m_offset_pos_handle,
            utils::compute_edge_length);
    m_offset_edge_length_transfer->run_on_all();

    //////////////////////////////////
    // Lambdas for priority
    //////////////////////////////////
    m_offset_prio_long_edges_first = [this](const simplex::Simplex& s) -> double {
        assert(s.primitive_type() == PrimitiveType::Edge);
        auto acc_el = m_offset_mesh->create_accessor<double>(m_offset_edge_length_attribute);
        return -acc_el.scalar_attribute(s.tuple());
    };
    m_offset_prio_short_edges_first = [this](const simplex::Simplex& s) -> double {
        assert(s.primitive_type() == PrimitiveType::Edge);
        auto acc = m_offset_mesh->create_accessor<double>(m_offset_edge_length_attribute);
        return acc.scalar_attribute(s.tuple());
    };

    //////////////////////////////////
    // Convergence attributes
    m_offset_vertex_converged_attribute = m_offset_mesh->register_attribute<int64_t>(
        "offset_optimization_vertex_converged",
        PrimitiveType::Vertex,
        1);
    m_offset_face_converged_attribute = m_offset_mesh->register_attribute<int64_t>(
        "offset_optimization_face_converged",
        pt_face,
        1);


    //////////////////////////////////
    // Mean ratio metric
    m_offset_mean_ratio_metric_attribute = m_offset_mesh->register_attribute<double>(
        "offset_optimization_mean_ratio_metric",
        pt_face,
        1);

    //////////////////////////////////
    // face distance
    m_offset_face_distance_attribute = m_offset_mesh->register_attribute<double>(
        "offset_optimization_face_distance",
        pt_face,
        1,
        std::numeric_limits<double>::max());

    m_offset_vertex_distance_attribute = m_offset_mesh->register_attribute<double>(
        "offset_optimization_vertex_distance",
        PrimitiveType::Vertex,
        1,
        std::numeric_limits<double>::max());

    //////////////////////////////////
    // face target distance
    m_offset_face_target_distance_attribute = m_offset_mesh->register_attribute<double>(
        "offset_optimization_face_target_distance",
        pt_face,
        1,
        std::numeric_limits<double>::max());

    m_offset_vertex_target_distance_attribute = m_offset_mesh->register_attribute<double>(
        "offset_optimization_vertex_target_distance",
        PrimitiveType::Vertex,
        1,
        std::numeric_limits<double>::max());

    //////////////////////////////////
    // Normal deviation
    m_offset_face_normal_deviation_attribute = m_offset_mesh->register_attribute<double>(
        "offset_optimization_normal_deviation",
        pt_face,
        1);

    //////////////////////////////////
    // Offset normal samples
    if (m_mesh.top_simplex_type() == PrimitiveType::Tetrahedron) {
        m_offset_normal_samples_attribute = m_offset_mesh->register_attribute<double>(
            "offset_optimization_normal_samples",
            pt_face,
            15);
    } else if (m_mesh.top_simplex_type() == PrimitiveType::Triangle) {
        m_offset_normal_samples_attribute = m_offset_mesh->register_attribute<double>(
            "offset_optimization_normal_samples",
            pt_face,
            8);
    }
    // offset normal samples transfer function
    if (m_mesh.top_simplex_type() == PrimitiveType::Tetrahedron) {
        auto offset_normal_samples_transfer_function_wout_convergence =
            [this](const Eigen::MatrixXd& P, const std::vector<Tuple>& tuples) -> Eigen::VectorXd {
            assert(P.rows() == 3); // rows --> attribute dimension
            assert(P.cols() == 3); // cols --> number of vertices

            using wmtk::utils::wmtk_orient3d;

            if (!m_offset_mesh->is_ccw(tuples[0])) {
                log_and_throw_error("Tuple is not CCW");
            }

            const simplex::Simplex f(*m_offset_mesh, PrimitiveType::Triangle, tuples[0]);

            ////////////////////////////////////////
            // normal samples

            const Eigen::Vector3d p0 = P.block(0, 0, 3, 1);
            const Eigen::Vector3d p1 = P.block(0, 1, 3, 1);
            const Eigen::Vector3d p2 = P.block(0, 2, 3, 1);
            const Eigen::Vector3d p_mid = (p0 + p1 + p2) / 3.;

            Eigen::Vector3d n_face = ((p1 - p0).cross(p2 - p0)).normalized();
            // check face normal orientation
            {
                Tuple tet_tuple = m_offset_mesh->map_to_parent_tuple(f);


                const auto tet_tag_acc = m_mesh.create_const_accessor<int64_t>(m_tet_tag_handle);
                // the tet tuple must be the one inside the offset (with a tag)
                if (tet_tag_acc.const_scalar_attribute(tet_tuple) == 0) {
                    tet_tuple = m_mesh.switch_tuple(tet_tuple, PrimitiveType::Tetrahedron);
                    // if (tet_tag_acc.const_scalar_attribute(tet_tuple) == 0) {
                    //     log_and_throw_error("no tagged tet incident to offset");
                    // }
                }
                // get opposite vertex
                Tuple opposite_vertex = m_mesh.switch_tuples(
                    tet_tuple,
                    {PrimitiveType::Triangle, PrimitiveType::Edge, PrimitiveType::Vertex});

                const auto tet_pos_acc = m_mesh.create_const_accessor<double>(m_pos_handle);
                const Eigen::Vector3d p_opp = tet_pos_acc.const_vector_attribute(opposite_vertex);

                const Eigen::Vector3d p_normal = p_mid + n_face;
                if (wmtk_orient3d(p0, p1, p2, p_opp) != wmtk_orient3d(p0, p1, p2, p_normal)) {
                    n_face *= -1;
                }
                // n_face points now in the direction of the input
            }

            constexpr double u = 0.1;

            const Eigen::Vector3d p0_in = (1 - u) * p0 + u * p_mid;
            const Eigen::Vector3d p1_in = (1 - u) * p1 + u * p_mid;
            const Eigen::Vector3d p2_in = (1 - u) * p2 + u * p_mid;

            auto [n_mid, sq_dist_mid, off_dist_mid] = compute_offset_normal(p_mid);
            auto [n0, sq_dist0, off_dist_0] = compute_offset_normal(p0_in);
            auto [n1, sq_dist1, off_dist_1] = compute_offset_normal(p1_in);
            auto [n2, sq_dist2, off_dist_2] = compute_offset_normal(p2_in);

            Eigen::Vector<double, 15> samples;
            samples.block(0, 0, 3, 1) = n_face;
            samples.block(3, 0, 3, 1) = n_mid;
            samples.block(6, 0, 3, 1) = n0;
            samples.block(9, 0, 3, 1) = n1;
            samples.block(12, 0, 3, 1) = n2;

            ////////////////////////////////////////
            // normal deviation
            double normal_deviation;
            {
                const double nd0 = utils::normal_angle(n_mid, n0);
                const double nd1 = utils::normal_angle(n_mid, n1);
                const double nd2 = utils::normal_angle(n_mid, n2);

                normal_deviation = std::max({nd0, nd1, nd2});

                auto nd_acc = m_offset_mesh->create_accessor(
                    m_offset_face_normal_deviation_attribute.as<double>());
                nd_acc.scalar_attribute(f) = normal_deviation;
            }

            ////////////////////////////////////////
            // mean ratio metric
            const double mrm = utils::mean_ratio_metric(p0, p1, p2);
            {
                auto mrm_acc = m_offset_mesh->create_accessor(
                    m_offset_mean_ratio_metric_attribute.as<double>());
                mrm_acc.scalar_attribute(f) = mrm;
            }

            ////////////////////////////////////////
            // face distance
            {
                auto fd_acc =
                    m_offset_mesh->create_accessor(m_offset_face_distance_attribute.as<double>());
                const double sq_dist_max = std::max({sq_dist_mid, sq_dist0, sq_dist1, sq_dist2});
                fd_acc.scalar_attribute(f) = std::sqrt(sq_dist_max);
            }
            ////////////////////////////////////////
            // face target distance
            {
                auto fd_acc = m_offset_mesh->create_accessor(
                    m_offset_face_target_distance_attribute.as<double>());
                const double dist_max =
                    std::max({off_dist_mid, off_dist_0, off_dist_1, off_dist_2});
                fd_acc.scalar_attribute(f) = dist_max;
            }
            // vertex distance
            {
                auto [vn0, vsq_dist0, voff_dist_0] = compute_offset_normal(p0);
                auto [vn1, vsq_dist1, voff_dist_1] = compute_offset_normal(p1);
                auto [vn2, vsq_dist2, voff_dist_2] = compute_offset_normal(p2);

                auto vd_target_acc =
                    m_offset_mesh->create_accessor(m_offset_vertex_distance_attribute.as<double>());
                vd_target_acc.scalar_attribute(tuples[0]) = std::sqrt(vsq_dist0);
                vd_target_acc.scalar_attribute(tuples[1]) = std::sqrt(vsq_dist1);
                vd_target_acc.scalar_attribute(tuples[2]) = std::sqrt(vsq_dist2);

                auto vd_acc = m_offset_mesh->create_accessor(
                    m_offset_vertex_target_distance_attribute.as<double>());
                vd_acc.scalar_attribute(tuples[0]) = voff_dist_0;
                vd_acc.scalar_attribute(tuples[1]) = voff_dist_1;
                vd_acc.scalar_attribute(tuples[2]) = voff_dist_2;
            }

            return samples;
        };

        auto offset_normal_samples_transfer_function =
            [this](const Eigen::MatrixXd& P, const std::vector<Tuple>& tuples) -> Eigen::VectorXd {
            assert(P.rows() == 3); // rows --> attribute dimension
            assert(P.cols() == 3); // cols --> number of vertices

            using wmtk::utils::wmtk_orient3d;

            if (!m_offset_mesh->is_ccw(tuples[0])) {
                log_and_throw_error("Tuple is not CCW");
            }

            const simplex::Simplex f(*m_offset_mesh, PrimitiveType::Triangle, tuples[0]);

            ////////////////////////////////////////
            // normal samples

            Eigen::Vector3d p0 = P.block(0, 0, 3, 1);
            Eigen::Vector3d p1 = P.block(0, 1, 3, 1);
            Eigen::Vector3d p2 = P.block(0, 2, 3, 1);
            const Eigen::Vector3d p_mid = (p0 + p1 + p2) / 3.;

            Eigen::Vector3d n_face = ((p1 - p0).cross(p2 - p0)).normalized();
            // check face normal orientation
            {
                Tuple tet_tuple = m_offset_mesh->map_to_parent_tuple(f);


                const auto tet_tag_acc = m_mesh.create_const_accessor<int64_t>(m_tet_tag_handle);
                // the tet tuple must be the one inside the offset (with a tag)
                if (tet_tag_acc.const_scalar_attribute(tet_tuple) == 0) {
                    tet_tuple = m_mesh.switch_tuple(tet_tuple, PrimitiveType::Tetrahedron);
                    // if (tet_tag_acc.const_scalar_attribute(tet_tuple) == 0) {
                    //     log_and_throw_error("no tagged tet incident to offset");
                    // }
                }
                // get opposite vertex
                Tuple opposite_vertex = m_mesh.switch_tuples(
                    tet_tuple,
                    {PrimitiveType::Triangle, PrimitiveType::Edge, PrimitiveType::Vertex});

                const auto tet_pos_acc = m_mesh.create_const_accessor<double>(m_pos_handle);
                const Eigen::Vector3d p_opp = tet_pos_acc.const_vector_attribute(opposite_vertex);

                const Eigen::Vector3d p_normal = p_mid + n_face;
                if (wmtk_orient3d(p0, p1, p2, p_opp) != wmtk_orient3d(p0, p1, p2, p_normal)) {
                    n_face *= -1;
                }
                // n_face points now in the direction of the input
            }

            constexpr double u = 0.1;

            const Eigen::Vector3d p0_in = (1 - u) * p0 + u * p_mid;
            const Eigen::Vector3d p1_in = (1 - u) * p1 + u * p_mid;
            const Eigen::Vector3d p2_in = (1 - u) * p2 + u * p_mid;

            auto [n_mid, sq_dist_mid, off_dist_mid] = compute_offset_normal(p_mid);
            auto [n0, sq_dist0, off_dist_0] = compute_offset_normal(p0_in);
            auto [n1, sq_dist1, off_dist_1] = compute_offset_normal(p1_in);
            auto [n2, sq_dist2, off_dist_2] = compute_offset_normal(p2_in);

            Eigen::Vector<double, 15> samples;
            samples.block(0, 0, 3, 1) = n_face;
            samples.block(3, 0, 3, 1) = n_mid;
            samples.block(6, 0, 3, 1) = n0;
            samples.block(9, 0, 3, 1) = n1;
            samples.block(12, 0, 3, 1) = n2;

            ////////////////////////////////////////
            // normal deviation
            double normal_deviation;
            {
                const double nd0 = utils::normal_angle(n_mid, n0);
                const double nd1 = utils::normal_angle(n_mid, n1);
                const double nd2 = utils::normal_angle(n_mid, n2);

                normal_deviation = std::max({nd0, nd1, nd2});

                auto nd_acc = m_offset_mesh->create_accessor(
                    m_offset_face_normal_deviation_attribute.as<double>());
                nd_acc.scalar_attribute(f) = normal_deviation;
            }

            ////////////////////////////////////////
            // mean ratio metric
            const double mrm = utils::mean_ratio_metric(p0, p1, p2);
            {
                auto mrm_acc = m_offset_mesh->create_accessor(
                    m_offset_mean_ratio_metric_attribute.as<double>());
                mrm_acc.scalar_attribute(f) = mrm;
            }

            ////////////////////////////////////////
            // face distance
            {
                auto fd_acc =
                    m_offset_mesh->create_accessor(m_offset_face_distance_attribute.as<double>());
                const double sq_dist_max = std::max({sq_dist_mid, sq_dist0, sq_dist1, sq_dist2});
                fd_acc.scalar_attribute(f) = std::sqrt(sq_dist_max);
            }
            ////////////////////////////////////////
            // face target distance
            {
                auto fd_acc = m_offset_mesh->create_accessor(
                    m_offset_face_target_distance_attribute.as<double>());
                const double dist_max =
                    std::max({off_dist_mid, off_dist_0, off_dist_1, off_dist_2});
                fd_acc.scalar_attribute(f) = dist_max;
            }
            // vertex distance
            {
                auto [vn0, vsq_dist0, voff_dist_0] = compute_offset_normal(p0);
                auto [vn1, vsq_dist1, voff_dist_1] = compute_offset_normal(p1);
                auto [vn2, vsq_dist2, voff_dist_2] = compute_offset_normal(p2);

                auto vd_target_acc =
                    m_offset_mesh->create_accessor(m_offset_vertex_distance_attribute.as<double>());
                vd_target_acc.scalar_attribute(tuples[0]) = std::sqrt(vsq_dist0);
                vd_target_acc.scalar_attribute(tuples[1]) = std::sqrt(vsq_dist1);
                vd_target_acc.scalar_attribute(tuples[2]) = std::sqrt(vsq_dist2);

                auto vd_acc = m_offset_mesh->create_accessor(
                    m_offset_vertex_target_distance_attribute.as<double>());
                vd_acc.scalar_attribute(tuples[0]) = voff_dist_0;
                vd_acc.scalar_attribute(tuples[1]) = voff_dist_1;
                vd_acc.scalar_attribute(tuples[2]) = voff_dist_2;
            }

            ///////////////////////////////////////
            // convergence
            {
                auto face_converged_acc =
                    m_offset_mesh->create_accessor(m_offset_face_converged_attribute.as<int64_t>());

                int64_t is_converged = 0;

                if (normal_deviation < m_max_normal_deviation && mrm > 0.5) {
                    const auto vertex_converged_acc = m_offset_mesh->create_const_accessor(
                        m_offset_vertex_converged_attribute.as<int64_t>());
                    const bool v0_converged =
                        (vertex_converged_acc.const_scalar_attribute(tuples[0]) == 1);
                    const bool v1_converged =
                        (vertex_converged_acc.const_scalar_attribute(tuples[1]) == 1);
                    const bool v2_converged =
                        (vertex_converged_acc.const_scalar_attribute(tuples[2]) == 1);

                    if (v0_converged && v1_converged && v2_converged) {
                        is_converged = 1;
                    }
                }

                face_converged_acc.scalar_attribute(f) = is_converged;
            }

            return samples;
        };

        m_offset_point_to_face_transfer =
            std::make_shared<wmtk::operations::SingleAttributeTransferStrategy<double, double>>(
                m_offset_normal_samples_attribute,
                m_offset_pos_handle,
                offset_normal_samples_transfer_function);

        m_offset_point_to_face_transfer_wout_convergence =
            std::make_shared<wmtk::operations::SingleAttributeTransferStrategy<double, double>>(
                m_offset_normal_samples_attribute,
                m_offset_pos_handle,
                offset_normal_samples_transfer_function_wout_convergence);
    } else if (m_mesh.top_simplex_type() == PrimitiveType::Triangle) {
        auto offset_normal_samples_transfer_function_wout_convergence =
            [this](const Eigen::MatrixXd& P, const std::vector<Tuple>& tuples) -> Eigen::VectorXd {
            assert(P.rows() == 2); // rows --> attribute dimension
            assert(P.cols() == 2); // cols --> number of vertices

            if (!m_offset_mesh->is_ccw(tuples[0])) {
                log_and_throw_error("Tuple is not CCW");
            }

            const simplex::Simplex f(*m_offset_mesh, PrimitiveType::Edge, tuples[0]);

            ////////////////////////////////////////
            // normal samples

            Eigen::Vector2d p0 = P.block(0, 0, 2, 1);
            Eigen::Vector2d p1 = P.block(0, 1, 2, 1);
            const Eigen::Vector2d p_mid = (p0 + p1) / 2.;

            const Eigen::Vector2d edge_vector_normalized = (p1 - p0).normalized();
            Eigen::Vector2d n_face(-edge_vector_normalized[1], edge_vector_normalized[0]);

            constexpr double u = 0.1;

            const Eigen::Vector2d p0_in = (1 - u) * p0 + u * p_mid;
            const Eigen::Vector2d p1_in = (1 - u) * p1 + u * p_mid;

            auto [n_mid, sq_dist_mid, off_dist_mid] = compute_offset_normal(p_mid);
            auto [n0, sq_dist0, off_dist_0] = compute_offset_normal(p0_in);
            auto [n1, sq_dist1, off_dist_1] = compute_offset_normal(p1_in);

            Eigen::Vector<double, 8> samples;
            samples.block(0, 0, 2, 1) = n_face;
            samples.block(2, 0, 2, 1) = n_mid;
            samples.block(4, 0, 2, 1) = n0;
            samples.block(6, 0, 2, 1) = n1;

            ////////////////////////////////////////
            // normal deviation
            double normal_deviation;
            {
                const double nd0 = utils::normal_angle(n_mid, n0);
                const double nd1 = utils::normal_angle(n_mid, n1);

                normal_deviation = std::max(nd0, nd1);

                auto nd_acc = m_offset_mesh->create_accessor(
                    m_offset_face_normal_deviation_attribute.as<double>());
                nd_acc.scalar_attribute(f) = normal_deviation;
            }

            ////////////////////////////////////////
            // mean ratio metric
            {
                // this metric has no meaning for edges
                auto mrm_acc = m_offset_mesh->create_accessor(
                    m_offset_mean_ratio_metric_attribute.as<double>());
                mrm_acc.scalar_attribute(f) = 1;
            }

            ////////////////////////////////////////
            // face distance
            {
                auto fd_acc =
                    m_offset_mesh->create_accessor(m_offset_face_distance_attribute.as<double>());
                const double sq_dist_max = std::max({sq_dist_mid, sq_dist0, sq_dist1});
                fd_acc.scalar_attribute(f) = std::sqrt(sq_dist_max);
            }
            ////////////////////////////////////////
            // face target distance
            {
                auto fd_acc = m_offset_mesh->create_accessor(
                    m_offset_face_target_distance_attribute.as<double>());
                const double dist_max = std::max({off_dist_mid, off_dist_0, off_dist_1});
                fd_acc.scalar_attribute(f) = dist_max;
            }
            // vertex distance
            {
                auto [vn0, vsq_dist0, voff_dist_0] = compute_offset_normal(p0);
                auto [vn1, vsq_dist1, voff_dist_1] = compute_offset_normal(p1);

                auto vd_target_acc =
                    m_offset_mesh->create_accessor(m_offset_vertex_distance_attribute.as<double>());
                vd_target_acc.scalar_attribute(tuples[0]) = std::sqrt(vsq_dist0);
                vd_target_acc.scalar_attribute(tuples[1]) = std::sqrt(vsq_dist1);

                auto vd_acc = m_offset_mesh->create_accessor(
                    m_offset_vertex_target_distance_attribute.as<double>());
                vd_acc.scalar_attribute(tuples[0]) = voff_dist_0;
                vd_acc.scalar_attribute(tuples[1]) = voff_dist_1;
            }

            return samples;
        };

        auto offset_normal_samples_transfer_function =
            [this](const Eigen::MatrixXd& P, const std::vector<Tuple>& tuples) -> Eigen::VectorXd {
            assert(P.rows() == 2); // rows --> attribute dimension
            assert(P.cols() == 2); // cols --> number of vertices

            if (!m_offset_mesh->is_ccw(tuples[0])) {
                log_and_throw_error("Tuple is not CCW");
            }

            const simplex::Simplex f(*m_offset_mesh, PrimitiveType::Edge, tuples[0]);

            ////////////////////////////////////////
            // normal samples

            Eigen::Vector2d p0 = P.block(0, 0, 2, 1);
            Eigen::Vector2d p1 = P.block(0, 1, 2, 1);
            const Eigen::Vector2d p_mid = (p0 + p1) / 2.;

            const Eigen::Vector2d edge_vector_normalized = (p1 - p0).normalized();
            Eigen::Vector2d n_face(-edge_vector_normalized[1], edge_vector_normalized[0]);

            constexpr double u = 0.1;

            const Eigen::Vector2d p0_in = (1 - u) * p0 + u * p_mid;
            const Eigen::Vector2d p1_in = (1 - u) * p1 + u * p_mid;

            auto [n_mid, sq_dist_mid, off_dist_mid] = compute_offset_normal(p_mid);
            auto [n0, sq_dist0, off_dist_0] = compute_offset_normal(p0_in);
            auto [n1, sq_dist1, off_dist_1] = compute_offset_normal(p1_in);

            Eigen::Vector<double, 8> samples;
            samples.block(0, 0, 2, 1) = n_face;
            samples.block(2, 0, 2, 1) = n_mid;
            samples.block(4, 0, 2, 1) = n0;
            samples.block(6, 0, 2, 1) = n1;

            ////////////////////////////////////////
            // normal deviation
            double normal_deviation;
            {
                const double nd0 = utils::normal_angle(n_mid, n0);
                const double nd1 = utils::normal_angle(n_mid, n1);

                normal_deviation = std::max(nd0, nd1);

                auto nd_acc = m_offset_mesh->create_accessor(
                    m_offset_face_normal_deviation_attribute.as<double>());
                nd_acc.scalar_attribute(f) = normal_deviation;
            }

            ////////////////////////////////////////
            // mean ratio metric
            {
                // this metric has no meaning for edges
                auto mrm_acc = m_offset_mesh->create_accessor(
                    m_offset_mean_ratio_metric_attribute.as<double>());
                mrm_acc.scalar_attribute(f) = 1;
            }

            ////////////////////////////////////////
            // face distance
            {
                auto fd_acc =
                    m_offset_mesh->create_accessor(m_offset_face_distance_attribute.as<double>());
                const double sq_dist_max = std::max({sq_dist_mid, sq_dist0, sq_dist1});
                fd_acc.scalar_attribute(f) = std::sqrt(sq_dist_max);
            }
            ////////////////////////////////////////
            // face target distance
            {
                auto fd_acc = m_offset_mesh->create_accessor(
                    m_offset_face_target_distance_attribute.as<double>());
                const double dist_max = std::max({off_dist_mid, off_dist_0, off_dist_1});
                fd_acc.scalar_attribute(f) = dist_max;
            }
            // vertex distance
            {
                auto [vn0, vsq_dist0, voff_dist_0] = compute_offset_normal(p0);
                auto [vn1, vsq_dist1, voff_dist_1] = compute_offset_normal(p1);

                auto vd_target_acc =
                    m_offset_mesh->create_accessor(m_offset_vertex_distance_attribute.as<double>());
                vd_target_acc.scalar_attribute(tuples[0]) = std::sqrt(vsq_dist0);
                vd_target_acc.scalar_attribute(tuples[1]) = std::sqrt(vsq_dist1);

                auto vd_acc = m_offset_mesh->create_accessor(
                    m_offset_vertex_target_distance_attribute.as<double>());
                vd_acc.scalar_attribute(tuples[0]) = voff_dist_0;
                vd_acc.scalar_attribute(tuples[1]) = voff_dist_1;
            }

            ///////////////////////////////////////
            // convergence
            {
                auto face_converged_acc =
                    m_offset_mesh->create_accessor(m_offset_face_converged_attribute.as<int64_t>());

                int64_t is_converged = 0;

                if (normal_deviation < m_max_normal_deviation) {
                    const auto vertex_converged_acc = m_offset_mesh->create_const_accessor(
                        m_offset_vertex_converged_attribute.as<int64_t>());
                    const bool v0_converged =
                        (vertex_converged_acc.const_scalar_attribute(tuples[0]) == 1);
                    const bool v1_converged =
                        (vertex_converged_acc.const_scalar_attribute(tuples[1]) == 1);

                    if (v0_converged && v1_converged) {
                        is_converged = 1;
                    }
                }

                face_converged_acc.scalar_attribute(f) = is_converged;
            }

            return samples;
        };

        m_offset_point_to_face_transfer =
            std::make_shared<wmtk::operations::SingleAttributeTransferStrategy<double, double>>(
                m_offset_normal_samples_attribute,
                m_offset_pos_handle,
                offset_normal_samples_transfer_function);

        m_offset_point_to_face_transfer_wout_convergence =
            std::make_shared<wmtk::operations::SingleAttributeTransferStrategy<double, double>>(
                m_offset_normal_samples_attribute,
                m_offset_pos_handle,
                offset_normal_samples_transfer_function_wout_convergence);
    }
    m_offset_point_to_face_transfer->run_on_all();

    // position transfer to embedding
    auto propagate_double = [](const Eigen::MatrixXd& P) -> Eigen::VectorXd { return P; };

    m_offset_position_to_embedding_transfer =
        std::make_shared<operations::SingleAttributeTransferStrategy<double, double>>(
            m_pos_handle,
            m_offset_pos_handle,
            propagate_double);

    //////////////////////////////////
    // Target edge length
    m_offset_target_edge_length_attribute = m_offset_mesh->register_attribute<double>(
        "offset_optimization_offset_target_edge_length",
        PrimitiveType::Edge,
        1,
        false,
        m_max_edge_length); // defaults to max edge length

    auto compute_target_edge_length =
        [this](const Eigen::MatrixXd& P, const std::vector<Tuple>& tuples) -> Eigen::VectorXd {
        assert(P.rows() == 1); // rows --> attribute dimension
        assert(P.cols() == 1 || P.cols() == 2);
        assert(tuples.size() == 1 || tuples.size() == 2);

        auto target_edge_length_accessor =
            m_offset_mesh->create_accessor(m_offset_target_edge_length_attribute.as<double>());

        const auto mrm_acc =
            m_offset_mesh->create_const_accessor(m_offset_mean_ratio_metric_attribute.as<double>());

        double mrm_min = std::numeric_limits<double>::max();
        for (const Tuple& t : tuples) {
            mrm_min = std::min(mrm_min, mrm_acc.const_scalar_attribute(t));
        }

        const double current_target_edge_length =
            target_edge_length_accessor.const_scalar_attribute(tuples[0]);

        const double normal_deviation = P.maxCoeff();

        // get all vertices and compute mrm

        double new_target_edge_length = current_target_edge_length;
        if (normal_deviation > m_max_normal_deviation || mrm_min < 0.5) {
            new_target_edge_length *= 0.5;
        } else if (normal_deviation < m_min_normal_deviation && mrm_min > 0.5) {
            new_target_edge_length *= 1.5;
        }

        // at maximum twice the incident target edge length
        {
            const auto vertices = simplex::faces_single_dimension(
                *m_offset_mesh,
                simplex::Simplex::edge(*m_offset_mesh, tuples[0]),
                PrimitiveType::Vertex);

            for (const simplex::Simplex& v : vertices) {
                auto edges = simplex::cofaces_single_dimension_iterable(
                    *m_offset_mesh,
                    v,
                    PrimitiveType::Edge);

                double min_tel = std::numeric_limits<double>::max();
                for (const Tuple& e : edges) {
                    min_tel =
                        std::min(min_tel, target_edge_length_accessor.const_scalar_attribute(e));
                }
                new_target_edge_length = std::min(new_target_edge_length, 1.5 * min_tel);
            }
        }

        new_target_edge_length = std::min(new_target_edge_length, m_max_edge_length); // upper bound
        new_target_edge_length = std::max(new_target_edge_length, m_min_edge_length); // lower bound

        return Eigen::VectorXd::Constant(1, new_target_edge_length);

        // return Eigen::VectorXd::Constant(1, m_max_edge_length);
    };
    m_offset_target_edge_length_transfer =
        std::make_shared<wmtk::operations::SingleAttributeTransferStrategy<double, double>>(
            m_offset_target_edge_length_attribute,
            m_offset_face_normal_deviation_attribute,
            compute_target_edge_length);

    // set target length to current length
    {
        auto target_edge_length_acc =
            m_offset_mesh->create_accessor(m_offset_target_edge_length_attribute.as<double>());
        const auto edge_length_acc =
            m_offset_mesh->create_const_accessor(m_offset_edge_length_attribute.as<double>());

        for (const Tuple& t : m_offset_mesh->get_all(PrimitiveType::Edge)) {
            target_edge_length_acc.scalar_attribute(t) = edge_length_acc.const_scalar_attribute(t);
        }
    }

    m_offset_target_edge_length_transfer->run_on_all();

    //////////////////////////////////
    // Invariants
    //////////////////////////////////
    m_offset_todo_larger = std::make_shared<TodoLargerInvariant>(
        *m_offset_mesh,
        m_offset_edge_length_attribute.as<double>(),
        m_offset_target_edge_length_attribute.as<double>(),
        4.0 / 3.0);

    m_offset_todo_smaller = std::make_shared<TodoSmallerInvariant>(
        *m_offset_mesh,
        m_offset_edge_length_attribute.as<double>(),
        m_offset_target_edge_length_attribute.as<double>(),
        3.0 / 4.0);

    if (m_mesh.top_simplex_type() == PrimitiveType::Tetrahedron) {
        m_offset_swap_invariant = std::make_shared<invariants::OffsetSwapInvariant>(
            m_offset_normal_samples_attribute,
            m_offset_pos_handle,
            m_max_normal_deviation);
    }


    m_normal_deviation_invariant = std::make_shared<invariants::NormalDeviationAfterInvariant>(
        m_offset_face_normal_deviation_attribute,
        m_max_normal_deviation,
        false);

    m_offset_convergence_invariant =
        std::make_shared<invariants::ConvergenceInvariant>(m_offset_face_converged_attribute);

    //////////////////////////////////
    // Smoothing function
    if (m_mesh.top_simplex_type() == PrimitiveType::Tetrahedron) {
        m_offset_smoothing_function = [this](Mesh& m, const simplex::Simplex& s) -> bool {
            assert(m == *m_offset_mesh);
            assert(s.primitive_type() == PrimitiveType::Vertex);

            auto pos_offset_acc = m.create_accessor<double>(m_offset_pos_handle);
            const Eigen::Vector3d p0 = pos_offset_acc.const_vector_attribute(s);

            // laplacian smoothing
            Eigen::Vector3d p_laplace = Eigen::Vector3d::Zero();
            {
                auto neighs = simplex::link_single_dimension_iterable(m, s, PrimitiveType::Vertex);
                int64_t n_neighs = 0;
                for (const Tuple& neigh : neighs) {
                    ++n_neighs;
                    p_laplace += pos_offset_acc.const_vector_attribute(neigh);
                }
                p_laplace /= n_neighs;
            }

            // project with quadrics
            Eigen::Vector3d p_optimal;
            {
                utils::Quadrics q(
                    *m_input_bvh,
                    m_offset_pos_handle,
                    m_tet_tag_handle,
                    m_input_pos_handle,
                    s);
                p_optimal = q.solve(p_laplace);
            }

            {
                double w = 0.5;
                p_optimal = (1 - w) * p0 + w * p_optimal;
            }

            auto vertex_converged_acc =
                m_offset_mesh->create_accessor(m_offset_vertex_converged_attribute.as<int64_t>());

            // bisection search
            {
                auto pos_emb_acc = m_mesh.create_accessor<double>(m_pos_handle);
                const simplex::Simplex s_emb = m.map_to_parent(s);

                const auto tet_tuples = simplex::top_dimension_cofaces_tuples(m_mesh, s_emb);
                auto p_emb = pos_emb_acc.vector_attribute(s_emb);
                assert((p_emb - p0).squaredNorm() == 0);

                p_emb = p_optimal;
                if (m_embedding_inversion_invariant->after({}, tet_tuples)) {
                    pos_offset_acc.vector_attribute(s) = p_emb;
                    // vertex_converged_acc.scalar_attribute(s) = 1;
                    return true;
                }
                vertex_converged_acc.scalar_attribute(s) = 0; // result is not optimal

                //// bisection search from laplacian
                // for (size_t i = 0; i < 10; ++i) {
                //    const double u = 1.0 / (2 << i);
                //    p_emb = (1 - u) * p_laplace + u * p_optimal;
                //
                //    if (m_embedding_inversion_invariant->after({}, tet_tuples)) {
                //        pos_offset_acc.vector_attribute(s.tuple()) = p_emb;
                //        return true;
                //    }
                //}

                // bisection search from p0
                for (size_t i = 0; i < 10; ++i) {
                    const double u = 1.0 / (2 << i);
                    p_emb = (1 - u) * p0 + u * p_optimal;

                    if (m_embedding_inversion_invariant->after({}, tet_tuples)) {
                        pos_offset_acc.vector_attribute(s) = p_emb;
                        return true;
                    }
                }

                p_emb = p0;
            }

            // Vertex position is not optimal but it cannot be moved either --> it is converged.
            // vertex_converged_acc.scalar_attribute(s) = 1;
            return false;
        };
    } else if (m_mesh.top_simplex_type() == PrimitiveType::Triangle) {
        m_offset_smoothing_function = [this](Mesh& m, const simplex::Simplex& s) -> bool {
            assert(&m == &*m_offset_mesh);
            assert(s.primitive_type() == PrimitiveType::Vertex);

            auto pos_offset_acc = m.create_accessor<double>(m_offset_pos_handle);
            const Eigen::Vector2d p0 = pos_offset_acc.const_vector_attribute(s);

            // laplacian smoothing
            Eigen::Vector2d p_laplace = Eigen::Vector2d::Zero();
            {
                auto neighs = simplex::link_single_dimension_iterable(m, s, PrimitiveType::Vertex);
                int64_t n_neighs = 0;
                for (const Tuple& neigh : neighs) {
                    ++n_neighs;
                    p_laplace += pos_offset_acc.const_vector_attribute(neigh);
                }
                p_laplace /= n_neighs;
            }

            // project with quadrics
            Eigen::Vector2d p_optimal;
            {
                utils::Quadrics q(
                    *m_input_bvh,
                    m_offset_pos_handle,
                    m_tet_tag_handle,
                    m_input_pos_handle,
                    s);
                const Eigen::Vector3d p_laplace_3(p_laplace[0], p_laplace[1], 0);
                const Eigen::Vector3d p_optimal_3 = q.solve(p_laplace_3);
                p_optimal = p_optimal_3.block(0, 0, 2, 1).eval();
            }

            {
                double w = 0.5;
                p_optimal = (1 - w) * p0 + w * p_optimal;
            }

            auto vertex_converged_acc =
                m_offset_mesh->create_accessor(m_offset_vertex_converged_attribute.as<int64_t>());

            // bisection search
            {
                auto pos_emb_acc = m_mesh.create_accessor<double>(m_pos_handle);
                const simplex::Simplex s_emb = m.map_to_parent(s);

                const auto tet_tuples = simplex::top_dimension_cofaces_tuples(m_mesh, s_emb);
                auto p_emb = pos_emb_acc.vector_attribute(s_emb);
                assert((p_emb - p0).squaredNorm() == 0);

                p_emb = p_optimal;
                if (m_embedding_inversion_invariant->after({}, tet_tuples)) {
                    pos_offset_acc.vector_attribute(s) = p_emb;
                    // vertex_converged_acc.scalar_attribute(s) = 1;
                    return true;
                }
                vertex_converged_acc.scalar_attribute(s) = 0; // result is not optimal

                // bisection search from p0
                for (size_t i = 0; i < 10; ++i) {
                    const double u = 1.0 / (2 << i);
                    p_emb = (1 - u) * p0 + u * p_optimal;

                    if (m_embedding_inversion_invariant->after({}, tet_tuples)) {
                        pos_offset_acc.vector_attribute(s) = p_emb;
                        return true;
                    }
                }

                p_emb = p0;
            }

            // Vertex position is not optimal but it cannot be moved either --> it is converged.
            // vertex_converged_acc.scalar_attribute(s) = 1;
            return false;
        };
    }


    //////////////////////////////////
    // Collapse strategies
    m_offset_pos_collapse_strat_1 =
        std::make_shared<operations::CollapseNewAttributeStrategy<double>>(m_offset_pos_handle);
    m_offset_pos_collapse_strat_1->set_strategy(operations::CollapseBasicStrategy::CopyOther);
    m_offset_pos_collapse_strat_2 =
        std::make_shared<operations::CollapseNewAttributeStrategy<double>>(m_offset_pos_handle);
    m_offset_pos_collapse_strat_2->set_strategy(operations::CollapseBasicStrategy::CopyTuple);

    m_offset_vertex_converged_collapse_strat_1 =
        std::make_shared<operations::CollapseNewAttributeStrategy<int64_t>>(
            m_offset_vertex_converged_attribute);
    m_offset_vertex_converged_collapse_strat_1->set_strategy(
        operations::CollapseBasicStrategy::CopyOther);
    m_offset_vertex_converged_collapse_strat_2 =
        std::make_shared<operations::CollapseNewAttributeStrategy<int64_t>>(
            m_offset_vertex_converged_attribute);
    m_offset_vertex_converged_collapse_strat_2->set_strategy(
        operations::CollapseBasicStrategy::CopyTuple);

    // report input metrics
    m_metrics_report["input"]["num_triangles"] =
        m_input_mesh->get_all(PrimitiveType::Triangle).size();
    m_metrics_report["input"]["num_vertices"] = m_input_mesh->get_all(PrimitiveType::Vertex).size();
    m_metrics_report["input"]["offset_distance"] = m_offset_distance;
    m_metrics_report["input"]["min_edge_length"] = m_min_edge_length;
    m_metrics_report["input"]["max_edge_length"] = m_max_edge_length;
    m_metrics_report["input"]["min_normal_deviation"] = m_min_normal_deviation;
    m_metrics_report["input"]["max_normal_deviation"] = m_max_normal_deviation;
}

void OffsetOptimization::optimize_embedding(const int64_t n_iterations)
{
    std::vector<attribute::MeshAttributeHandle> none_attributes(
        {m_embedding_edge_length_attribute, m_embedding_amips_attribute});

    auto add_transfers = [this](operations::Operation& op) {
        op.add_transfer_strategy(m_embedding_edge_length_transfer);
        op.add_transfer_strategy(m_embedding_amips_transfer);
        // op.add_transfer_strategy(m_embedding_target_edge_length_transfer);
    };

    auto tel_update = [this]() {
        // set all TEL to EL
        auto tel_acc =
            m_mesh.create_accessor(m_embedding_target_edge_length_attribute.as<double>());
        const auto el_acc =
            m_mesh.create_const_accessor(m_embedding_edge_length_attribute.as<double>());

        for (const Tuple& t : m_mesh.get_all(PrimitiveType::Edge)) {
            tel_acc.scalar_attribute(t) = el_acc.const_scalar_attribute(t);
        }

        const auto amips_acc =
            m_mesh.create_const_accessor(m_embedding_amips_attribute.as<double>());

        for (const Tuple& t : m_mesh.get_all(PrimitiveType::Tetrahedron)) {
            if (amips_acc.const_scalar_attribute(t) < 100) {
                continue;
            }
            const auto edges = simplex::faces_single_dimension(
                m_mesh,
                simplex::Simplex::tetrahedron(m_mesh, t),
                PrimitiveType::Edge);

            const auto vertices = simplex::faces_single_dimension(
                m_mesh,
                simplex::Simplex::tetrahedron(m_mesh, t),
                PrimitiveType::Vertex);

            double el_min = std::numeric_limits<double>::max();
            double el_max = std::numeric_limits<double>::lowest();

            for (const simplex::Simplex& e : edges) {
                const double el = el_acc.const_scalar_attribute(e);
                el_min = std::min(el_min, el);
                el_max = std::max(el_max, el);
            }

            const double el_med = 0.5 * (el_max - el_min);
            for (const simplex::Simplex& v : vertices) {
                auto incident_edges =
                    simplex::cofaces_single_dimension_iterable(m_mesh, v, PrimitiveType::Edge);
                for (const Tuple& e : incident_edges) {
                    tel_acc.scalar_attribute(e) = el_med;
                }
            }
        }
    };


    //////////////////////////////////
    // 1) EdgeSplit
    //////////////////////////////////
    auto split = std::make_shared<operations::EdgeSplit>(m_mesh);
    split->set_priority(m_embedding_prio_long_edges_first);

    split->add_invariant(m_embedding_todo_larger);
    split->add_invariant(m_embedding_roi_invariant);
    split->add_invariant(m_embedding_no_child_invariants);
    split->add_invariant(m_embedding_inversion_invariant); // we could reach floating point
                                                           // precision so it's better to have that

    split->set_new_attribute_strategy(m_pos_handle);
    split->set_new_attribute_strategy(
        m_embedding_target_edge_length_attribute,
        operations::SplitBasicStrategy::Copy,
        operations::SplitRibBasicStrategy::Mean);
    split->set_new_attribute_strategy(m_tet_tag_handle);
    split->set_new_attribute_strategy(m_embedding_roi_attribute);
    for (const auto& attr : none_attributes) {
        split->set_new_attribute_strategy(
            attr,
            operations::SplitBasicStrategy::None,
            operations::SplitRibBasicStrategy::None);
    }
    add_transfers(*split);

    //////////////////////////////////
    // 2) EdgeCollapse
    //////////////////////////////////

    // auto setup_halfedge_collapse =
    //     [this, &none_attributes, &add_transfers](
    //         const int64_t collapse_type) -> std::shared_ptr<operations::EdgeCollapse> {
    //     auto _collapse = std::make_shared<operations::EdgeCollapse>(m_mesh);
    //     _collapse->set_priority(m_embedding_prio_short_edges_first);
    //
    //    // use different before energy depending on the type
    //    _collapse->add_invariant(
    //        std::make_shared<wmtk::invariants::CollapseEnergyBeforeInvariant<double>>(
    //            m_mesh,
    //            m_pos_handle.as<double>(),
    //            m_embedding_amips_attribute.as<double>(),
    //            collapse_type - 1));
    //
    //    _collapse->add_invariant(m_embedding_inversion_invariant);
    //
    //    // use different strategies for position depending on the type
    //    switch (collapse_type) {
    //    case 1:
    //        _collapse->set_new_attribute_strategy(m_pos_handle, m_embedding_pos_collapse_strat_1);
    //        break;
    //    case 2:
    //        _collapse->set_new_attribute_strategy(m_pos_handle, m_embedding_pos_collapse_strat_2);
    //        break;
    //    default: log_and_throw_error("Unknown halfedge collapse type: {}", collapse_type); break;
    //    }
    //    _collapse->set_new_attribute_strategy(
    //        m_embedding_target_edge_length_attribute,
    //        operations::CollapseBasicStrategy::Mean);
    //    _collapse->set_new_attribute_strategy(m_tet_tag_handle);
    //    for (const auto& attr : none_attributes) {
    //        _collapse->set_new_attribute_strategy(attr, operations::CollapseBasicStrategy::None);
    //    }
    //    add_transfers(*_collapse);
    //
    //    return _collapse;
    //};
    //
    // auto collapse_1 = setup_halfedge_collapse(1);
    // auto collapse_2 = setup_halfedge_collapse(2);
    // auto collapse = std::make_shared<operations::OrOperationSequence>(m_mesh);
    // collapse->add_operation(collapse_1);
    // collapse->add_operation(collapse_2);
    // collapse->add_invariant(m_embedding_todo_smaller);
    // collapse->add_invariant(m_embedding_no_child_invariants);
    // collapse->add_invariant(m_embedding_separate_substructures_invariant);
    // collapse->add_invariant(m_link_conditions);

    auto collapse = std::make_shared<operations::EdgeCollapse>(m_mesh);
    collapse->set_priority(m_embedding_prio_short_edges_first);

    collapse->add_invariant(m_embedding_todo_smaller);
    collapse->add_invariant(m_embedding_roi_invariant);
    collapse->add_invariant(m_embedding_no_child_invariants);
    collapse->add_invariant(m_embedding_separate_substructures_invariant);
    collapse->add_invariant(m_link_conditions);
    collapse->add_invariant(m_embedding_inversion_invariant);
    collapse->add_invariant(std::make_shared<wmtk::invariants::MaxFunctionInvariant>(
        PrimitiveType::Tetrahedron,
        m_embedding_amips));

    collapse->set_new_attribute_strategy(m_pos_handle, m_embedding_pos_collapse_strat);
    collapse->set_new_attribute_strategy(m_embedding_target_edge_length_attribute);
    collapse->set_new_attribute_strategy(m_tet_tag_handle);
    collapse->set_new_attribute_strategy(m_embedding_roi_attribute);
    for (const auto& attr : none_attributes) {
        collapse->set_new_attribute_strategy(attr, operations::CollapseBasicStrategy::None);
    }
    add_transfers(*collapse);

    //////////////////////////////////
    // 3) Swap
    //////////////////////////////////
    auto setup_swap = [this, &none_attributes, &add_transfers](
                          operations::Operation& op,
                          operations::EdgeCollapse& collapse,
                          operations::EdgeSplit& split) {
        // op.set_priority(m_embedding_prio_long_edges_first);

        op.add_invariant(m_embedding_inversion_invariant);

        add_transfers(op);

        collapse.add_invariant(m_link_conditions);
        collapse.add_invariant(
            m_embedding_separate_substructures_invariant); // TODO pretty sure this invariant is not
                                                           // necessary

        collapse.set_new_attribute_strategy(
            m_pos_handle,
            operations::CollapseBasicStrategy::CopyOther);
        collapse.set_new_attribute_strategy(
            m_embedding_target_edge_length_attribute,
            operations::CollapseBasicStrategy::Mean);
        collapse.set_new_attribute_strategy(m_tet_tag_handle);
        collapse.set_new_attribute_strategy(m_embedding_roi_attribute);
        for (const auto& attr : none_attributes) {
            collapse.set_new_attribute_strategy(attr, operations::CollapseBasicStrategy::None);
        }

        split.set_new_attribute_strategy(m_pos_handle);
        split.set_new_attribute_strategy(
            m_embedding_target_edge_length_attribute,
            operations::SplitBasicStrategy::Copy,
            operations::SplitRibBasicStrategy::Mean);
        split.set_new_attribute_strategy(m_tet_tag_handle);
        split.set_new_attribute_strategy(m_embedding_roi_attribute);
        for (const auto& attr : none_attributes) {
            split.set_new_attribute_strategy(
                attr,
                operations::SplitBasicStrategy::None,
                operations::SplitRibBasicStrategy::None);
        }
    };

    //// 3 - 1 - 1) TetEdgeSwap 4-4 1
    // auto swap44 = std::make_shared<operations::composite::TetEdgeSwap>(m_mesh, 0);
    // setup_swap(*swap44, swap44->collapse(), swap44->split());
    //swap44->add_invariant(m_embedding_valence_4_invariant); // extra edge valance invariant
    //swap44->add_invariant(m_embedding_swap44_invariant); // check energy before swap
    //
    //// 3 - 1 - 2) TetEdgeSwap 4-4 2
    // auto swap44_2 = std::make_shared<operations::composite::TetEdgeSwap>(m_mesh, 1);
    // setup_swap(*swap44_2, swap44_2->collapse(), swap44_2->split());
    //swap44_2->add_invariant(m_embedding_valence_4_invariant); // extra edge valance invariant
    //swap44_2->add_invariant(m_embedding_swap44_2_invariant); // check energy before swap
    //
    //// 3 - 2) TetEdgeSwap 3-2
    // auto swap32 = std::make_shared<operations::composite::TetEdgeSwap>(m_mesh, 0);
    // setup_swap(*swap32, swap32->collapse(), swap32->split());
    //swap32->add_invariant(m_embedding_valence_3_invariant); // extra edge valance invariant
    //swap32->add_invariant(m_embedding_swap32_invariant); // check energy before swap
    //
    //// 3 - 3) TetFaceSwap 2-3
    // auto swap23 = std::make_shared<operations::composite::TetFaceSwap>(m_mesh);
    // setup_swap(*swap23, swap23->collapse(), swap23->split());
    //swap23->add_invariant(m_embedding_swap23_invariant); // check energy before swap

    auto swap_all = std::make_shared<operations::OrOperationSequence>(m_mesh);
    if (m_mesh.top_simplex_type() == PrimitiveType::Tetrahedron) {
        auto swap56 = std::make_shared<operations::MinOperationSequence>(m_mesh);
        for (int i = 0; i < 5; ++i) {
            auto swap = std::make_shared<operations::composite::TetEdgeSwap>(m_mesh, i);

            swap->add_invariant(std::make_shared<Swap56EnergyBeforeInvariant<double>>(
                m_mesh,
                m_pos_handle.as<double>(),
                i));

            setup_swap(*swap, swap->collapse(), swap->split());

            swap56->add_operation(swap);
        }

        auto swap56_energy_check = [&](int64_t idx, const simplex::Simplex& t) -> double {
            constexpr static PrimitiveType PV = PrimitiveType::Vertex;
            constexpr static PrimitiveType PE = PrimitiveType::Edge;
            constexpr static PrimitiveType PF = PrimitiveType::Triangle;
            constexpr static PrimitiveType PT = PrimitiveType::Tetrahedron;

            auto accessor = m_mesh.create_const_accessor(m_pos_handle.as<double>());

            const Tuple e0 = t.tuple();
            const Tuple e1 = m_mesh.switch_tuple(e0, PV);

            std::array<Tuple, 5> v;
            auto iter_tuple = e0;
            for (int64_t i = 0; i < 5; ++i) {
                v[i] = m_mesh.switch_tuples(iter_tuple, {PE, PV});
                iter_tuple = m_mesh.switch_tuples(iter_tuple, {PF, PT});
            }
            if (iter_tuple != e0) return 0;
            assert(iter_tuple == e0);

            // five iterable vertices remap to 0-4 by m_collapse_index, 0: m_collapse_index, 5:
            // e0, 6: e1
            std::array<Eigen::Vector3<double>, 7> positions = {
                {accessor.const_vector_attribute(v[(idx + 0) % 5]),
                 accessor.const_vector_attribute(v[(idx + 1) % 5]),
                 accessor.const_vector_attribute(v[(idx + 2) % 5]),
                 accessor.const_vector_attribute(v[(idx + 3) % 5]),
                 accessor.const_vector_attribute(v[(idx + 4) % 5]),
                 accessor.const_vector_attribute(e0),
                 accessor.const_vector_attribute(e1)}};

            std::array<std::array<int, 4>, 6> new_tets = {
                {{{0, 1, 2, 5}},
                 {{0, 2, 3, 5}},
                 {{0, 3, 4, 5}},
                 {{0, 1, 2, 6}},
                 {{0, 2, 3, 6}},
                 {{0, 3, 4, 6}}}};

            double new_energy_max = std::numeric_limits<double>::lowest();

            for (int i = 0; i < 6; ++i) {
                if (wmtk::utils::wmtk_orient3d(
                        positions[new_tets[i][0]],
                        positions[new_tets[i][1]],
                        positions[new_tets[i][2]],
                        positions[new_tets[i][3]]) > 0) {
                    auto energy = wmtk::function::utils::Tet_AMIPS_energy({{
                        positions[new_tets[i][0]][0],
                        positions[new_tets[i][0]][1],
                        positions[new_tets[i][0]][2],
                        positions[new_tets[i][1]][0],
                        positions[new_tets[i][1]][1],
                        positions[new_tets[i][1]][2],
                        positions[new_tets[i][2]][0],
                        positions[new_tets[i][2]][1],
                        positions[new_tets[i][2]][2],
                        positions[new_tets[i][3]][0],
                        positions[new_tets[i][3]][1],
                        positions[new_tets[i][3]][2],
                    }});


                    if (energy > new_energy_max) new_energy_max = energy;
                } else {
                    auto energy = wmtk::function::utils::Tet_AMIPS_energy({{
                        positions[new_tets[i][1]][0],
                        positions[new_tets[i][1]][1],
                        positions[new_tets[i][1]][2],
                        positions[new_tets[i][0]][0],
                        positions[new_tets[i][0]][1],
                        positions[new_tets[i][0]][2],
                        positions[new_tets[i][2]][0],
                        positions[new_tets[i][2]][1],
                        positions[new_tets[i][2]][2],
                        positions[new_tets[i][3]][0],
                        positions[new_tets[i][3]][1],
                        positions[new_tets[i][3]][2],
                    }});


                    if (energy > new_energy_max) new_energy_max = energy;
                }
            }

            return new_energy_max;
        };

        swap56->set_value_function(swap56_energy_check);
        swap56->add_invariant(std::make_shared<wmtk::invariants::EdgeValenceInvariant>(m_mesh, 5));

        auto swap44 = std::make_shared<operations::MinOperationSequence>(m_mesh);
        for (int i = 0; i < 2; ++i) {
            auto swap = std::make_shared<operations::composite::TetEdgeSwap>(m_mesh, i);

            swap->add_invariant(std::make_shared<Swap44EnergyBeforeInvariant<double>>(
                m_mesh,
                m_pos_handle.as<double>(),
                i));

            setup_swap(*swap, swap->collapse(), swap->split());

            swap44->add_operation(swap);
        }

        auto swap44_energy_check = [&](int64_t idx, const simplex::Simplex& t) -> double {
            constexpr static PrimitiveType PV = PrimitiveType::Vertex;
            constexpr static PrimitiveType PE = PrimitiveType::Edge;
            constexpr static PrimitiveType PF = PrimitiveType::Triangle;
            constexpr static PrimitiveType PT = PrimitiveType::Tetrahedron;

            auto accessor = m_mesh.create_const_accessor(m_pos_handle.as<double>());

            // get the coords of the vertices
            // input edge end points
            const Tuple e0 = t.tuple();
            const Tuple e1 = m_mesh.switch_tuple(e0, PV);
            // other four vertices
            std::array<Tuple, 4> v;
            auto iter_tuple = e0;
            for (int64_t i = 0; i < 4; ++i) {
                v[i] = m_mesh.switch_tuples(iter_tuple, {PE, PV});
                iter_tuple = m_mesh.switch_tuples(iter_tuple, {PF, PT});
            }

            if (iter_tuple != e0) return 0;
            assert(iter_tuple == e0);

            std::array<Eigen::Vector3<double>, 6> positions = {
                {accessor.const_vector_attribute(v[(idx + 0) % 4]),
                 accessor.const_vector_attribute(v[(idx + 1) % 4]),
                 accessor.const_vector_attribute(v[(idx + 2) % 4]),
                 accessor.const_vector_attribute(v[(idx + 3) % 4]),
                 accessor.const_vector_attribute(e0),
                 accessor.const_vector_attribute(e1)}};

            std::array<std::array<int, 4>, 4> new_tets = {
                {{{0, 1, 2, 4}}, {{0, 2, 3, 4}}, {{0, 1, 2, 5}}, {{0, 2, 3, 5}}}};

            double new_energy_max = std::numeric_limits<double>::lowest();

            for (int i = 0; i < 4; ++i) {
                if (wmtk::utils::wmtk_orient3d(
                        positions[new_tets[i][0]],
                        positions[new_tets[i][1]],
                        positions[new_tets[i][2]],
                        positions[new_tets[i][3]]) > 0) {
                    auto energy = wmtk::function::utils::Tet_AMIPS_energy({{
                        positions[new_tets[i][0]][0],
                        positions[new_tets[i][0]][1],
                        positions[new_tets[i][0]][2],
                        positions[new_tets[i][1]][0],
                        positions[new_tets[i][1]][1],
                        positions[new_tets[i][1]][2],
                        positions[new_tets[i][2]][0],
                        positions[new_tets[i][2]][1],
                        positions[new_tets[i][2]][2],
                        positions[new_tets[i][3]][0],
                        positions[new_tets[i][3]][1],
                        positions[new_tets[i][3]][2],
                    }});

                    if (energy > new_energy_max) new_energy_max = energy;
                } else {
                    auto energy = wmtk::function::utils::Tet_AMIPS_energy({{
                        positions[new_tets[i][1]][0],
                        positions[new_tets[i][1]][1],
                        positions[new_tets[i][1]][2],
                        positions[new_tets[i][0]][0],
                        positions[new_tets[i][0]][1],
                        positions[new_tets[i][0]][2],
                        positions[new_tets[i][2]][0],
                        positions[new_tets[i][2]][1],
                        positions[new_tets[i][2]][2],
                        positions[new_tets[i][3]][0],
                        positions[new_tets[i][3]][1],
                        positions[new_tets[i][3]][2],
                    }});

                    if (energy > new_energy_max) new_energy_max = energy;
                }
            }

            return new_energy_max;
        };

        swap44->set_value_function(swap44_energy_check);
        swap44->add_invariant(std::make_shared<wmtk::invariants::EdgeValenceInvariant>(m_mesh, 4));

        auto swap32 = std::make_shared<operations::composite::TetEdgeSwap>(m_mesh, 0);
        setup_swap(*swap32, swap32->collapse(), swap32->split());
        swap32->add_invariant(m_embedding_valence_3_invariant); // extra edge valance invariant
        swap32->add_invariant(m_embedding_swap32_invariant); // check energy before swap

        swap_all->add_operation(swap32);
        swap_all->add_operation(swap44);
        swap_all->add_operation(swap56);
    } else if (m_mesh.top_simplex_type() == PrimitiveType::Triangle) {
        auto swap = std::make_shared<operations::composite::TriEdgeSwap>(m_mesh);

        setup_swap(*swap, swap->collapse(), swap->split());

        swap_all->add_operation(swap);
    }
    // no need for the ROI invariant as it has the AmipsBeforeInvariant
    swap_all->add_invariant(
        std::make_shared<invariants::AmipsBeforeInvariant>(m_embedding_amips_attribute, 100));
    swap_all->add_invariant(m_embedding_no_child_invariants);

    swap_all->set_priority(m_embedding_prio_long_edges_first);

    //////////////////////////////////
    // 4) Smoothing
    //////////////////////////////////
    // auto smooth = std::make_shared<operations::OptimizationSmoothing>(m_amips_energy);
    auto smooth = std::make_shared<operations::AmipsOptimization>(m_pos_handle);
    auto solver_settings = smooth->nonlinear_solver_params();
    solver_settings["max_iterations"] = 2;
    smooth->set_nonlinear_solver_params(solver_settings);

    smooth->add_invariant(m_embedding_roi_invariant);
    smooth->add_invariant(m_embedding_no_child_invariants);
    smooth->add_invariant(m_embedding_inversion_invariant);

    std::vector<std::shared_ptr<operations::Operation>> ops{split, collapse, swap_all};
    std::vector<std::string> ops_name{"split", "collapse", "swap_all"};

    m_embedding_edge_length_transfer->run_on_all();
    m_embedding_amips_transfer->run_on_all();

    static int64_t write_counter = 0;
    if (write_counter == 0) {
        components::utils::write_mesh(
            m_pos_handle,
            fmt::format("embedding_opt_{}", write_counter++),
            m_debug_print_embedding);
    }

    Scheduler scheduler;
    for (int64_t i = 0; i < n_iterations; ++i) {
        logger().info("Embedding Pass {}/{}", i + 1, n_iterations);
        SchedulerStats pass_stats;

        logger().info("Update target edge lengths and ROI.");
        m_embedding_target_edge_length_transfer->run_on_all();
        // tel_update();
        roi_update();

        logger().info("Perform operations.");
        int jj = 0;
        for (auto& op : ops) {
            auto stats = scheduler.run_operation_on_all_parallel_prefilter(*op);
            pass_stats += stats;
            logger().info(
                "Executed {}, {} ops (S/F) {}/{}. Time: executing: {}",
                ops_name[jj],
                stats.number_of_performed_operations(),
                stats.number_of_successful_operations(),
                stats.number_of_failed_operations(),
                stats.executing_time);
            ++jj;

            components::utils::write_mesh(
                m_pos_handle,
                fmt::format("embedding_opt_{}", write_counter++),
                m_debug_print_embedding);

            check_tet_tag_validity();
        }
        // smoothing
        for (int64_t i = 0; i < 5; ++i) {
            auto stats = scheduler.run_operation_on_all_parallel_prefilter(*smooth);
            pass_stats += stats;
            logger().info(
                "Executed smooth {}, {} ops (S/F) {}/{}. Time: executing: {}",
                i,
                stats.number_of_performed_operations(),
                stats.number_of_successful_operations(),
                stats.number_of_failed_operations(),
                stats.executing_time);


            components::utils::write_mesh(
                m_pos_handle,
                fmt::format("embedding_opt_{}", write_counter++),
                m_debug_print_embedding);
        }

        m_embedding_edge_length_transfer->run_on_all();
        m_embedding_amips_transfer->run_on_all();

        logger().info(
            "Executed {} ops (S/F) {}/{}. Time: executing: {}",
            pass_stats.number_of_performed_operations(),
            pass_stats.number_of_successful_operations(),
            pass_stats.number_of_failed_operations(),
            pass_stats.executing_time);

        multimesh::consolidate(m_mesh);
    }
}

void OffsetOptimization::optimize_offset(const int64_t n_iterations)
{
    std::vector<attribute::MeshAttributeHandle> none_attributes(
        {m_offset_edge_length_attribute,
         m_offset_normal_samples_attribute,
         m_offset_mean_ratio_metric_attribute,
         m_offset_face_distance_attribute,
         m_offset_face_target_distance_attribute,
         m_offset_vertex_distance_attribute,
         m_offset_vertex_target_distance_attribute,
         m_offset_face_normal_deviation_attribute,
         m_offset_face_converged_attribute,
         // m_pos_handle,
         m_embedding_edge_length_attribute,
         m_embedding_target_edge_length_attribute,
         m_embedding_amips_attribute,
         m_embedding_roi_attribute});

    auto add_transfers = [this](operations::Operation& op) {
        // op.add_transfer_strategy(m_offset_position_to_embedding_transfer);
        op.add_transfer_strategy(m_offset_edge_length_transfer);
        op.add_transfer_strategy(m_offset_point_to_face_transfer);
        // op.add_transfer_strategy(m_offset_mean_ratio_metric_transfer);
        // op.add_transfer_strategy(m_offset_face_normal_deviation_transfer);
        //  op.add_transfer_strategy(m_embedding_edge_length_transfer);
        //  op.add_transfer_strategy(m_embedding_amips_transfer);
        //  op.add_transfer_strategy(m_embedding_target_edge_length_transfer);
    };

    //////////////////////////////////
    // 1) EdgeSplit
    //////////////////////////////////
    auto split = std::make_shared<operations::EdgeSplit>(*m_offset_mesh);
    split->set_priority(m_offset_prio_long_edges_first);

    split->add_invariant(m_offset_todo_larger);
    split->add_invariant(m_offset_convergence_invariant);
    split->add_invariant(m_embedding_inversion_invariant); // we could reach floating point
                                                           // precision so it's better to have that

    if (m_input_mesh->top_simplex_type() == PrimitiveType::Triangle) {
        split->add_invariant(
            std::make_shared<invariants::LongestEdgeInvariant>(m_offset_pos_handle));
    }

    {
        split->set_new_attribute_strategy(m_offset_pos_handle);
        split->set_new_attribute_strategy(m_pos_handle);
        split->set_new_attribute_strategy(
            m_offset_vertex_converged_attribute,
            operations::SplitBasicStrategy::None,
            operations::SplitRibBasicStrategy::None);
        split->set_new_attribute_strategy(
            m_offset_target_edge_length_attribute,
            operations::SplitBasicStrategy::Copy,
            operations::SplitRibBasicStrategy::Mean);
        split->set_new_attribute_strategy(m_tet_tag_handle);
        for (const auto& attr : none_attributes) {
            split->set_new_attribute_strategy(
                attr,
                operations::SplitBasicStrategy::None,
                operations::SplitRibBasicStrategy::None);
        }
    }

    add_transfers(*split);

    //////////////////////////////////
    // 2) EdgeCollapse
    //////////////////////////////////

    auto setup_halfedge_collapse =
        [this, &none_attributes, &add_transfers](
            const int64_t collapse_type) -> std::shared_ptr<operations::EdgeCollapse> {
        auto _collapse = std::make_shared<operations::EdgeCollapse>(*m_offset_mesh);
        _collapse->set_priority(m_offset_prio_short_edges_first);

        _collapse->add_invariant(m_normal_deviation_invariant);
        _collapse->add_invariant(m_embedding_inversion_invariant);

        //// use different before energy depending on the type
        //_collapse->add_invariant(
        //    std::make_shared<wmtk::invariants::CollapseEnergyBeforeInvariant<double>>(
        //        m_mesh,
        //        m_pos_handle.as<double>(),
        //        m_embedding_amips_attribute.as<double>(),
        //        collapse_type - 1));

        if (m_mesh.top_simplex_type() == PrimitiveType::Tetrahedron) {
            _collapse->add_invariant(std::make_shared<invariants::OffsetCollapseBeforeInvariant>(
                m_offset_normal_samples_attribute,
                m_offset_pos_handle,
                m_max_normal_deviation,
                collapse_type));
        }


        // use different strategies for position depending on the type
        switch (collapse_type) {
        case 1:
            _collapse->set_new_attribute_strategy(
                m_offset_pos_handle,
                m_offset_pos_collapse_strat_1);
            _collapse->set_new_attribute_strategy(
                m_offset_vertex_converged_attribute,
                m_offset_vertex_converged_collapse_strat_1);
            _collapse->set_new_attribute_strategy(m_pos_handle, m_embedding_pos_collapse_strat_1);
            break;
        case 2:
            _collapse->set_new_attribute_strategy(
                m_offset_pos_handle,
                m_offset_pos_collapse_strat_2);
            _collapse->set_new_attribute_strategy(
                m_offset_vertex_converged_attribute,
                m_offset_vertex_converged_collapse_strat_2);
            _collapse->set_new_attribute_strategy(m_pos_handle, m_embedding_pos_collapse_strat_2);
            break;
        default: log_and_throw_error("Unknown halfedge collapse type: {}", collapse_type); break;
        }

        // _collapse->set_new_attribute_strategy(m_offset_pos_handle);

        _collapse->set_new_attribute_strategy(
            m_offset_target_edge_length_attribute,
            operations::CollapseBasicStrategy::Mean);
        _collapse->set_new_attribute_strategy(m_tet_tag_handle);
        for (const auto& attr : none_attributes) {
            _collapse->set_new_attribute_strategy(attr, operations::CollapseBasicStrategy::None);
        }

        add_transfers(*_collapse);

        return _collapse;
    };

    auto collapse_1 = setup_halfedge_collapse(1);
    auto collapse_2 = setup_halfedge_collapse(2);
    auto collapse = std::make_shared<operations::OrOperationSequence>(*m_offset_mesh);
    collapse->add_operation(collapse_1);
    collapse->add_operation(collapse_2);
    collapse->add_invariant(m_offset_todo_smaller);
    collapse->add_invariant(m_offset_convergence_invariant);
    collapse->add_invariant(m_link_conditions);

    //////////////////////////////////
    // 3) Swap
    //////////////////////////////////
    std::shared_ptr<operations::composite::TriEdgeSwap> swap;
    if (m_mesh.top_simplex_type() == PrimitiveType::Tetrahedron) {
        swap = std::make_shared<operations::composite::TriEdgeSwap>(*m_offset_mesh);
        swap->set_priority(m_offset_prio_long_edges_first);

        swap->add_invariant(m_offset_convergence_invariant);
        swap->add_invariant(m_offset_swap_invariant);
        // swap->add_invariant(m_normal_deviation_invariant);

        swap->collapse().add_invariant(m_link_conditions);
        swap->collapse().add_invariant(m_embedding_inversion_invariant);

        {
            swap->split().set_new_attribute_strategy(
                m_offset_pos_handle,
                operations::SplitBasicStrategy::None,
                operations::SplitRibBasicStrategy::None);
            swap->split().set_new_attribute_strategy(
                m_pos_handle,
                operations::SplitBasicStrategy::None,
                operations::SplitRibBasicStrategy::None);
            swap->split().set_new_attribute_strategy(
                m_offset_vertex_converged_attribute,
                operations::SplitBasicStrategy::None,
                operations::SplitRibBasicStrategy::None);
            swap->split().set_new_attribute_strategy(
                m_offset_target_edge_length_attribute,
                operations::SplitBasicStrategy::Copy,
                operations::SplitRibBasicStrategy::Mean);
            swap->split().set_new_attribute_strategy(m_tet_tag_handle);
            for (const auto& attr : none_attributes) {
                swap->split().set_new_attribute_strategy(
                    attr,
                    operations::SplitBasicStrategy::None,
                    operations::SplitRibBasicStrategy::None);
            }

            swap->collapse().set_new_attribute_strategy(
                m_offset_pos_handle,
                operations::CollapseBasicStrategy::CopyOther);
            swap->collapse().set_new_attribute_strategy(
                m_pos_handle,
                operations::CollapseBasicStrategy::CopyOther);
            swap->collapse().set_new_attribute_strategy(
                m_offset_vertex_converged_attribute,
                operations::CollapseBasicStrategy::CopyOther);
            swap->collapse().set_new_attribute_strategy(
                m_offset_target_edge_length_attribute,
                operations::CollapseBasicStrategy::Mean);
            swap->collapse().set_new_attribute_strategy(m_tet_tag_handle);
            for (const auto& attr : none_attributes) {
                swap->collapse().set_new_attribute_strategy(
                    attr,
                    operations::CollapseBasicStrategy::None);
            }
        }

        add_transfers(*swap);
    }

    //////////////////////////////////
    // 4) Smoothing
    //////////////////////////////////

    // Smoothing with binary search for valid position
    auto smooth = std::make_shared<operations::AttributesUpdateWithFunction>(*m_offset_mesh);
    smooth->set_function(m_offset_smoothing_function);

    // Smoothing with optimization
    // auto smooth = std::make_shared<operations::AmipsOffsetOptimization>(
    //     m_pos_handle,
    //     m_offset_pos_handle,
    //     *m_input_bvh,
    //     m_offset_distance);

    smooth->add_invariant(m_offset_convergence_invariant);
    smooth->add_invariant(m_embedding_inversion_invariant);
    smooth->add_invariant(std::make_shared<invariants::NormalDeviationAfterInvariant>(
        m_offset_face_normal_deviation_attribute,
        60,
        true)); // do not mess up normal deviation completely

    std::vector<std::shared_ptr<operations::Operation>> ops;
    std::vector<std::string> ops_name;
    if (m_mesh.top_simplex_type() == PrimitiveType::Tetrahedron) {
        ops = std::vector<std::shared_ptr<operations::Operation>>{split, collapse, swap};
        ops_name = std::vector<std::string>{"split", "collapse", "swap"};
    } else if (m_mesh.top_simplex_type() == PrimitiveType::Triangle) {
        ops = std::vector<std::shared_ptr<operations::Operation>>{split, collapse};
        ops_name = std::vector<std::string>{"split", "collapse"};
    }

    static int64_t write_counter = 0;
    if (write_counter == 0) {
        components::utils::write_mesh(
            m_offset_pos_handle,
            fmt::format("offset_opt_{}", write_counter++),
            m_debug_print_offset);
    }

    Scheduler scheduler;
    for (int64_t i = 0; i < n_iterations; ++i) {
        logger().info("Offset Pass {}/{}", i + 1, n_iterations);
        SchedulerStats pass_stats;

        logger().info("Update target edge lengths.");
        m_offset_target_edge_length_transfer->run_on_all();

        logger().info("Perform operations.");
        int jj = 0;
        for (auto& op : ops) {
            auto stats = scheduler.run_operation_on_all_parallel_prefilter(*op);
            pass_stats += stats;
            logger().info(
                "Executed {}, {} ops (S/F) {}/{}. Time: executing: {}",
                ops_name[jj],
                stats.number_of_performed_operations(),
                stats.number_of_successful_operations(),
                stats.number_of_failed_operations(),
                stats.executing_time);

            ++jj;

            components::utils::write_mesh(
                m_offset_pos_handle,
                fmt::format("offset_opt_{}", write_counter++),
                m_debug_print_offset);

            check_tet_tag_validity();
        }

        // smoothing
        for (int64_t i = 0; i < 5; ++i) {
            auto stats = scheduler.run_operation_on_all_parallel_prefilter(*smooth);
            pass_stats += stats;
            logger().info(
                "Executed smooth {}, {} ops (S/F) {}/{}. Time: executing: {}",
                i,
                stats.number_of_performed_operations(),
                stats.number_of_successful_operations(),
                stats.number_of_failed_operations(),
                stats.executing_time);

            components::utils::write_mesh(
                m_offset_pos_handle,
                fmt::format("offset_opt_{}", write_counter++),
                m_debug_print_offset);
        }

        m_offset_edge_length_transfer->run_on_all();
        m_offset_point_to_face_transfer->run_on_all();

        // m_offset_target_edge_length_transfer->run_on_all();
        //// limit TEL difference between two edges to a factor of two
        //{
        //    auto tel_acc =
        //        m_offset_mesh->create_accessor(m_offset_target_edge_length_attribute.as<double>());
        //    for (const Tuple& v : m_offset_mesh->get_all(PrimitiveType::Vertex)) {
        //        const auto edges = simplex::cofaces_single_dimension_tuples(
        //            *m_offset_mesh,
        //            simplex::Simplex::vertex(*m_offset_mesh, v),
        //            PrimitiveType::Edge);
        //
        //        double min_tel = std::numeric_limits<double>::max();
        //        for (const Tuple& e : edges) {
        //            min_tel = std::min(min_tel, tel_acc.const_scalar_attribute(e));
        //        }
        //        for (const Tuple& e : edges) {
        //            double new_tel = std::min(tel_acc.const_scalar_attribute(e), 2 * min_tel);
        //            tel_acc.scalar_attribute(e) = new_tel;
        //        }
        //    }
        //}

        logger().info(
            "Pass summary: Executed {} ops (S/F) {}/{}. Time: executing: {}",
            pass_stats.number_of_performed_operations(),
            pass_stats.number_of_successful_operations(),
            pass_stats.number_of_failed_operations(),
            pass_stats.executing_time);

        multimesh::consolidate(*m_offset_mesh);

        components::utils::write_mesh(
            m_offset_pos_handle,
            fmt::format("offset_opt_{}", write_counter++),
            m_debug_print_offset);
    }
}

void OffsetOptimization::smooth_all(const int64_t n_iterations)
{
    auto smooth_offset = std::make_shared<operations::AttributesUpdateWithFunction>(*m_offset_mesh);
    smooth_offset->set_function(m_offset_smoothing_function);

    smooth_offset->add_invariant(m_offset_convergence_invariant);
    smooth_offset->add_invariant(m_embedding_inversion_invariant);
    smooth_offset->add_invariant(std::make_shared<invariants::NormalDeviationAfterInvariant>(
        m_offset_face_normal_deviation_attribute,
        60,
        true)); // do not mess up normal deviation completely

    auto smooth_embedding = std::make_shared<operations::AmipsOptimization>(m_pos_handle);
    auto solver_settings = smooth_embedding->nonlinear_solver_params();
    solver_settings["max_iterations"] = 2;
    smooth_embedding->set_nonlinear_solver_params(solver_settings);

    smooth_embedding->add_invariant(m_embedding_no_child_invariants);
    smooth_embedding->add_invariant(m_embedding_inversion_invariant);

    std::vector<std::shared_ptr<operations::Operation>> ops{smooth_offset, smooth_embedding};
    std::vector<std::string> ops_name{"smooth_offset", "smooth_embedding"};

    static int64_t write_counter = 0;
    components::utils::write_mesh(
        m_offset_pos_handle,
        fmt::format("smooth_opt_{}", write_counter++),
        m_debug_print_smooth_all);

    Scheduler scheduler;
    for (int64_t i = 0; i < n_iterations; ++i) {
        logger().info("Smooth All Pass {}/{}", i + 1, n_iterations);
        SchedulerStats pass_stats;

        logger().info("Perform operations.");
        int jj = 0;
        for (auto& op : ops) {
            auto stats = scheduler.run_operation_on_all_parallel_prefilter(*op);
            pass_stats += stats;
            logger().info(
                "Executed {}, {} ops (S/F) {}/{}. Time: executing: {}",
                ops_name[jj],
                stats.number_of_performed_operations(),
                stats.number_of_successful_operations(),
                stats.number_of_failed_operations(),
                stats.executing_time);
            ++jj;
        }

        m_offset_point_to_face_transfer_wout_convergence->run_on_all();

        components::utils::write_mesh(
            m_offset_pos_handle,
            fmt::format("smooth_opt_{}", write_counter++),
            m_debug_print_smooth_all);
    }

    m_embedding_amips_transfer->run_on_all();
}

std::vector<double> OffsetOptimization::get_error_metrics()
{
    const auto mrm_acc =
        m_offset_mesh->create_const_accessor(m_offset_mean_ratio_metric_attribute.as<double>());
    const auto dist_acc =
        m_offset_mesh->create_const_accessor(m_offset_face_distance_attribute.as<double>());
    const auto nd_acc =
        m_offset_mesh->create_const_accessor(m_offset_face_normal_deviation_attribute.as<double>());

    const auto pos_acc = m_offset_mesh->create_const_accessor(m_offset_pos_handle.as<double>());

    const std::vector<Tuple> triangles = m_offset_mesh->get_all(m_offset_mesh->top_simplex_type());

    double mrm_avg = 0;
    double mrm_min = std::numeric_limits<double>::max();

    double face_dist_err_avg = 0;
    double face_dist_err_max = std::numeric_limits<double>::lowest();

    double nd_avg = 0;
    double nd_max = std::numeric_limits<double>::lowest();

    for (const Tuple& t : triangles) {
        const double mrm = mrm_acc.const_scalar_attribute(t);
        mrm_avg += mrm;
        mrm_min = std::min(mrm_min, mrm);

        const double dist_err = std::abs(dist_acc.const_scalar_attribute(t) - m_offset_distance);
        face_dist_err_avg += dist_err;
        face_dist_err_max = std::max(face_dist_err_max, dist_err);

        const double nd = nd_acc.const_scalar_attribute(t);
        nd_avg += nd;
        nd_max = std::max(nd_max, nd);
    }

    mrm_avg /= triangles.size();
    face_dist_err_avg /= triangles.size() / m_offset_distance;
    face_dist_err_max /= m_offset_distance;
    nd_avg /= triangles.size();

    double v_rel_adapted_dist_err_avg = 0;
    double v_rel_adapted_dist_err_max = std::numeric_limits<double>::lowest();
    const std::vector<Tuple> vertices = m_offset_mesh->get_all(PrimitiveType::Vertex);
    for (const Tuple& t : vertices) {
        const auto p = pos_acc.const_vector_attribute(t);
        const auto [sq_dist, od] = m_input_bvh->sq_dist_and_offset_distance(p);
        const double dist = std::sqrt(sq_dist);
        const double rel_dist = std::abs(dist - od) / m_offset_distance;
        v_rel_adapted_dist_err_avg += rel_dist;
        v_rel_adapted_dist_err_max = std::max(v_rel_adapted_dist_err_max, rel_dist);
    }

    v_rel_adapted_dist_err_avg /= vertices.size();

    nlohmann::json report;

    report["mean_ratio_metric_average"] = mrm_avg;
    report["mean_ratio_metric_min"] = mrm_min;
    report["relative_distance_average"] = face_dist_err_avg;
    report["relative_distance_max"] = face_dist_err_max;
    report["v_relative_adapted_distance_average"] = v_rel_adapted_dist_err_avg;
    report["v_relative_adapted_distance_max"] = v_rel_adapted_dist_err_max;
    report["normal_deviation_average"] = nd_avg;
    report["normal_deviation_max"] = nd_max;
    report["num_offset_triangles"] = triangles.size();
    report["num_offset_vertices"] = m_offset_mesh->get_all(PrimitiveType::Vertex).size();
    report["num_embedding_tets"] = m_mesh.get_all(PrimitiveType::Tetrahedron).size();
    report["num_embedding_vertices"] = m_mesh.get_all(PrimitiveType::Vertex).size();
    m_metrics_report["metrics"].push_back(report);

    logger().info("Relative adapted distance average = {:.4}%", 100 * v_rel_adapted_dist_err_avg);
    logger().info("Relative adapted distance max = {:.4}%", 100 * v_rel_adapted_dist_err_max);
    logger().info("Normal deviation average = {:.4}", nd_avg);

    return {100 * v_rel_adapted_dist_err_avg, 100 * v_rel_adapted_dist_err_max, nd_avg};
}

nlohmann::json& OffsetOptimization::metrics_json()
{
    return m_metrics_report;
}

void OffsetOptimization::set_debug_prints(bool embedding, bool offset, bool smooth)
{
    m_debug_print_embedding = embedding;
    m_debug_print_offset = offset;
    m_debug_print_smooth_all = smooth;
}

void OffsetOptimization::roi_update()
{
    auto roi_acc = m_mesh.create_accessor(m_embedding_roi_attribute.as<int64_t>());

    const auto amips_acc = m_mesh.create_const_accessor(m_embedding_amips_attribute.as<double>());

    const auto tets = m_mesh.get_all(PrimitiveType::Tetrahedron);

    // set all tags to 0
    // tag all tets with high energy
    for (const Tuple& t : tets) {
        if (amips_acc.const_scalar_attribute(t) > 100) {
            roi_acc.scalar_attribute(t) = 1;
        } else {
            roi_acc.scalar_attribute(t) = 0;
        }
    }

    // spread tag to neighbors
    for (int64_t i = 1; i < 3; ++i) {
        for (const Tuple& t : tets) {
            if (roi_acc.const_scalar_attribute(t) != i) {
                continue;
            }

            const auto vertices = simplex::faces_single_dimension(
                m_mesh,
                simplex::Simplex::tetrahedron(m_mesh, t),
                PrimitiveType::Vertex);

            for (const simplex::Simplex& v : vertices) {
                // do not spread through substructures
                bool is_child = false;
                for (const auto& child : m_mesh.get_child_meshes()) {
                    if (m_mesh.simplex_is_in_child(*child, v)) {
                        is_child = true;
                        break;
                    }
                }
                if (is_child) {
                    continue;
                }

                for (const Tuple& vtt : simplex::top_dimension_cofaces_iterable(m_mesh, v)) {
                    if (roi_acc.const_scalar_attribute(vtt) < i) {
                        roi_acc.scalar_attribute(vtt) = i + 1;
                    }
                }
            }
        }
    }

    for (const Tuple& t : tets) {
        if (roi_acc.const_scalar_attribute(t) != 0) {
            roi_acc.scalar_attribute(t) = 1;
        }
    }
}

std::tuple<Eigen::Vector3d, double, double> OffsetOptimization::compute_offset_normal(
    const Eigen::Vector3d& p)
{
    const auto [sq_dist, offset_distance, nearest_point] =
        m_input_bvh->sq_dist_offset_distance_nearest_point(p);

    Eigen::Vector3d n = (nearest_point - p).normalized();

    return std::make_tuple(n, sq_dist, offset_distance);
}

std::tuple<Eigen::Vector2d, double, double> OffsetOptimization::compute_offset_normal(
    const Eigen::Vector2d& p)
{
    const auto [sq_dist, offset_distance, nearest_point] =
        m_input_bvh->sq_dist_offset_distance_nearest_point(p);

    Eigen::Vector2d n = (nearest_point - p).normalized();

    return std::make_tuple(n, sq_dist, offset_distance);
}

void OffsetOptimization::check_tet_tag_validity()
{
    const auto tet_tag_acc = m_mesh.create_const_accessor<int64_t>(m_tet_tag_handle);
    for (const Tuple& t : m_offset_mesh->get_all(m_offset_mesh->top_simplex_type())) {
        const simplex::Simplex f(*m_offset_mesh, m_offset_mesh->top_simplex_type(), t);
        Tuple tet_tuple = m_offset_mesh->map_to_parent_tuple(f);
        if (tet_tag_acc.const_scalar_attribute(tet_tuple) == 0) {
            const int64_t gid0 = wmtk::utils::TupleInspector::global_cid(tet_tuple);
            tet_tuple = m_mesh.switch_tuple(tet_tuple, m_mesh.top_simplex_type());
            if (tet_tag_acc.const_scalar_attribute(tet_tuple) == 0) {
                // components::utils::write_mesh(m_pos_handle, "embedding_tag_missing", true);
                // components::utils::write_mesh(m_offset_pos_handle, "offset_tag_missing", true);
                // components::utils::write_mesh(m_input_pos_handle, "input_tag_missing", true);
                const int64_t gid1 = wmtk::utils::TupleInspector::global_cid(tet_tuple);
                logger().warn("no tagged tet incident to offset. Tuples: {}, {}", gid0, gid1);
                return;
            }
        }
    }
}

} // namespace wmtk::components::internal