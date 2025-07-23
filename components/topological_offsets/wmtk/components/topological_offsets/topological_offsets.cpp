#include "topological_offsets.hpp"

#include <wmtk/Mesh.hpp>
#include <wmtk/Scheduler.hpp>
#include <wmtk/TetMesh.hpp>
#include <wmtk/TriMesh.hpp>
#include <wmtk/io/MeshReader.hpp>
#include <wmtk/io/ParaviewWriter.hpp>
#include <wmtk/multimesh/consolidate.hpp>
#include <wmtk/operations/attribute_update/AttributeTransferStrategy.hpp>
#include <wmtk/simplex/faces_single_dimension.hpp>
#include <wmtk/utils/orient.hpp>

#include <wmtk/utils/Logger.hpp>
#include <wmtk/utils/Stopwatch.hpp>
#include <wmtk/utils/bbox_from_mesh.hpp>
#include <wmtk/utils/primitive_range.hpp>

#include <wmtk/components/topological_offsets/internal/utils/bvh_from_mesh.hpp>
#include <wmtk/components/utils/get_attributes.hpp>

#include "internal/OffsetOptimization.hpp"
#include "internal/TopologicalOffsetGenerator.hpp"
#include "internal/utils/ConvergenceRate.hpp"
#include "internal/utils/OffsetBvh.hpp"
#include "internal/utils/apply_distance_function.hpp"
#include "internal/utils/write_mesh.hpp"

namespace wmtk::components {

namespace {

PrimitiveType find_tag_primitive_type(
    std::map<PrimitiveType, attribute::MeshAttributeHandle> tag_handles,
    const int64_t tag_value)
{
    const Mesh& m = tag_handles[PrimitiveType::Vertex].mesh();

    for (const PrimitiveType pt : wmtk::utils::primitive_below(m.top_simplex_type())) {
        auto tag_acc = m.create_const_accessor(tag_handles[pt].as<int64_t>());
        for (const Tuple& t : m.get_all(pt)) {
            if (tag_acc.const_scalar_attribute(t) == tag_value) {
                return pt;
            }
        }
    }

    log_and_throw_error("Cannot generate offset. Tag {} not found.", tag_value);
}

attribute::MeshAttributeHandle register_and_transfer_positions(
    attribute::MeshAttributeHandle parent_pos_handle,
    Mesh& child,
    const std::string& pos_attr_name = "vertices")
{
    auto child_pos_handle = child.register_attribute<double>(
        pos_attr_name,
        PrimitiveType::Vertex,
        parent_pos_handle.dimension());

    auto propagate_to_child_position = [](const Eigen::MatrixXd& P) -> Eigen::VectorXd {
        return P;
    };

    auto offset_pos_transfer =
        std::make_shared<wmtk::operations::SingleAttributeTransferStrategy<double, double>>(
            child_pos_handle,
            parent_pos_handle,
            propagate_to_child_position);
    offset_pos_transfer->run_on_all();

    return child_pos_handle;
}

} // namespace

TopologicalOffsetsOutput topological_offsets(TopologicalOffsetsOptions options)
{
    attribute::MeshAttributeHandle pos_handle = options.embedding_position_handle;
    Mesh& mesh = pos_handle.mesh();
    Mesh& input_mesh = *options.input_mesh;

    if (!mesh.is_connectivity_valid()) {
        log_and_throw_error("Input mesh for topological_offsets connectivity invalid.");
    }
    if (mesh.top_simplex_type() != PrimitiveType::Tetrahedron &&
        mesh.top_simplex_type() != PrimitiveType::Triangle) {
        log_and_throw_error(
            "Input mesh for topological_offsets must be a tetrahedral or triangle mesh.");
    }
    if (input_mesh.top_simplex_type() != PrimitiveType::Triangle &&
        input_mesh.top_simplex_type() != PrimitiveType::Edge) {
        log_and_throw_error(
            "topological_offsets was not tested for anything else but triangles or edges.");
    }

    if (options.distance < 0) {
        log_and_throw_error("Negative distance is not supported. Use child-meshes to tag the "
                            "exterior or interior region, if a one-sided offset is desired.");
    }

    if (!options.min_edge_length) {
        const double r = options.distance;
        const double nd_max = options.max_normal_deviation;
        const double l_min = 2 * r * std::sin(nd_max * M_PI / 180.);
        logger().warn(
            "No min edge length given. Computing it from distance (r) and max normal "
            "deviation (nd): 2 * r * sin(nd) --> min_edge_length = {}.",
            l_min);
        options.min_edge_length = l_min;
    }
    if (!options.max_edge_length) {
        constexpr double l_max = std::numeric_limits<double>::infinity();
        logger().warn(
            "No max edge length given. Setting it to infinity. --> max_edge_length = {}",
            l_max);
        options.max_edge_length = l_max;
    }

    TopologicalOffsetsOutput output;

    std::shared_ptr<internal::utils::OffsetBvh> offset_bvh;
    // init bvh + make distance absolute
    {
        wmtk::utils::StopWatch sw("bvh_initialization");

        auto input_pos_handle = register_and_transfer_positions(pos_handle, input_mesh);

        if (options.relative_distance_and_length) {
            const double diag = wmtk::utils::bbox_diagonal_from_mesh(input_pos_handle);
            options.distance *= diag;
            options.min_edge_length.value() *= diag;
            options.max_edge_length.value() *= diag;

            logger().info(
                "Input distance was relative to bbox diagonal = {}. \ndistance = "
                "{}\nmin_edge_length = "
                "{}\nmax_edge_length = {}",
                diag,
                options.distance,
                options.min_edge_length.value(),
                options.max_edge_length.value());
        }

        offset_bvh = std::make_shared<internal::utils::OffsetBvh>(input_pos_handle);

        if (options.restrict_min_edge_length_to_input_avg) {
            // compute avg input mesh edge length
            const auto p_acc = input_mesh.create_const_accessor<double>(input_pos_handle);

            double l_avg = 0;
            const auto edges = input_mesh.get_all(PrimitiveType::Edge);
            for (const Tuple& e_tuple : edges) {
                const simplex::Simplex v0(input_mesh, PrimitiveType::Vertex, e_tuple);
                const simplex::Simplex v1(
                    input_mesh,
                    PrimitiveType::Vertex,
                    input_mesh.switch_tuple(e_tuple, PrimitiveType::Vertex));
                const auto p0 = p_acc.const_vector_attribute(v0);
                const auto p1 = p_acc.const_vector_attribute(v1);
                const double l = (p1 - p0).norm();
                l_avg += l;
            }
            l_avg /= edges.size();
            options.min_edge_length.value() = std::max(options.min_edge_length.value(), l_avg);
            logger().info(
                "Input average edge length = {} --> min_edge_length = {}",
                l_avg,
                options.min_edge_length.value());
        }

        input_mesh.clear_attributes();

        sw.stop();
        output.report["timings_seconds"]["bvh_init"] = sw.getElapsedTime();
    }

    attribute::MeshAttributeHandle top_simplex_tag_handle;

    // generate topological offset
    {
        wmtk::utils::StopWatch sw("top_init");

        TopologicalOffsetGenerator tog(
            mesh,
            pos_handle,
            options.distance,
            options.intermediate_output);

        tog.add_input_substructure(input_mesh);
        if (options.inside_mesh) {
            tog.add_input_substructure(*options.inside_mesh);
        }
        // tog.init_offset_distance_field(input_mesh, options.distance);
        {
            attribute::MeshAttributeHandle input_pos_handle =
                register_and_transfer_positions(pos_handle, input_mesh);

            offset_bvh->update_position_handle();
            offset_bvh->register_offset_distance_attribute(options.distance);
            tog.register_bvh(offset_bvh);

            auto v_offset_distance_handle = offset_bvh->v_offset_distance_handle();
            if (!options.distance_function.empty()) {
                internal::utils::apply_distance_function(
                    input_pos_handle,
                    v_offset_distance_handle,
                    options.distance,
                    options.distance_function);
            }
        }

        logger().info("Tag d-simplex as inside");
        {
            wmtk::utils::StopWatch sw_in("Tag inside tets");

            if (options.finite_offset) {
                tog.expand_finite_offset();
            } else {
                if (options.adapt_offset_distance) {
                    tog.adapt_offset_distance();
                }
                if (options.use_warm_start) {
                    tog.expand_topological_offset_input();
                }
                tog.regularize(false);
            }
        }

        {
            offset_bvh->update_offset_distance_handle();
            const auto f_offset_distance_handle = offset_bvh->f_offset_distance_handle();
            const auto v_offset_distance_handle = offset_bvh->v_offset_distance_handle();
            input_mesh.clear_attributes({f_offset_distance_handle, v_offset_distance_handle});
            offset_bvh->update_offset_distance_handle();
        }

        logger().info("Generate simplicial embedding of the input.");
        tog.regularize(options.use_simplicial_embedding);

        logger().info("Generate topological offset.");

        tog.marching(offset_bvh->bvh(), 0);
        if (options.expand_topo && !options.finite_offset) {
            top_simplex_tag_handle = tog.tag_offset_tets();
            input_mesh.consolidate();
            offset_bvh->update_offset_distance_handle();
            const auto f_offset_distance_handle = offset_bvh->f_offset_distance_handle();
            const auto v_offset_distance_handle = offset_bvh->v_offset_distance_handle();
            input_mesh.clear_attributes({f_offset_distance_handle, v_offset_distance_handle});
            register_and_transfer_positions(pos_handle, input_mesh);
            offset_bvh->update_position_handle();
            offset_bvh->update_offset_distance_handle();
            offset_bvh->update_face_id_map();
            tog.expand_topological_offset_post();
            output.offset_mesh = tog.generate_substructure_from_offset_tag();
        } else {
            // tagging of offsets tets needs to be done afterwards for that case
            output.offset_mesh = tog.generate_substructure_from_offset_tag();
            top_simplex_tag_handle = tog.tag_offset_tets();
        }

        if (mesh.top_simplex_type() == PrimitiveType::Tetrahedron &&
            mesh.has_attribute<int64_t>("img_tag", PrimitiveType::Tetrahedron)) {
            mesh.clear_attributes(
                {pos_handle,
                 top_simplex_tag_handle,
                 mesh.get_attribute_handle<int64_t>("img_tag", PrimitiveType::Tetrahedron)});
        } else {
            mesh.clear_attributes({pos_handle, top_simplex_tag_handle});
        }

        top_simplex_tag_handle = mesh.get_attribute_handle<int64_t>(
            "topological_offsets_top_simplex_tag",
            mesh.top_simplex_type());

        output.offset_position_handle =
            register_and_transfer_positions(pos_handle, *output.offset_mesh);

        sw.stop();
        output.report["timings_seconds"]["tog_init"] = sw.getElapsedTime();
    }

    Mesh& offset_mesh = *output.offset_mesh;
    attribute::MeshAttributeHandle offset_pos_handle = output.offset_position_handle;

    attribute::MeshAttributeHandle input_pos_handle =
        register_and_transfer_positions(pos_handle, input_mesh);

    multimesh::consolidate(mesh);
    offset_bvh->update_position_handle();
    offset_bvh->update_offset_distance_handle();
    offset_bvh->update_face_id_map();

    // check offset validity
    if (mesh.top_simplex_type() == PrimitiveType::Tetrahedron) {
        using wmtk::utils::wmtk_orient3d;

        for (const Tuple& tuple_face : offset_mesh.get_all(offset_mesh.top_simplex_type())) {
            const simplex::Simplex s(offset_mesh, offset_mesh.top_simplex_type(), tuple_face);

            const auto vertices =
                simplex::faces_single_dimension_tuples(offset_mesh, s, PrimitiveType::Vertex);
            assert(vertices.size() == 3);

            const auto offset_pos_acc =
                offset_mesh.create_const_accessor(offset_pos_handle.as<double>());

            const Eigen::Vector3d p0 = offset_pos_acc.const_vector_attribute(vertices[0]);
            const Eigen::Vector3d p1 = offset_pos_acc.const_vector_attribute(vertices[1]);
            const Eigen::Vector3d p2 = offset_pos_acc.const_vector_attribute(vertices[2]);

            const Eigen::Vector3d p_mid = (p0 + p1 + p2) / 3.;

            Eigen::Vector3d n_face = ((p1 - p0).cross(p2 - p0)).normalized();

            // check face normal orientation
            {
                Tuple tet_tuple = offset_mesh.map_to_parent_tuple(s);

                const auto tet_tag_acc =
                    mesh.create_const_accessor<int64_t>(top_simplex_tag_handle);
                // the tet tuple must be the one inside the offset (with a tag)
                if (tet_tag_acc.const_scalar_attribute(tet_tuple) == 0) {
                    tet_tuple = mesh.switch_tuple(tet_tuple, mesh.top_simplex_type());
                    if (tet_tag_acc.const_scalar_attribute(tet_tuple) == 0) {
                        log_and_throw_error("no tagged tet incident to offset");
                    }
                }
                // get opposite vertex
                Tuple opposite_vertex = mesh.switch_tuples(
                    tet_tuple,
                    {PrimitiveType::Triangle, PrimitiveType::Edge, PrimitiveType::Vertex});

                const auto tet_pos_acc = mesh.create_const_accessor<double>(pos_handle);
                const Eigen::Vector3d p3 = tet_pos_acc.const_vector_attribute(opposite_vertex);

                const Eigen::Vector3d p_normal = p_mid + n_face;
                if (wmtk_orient3d(p0, p1, p2, p3) != wmtk_orient3d(p0, p1, p2, p_normal)) {
                    n_face *= -1;
                }
                // n_face points now in the direction of the input
            }

            auto [sq_dist, offset_distance, nearest_point] =
                offset_bvh->sq_dist_offset_distance_nearest_point(p_mid);

            const Eigen::Vector3d nearest_normal = (nearest_point - p_mid).normalized();

            if (n_face.dot(nearest_normal) < 0) {
                logger().warn(
                    "Bad face orientation in initial offset (n * n_proj = {}). Consider refining "
                    "the embedding.",
                    n_face.dot(nearest_normal));
                break;
            }
        }
    }

    // optimize offset
    {
        logger().info("Optimize infinitesimal offset.");

        // log input distance
        {
            double sq_dist;
            SimpleBVH::VectorMax3d nearest_point;

            auto p_acc = offset_mesh.create_const_accessor<double>(offset_pos_handle);

            double dist_avg = 0;
            const auto vertices = offset_mesh.get_all(PrimitiveType::Vertex);

            for (const Tuple& t : vertices) {
                const simplex::Simplex v(offset_mesh, PrimitiveType::Vertex, t);
                const auto p = p_acc.const_vector_attribute(v);

                offset_bvh->bvh().nearest_facet(p, nearest_point, sq_dist);
                dist_avg += std::sqrt(sq_dist);
            }
            dist_avg /= vertices.size();

            wmtk::logger().info("Initial offset distance: {}", dist_avg);
            output.report["optimization"]["initial_offset_distance"] = dist_avg;
        }

        internal::OffsetOptimization oo(mesh, pos_handle);
        oo.set_offset_mesh(
            output.offset_mesh,
            options.input_mesh,
            offset_pos_handle,
            input_pos_handle);
        oo.set_input_bvh(offset_bvh);
        oo.set_offset_distance(options.distance);
        oo.set_edge_length_constraints(
            options.min_edge_length.value(),
            options.max_edge_length.value());
        oo.set_normal_deviation(options.min_normal_deviation, options.max_normal_deviation);
        oo.set_tetrahedron_tag_handle(top_simplex_tag_handle);
        oo.set_debug_prints(
            options.DEBUG_print_embedding,
            options.DEBUG_print_offset,
            options.DEBUG_print_smooth);

        oo.init_embedding_optimization();
        oo.init_offset_optimization();

        utils::write_mesh(pos_handle, fmt::format("embedding_{}", 0), options.intermediate_output);
        utils::write_mesh(
            input_pos_handle,
            fmt::format("input_{}", 0),
            options.intermediate_output);
        utils::write_mesh(
            offset_pos_handle,
            fmt::format("offset_{}", 0),
            options.intermediate_output);
        wmtk::utils::StopWatch sw("offset optimization");

        auto errs = oo.get_error_metrics(); // dist, nd
        components::internal::utils::ConvergenceRate convergence_rates(
            errs,
            options.convergence_max);

        for (size_t i = 0; i < options.passes; ++i) {
            logger().info(
                "==================== Outer loop {}/{} ====================",
                i + 1,
                options.passes);

            if (convergence_rates.is_converged()) {
                logger().info(">>>>>>>>>>>>>>>>>>>> offset converged <<<<<<<<<<<<<<<<<<<<");
                break;
            }

            oo.smooth_all(5);

            oo.optimize_offset(options.offset_passes);
            utils::write_mesh(
                offset_pos_handle,
                fmt::format("offset_{}", i + 1),
                options.intermediate_output);
            oo.optimize_embedding(options.embedding_passes);
            // utils::write_mesh(
            //     pos_handle,
            //     fmt::format("embedding_{}", i + 1),
            //     options.intermediate_output);

            errs = oo.get_error_metrics(); // dist, nd
            convergence_rates.update_vals(errs);
        }
        logger().info("Finished infinitesimal offset optimization.");

        sw.stop();

        output.report["timings_seconds"]["bvh_init"] = sw.getElapsedTime();
        output.report["optimization_metrics"] = oo.metrics_json();
    }

    //// reconstruct tags from substructures
    //{
    //    tag_handles = get_tag_attributes(mesh, options);
    //    for (const auto& [tag, substr] : substructures) {
    //        if (tag == -1) {
    //            continue;
    //        }
    //        const PrimitiveType pt = substr->top_simplex_type();
    //        auto acc = mesh.create_accessor(tag_handles[pt].as<int64_t>());
    //        for (const Tuple& t : substr->get_all(pt)) {
    //            const Tuple tt = substr->map_to_parent_tuple(simplex::Simplex(*substr, pt, t));
    //            acc.scalar_attribute(tt) = tag;
    //        }
    //    }

    //    const auto children = mesh.get_all_child_meshes();
    //    for (const auto& sub : children) {
    //        mesh.deregister_child_mesh(sub);
    //    }
    //    utils::write_mesh(
    //        pos_handle,
    //        fmt::format("output_with_tags_{}", 0),
    //        options.intermediate_output);
    //}

    //// delete all attributes besides tags and position
    //{
    //    std::vector<attribute::MeshAttributeHandle> keeps{
    //        {pos_handle,
    //         tag_handles[PrimitiveType::Tetrahedron],
    //         tag_handles[PrimitiveType::Triangle],
    //         tag_handles[PrimitiveType::Edge],
    //         tag_handles[PrimitiveType::Vertex]}};
    //    mesh.clear_attributes(keeps);
    //}

    //// output
    // cache.write_mesh(mesh, options.output);

    return output;
}
} // namespace wmtk::components
