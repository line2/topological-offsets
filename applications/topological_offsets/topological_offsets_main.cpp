#include <jse/jse.h>
#include <CLI/CLI.hpp>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <vector>

#include <wmtk/Mesh.hpp>
#include <wmtk/operations/attribute_update/AttributeTransferStrategy.hpp>
#include <wmtk/simplex/faces_single_dimension.hpp>
#include <wmtk/utils/Logger.hpp>

#include <wmtk/components/input/input.hpp>
#include <wmtk/components/multimesh_from_tag/multimesh_from_tag.hpp>
#include <wmtk/components/output/output.hpp>
#include <wmtk/components/topological_offsets/internal/utils/OffsetBvh.hpp>
#include <wmtk/components/topological_offsets/topological_offsets.hpp>
#include <wmtk/components/utils/json_utils.hpp>
#include <wmtk/components/utils/resolve_path.hpp>
#include <wmtk/utils/Stopwatch.hpp>
#include <wmtk/utils/primitive_range.hpp>

#include "utils/igl_input.hpp"

#include "topological_offsets_spec.hpp"

using namespace wmtk;
namespace fs = std::filesystem;

using wmtk::components::utils::resolve_paths;

constexpr int WMTK_FAIL = 1;

PrimitiveType find_tag_primitive_type(
    std::map<PrimitiveType, attribute::MeshAttributeHandle> tag_handles,
    const int64_t tag_value)
{
    const Mesh& m = tag_handles.begin()->second.mesh();

    for (const PrimitiveType pt : wmtk::utils::primitive_below(m.top_simplex_type())) {
        if (tag_handles.count(pt) == 0) {
            continue;
        }
        auto tag_acc = m.create_const_accessor(tag_handles[pt].as<int64_t>());
        for (const Tuple& t : m.get_all(pt)) {
            if (tag_acc.const_scalar_attribute(t) == tag_value) {
                return pt;
            }
        }
    }

    log_and_throw_error("Cannot generate offset. Tag {} not found.", tag_value);
}

/**
 * @brief Create a child mesh from the boundary.
 * This function registers a new attribute. Call `mesh.clear_attributes()` afterwards.
 */
void make_boundary_child_mesh(std::shared_ptr<Mesh>& mesh)
{
    auto tag_attr = mesh->register_attribute<int64_t>(
        "topological_offsets_boundary_tag",
        get_primitive_type_from_id(mesh->top_cell_dimension() - 1),
        1);
    auto acc = mesh->create_accessor<int64_t>(tag_attr);
    for (const Tuple& t : mesh->get_all(tag_attr.primitive_type())) {
        if (mesh->is_boundary(tag_attr.primitive_type(), t)) {
            acc.scalar_attribute(t) = 1;
        }
    }

    components::multimesh_from_tag(mesh, tag_attr, 1);
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

int main(int argc, char* argv[])
{
    CLI::App app{argv[0]};

    app.ignore_case();

    fs::path json_input_file;
    app.add_option("-j, --json", json_input_file, "json specification file")
        ->required(true)
        ->check(CLI::ExistingFile);
    CLI11_PARSE(app, argc, argv);

    nlohmann::json j;
    {
        std::ifstream ifs(json_input_file);
        j = nlohmann::json::parse(ifs);

        jse::JSE spec_engine;
        bool r = spec_engine.verify_json(j, topological_offsets_spec);
        if (!r) {
            logger().error("{}", spec_engine.log2str());
            return WMTK_FAIL;
        } else {
            j = spec_engine.inject_defaults(j, topological_offsets_spec);
        }
    }

    // logger settings
    {
        std::string log_file_name = j["log_file"];
        if (!log_file_name.empty()) {
            wmtk::set_file_logger(log_file_name);
        }

        if (!wmtk::has_user_overloaded_logger_level()) {
            wmtk::opt_logger().set_level(j["opt_log_level"]);
        }
    }

    utils::StopWatch sw_total("Topological Offsets Total");

    const std::string pos_attr_name = "vertices";

    std::map<PrimitiveType, std::string> tag_attr_names;
    tag_attr_names[PrimitiveType::Tetrahedron] = "tet_tag";
    tag_attr_names[PrimitiveType::Triangle] = "face_tag";
    tag_attr_names[PrimitiveType::Edge] = "edge_tag";
    tag_attr_names[PrimitiveType::Vertex] = "vertex_tag";

    const int64_t input_tag = 1;
    const int64_t embedding_tag = 0;


    const fs::path input_file = resolve_paths(json_input_file, {j["input_path"], j["input"]});
    const fs::path output_file = std::string(j["output"]).empty()
                                     ? (json_input_file.stem().string() + "_to")
                                     : std::string(j["output"]);

    std::map<PrimitiveType, attribute::MeshAttributeHandle> tag_handles;

    std::shared_ptr<Mesh> mesh_in;
    if (input_file.extension() == ".msh") {
        if (j["tag_input"] == "") {
            logger().info("Expecting file to be from tetwild and containing the in_out "
                          "attribute. Converting it into a tag");
            mesh_in = components::input::input(input_file, false, {"in_out"});
            Mesh& mesh = *mesh_in;

            auto in_out_handle =
                mesh.get_attribute_handle<double>("in_out", PrimitiveType::Tetrahedron);
            auto in_out_acc = mesh.create_accessor<double>(in_out_handle);

            tag_handles[mesh.top_simplex_type()] = mesh.register_attribute<int64_t>(
                tag_attr_names[mesh.top_simplex_type()],
                mesh.top_simplex_type(),
                1,
                false,
                embedding_tag);
            auto tag_acc = mesh.create_accessor<int64_t>(tag_handles[mesh.top_simplex_type()]);

            bool has_tag = false;

            for (const Tuple& t : mesh.get_all(mesh.top_simplex_type())) {
                if (in_out_acc.const_scalar_attribute(t) < 0.5) {
                    tag_acc.scalar_attribute(t) = input_tag;
                    has_tag = true;
                }
            }

            if (!has_tag) {
                logger().error("Cannot generate offset. Input tag not found.");
                return WMTK_FAIL;
            }

            // check if boundary tets are tagged
            for (const Tuple& t : mesh.get_all(PrimitiveType::Tetrahedron)) {
                const simplex::Simplex tet(mesh, PrimitiveType::Tetrahedron, t);

                if (tag_acc.const_scalar_attribute(tet) != 1) {
                    continue;
                }

                bool tet_is_boundary = false;
                const auto vertices =
                    simplex::faces_single_dimension(mesh, tet, PrimitiveType::Vertex);
                for (const simplex::Simplex& v : vertices) {
                    if (mesh.is_boundary(v)) {
                        logger().error(
                            "Cannot generate offset. Input contains tag at boundary tets.");
                        return WMTK_FAIL;
                    }
                }
            }


            auto pos_handle =
                mesh_in->get_attribute_handle<double>(pos_attr_name, PrimitiveType::Vertex);
            mesh.clear_attributes({tag_handles[mesh.top_simplex_type()], pos_handle});
        } else {
            logger().info("A tag input was given.");

            mesh_in = components::input::input(input_file, j["2d"]);
            Mesh& mesh = *mesh_in;

            auto pos_handle =
                mesh.get_attribute_handle<double>(pos_attr_name, PrimitiveType::Vertex);

            const fs::path tag_input_file =
                resolve_paths(json_input_file, {j["input_path"], j["tag_input"]});
            std::shared_ptr<Mesh> tag_mesh_in;
            if (tag_input_file.extension() == ".msh") {
                tag_mesh_in = components::input::input(tag_input_file, j["2d"]);
            } else {
                logger().info("Use igl input for tag mesh");
                tag_mesh_in = utils::igl_input(tag_input_file);
                // components::output::output(*tag_mesh_in, "tag_mesh", pos_attr_name);
            }
            Mesh& tag_mesh = *tag_mesh_in;
            auto tag_pos_handle =
                tag_mesh.get_attribute_handle<double>(pos_attr_name, PrimitiveType::Vertex);


            tag_handles[tag_mesh.top_simplex_type()] = mesh.register_attribute<int64_t>(
                tag_attr_names[tag_mesh.top_simplex_type()],
                tag_mesh.top_simplex_type(),
                1,
                false,
                embedding_tag);
            auto tag_acc = mesh.create_accessor<int64_t>(tag_handles[tag_mesh.top_simplex_type()]);


            auto get_tet_center = [&tag_mesh, &tag_pos_handle](const Tuple& _t) -> Eigen::VectorXd {
                auto pos_acc = tag_mesh.create_const_accessor<double>(tag_pos_handle);

                Eigen::VectorXd mid_point;
                mid_point.setZero(pos_acc.dimension());

                const auto vertices = tag_mesh.orient_vertices(_t);
                for (const Tuple& v : vertices) {
                    mid_point += pos_acc.const_vector_attribute(v);
                }
                mid_point /= vertices.size();

                return mid_point;
            };

            // use geometry to find for each top simplex the corresponding simplex in the actual
            // mesh
            components::internal::utils::OffsetBvh bvh(pos_handle, tag_mesh.top_simplex_type());

            const auto tag_tets = tag_mesh.get_all(tag_mesh.top_simplex_type());
            if (tag_tets.empty()) {
                log_and_throw_error("Cannot generate offset. Tag mesh is empty.");
            }
            for (const Tuple& t : tag_tets) {
                Eigen::VectorXd p = get_tet_center(t);

                auto [sq_dist, tt] = bvh.sq_dist_and_tuple(p);
                tag_acc.scalar_attribute(tt) = input_tag;
            }
            // components::output::output(mesh, "mesh_tagged", pos_attr_name);
        }
    } else {
        mesh_in = components::input::input(input_file);

        Mesh& m = *mesh_in;

        for (const PrimitiveType pt : wmtk::utils::primitive_below(m.top_simplex_type())) {
            if (m.has_attribute<int64_t>(tag_attr_names[pt], pt)) {
                tag_handles[pt] = m.get_attribute_handle<int64_t>(tag_attr_names[pt], pt);
            }
        }

        auto pos_handle = m.get_attribute_handle<double>(pos_attr_name, PrimitiveType::Vertex);
        m.clear_attributes({tag_handles[m.top_simplex_type()], pos_handle});

        if (tag_handles.empty()) {
            log_and_throw_error("Tag attributes could not be found in input mesh.");
        }
    }

    Mesh& mesh = *mesh_in;

    if (j["invert_tags"]) {
        auto tag_acc = mesh.create_accessor<int64_t>(tag_handles[mesh.top_simplex_type()]);

        for (const Tuple& t : mesh.get_all(mesh.top_simplex_type())) {
            if (tag_acc.const_scalar_attribute(t) == input_tag) {
                tag_acc.scalar_attribute(t) = embedding_tag;
            } else if (tag_acc.const_scalar_attribute(t) == embedding_tag) {
                tag_acc.scalar_attribute(t) = input_tag;
            }
        }
    }

    if (mesh.top_simplex_type() == PrimitiveType::Triangle) {
        tag_attr_names.erase(PrimitiveType::Tetrahedron);
    }

    attribute::MeshAttributeHandle embedding_pos_handle =
        mesh.get_attribute_handle<double>(pos_attr_name, PrimitiveType::Vertex);

    std::map<int64_t, std::shared_ptr<Mesh>> substructures;

    // convert tags into child-meshes
    {
        // collect all tags
        std::set<int64_t> tag_values;
        for (const auto& [pt, tag_handle] : tag_handles) {
            const auto acc = mesh.create_const_accessor<int64_t>(tag_handle);
            for (const Tuple t : mesh.get_all(pt)) {
                tag_values.insert(acc.const_scalar_attribute(t));
            }
        }

        // create child-mesh for each tag
        for (const int64_t tag_value : tag_values) {
            logger().info("Generate substructures from tag {}.", tag_value);

            const PrimitiveType pt = find_tag_primitive_type(tag_handles, tag_value);

            // create substructure
            substructures[tag_value] =
                components::multimesh_from_tag(mesh_in, tag_handles[pt], tag_value);
        }
    }

    mesh.clear_attributes({embedding_pos_handle});

    // make child-mesh from boundary
    make_boundary_child_mesh(mesh_in);
    mesh.clear_attributes({embedding_pos_handle});


    std::shared_ptr<Mesh> input_mesh;
    std::shared_ptr<Mesh> inside_mesh;

    // use boundary only if input consists of top-dimension simplices
    if (substructures[input_tag]->top_simplex_type() == mesh.top_simplex_type()) {
        Mesh& mm = *substructures[input_tag];
        auto tag_attr = mesh.register_attribute<int64_t>(
            "topological_offsets_input_boundary_tag",
            get_primitive_type_from_id(mesh.top_cell_dimension() - 1),
            1);
        auto acc = mesh.create_accessor<int64_t>(tag_attr);
        for (const Tuple& t : mm.get_all(tag_attr.primitive_type())) {
            if (mm.is_boundary(tag_attr.primitive_type(), t)) {
                const Tuple t_parent =
                    mm.map_to_parent_tuple(simplex::Simplex(mm, tag_attr.primitive_type(), t));
                acc.scalar_attribute(t_parent) = 1;
            }
        }

        input_mesh = components::multimesh_from_tag(mesh_in, tag_attr, 1);
        mesh.clear_attributes({embedding_pos_handle});

        inside_mesh = substructures[input_tag];
    } else {
        input_mesh = substructures[input_tag];
    }

    const std::string side = j["side"];
    if (side != "out" && !inside_mesh) {
        log_and_throw_error(
            "Cannot perform offset on side '{}' if the input mesh is not top dimension.",
            side);
    }
    if (side == "out") {
        mesh.deregister_child_mesh(substructures[embedding_tag]);
    }
    if (side == "in") {
        mesh.deregister_child_mesh(inside_mesh);
        inside_mesh = substructures[embedding_tag];
    }
    if (side == "double" && inside_mesh) {
        mesh.deregister_child_mesh(inside_mesh);
        inside_mesh.reset();
        mesh.deregister_child_mesh(substructures[embedding_tag]);
    }

    mesh.consolidate();

    nlohmann::json out_json;

    // topological offsets
    if (j["distance"].is_number()) {
        const double distance = j["distance"];

        components::TopologicalOffsetsOutput tog_output;

        {
            components::TopologicalOffsetsOptions options;
            options.embedding_position_handle = embedding_pos_handle;
            options.input_mesh = input_mesh;
            options.inside_mesh = inside_mesh;
            options.distance = std::abs(distance);

            options.passes = j["passes"];
            options.embedding_passes = j["embedding_passes"];
            options.offset_passes = j["offset_passes"];

            if (!j.at("min_edge_length").is_null()) {
                options.min_edge_length = j.at("min_edge_length");
            }
            if (!j.at("max_edge_length").is_null()) {
                options.max_edge_length = j.at("max_edge_length");
            }
            options.min_normal_deviation = j["min_normal_deviation"];
            options.max_normal_deviation = j["max_normal_deviation"];
            options.convergence_max = j["convergence_max"];
            options.distance_function = j["distance_function"];
            options.finite_offset = j["finite_offset"];
            options.relative_distance_and_length = j["relative_distance_and_length"];
            options.intermediate_output = j["intermediate_output"];
            options.use_simplicial_embedding = j["use_simplicial_embedding"];
            options.use_warm_start = j["use_warm_start"];
            options.expand_topo = j["expand_topo"];
            options.adapt_offset_distance = j["adapt_offset_distance"];
            options.restrict_min_edge_length_to_input_avg =
                j["restrict_min_edge_length_to_input_avg"];
            options.DEBUG_print_embedding = j["DEBUG_print_embedding"];
            options.DEBUG_print_offset = j["DEBUG_print_offset"];
            options.DEBUG_print_smooth = j["DEBUG_print_smooth"];

            utils::StopWatch sw("Topological Offsets");
            tog_output = components::topological_offsets(options);
        }

        wmtk::components::output::output(
            *tog_output.offset_mesh,
            output_file,
            tog_output.offset_position_handle);

        out_json["stats"]["vertices"] = mesh.get_all(PrimitiveType::Vertex).size();
        out_json["stats"]["edges"] = mesh.get_all(PrimitiveType::Edge).size();
        out_json["stats"]["triangles"] = mesh.get_all(PrimitiveType::Triangle).size();
        out_json["stats"]["tets"] = mesh.get_all(PrimitiveType::Tetrahedron).size();
        out_json["stats"]["tog"] = tog_output.report;

    } else {
        if (!j["distance"].is_array()) {
            log_and_throw_error("Distance is neither a number nor an array.");
        }
        const std::vector<double> distances = j["distance"];

        if (!std::is_sorted(distances.begin(), distances.end(), std::greater<>())) {
            log_and_throw_error("Distances must be sorted in ascending order.");
        }

        for (size_t i = 0; i < distances.size(); ++i) {
            components::TopologicalOffsetsOutput tog_output;

            mesh.clear_attributes({embedding_pos_handle});
            embedding_pos_handle =
                mesh.get_attribute_handle<double>("vertices", PrimitiveType::Vertex);
            for (auto& child : mesh.get_all_child_meshes()) {
                child->clear_attributes();
            }

            {
                components::TopologicalOffsetsOptions options;
                options.embedding_position_handle = embedding_pos_handle;
                options.input_mesh = input_mesh;
                options.inside_mesh = inside_mesh;
                options.distance = std::abs(distances[i]);

                options.passes = j["passes"];
                options.embedding_passes = j["embedding_passes"];
                options.offset_passes = j["offset_passes"];

                if (!j.at("min_edge_length").is_null()) {
                    options.min_edge_length = j.at("min_edge_length");
                }
                if (!j.at("max_edge_length").is_null()) {
                    options.max_edge_length = j.at("max_edge_length");
                }
                options.finite_offset = j["finite_offset"];
                options.relative_distance_and_length = j["relative_distance_and_length"];
                options.intermediate_output = j["intermediate_output"];
                options.use_simplicial_embedding = j["use_simplicial_embedding"];
                options.use_warm_start = j["use_warm_start"];
                options.adapt_offset_distance = j["adapt_offset_distance"];
                options.restrict_min_edge_length_to_input_avg =
                    j["restrict_min_edge_length_to_input_avg"];

                utils::StopWatch sw("Topological Offsets");
                tog_output = components::topological_offsets(options);
            }

            wmtk::components::output::output(
                *tog_output.offset_mesh,
                fmt::format("{}_{}", output_file.string(), i),
                tog_output.offset_position_handle);

            out_json["stats"][fmt::format("vertices_{}", i)] =
                mesh.get_all(PrimitiveType::Vertex).size();
            out_json["stats"][fmt::format("edges_{}", i)] =
                mesh.get_all(PrimitiveType::Edge).size();
            out_json["stats"][fmt::format("triangles_{}", i)] =
                mesh.get_all(PrimitiveType::Triangle).size();
            out_json["stats"][fmt::format("tets_{}", i)] =
                mesh.get_all(PrimitiveType::Tetrahedron).size();
            out_json["stats"][fmt::format("tog_{}", i)] = tog_output.report;
        }
    }

    // input output
    wmtk::components::output::output(
        *input_mesh,
        fmt::format("{}_input", output_file.string()),
        "vertices");


    // embedding output
    {
        // add tet tags
        logger().info("Add tet tags");
        const PrimitiveType pt_top = mesh.top_simplex_type();
        const PrimitiveType pt_face = get_primitive_type_from_id(get_primitive_type_id(pt_top) - 1);

        const auto tets = mesh.get_all(pt_top);
        logger().info("#tets = {}", tets.size());
        auto offset_tag_handle = mesh.register_attribute<int64_t>("offset_tags", pt_top, 1);
        auto tag_acc = mesh.create_accessor<int64_t>(offset_tag_handle);

        if (inside_mesh) {
            logger().info("#inside tets = {}", inside_mesh->get_all(pt_top).size());
            for (const Tuple& tet_tuple : inside_mesh->get_all(pt_top)) {
                const simplex::Simplex tet(*inside_mesh, pt_top, tet_tuple);
                const Tuple parent_tuple = inside_mesh->map_to_parent_tuple(tet);
                tag_acc.scalar_attribute(parent_tuple) = -1;
            }
        }

        int64_t tag_counter = 0;
        for (const Tuple& tet_tuple : tets) {
            if (tag_acc.const_scalar_attribute(tet_tuple) != 0) {
                continue;
            }
            const simplex::Simplex tet(mesh, pt_top, tet_tuple);

            // tet is neither in a child nor tagged --> add tag and flood fill
            ++tag_counter;

            logger().info("Adding tet tag {}", tag_counter);

            std::vector<Tuple> q(200);
            size_t q_front = 0;
            size_t q_back = 0;

            q[q_back++] = tet_tuple;

            while (q_front != q_back) {
                const Tuple q_tuple = q[q_front++];
                const simplex::Simplex q_simplex(mesh, pt_top, q_tuple);

                const auto faces = simplex::faces_single_dimension(mesh, q_simplex, pt_face);
                for (const simplex::Simplex& face : faces) {
                    if (mesh.is_boundary(face)) {
                        continue;
                    }
                    bool is_in_child = false;
                    for (auto& child : mesh.get_all_child_meshes()) {
                        if (mesh.simplex_is_in_child(*child, face)) {
                            is_in_child = true;
                            break;
                        }
                    }
                    if (is_in_child) {
                        continue;
                    }

                    const Tuple ts = mesh.switch_tuple(face.tuple(), pt_top);
                    const simplex::Simplex ts_simplex(mesh, pt_top, ts);
                    if (tag_acc.const_scalar_attribute(ts) != 0) {
                        continue;
                    }
                    tag_acc.scalar_attribute(ts) = tag_counter;
                    if (q_back + 1 == q.size()) {
                        q.resize(1.5 * q.size());
                    }
                    q[q_back++] = ts;
                }
            }
        }

        wmtk::components::output::output(
            mesh,
            fmt::format("{}_embedding", output_file.string()),
            embedding_pos_handle);
    }

    sw_total.stop();

    out_json["input"] = j;
    out_json["stats"]["total_runtime_seconds"] = sw_total.getElapsedTime();

    const std::string report = j["report"];
    if (!report.empty()) {
        std::ofstream ofs(report);
        ofs << std::setw(4) << out_json;
    }


    return 0;
}
