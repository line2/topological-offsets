#include "TopologicalOffsetGenerator.hpp"

#include <igl/harmonic.h>
#include <set>
#include <wmtk/io/ParaviewWriter.hpp>
#include <wmtk/simplex/IdSimplexCollection.hpp>
#include <wmtk/simplex/faces.hpp>
#include <wmtk/simplex/faces_single_dimension.hpp>
#include <wmtk/simplex/internal/VisitedArray.hpp>
#include <wmtk/simplex/link_single_dimension_iterable.hpp>
#include <wmtk/simplex/top_dimension_cofaces.hpp>
#include <wmtk/simplex/top_dimension_cofaces_iterable.hpp>
#include <wmtk/utils/Logger.hpp>
#include <wmtk/utils/TupleInspector.hpp>
#include <wmtk/utils/bbox_from_mesh.hpp>
#include <wmtk/utils/primitive_range.hpp>

#include <wmtk/components/marching/marching.hpp>
#include <wmtk/components/multimesh_from_tag/multimesh_from_tag.hpp>
#include <wmtk/components/regular_space/regular_space.hpp>

#include "utils/bvh_to_mesh_index_map.hpp"
#include "utils/is_tet_in_offset.hpp"
#include "utils/tet_conains_connected_component.hpp"
#include "utils/write_mesh.hpp"

namespace wmtk::components {

TopologicalOffsetGenerator::TopologicalOffsetGenerator(
    Mesh& mesh,
    attribute::MeshAttributeHandle& pos_handle,
    const double offset_distance,
    bool debug_print)
    : m_mesh(mesh)
    , m_pos_handle(pos_handle)
    , m_offset_distance(offset_distance)
    , m_debug_print(debug_print)
    , m_pt_top(mesh.top_simplex_type())
    , m_pt_face(get_primitive_type_from_id(get_primitive_type_id(mesh.top_simplex_type()) - 1))
{
    assert(
        m_mesh.top_simplex_type() == PrimitiveType::Tetrahedron ||
        m_mesh.top_simplex_type() == PrimitiveType::Triangle);

    if (mesh.top_simplex_type() == PrimitiveType::Tetrahedron) {
        m_tag_handles[PrimitiveType::Tetrahedron] =
            m_mesh.register_attribute<int64_t>("tog_tet_tags", PrimitiveType::Tetrahedron, 1);
    }
    m_tag_handles[PrimitiveType::Triangle] =
        m_mesh.register_attribute<int64_t>("tog_tri_tags", PrimitiveType::Triangle, 1);
    m_tag_handles[PrimitiveType::Edge] =
        m_mesh.register_attribute<int64_t>("tog_edge_tags", PrimitiveType::Edge, 1);
    m_tag_handles[PrimitiveType::Vertex] =
        m_mesh.register_attribute<int64_t>("tog_vertex_tags", PrimitiveType::Vertex, 1);

    m_offset_tag_handle = m_mesh.register_attribute<int64_t>(
        "topological_offsets_top_simplex_tag",
        m_mesh.top_simplex_type(),
        1);

    if (m_debug_print) {
        logger().info("Constructing TopologicalOffsetGenerator");
        const std::string name = fmt::format("tog_debug_{}", m_debug_print_counter++);
        utils::write_mesh(m_pos_handle, name, m_debug_print);
    }
}

void TopologicalOffsetGenerator::add_tags(attribute::MeshAttributeHandle& tag_handle)
{
    log_and_throw_error("Function add_tags() is deprecated");
    const PrimitiveType pt = tag_handle.primitive_type();
    const auto tag_acc = m_mesh.create_const_accessor(tag_handle.as<int64_t>());
    auto inside_tag_acc = m_mesh.create_accessor(m_tag_handles[pt].as<int64_t>());

    for (const Tuple& t : m_mesh.get_all(pt)) {
        if (tag_acc.const_scalar_attribute(t) > 0) {
            inside_tag_acc.scalar_attribute(t) = m_inside_tag;
        }
    }
}

void TopologicalOffsetGenerator::add_input_substructure(Mesh& child)
{
    const PrimitiveType pt_top = child.top_simplex_type();

    for (const PrimitiveType pt : wmtk::utils::primitive_below(pt_top)) {
        auto inside_tag_acc = m_mesh.create_accessor(m_tag_handles[pt].as<int64_t>());

        for (const Tuple& t_child : child.get_all(pt)) {
            const Tuple t_parent = child.map_to_parent_tuple(simplex::Simplex(child, pt, t_child));
            inside_tag_acc.scalar_attribute(t_parent) = m_inside_tag;
        }
    }
}

void TopologicalOffsetGenerator::regularize(bool use_simplicial_embedding)
{
    std::vector<attribute::MeshAttributeHandle> tag_vec;
    tag_vec.reserve(m_tag_handles.size());
    for (const PrimitiveType pt : wmtk::utils::primitive_below(m_mesh.top_simplex_type())) {
        tag_vec.emplace_back(m_tag_handles[pt]);
    }

    RegularSpaceOptions options;
    options.pass_through_attributes.emplace_back(m_pos_handle);
    options.pass_through_attributes.emplace_back(m_offset_tag_handle);
    options.pass_through_attributes.emplace_back(m_bvh->f_offset_distance_handle());
    options.pass_through_attributes.emplace_back(m_bvh->v_offset_distance_handle());
    assert(
        m_mesh.top_simplex_type() == PrimitiveType::Tetrahedron &&
        m_mesh.has_attribute<int64_t>("img_tag", PrimitiveType::Tetrahedron));
    options.pass_through_attributes.emplace_back(
        m_mesh.get_attribute_handle<int64_t>("img_tag", PrimitiveType::Tetrahedron));

    options.tag_attributes = m_tag_handles;
    options.value = m_inside_tag;
    options.generate_simplicial_embedding = use_simplicial_embedding;

    regular_space(m_mesh, options);

    if (m_debug_print && use_simplicial_embedding) {
        logger().info("Performed simplicial embedding.");
        const std::string name = fmt::format("tog_debug_{}", m_debug_print_counter++);
        utils::write_mesh(m_pos_handle, name, m_debug_print);
    }
}

void TopologicalOffsetGenerator::transfer_outside(
    const std::map<PrimitiveType, attribute::MeshAttributeHandle>& tag_handles)
{
    for (const auto& [pt, tag_handle] : tag_handles) {
        auto tag_acc = m_mesh.create_accessor<int64_t>(tag_handle);
        auto tog_tag_acc = m_mesh.create_const_accessor<int64_t>(m_tag_handles[pt]);

        for (const Tuple& t : m_mesh.get_all(pt)) {
            if (tog_tag_acc.const_scalar_attribute(t) == m_outside_tag) {
                tag_acc.scalar_attribute(t) = m_outside_tag;
            }
        }
    }
}

void TopologicalOffsetGenerator::register_bvh(
    const std::shared_ptr<internal::utils::OffsetBvh>& bvh)
{
    m_bvh = bvh;
}

void TopologicalOffsetGenerator::tag_tets_as_inside(
    bool keep_topology,
    bool use_warm_start,
    bool adapt_offset_distance)
{
    log_and_throw_error("deprecated");

    internal::utils::OffsetBvh& bvh = *m_bvh;
    const Eigen::MatrixXd bbox = wmtk::utils::bbox_from_mesh(m_pos_handle);

    auto within_od_handle =
        m_mesh.register_attribute<int64_t>("within_offset_distance", m_pt_top, 1);
    auto tet_touches_boundary_handle =
        m_mesh.register_attribute<int64_t>("tet_touches_boundary", m_pt_top, 1);
    auto visited_vertex_handle =
        m_mesh.register_attribute<int64_t>("tog_visited_vertex", PrimitiveType::Vertex, 1);
    auto v_tet_offset_distance_handle = m_mesh.register_attribute<double>(
        "tog_vertex_offset_distance",
        PrimitiveType::Vertex,
        1,
        false,
        m_offset_distance);
    auto laplacian_boundary_handle =
        m_mesh.register_attribute<int64_t>("laplacian_boundary", PrimitiveType::Vertex, 1);
    auto v_dist = m_mesh.register_attribute<double>("v_dist", PrimitiveType::Vertex, 1);
    auto v_inside = m_mesh.register_attribute<int64_t>("v_inside", PrimitiveType::Vertex, 1);

    Mesh& input_mesh = bvh.mesh();
    auto f_offset_distance_acc =
        input_mesh.create_accessor<double>(m_bvh->f_offset_distance_handle());
    auto v_offset_distance_acc =
        input_mesh.create_accessor<double>(m_bvh->v_offset_distance_handle());

    const auto pos_acc = m_mesh.create_const_accessor<double>(m_pos_handle);

    std::map<PrimitiveType, attribute::Accessor<int64_t>> tag_accs;
    for (const auto& [pt, tag_handle] : m_tag_handles) {
        tag_accs.emplace(pt, m_mesh.create_accessor<int64_t>(tag_handle));
    }
    auto offset_tag_acc = m_mesh.create_accessor<int64_t>(m_offset_tag_handle);

    auto within_od_acc = m_mesh.create_accessor<int64_t>(within_od_handle);
    auto tet_touches_boundary_acc = m_mesh.create_accessor<int64_t>(tet_touches_boundary_handle);

    for (const Tuple& t : m_mesh.get_all(m_pt_top)) {
        const simplex::Simplex tet(m_mesh, m_pt_top, t);
        if (tet_touches_boundary(tet)) {
            tet_touches_boundary_acc.scalar_attribute(t) = 1;
        } else {
            tet_touches_boundary_acc.scalar_attribute(t) = 0;
        }
    }

    auto visited_vertex_acc = m_mesh.create_accessor<int64_t>(visited_vertex_handle);
    auto v_tet_offset_distance_acc = m_mesh.create_accessor<double>(v_tet_offset_distance_handle);

    for (const Tuple& t : m_mesh.get_all(PrimitiveType::Vertex)) {
        const auto p = pos_acc.const_vector_attribute(t);
        auto [sq_dist, od] = bvh.sq_dist_and_offset_distance(p);
        v_tet_offset_distance_acc.scalar_attribute(t) = od;
    }

    auto update_offset_distance_jacobian = [&]() -> bool {
        const auto pos_acc = m_mesh.create_const_accessor<double>(m_pos_handle);
        const auto vertices = m_mesh.get_all(PrimitiveType::Vertex);

        auto laplacian_boundary_acc = m_mesh.create_accessor<int64_t>(laplacian_boundary_handle);

        // reset visisted tags
        for (const Tuple& t : vertices) {
            const auto p = pos_acc.const_vector_attribute(t);
            auto [sq_dist, od] = bvh.sq_dist_and_offset_distance(p);
            v_tet_offset_distance_acc.scalar_attribute(t) = od;
            visited_vertex_acc.scalar_attribute(t) = 0;
        }

        logger().info("Check for offset intersections with advancing front");
        bool found_bad_tet = false;

        // tag vertices that are outside
        for (const Tuple& tet_tuple : m_mesh.get_all(m_pt_top)) {
            const simplex::Simplex tet(m_mesh, m_pt_top, tet_tuple);

            if (tag_accs.at(m_pt_top).const_scalar_attribute(tet) != m_outside_tag) {
                continue;
            }

            double tet_dist = m_offset_distance;

            if (within_od_acc.const_scalar_attribute(tet) == 1) {
                // tet should be in offset but is not -> adapt offset distance

                //tag_accs.at(pt_top).scalar_attribute(tet) = 3; // DEBUG tag bad tets
                tet_dist = get_tet_small_dist(tet);
                found_bad_tet = true;
            }

            const auto tet_vertices = m_mesh.orient_vertices(tet_tuple);

            for (const Tuple& v_tuple : tet_vertices) {
                const simplex::IdSimplex v_ids =
                    m_mesh.get_id_simplex(v_tuple, PrimitiveType::Vertex);

                const double d = v_tet_offset_distance_acc.const_scalar_attribute(v_ids);
                v_tet_offset_distance_acc.scalar_attribute(v_ids) = std::min(d, tet_dist);
                visited_vertex_acc.scalar_attribute(v_ids) = 1;
            }
        }

        if (!found_bad_tet) {
            return false;
        }

        // bound difference on boundary vertices
        for (int64_t i = 0; i < 10; i++) {
            for (const Tuple& v_tuple : vertices) {
                const simplex::Simplex v(m_mesh, PrimitiveType::Vertex, v_tuple);
                if (visited_vertex_acc.const_scalar_attribute(v) != 1) {
                    continue;
                }

                auto neighs =
                    simplex::link_single_dimension_iterable(m_mesh, v, PrimitiveType::Vertex);
                double d_min = v_tet_offset_distance_acc.const_scalar_attribute(v);
                for (const Tuple& n : neighs) {
                    const double d = v_tet_offset_distance_acc.const_scalar_attribute(n);
                    d_min = std::min(d_min, 1.5 * d);
                }
                v_tet_offset_distance_acc.scalar_attribute(v) = d_min;
            }
        }

        // use igl::harmonic
        std::map<simplex::IdSimplex, int64_t> vids;
        const auto tets = m_mesh.get_all(m_pt_top);

        for (int64_t i = 0; i < vertices.size(); ++i) {
            const simplex::IdSimplex v = m_mesh.get_id_simplex(vertices[i], PrimitiveType::Vertex);
            vids[v] = i;
        }

        if (m_debug_print) {
            logger().info("Laplacian init");
            const std::string name = fmt::format("tog_debug_{}", m_debug_print_counter++);
            utils::write_mesh(m_pos_handle, name, m_debug_print);
            utils::write_mesh(m_pos_handle, name, PrimitiveType::Vertex, m_debug_print);
        }

        Eigen::MatrixXd V; // #V by dim vertex positions
        Eigen::MatrixXi F; // #F by simplex-size list of element indices
        Eigen::VectorXi b; // #b boundary indices into V
        Eigen::MatrixXd bc; // #b by #W list of boundary values
        const int64_t k = 1; // power of harmonic operation (1: harmonic, 2: biharmonic, etc)
        Eigen::MatrixXd W; // #V by #W list of weights

        // V
        V.resize(vertices.size(), pos_acc.dimension());
        std::set<int> boundary_ids;
        for (int64_t i = 0; i < vertices.size(); ++i) {
            V.row(i) = pos_acc.const_vector_attribute(vertices[i]);
            if (visited_vertex_acc.const_scalar_attribute(vertices[i]) == 1) {
                laplacian_boundary_acc.scalar_attribute(vertices[i]) = 1;
                boundary_ids.insert(i);
                continue;
            }

            // check if vertex is inside
            bool is_inside = false;
            auto neighs = simplex::top_dimension_cofaces_iterable(
                m_mesh,
                simplex::Simplex(m_mesh, PrimitiveType::Vertex, vertices[i]));
            for (const Tuple& t : neighs) {
                if (offset_tag_acc.const_scalar_attribute(t) == m_offset_tag) {
                    is_inside = true;
                    break;
                }
            }
            if (!is_inside) {
                laplacian_boundary_acc.scalar_attribute(vertices[i]) = 2;
                boundary_ids.insert(i);
            }
        }
        // b & bc
        {
            std::vector<int> bvec(boundary_ids.begin(), boundary_ids.end());
            b = Eigen::Map<Eigen::VectorXi>(bvec.data(), bvec.size());
            bc.resize(bvec.size(), 1);
            for (int64_t i = 0; i < bvec.size(); ++i) {
                bc(i, 0) = v_tet_offset_distance_acc.const_scalar_attribute(vertices[bvec[i]]);
            }
        }

        int64_t nF = 0;
        for (int64_t i = 0; i < tets.size(); ++i) {
            if (offset_tag_acc.const_scalar_attribute(tets[i]) == m_offset_tag) {
                ++nF;
            }
        }

        // F
        F.resize(nF, pos_acc.dimension() + 1);
        int64_t f_counter = 0;
        for (int64_t i = 0; i < tets.size(); ++i) {
            if (offset_tag_acc.const_scalar_attribute(tets[i]) != m_offset_tag) {
                continue;
            }
            const auto tet_vertices = m_mesh.orient_vertices(tets[i]);
            for (int64_t j = 0; j < tet_vertices.size(); ++j) {
                const simplex::IdSimplex v =
                    m_mesh.get_id_simplex(tet_vertices[j], PrimitiveType::Vertex);
                F(f_counter, j) = vids[v];
            }
            ++f_counter;
        }

        logger().info("Solve igl::harmonic...");
        if (!igl::harmonic(V, F, b, bc, k, W)) {
            log_and_throw_error("Could not solve igl::harmonic!");
        }
        logger().info("... done");
        assert(W.rows() == V.rows());
        assert(W.cols() == bc.cols());

        for (int64_t i = 0; i < vertices.size(); ++i) {
            v_tet_offset_distance_acc.scalar_attribute(vertices[i]) = W(i, 0);
        }


        // transfer values from tet vertices to tri faces and vertices
        for (const Tuple& face_tuple : input_mesh.get_all(m_pt_face)) {
            const auto child_vertices = input_mesh.orient_vertices(face_tuple);
            double min_dist = f_offset_distance_acc.const_scalar_attribute(face_tuple);
            for (const Tuple& child_vertex : child_vertices) {
                const simplex::Simplex child_vertex_simplex(
                    input_mesh,
                    PrimitiveType::Vertex,
                    child_vertex);
                const simplex::Simplex parent_vertex =
                    input_mesh.map_to_parent(child_vertex_simplex);
                const double d = v_tet_offset_distance_acc.const_scalar_attribute(parent_vertex);
                min_dist = std::min(d, min_dist);
            }
            f_offset_distance_acc.scalar_attribute(face_tuple) = min_dist;
        }

        return true;
    };

    if (keep_topology) {
        // topological offset
        if (use_warm_start) {
            logger().info("Use warm start");
            if (adapt_offset_distance) {
                logger().info("Adapt offset distance");
                update_within_od_aggressive(within_od_handle, bbox);

                tag_topological_offset_tets(within_od_handle, tet_touches_boundary_handle, true);

                bool found_bad_tet = update_offset_distance_jacobian();

                auto v_dist_acc = m_mesh.create_accessor<double>(v_dist);
                auto v_inside_acc = m_mesh.create_accessor<int64_t>(v_inside);
                const auto p_acc = m_mesh.create_const_accessor<double>(m_pos_handle);
                for (const Tuple& t : m_mesh.get_all(PrimitiveType::Vertex)) {
                    const auto p = p_acc.const_vector_attribute(t);
                    const auto [sq_dist, od] = bvh.sq_dist_and_offset_distance(p);
                    const double d = std::sqrt(sq_dist);
                    v_dist_acc.scalar_attribute(t) = d;
                    if (d < od) {
                        v_inside_acc.scalar_attribute(t) = 1;
                    }
                }

                if (m_debug_print) {
                    logger().info("Performed aggressive tagging");
                    const std::string name = fmt::format("tog_debug_{}", m_debug_print_counter++);
                    utils::write_mesh(m_pos_handle, name, m_debug_print);
                    utils::write_mesh(m_pos_handle, name, PrimitiveType::Vertex, m_debug_print);
                }
                if (!found_bad_tet) {
                    logger().info("No colliding offsets detected.");
                } else {
                    logger().info("Colliding offsets detected. Adapting offset distance.");
                }

                f_to_v_offset_distance();
            }

            // update "within offset" attribute
            update_within_od_conservative(within_od_handle, bbox);

            reset_tet_tags();
            tag_topological_offset_tets(within_od_handle, tet_touches_boundary_handle, false);
        } else {
            logger().info("Use cold start");
        }
    } else {
        logger().info("Finite offset tagging");

        update_within_od_conservative(within_od_handle, bbox);

        auto tet_tag_acc = m_mesh.create_accessor<int64_t>(m_tag_handles[m_pt_top]);
        auto offset_tag_acc = m_mesh.create_accessor<int64_t>(m_offset_tag_handle);

        // finite offset
        for (const Tuple& t : m_mesh.get_all(m_pt_top)) {
            const simplex::Simplex tet(m_mesh, m_pt_top, t);

            if (tet_tag_acc.const_scalar_attribute(tet) == m_inside_tag) {
                continue;
            }

            if (tet_touches_boundary_acc.scalar_attribute(tet) == 1) {
                continue;
            }

            if (within_od_acc.const_scalar_attribute(tet) == 1) {
                tet_tag_acc.scalar_attribute(tet) = m_inside_tag;
                offset_tag_acc.scalar_attribute(tet) = m_offset_tag;
            }
        }
    }

    // transfer offset distance from faces to vertices
    f_to_v_offset_distance();

    // use regular space component to copy tags into the tets' faces
    regularize(false);

    if (m_debug_print) {
        logger().info("Finished tagging tets as inside (and adapting offset distance)");
        const std::string name = fmt::format("tog_debug_{}", m_debug_print_counter++);
        utils::write_mesh(m_pos_handle, name, m_debug_print);
    }

    // m_mesh.delete_attribute(v_inside);
    // m_mesh.delete_attribute(v_dist);
    // m_mesh.delete_attribute(laplacian_boundary_handle);
    // m_mesh.delete_attribute(v_tet_offset_distance_handle);
    // m_mesh.delete_attribute(visited_vertex_handle);
    // m_mesh.delete_attribute(tet_touches_boundary_handle);
    // m_mesh.delete_attribute(within_od_handle);
}

void TopologicalOffsetGenerator::adapt_offset_distance()
{
    const auto pos_acc = m_mesh.create_const_accessor<double>(m_pos_handle);

    const Eigen::MatrixXd bbox = wmtk::utils::bbox_from_mesh(m_pos_handle);

    auto within_od_handle =
        m_mesh.register_attribute<int64_t>("within_offset_distance", m_pt_top, 1);
    auto tet_touches_boundary_handle =
        m_mesh.register_attribute<int64_t>("tet_touches_boundary", m_pt_top, 1);
    auto visited_vertex_handle =
        m_mesh.register_attribute<int64_t>("tog_visited_vertex", PrimitiveType::Vertex, 1);
    auto v_tet_offset_distance_handle = m_mesh.register_attribute<double>(
        "tog_vertex_offset_distance",
        PrimitiveType::Vertex,
        1,
        false,
        m_offset_distance);
    auto laplacian_boundary_handle =
        m_mesh.register_attribute<int64_t>("laplacian_boundary", PrimitiveType::Vertex, 1);
    auto v_dist = m_mesh.register_attribute<double>("v_dist", PrimitiveType::Vertex, 1);
    auto v_inside = m_mesh.register_attribute<int64_t>("v_inside", PrimitiveType::Vertex, 1);

    std::map<PrimitiveType, attribute::Accessor<int64_t>> tag_accs;
    for (const auto& [pt, tag_handle] : m_tag_handles) {
        tag_accs.emplace(pt, m_mesh.create_accessor<int64_t>(tag_handle));
    }
    auto offset_tag_acc = m_mesh.create_accessor<int64_t>(m_offset_tag_handle);
    auto within_od_acc = m_mesh.create_accessor<int64_t>(within_od_handle);
    auto tet_touches_boundary_acc = m_mesh.create_accessor<int64_t>(tet_touches_boundary_handle);

    for (const Tuple& t : m_mesh.get_all(m_pt_top)) {
        const simplex::Simplex tet(m_mesh, m_pt_top, t);
        if (tet_touches_boundary(tet)) {
            tet_touches_boundary_acc.scalar_attribute(t) = 1;
        } else {
            tet_touches_boundary_acc.scalar_attribute(t) = 0;
        }
    }

    auto visited_vertex_acc = m_mesh.create_accessor<int64_t>(visited_vertex_handle);
    auto v_tet_offset_distance_acc = m_mesh.create_accessor<double>(v_tet_offset_distance_handle);
    auto laplacian_boundary_acc = m_mesh.create_accessor<int64_t>(laplacian_boundary_handle);
    auto v_dist_acc = m_mesh.create_accessor<double>(v_dist);
    auto v_inside_acc = m_mesh.create_accessor<int64_t>(v_inside);

    for (const Tuple& t : m_mesh.get_all(PrimitiveType::Vertex)) {
        const auto p = pos_acc.const_vector_attribute(t);
        auto [sq_dist, od] = m_bvh->sq_dist_and_offset_distance(p);
        v_tet_offset_distance_acc.scalar_attribute(t) = od;
    }

    auto update_offset_distance_jacobian = [&]() -> bool {
        const auto pos_acc = m_mesh.create_const_accessor<double>(m_pos_handle);
        const auto vertices = m_mesh.get_all(PrimitiveType::Vertex);

        // reset visisted tags
        for (const Tuple& t : vertices) {
            const auto p = pos_acc.const_vector_attribute(t);
            auto [sq_dist, od] = m_bvh->sq_dist_and_offset_distance(p);
            v_tet_offset_distance_acc.scalar_attribute(t) = od;
            visited_vertex_acc.scalar_attribute(t) = 0;
        }

        logger().info("Check for offset intersections with advancing front");
        bool found_bad_tet = false;

        // tag vertices that are outside
        for (const Tuple& tet_tuple : m_mesh.get_all(m_pt_top)) {
            const simplex::Simplex tet(m_mesh, m_pt_top, tet_tuple);

            if (tag_accs.at(m_pt_top).const_scalar_attribute(tet) != m_outside_tag) {
                continue;
            }

            double tet_dist = m_offset_distance;

            if (within_od_acc.const_scalar_attribute(tet) == 1) {
                // tet should be in offset but is not -> adapt offset distance

                //tag_accs.at(pt_top).scalar_attribute(tet) = 3; // DEBUG tag bad tets
                tet_dist = get_tet_small_dist(tet);
                found_bad_tet = true;
            }

            const auto tet_vertices = m_mesh.orient_vertices(tet_tuple);

            for (const Tuple& v_tuple : tet_vertices) {
                const simplex::IdSimplex v_ids =
                    m_mesh.get_id_simplex(v_tuple, PrimitiveType::Vertex);

                const double d = v_tet_offset_distance_acc.const_scalar_attribute(v_ids);
                v_tet_offset_distance_acc.scalar_attribute(v_ids) = std::min(d, tet_dist);
                visited_vertex_acc.scalar_attribute(v_ids) = 1;
            }
        }

        if (!found_bad_tet) {
            return false;
        }

        // bound difference on boundary vertices
        for (int64_t i = 0; i < 10; i++) {
            for (const Tuple& v_tuple : vertices) {
                const simplex::Simplex v(m_mesh, PrimitiveType::Vertex, v_tuple);
                if (visited_vertex_acc.const_scalar_attribute(v) != 1) {
                    continue;
                }

                auto neighs =
                    simplex::link_single_dimension_iterable(m_mesh, v, PrimitiveType::Vertex);
                double d_min = v_tet_offset_distance_acc.const_scalar_attribute(v);
                for (const Tuple& n : neighs) {
                    const double d = v_tet_offset_distance_acc.const_scalar_attribute(n);
                    d_min = std::min(d_min, 1.5 * d);
                }
                v_tet_offset_distance_acc.scalar_attribute(v) = d_min;
            }
        }

        for (int64_t i = 0; i < vertices.size(); ++i) {
            const simplex::IdSimplex v = m_mesh.get_id_simplex(vertices[i], PrimitiveType::Vertex);
            const double d = v_tet_offset_distance_acc.const_scalar_attribute(v);
            if (d < 0) {
                log_and_throw_error("Negative distance before igl::harmonic: {}", d);
            }
        }

        // use igl::harmonic
        std::map<simplex::IdSimplex, int64_t> vids;
        const auto tets = m_mesh.get_all(m_pt_top);

        for (int64_t i = 0; i < vertices.size(); ++i) {
            const simplex::IdSimplex v = m_mesh.get_id_simplex(vertices[i], PrimitiveType::Vertex);
            vids[v] = i;
        }

        if (m_debug_print) {
            logger().info("Laplacian init");
            const std::string name = fmt::format("tog_debug_{}", m_debug_print_counter++);
            utils::write_mesh(m_pos_handle, name, m_debug_print);
            utils::write_mesh(m_pos_handle, name, PrimitiveType::Vertex, m_debug_print);
        }

        Eigen::MatrixXd V; // #V by dim vertex positions
        Eigen::MatrixXi F; // #F by simplex-size list of element indices
        Eigen::VectorXi b; // #b boundary indices into V
        Eigen::MatrixXd bc; // #b by #W list of boundary values
        const int64_t k = 1; // power of harmonic operation (1: harmonic, 2: biharmonic, etc)
        Eigen::MatrixXd W; // #V by #W list of weights

        // V
        V.resize(vertices.size(), pos_acc.dimension());
        std::set<int> boundary_ids;
        for (int64_t i = 0; i < vertices.size(); ++i) {
            V.row(i) = pos_acc.const_vector_attribute(vertices[i]);
            if (visited_vertex_acc.const_scalar_attribute(vertices[i]) == 1) {
                laplacian_boundary_acc.scalar_attribute(vertices[i]) = 1;
                boundary_ids.insert(i);
                continue;
            }

            // check if vertex is inside
            bool is_inside = false;
            auto neighs = simplex::top_dimension_cofaces_iterable(
                m_mesh,
                simplex::Simplex(m_mesh, PrimitiveType::Vertex, vertices[i]));
            for (const Tuple& t : neighs) {
                if (offset_tag_acc.const_scalar_attribute(t) == m_offset_tag) {
                    is_inside = true;
                    break;
                }
            }
            if (!is_inside) {
                laplacian_boundary_acc.scalar_attribute(vertices[i]) = 2;
                boundary_ids.insert(i);
            }
        }
        // b & bc
        {
            std::vector<int> bvec(boundary_ids.begin(), boundary_ids.end());
            b = Eigen::Map<Eigen::VectorXi>(bvec.data(), bvec.size());
            bc.resize(bvec.size(), 1);
            for (int64_t i = 0; i < bvec.size(); ++i) {
                bc(i, 0) = v_tet_offset_distance_acc.const_scalar_attribute(vertices[bvec[i]]);
            }
        }

        int64_t nF = 0;
        for (int64_t i = 0; i < tets.size(); ++i) {
            if (offset_tag_acc.const_scalar_attribute(tets[i]) == m_offset_tag) {
                ++nF;
            }
        }

        // F
        F.resize(nF, pos_acc.dimension() + 1);
        int64_t f_counter = 0;
        for (int64_t i = 0; i < tets.size(); ++i) {
            if (offset_tag_acc.const_scalar_attribute(tets[i]) != m_offset_tag) {
                continue;
            }
            const auto tet_vertices = m_mesh.orient_vertices(tets[i]);
            for (int64_t j = 0; j < tet_vertices.size(); ++j) {
                const simplex::IdSimplex v =
                    m_mesh.get_id_simplex(tet_vertices[j], PrimitiveType::Vertex);
                F(f_counter, j) = vids[v];
            }
            ++f_counter;
        }

        logger().info("Solve igl::harmonic...");
        if (!igl::harmonic(V, F, b, bc, k, W)) {
            log_and_throw_error("Could not solve igl::harmonic!");
        }
        logger().info("... done");
        assert(W.rows() == V.rows());
        assert(W.cols() == bc.cols());

        for (int64_t i = 0; i < vertices.size(); ++i) {
            v_tet_offset_distance_acc.scalar_attribute(vertices[i]) = W(i, 0);
        }

        for (int64_t i = 0; i < vertices.size(); ++i) {
            const simplex::IdSimplex v = m_mesh.get_id_simplex(vertices[i], PrimitiveType::Vertex);
            const double d = v_tet_offset_distance_acc.const_scalar_attribute(v);
            if (d < 0) {
                double clamp_val = 1e-10;
                logger().warn(
                    "Negative distance after igl::harmonic: {}. Setting it to {}",
                    d,
                    clamp_val);
                v_tet_offset_distance_acc.scalar_attribute(v) = clamp_val;
            }
        }

        Mesh& input_mesh = m_bvh->mesh();
        auto f_offset_distance_acc =
            input_mesh.create_accessor<double>(m_bvh->f_offset_distance_handle());
        auto v_offset_distance_acc =
            input_mesh.create_accessor<double>(m_bvh->v_offset_distance_handle());

        // transfer values from tet vertices to tri faces and vertices
        for (const Tuple& face_tuple : input_mesh.get_all(m_pt_face)) {
            const auto child_vertices = input_mesh.orient_vertices(face_tuple);
            double min_dist = f_offset_distance_acc.const_scalar_attribute(face_tuple);
            for (const Tuple& child_vertex : child_vertices) {
                const simplex::Simplex child_vertex_simplex(
                    input_mesh,
                    PrimitiveType::Vertex,
                    child_vertex);
                const simplex::Simplex parent_vertex =
                    input_mesh.map_to_parent(child_vertex_simplex);
                const double d = v_tet_offset_distance_acc.const_scalar_attribute(parent_vertex);
                min_dist = std::min(d, min_dist);
            }
            f_offset_distance_acc.scalar_attribute(face_tuple) = min_dist;
        }

        return true;
    };

    logger().info("Adapt offset distance");
    update_within_od_aggressive(within_od_handle, bbox);

    tag_topological_offset_tets(within_od_handle, tet_touches_boundary_handle, true);

    bool found_bad_tet = update_offset_distance_jacobian();

    for (const Tuple& t : m_mesh.get_all(PrimitiveType::Vertex)) {
        const auto p = pos_acc.const_vector_attribute(t);
        const auto [sq_dist, od] = m_bvh->sq_dist_and_offset_distance(p);
        const double d = std::sqrt(sq_dist);
        v_dist_acc.scalar_attribute(t) = d;
        if (d < od) {
            v_inside_acc.scalar_attribute(t) = 1;
        }
    }

    if (m_debug_print) {
        logger().info("Performed aggressive tagging");
        const std::string name = fmt::format("tog_debug_{}", m_debug_print_counter++);
        utils::write_mesh(m_pos_handle, name, m_debug_print);
        utils::write_mesh(m_pos_handle, name, PrimitiveType::Vertex, m_debug_print);
    }
    if (!found_bad_tet) {
        logger().info("No colliding offsets detected.");
    } else {
        logger().info("Colliding offsets detected. Adapting offset distance.");
    }

    f_to_v_offset_distance();

    reset_tet_tags();

    if (m_debug_print) {
        logger().info("Finished adapting offset distance");
        const std::string name = fmt::format("tog_debug_{}", m_debug_print_counter++);
        utils::write_mesh(m_pos_handle, name, m_debug_print);
    }

    // m_mesh.delete_attribute(v_inside);
    // m_mesh.delete_attribute(v_dist);
    // m_mesh.delete_attribute(laplacian_boundary_handle);
    // m_mesh.delete_attribute(v_tet_offset_distance_handle);
    // m_mesh.delete_attribute(visited_vertex_handle);
    // m_mesh.delete_attribute(tet_touches_boundary_handle);
    // m_mesh.delete_attribute(within_od_handle);
}

void TopologicalOffsetGenerator::expand_topological_offset_input()
{
    const Eigen::MatrixXd bbox = wmtk::utils::bbox_from_mesh(m_pos_handle);

    auto within_od_handle =
        m_mesh.register_attribute<int64_t>("within_offset_distance", m_pt_top, 1, true);
    auto tet_touches_boundary_handle =
        m_mesh.register_attribute<int64_t>("tet_touches_boundary", m_pt_top, 1, true);

    std::map<PrimitiveType, attribute::Accessor<int64_t>> tag_accs;
    for (const auto& [pt, tag_handle] : m_tag_handles) {
        tag_accs.emplace(pt, m_mesh.create_accessor<int64_t>(tag_handle));
    }

    auto within_od_acc = m_mesh.create_accessor<int64_t>(within_od_handle);
    auto tet_touches_boundary_acc = m_mesh.create_accessor<int64_t>(tet_touches_boundary_handle);

    for (const Tuple& t : m_mesh.get_all(m_pt_top)) {
        const simplex::Simplex tet(m_mesh, m_pt_top, t);
        if (tet_touches_boundary(tet)) {
            tet_touches_boundary_acc.scalar_attribute(t) = 1;
        } else {
            tet_touches_boundary_acc.scalar_attribute(t) = 0;
        }
    }

    logger().info("Expand topological offset");

    update_within_od_conservative(within_od_handle, bbox);

    reset_tet_tags();
    tag_topological_offset_tets(within_od_handle, tet_touches_boundary_handle, false);

    // transfer offset distance from faces to vertices
    f_to_v_offset_distance();

    if (m_debug_print) {
        logger().info("Finished expanding the topological offset");
        const std::string name = fmt::format("tog_debug_{}", m_debug_print_counter++);
        utils::write_mesh(m_pos_handle, name, m_debug_print);
    }

    // m_mesh.delete_attribute(tet_touches_boundary_handle);
    // m_mesh.delete_attribute(within_od_handle);
}

void TopologicalOffsetGenerator::expand_finite_offset()
{
    const Eigen::MatrixXd bbox = wmtk::utils::bbox_from_mesh(m_pos_handle);

    auto within_od_handle =
        m_mesh.register_attribute<int64_t>("within_offset_distance", m_pt_top, 1, true);
    auto tet_touches_boundary_handle =
        m_mesh.register_attribute<int64_t>("tet_touches_boundary", m_pt_top, 1, true);

    std::map<PrimitiveType, attribute::Accessor<int64_t>> tag_accs;
    for (const auto& [pt, tag_handle] : m_tag_handles) {
        tag_accs.emplace(pt, m_mesh.create_accessor<int64_t>(tag_handle));
    }

    auto within_od_acc = m_mesh.create_accessor<int64_t>(within_od_handle);
    auto tet_touches_boundary_acc = m_mesh.create_accessor<int64_t>(tet_touches_boundary_handle);

    for (const Tuple& t : m_mesh.get_all(m_pt_top)) {
        const simplex::Simplex tet(m_mesh, m_pt_top, t);
        if (tet_touches_boundary(tet)) {
            tet_touches_boundary_acc.scalar_attribute(t) = 1;
        } else {
            tet_touches_boundary_acc.scalar_attribute(t) = 0;
        }
    }

    logger().info("Expand finite offset");

    update_within_od_conservative(within_od_handle, bbox);

    auto tet_tag_acc = m_mesh.create_accessor<int64_t>(m_tag_handles[m_pt_top]);
    auto offset_tag_acc = m_mesh.create_accessor<int64_t>(m_offset_tag_handle);

    // finite offset
    for (const Tuple& t : m_mesh.get_all(m_pt_top)) {
        const simplex::Simplex tet(m_mesh, m_pt_top, t);

        if (tet_tag_acc.const_scalar_attribute(tet) == m_inside_tag) {
            continue;
        }

        if (tet_touches_boundary_acc.scalar_attribute(tet) == 1) {
            continue;
        }

        if (within_od_acc.const_scalar_attribute(tet) == 1) {
            tet_tag_acc.scalar_attribute(tet) = m_inside_tag;
            offset_tag_acc.scalar_attribute(tet) = m_offset_tag;
        }
    }

    // transfer offset distance from faces to vertices
    f_to_v_offset_distance();

    // use regular space component to copy tags into the tets' faces
    regularize(false);

    if (m_debug_print) {
        logger().info("Finished expanding the finite offset");
        const std::string name = fmt::format("tog_debug_{}", m_debug_print_counter++);
        utils::write_mesh(m_pos_handle, name, m_debug_print);
    }

    // m_mesh.delete_attribute(tet_touches_boundary_handle);
    // m_mesh.delete_attribute(within_od_handle);
}

void TopologicalOffsetGenerator::expand_topological_offset_post()
{
    const Eigen::MatrixXd bbox = wmtk::utils::bbox_from_mesh(m_pos_handle);

    auto within_od_handle =
        m_mesh.register_attribute<int64_t>("within_offset_distance", m_pt_top, 1, true);
    auto tet_touches_boundary_handle =
        m_mesh.register_attribute<int64_t>("tet_touches_boundary", m_pt_top, 1, true);

    auto tet_offset_tag_acc = m_mesh.create_accessor<int64_t>(m_offset_tag_handle);

    std::map<PrimitiveType, attribute::Accessor<int64_t>> tag_accs;
    for (const auto& [pt, tag_handle] : m_tag_handles) {
        tag_accs.emplace(pt, m_mesh.create_accessor<int64_t>(tag_handle));
    }

    auto within_od_acc = m_mesh.create_accessor<int64_t>(within_od_handle);
    auto tet_touches_boundary_acc = m_mesh.create_accessor<int64_t>(tet_touches_boundary_handle);

    for (const Tuple& t : m_mesh.get_all(m_pt_top)) {
        const simplex::Simplex tet(m_mesh, m_pt_top, t);
        if (tet_touches_boundary(tet)) {
            tet_touches_boundary_acc.scalar_attribute(t) = 1;
        } else {
            tet_touches_boundary_acc.scalar_attribute(t) = 0;
        }
    }

    logger().info("Expand topological offset post");

    update_within_od_conservative(within_od_handle, bbox);

    tag_topological_offset_tets_through_face(within_od_handle, tet_touches_boundary_handle);

    regularize(false);

    // tag boundary of "input" as offset triangles
    for (const Tuple& face_tuple : m_mesh.get_all(m_pt_face)) {
        if (m_mesh.is_boundary(m_pt_face, face_tuple)) {
            continue;
        }

        if (tag_accs.at(m_pt_face).const_scalar_attribute(face_tuple) != m_inside_tag) {
            continue;
        }

        const simplex::Simplex face_simplex(m_mesh, m_pt_face, face_tuple);
        const auto top_tuples = simplex::top_dimension_cofaces_tuples(m_mesh, face_simplex);
        assert(top_tuples.size() == 2);
        for (const Tuple& tt : top_tuples) {
            if (tag_accs.at(m_pt_top).const_scalar_attribute(tt) == m_outside_tag) {
                tag_accs.at(m_pt_face).scalar_attribute(face_tuple) = m_offset_tag;
                break;
            }
        }
    }

    if (m_debug_print) {
        logger().info("Finished expanding the topological offset");
        const std::string name = fmt::format("tog_debug_{}", m_debug_print_counter++);
        utils::write_mesh(m_pos_handle, name, m_debug_print);
    }

    // m_mesh.delete_attribute(tet_touches_boundary_handle);
    // m_mesh.delete_attribute(within_od_handle);
}

void TopologicalOffsetGenerator::marching(const SimpleBVH::BVH& bvh, double offset_distance)
{
    /*
     * Using a distance field turned out to be not beneficial for the process as it might cause
     * vertices to snap to the wrong side and that messes with finding the nearest neighbor
     * later on.
     */
    // compute distance field
    auto df_handle = m_mesh.register_attribute<double>(
        "tog_distance_field",
        PrimitiveType::Vertex,
        1,
        false,
        std::numeric_limits<double>::max());

    const auto p_acc = m_mesh.create_const_accessor(m_pos_handle.as<double>());
    auto df_acc = m_mesh.create_accessor(df_handle.as<double>());

    if (offset_distance == 0) {
        auto v_tag_acc =
            m_mesh.create_const_accessor<int64_t>(m_tag_handles[PrimitiveType::Vertex]);
        for (const Tuple& t : m_mesh.get_all(PrimitiveType::Vertex)) {
            if (v_tag_acc.const_scalar_attribute(t) == m_inside_tag) {
                df_acc.scalar_attribute(t) = 0;
            } else {
                df_acc.scalar_attribute(t) = 1;
            }
        }
        offset_distance = 0.25;
    } else {
        double sq_dist;
        SimpleBVH::VectorMax3d nearest_point;

        for (const Tuple& t : m_mesh.get_all(PrimitiveType::Vertex)) {
            const simplex::Simplex v(m_mesh, PrimitiveType::Vertex, t);
            const Eigen::Vector3d p = p_acc.const_vector_attribute(v);

            bvh.nearest_facet(p, nearest_point, sq_dist);
            df_acc.scalar_attribute(v) = std::sqrt(sq_dist);
        }
    }

    MarchingOptions options;
    options.position_handle = m_pos_handle;
    options.label_handles[PrimitiveType::Vertex] = m_tag_handles[PrimitiveType::Vertex];
    options.label_handles[PrimitiveType::Edge] = m_tag_handles[PrimitiveType::Edge];

    if (m_mesh.top_simplex_type() == PrimitiveType::Tetrahedron) {
        options.label_handles[PrimitiveType::Triangle] = m_tag_handles[PrimitiveType::Triangle];
        options.pass_through_attributes.emplace_back(m_tag_handles[PrimitiveType::Tetrahedron]);
    } else if (m_mesh.top_simplex_type() == PrimitiveType::Triangle) {
        options.pass_through_attributes.emplace_back(m_tag_handles[PrimitiveType::Triangle]);
    }

    options.pass_through_attributes.emplace_back(m_offset_tag_handle);
    assert(
        m_mesh.top_simplex_type() == PrimitiveType::Tetrahedron &&
        m_mesh.has_attribute<int64_t>("img_tag", PrimitiveType::Tetrahedron));
    options.pass_through_attributes.emplace_back(
        m_mesh.get_attribute_handle<int64_t>("img_tag", PrimitiveType::Tetrahedron));
    options.input_values = {m_inside_tag};
    options.output_value = m_offset_tag;

    options.scalar_field = df_handle;
    options.isovalue = offset_distance;

    wmtk::components::marching(m_mesh, options);

    if (m_debug_print) {
        logger().info("Marching done");
        const std::string name = fmt::format("tog_debug_{}", m_debug_print_counter++);
        utils::write_mesh(m_pos_handle, name, m_debug_print);
    }
}

void TopologicalOffsetGenerator::marching_sampling(
    const SimpleBVH::BVH& bvh,
    double offset_distance)
{
    /*
     * Using a distance field turned out to be not beneficial for the process as it might cause
     * vertices to snap to the wrong side and that messes with finding the nearest neighbor
     * later on.
     */

    MarchingOptions options;
    options.position_handle = m_pos_handle;
    options.label_handles[PrimitiveType::Vertex] = m_tag_handles[PrimitiveType::Vertex];
    options.label_handles[PrimitiveType::Edge] = m_tag_handles[PrimitiveType::Edge];

    if (m_mesh.top_simplex_type() == PrimitiveType::Tetrahedron) {
        options.label_handles[PrimitiveType::Triangle] = m_tag_handles[PrimitiveType::Triangle];
        options.pass_through_attributes.emplace_back(m_tag_handles[PrimitiveType::Tetrahedron]);
    } else if (m_mesh.top_simplex_type() == PrimitiveType::Triangle) {
        options.pass_through_attributes.emplace_back(m_tag_handles[PrimitiveType::Triangle]);
    }

    options.pass_through_attributes.emplace_back(m_offset_tag_handle);
    assert(
        m_mesh.top_simplex_type() == PrimitiveType::Tetrahedron &&
        m_mesh.has_attribute<int64_t>("img_tag", PrimitiveType::Tetrahedron));
    options.pass_through_attributes.emplace_back(
        m_mesh.get_attribute_handle<int64_t>("img_tag", PrimitiveType::Tetrahedron));
    options.input_values = {m_inside_tag};
    options.output_value = m_offset_tag;

    auto oracle = [&bvh](const Eigen::VectorXd& p) -> double {
        double sq_dist;
        SimpleBVH::VectorMax3d nearest_point;
        bvh.nearest_facet(p, nearest_point, sq_dist);
        return std::sqrt(sq_dist);
    };
    options.oracle = oracle;
    options.isovalue = offset_distance;

    wmtk::components::marching(m_mesh, options);

    {
        auto df_handle = m_mesh.register_attribute<double>(
            "DEBUG_distance",
            PrimitiveType::Vertex,
            1,
            false,
            std::numeric_limits<double>::max());

        const auto p_acc = m_mesh.create_const_accessor(m_pos_handle.as<double>());
        auto df_acc = m_mesh.create_accessor(df_handle.as<double>());

        for (const Tuple& t : m_mesh.get_all(PrimitiveType::Vertex)) {
            const simplex::Simplex v(m_mesh, PrimitiveType::Vertex, t);
            const Eigen::Vector3d p = p_acc.const_vector_attribute(v);
            df_acc.scalar_attribute(v) = oracle(p);
        }
    }

    if (m_debug_print) {
        logger().info("Marching done");
        const std::string name = fmt::format("tog_debug_{}", m_debug_print_counter++);
        utils::write_mesh(m_pos_handle, name, m_debug_print);
    }
}

void TopologicalOffsetGenerator::marching()
{
    MarchingOptions options;
    options.position_handle = m_pos_handle;
    options.label_handles[PrimitiveType::Vertex] = m_tag_handles[PrimitiveType::Vertex];
    options.label_handles[PrimitiveType::Edge] = m_tag_handles[PrimitiveType::Edge];

    if (m_mesh.top_simplex_type() == PrimitiveType::Tetrahedron) {
        options.label_handles[PrimitiveType::Triangle] = m_tag_handles[PrimitiveType::Triangle];
        options.pass_through_attributes.emplace_back(m_tag_handles[PrimitiveType::Tetrahedron]);
    } else if (m_mesh.top_simplex_type() == PrimitiveType::Triangle) {
        options.pass_through_attributes.emplace_back(m_tag_handles[PrimitiveType::Triangle]);
    }

    options.pass_through_attributes.emplace_back(m_offset_tag_handle);
    assert(
        m_mesh.top_simplex_type() == PrimitiveType::Tetrahedron &&
        m_mesh.has_attribute<int64_t>("img_tag", PrimitiveType::Tetrahedron));
    options.pass_through_attributes.emplace_back(
        m_mesh.get_attribute_handle<int64_t>("img_tag", PrimitiveType::Tetrahedron));
    options.input_values = {m_inside_tag};
    options.output_value = m_offset_tag;

    wmtk::components::marching(m_mesh, options);

    if (m_debug_print) {
        logger().info("Marching done");
        const std::string name = fmt::format("tog_debug_{}", m_debug_print_counter++);
        utils::write_mesh(m_pos_handle, name, m_debug_print);
    }
}

attribute::MeshAttributeHandle TopologicalOffsetGenerator::tag_offset_tets()
{
    auto tog_v_acc = m_mesh.create_const_accessor<int64_t>(m_tag_handles.at(PrimitiveType::Vertex));
    auto tog_tet_acc = m_mesh.create_const_accessor<int64_t>(m_tag_handles.at(m_pt_top));
    auto tet_acc = m_mesh.create_accessor<int64_t>(m_offset_tag_handle);

    bool offset_tet_found = false;

    for (const Tuple& t : m_mesh.get_all(m_mesh.top_simplex_type())) {
        if (tet_acc.const_scalar_attribute(t) != 0) {
            continue;
        }

        if (tog_tet_acc.const_scalar_attribute(t) == m_inside_tag) {
            // tet_acc.scalar_attribute(t) = offset_tag_value;
            continue;
        }

        const auto vertices = simplex::faces_single_dimension_tuples(
            m_mesh,
            simplex::Simplex(m_mesh, m_mesh.top_simplex_type(), t),
            PrimitiveType::Vertex);

        int64_t offset_vertices = 0;
        int64_t input_vertices = 0;
        for (const Tuple& v : vertices) {
            const int64_t v_tag = tog_v_acc.const_scalar_attribute(v);
            if (v_tag == m_inside_tag) {
                input_vertices++;
            } else if (v_tag == m_offset_tag) {
                offset_vertices++;
            } else {
                break;
            }
        }

        if (offset_vertices > 0 && input_vertices > 0) {
            tet_acc.scalar_attribute(t) = m_offset_tag;
            tog_tet_acc.scalar_attribute(t) = m_inside_tag;
            offset_tet_found = true;
        }
    }

    regularize(false);

    if (m_debug_print) {
        logger().info("Tagging offset tets done");
        const std::string name = fmt::format("tog_debug_{}", m_debug_print_counter++);
        utils::write_mesh(m_pos_handle, name, m_debug_print);
    }

    return m_offset_tag_handle;
}

void TopologicalOffsetGenerator::tag_offset_triangles(
    const attribute::MeshAttributeHandle& tri_tag_handle,
    const int64_t offset_tag_value) const
{
    auto tog_acc = m_mesh.create_const_accessor<int64_t>(m_tag_handles.at(PrimitiveType::Triangle));
    auto tag_acc = m_mesh.create_accessor<int64_t>(tri_tag_handle);

    for (const Tuple& t : m_mesh.get_all(PrimitiveType::Triangle)) {
        if (tog_acc.const_scalar_attribute(t) == m_offset_tag) {
            tag_acc.scalar_attribute(t) = offset_tag_value;
        }
    }
}

std::shared_ptr<Mesh> TopologicalOffsetGenerator::generate_substructure_from_offset_tag()
{
    std::shared_ptr<Mesh> mesh_in = m_mesh.shared_from_this();
    auto sub = components::multimesh_from_tag(mesh_in, m_tag_handles.at(m_pt_face), m_offset_tag);

    Mesh& m = *sub;
    const PrimitiveType pt_face_face =
        get_primitive_type_from_id(get_primitive_type_id(m.top_simplex_type()) - 1);
    for (const Tuple& t : m.get_all(pt_face_face)) {
        if (m.is_boundary(pt_face_face, t)) {
            log_and_throw_error("Offset mesh is not closed.");
        }
    }

    return sub;
}

bool TopologicalOffsetGenerator::tet_touches_boundary(const simplex::Simplex& tet)
{
    Mesh& input_mesh = m_bvh->mesh();

    const auto children = m_mesh.get_child_meshes();

    for (const auto& child_mesh : children) {
        if (&*child_mesh == &input_mesh) {
            continue;
        }
        if (m_mesh.simplex_is_in_child(*child_mesh, tet)) {
            return true;
        }
    }

    bool tet_is_boundary = false;
    const auto vertices = simplex::faces_single_dimension(m_mesh, tet, PrimitiveType::Vertex);
    for (const simplex::Simplex& v : vertices) {
        if (m_mesh.is_boundary(v)) {
            return true;
        }

        // tet must not touch any child except for the input
        if (!m_mesh.simplex_is_in_child(input_mesh, v)) {
            for (const auto& child_mesh : children) {
                if (&*child_mesh == &input_mesh) {
                    continue;
                }
                if (m_mesh.simplex_is_in_child(*child_mesh, v)) {
                    return true;
                }
            }
        }
    }
    return false;
}

void TopologicalOffsetGenerator::update_within_od_conservative(
    const attribute::MeshAttributeHandle& within_od_handle,
    const Eigen::MatrixXd& bbox)
{
    using internal::utils::is_tet_in_offset_conservative_sampling;
    auto within_od_acc = m_mesh.create_accessor<int64_t>(within_od_handle);

    for (const Tuple& t : m_mesh.get_all(m_pt_top)) {
        const simplex::Simplex tet(m_mesh, m_pt_top, t);
        if (is_tet_in_offset_conservative_sampling(
                bbox,
                m_pos_handle,
                *m_bvh,
                tet,
                0.1 * m_offset_distance)) {
            // if (is_tet_in_offset(m_pos_handle, bvh, tet)) {
            within_od_acc.scalar_attribute(t) = 1;
        } else {
            within_od_acc.scalar_attribute(t) = 0;
        }
    }
}

void TopologicalOffsetGenerator::update_within_od_aggressive(
    const attribute::MeshAttributeHandle& within_od_handle,
    const Eigen::MatrixXd& bbox)
{
    using internal::utils::is_tet_in_offset_aggressive_sampling;
    auto within_od_acc = m_mesh.create_accessor<int64_t>(within_od_handle);

    for (const Tuple& t : m_mesh.get_all(m_pt_top)) {
        const simplex::Simplex tet(m_mesh, m_pt_top, t);
        if (is_tet_in_offset_aggressive_sampling(
                bbox,
                m_pos_handle,
                *m_bvh,
                tet,
                0.1 * m_offset_distance)) {
            // if (is_tet_in_offset(m_pos_handle, bvh, tet)) {
            within_od_acc.scalar_attribute(t) = 1;
        } else {
            within_od_acc.scalar_attribute(t) = 0;
        }
    }
}

void TopologicalOffsetGenerator::tag_topological_offset_tets(
    const attribute::MeshAttributeHandle& within_od_handle,
    const attribute::MeshAttributeHandle& tet_touches_boundary_handle,
    bool ignore_boundary)
{
    using internal::utils::tet_conains_connected_component;

    auto within_od_acc = m_mesh.create_accessor<int64_t>(within_od_handle);
    auto tet_touches_boundary_acc = m_mesh.create_accessor<int64_t>(tet_touches_boundary_handle);

    std::map<PrimitiveType, attribute::Accessor<int64_t>> tag_accs;
    for (const auto& [pt, tag_handle] : m_tag_handles) {
        tag_accs.emplace(pt, m_mesh.create_accessor<int64_t>(tag_handle));
    }
    auto offset_tag_acc = m_mesh.create_accessor<int64_t>(m_offset_tag_handle);

    auto is_tet_in_topological_offset = [&](const simplex::Simplex& tet,
                                            bool ignore_boundary = false) -> bool {
        if (tag_accs.at(m_pt_top).const_scalar_attribute(tet) != m_outside_tag) {
            return false;
        }

        // check if tet contains one connected component of inside tagged simplices
        if (!tet_conains_connected_component(m_tag_handles, m_inside_tag, tet)) {
            return false;
        }

        if (within_od_acc.const_scalar_attribute(tet) == 0) {
            return false;
        }

        if (!ignore_boundary && tet_touches_boundary_acc.scalar_attribute(tet) == 1) {
            return false;
        }

        return true;
    };

    logger().info("Tag tets that are inside offset distance.");
    std::vector<Tuple> q(200);
    size_t q_front = 0;
    size_t q_back = 0;

    for (const Tuple& tet_tuple : m_mesh.get_all(m_pt_top)) {
        const simplex::Simplex tet(m_mesh, m_pt_top, tet_tuple);

        if (!is_tet_in_topological_offset(tet, ignore_boundary)) {
            continue;
        }

        q[q_back++] = tet_tuple;
        if (q_back + 4 >= q.size()) {
            q.resize(q.size() * 1.5);
        }
    }

    // BFS on neighboring tets
    while (q_front != q_back) {
        const Tuple tet_tuple = q[q_front++];
        const simplex::Simplex tet(m_mesh, m_pt_top, tet_tuple);

        if (!is_tet_in_topological_offset(tet, ignore_boundary)) {
            continue;
        }

        tag_accs.at(m_pt_top).scalar_attribute(tet) = m_inside_tag;
        offset_tag_acc.scalar_attribute(tet) = m_offset_tag;
        for (const simplex::Simplex& f : simplex::faces(m_mesh, tet)) {
            tag_accs.at(f.primitive_type()).scalar_attribute(f) = m_inside_tag;
        }

        const auto vertices = m_mesh.orient_vertices(tet_tuple);
        simplex::IdSimplexCollection neighbors(m_mesh);
        for (const Tuple& v_tuple : vertices) {
            const simplex::Simplex v(m_mesh, PrimitiveType::Vertex, v_tuple);
            for (const Tuple& neighbor : simplex::top_dimension_cofaces_iterable(m_mesh, v)) {
                neighbors.add(m_pt_top, neighbor);
            }
        }
        neighbors.sort_and_clean();
        for (const simplex::IdSimplex neighbor : neighbors) {
            if (q_back + 4 >= q.size()) {
                q.resize(q.size() * 1.5);
            }
            q[q_back++] = m_mesh.get_tuple_from_id_simplex(neighbor);
        }
    }
}

void TopologicalOffsetGenerator::tag_topological_offset_tets_through_face(
    const attribute::MeshAttributeHandle& within_od_handle,
    const attribute::MeshAttributeHandle& tet_touches_boundary_handle)
{
    using internal::utils::tet_conains_face_connected_component;

    auto within_od_acc = m_mesh.create_accessor<int64_t>(within_od_handle);
    auto tet_touches_boundary_acc = m_mesh.create_accessor<int64_t>(tet_touches_boundary_handle);

    std::map<PrimitiveType, attribute::Accessor<int64_t>> tag_accs;
    for (const auto& [pt, tag_handle] : m_tag_handles) {
        tag_accs.emplace(pt, m_mesh.create_accessor<int64_t>(tag_handle));
    }
    auto offset_tag_acc = m_mesh.create_accessor<int64_t>(m_offset_tag_handle);

    auto is_tet_in_topological_offset = [&](const simplex::Simplex& tet,
                                            bool ignore_boundary = false) -> bool {
        if (tag_accs.at(m_pt_top).const_scalar_attribute(tet) != m_outside_tag) {
            return false;
        }

        // check if tet contains one connected component of inside tagged simplices
        if (!tet_conains_face_connected_component(m_tag_handles, m_inside_tag, tet)) {
            return false;
        }

        if (within_od_acc.const_scalar_attribute(tet) == 0) {
            return false;
        }

        if (!ignore_boundary && tet_touches_boundary_acc.scalar_attribute(tet) == 1) {
            return false;
        }

        return true;
    };

    logger().info("Tag tets that are inside offset distance.");
    std::vector<Tuple> q(200);
    size_t q_front = 0;
    size_t q_back = 0;

    for (const Tuple& tet_tuple : m_mesh.get_all(m_pt_top)) {
        const simplex::Simplex tet(m_mesh, m_pt_top, tet_tuple);

        if (tag_accs.at(m_pt_top).const_scalar_attribute(tet) == m_outside_tag) {
            continue;
        }

        if (q_back + 4 >= q.size()) {
            q.resize(q.size() * 1.5);
        }

        const auto faces = simplex::faces_single_dimension_tuples(m_mesh, tet, m_pt_face);
        for (Tuple t : faces) {
            if (m_mesh.is_boundary(m_pt_face, t)) {
                continue;
            }
            t = m_mesh.switch_tuple(t, m_pt_top);

            simplex::IdSimplex neighbor = m_mesh.get_id_simplex(t, m_pt_top);
            if (tag_accs.at(m_pt_top).const_scalar_attribute(neighbor) != m_outside_tag) {
                continue;
            }
            q[q_back++] = m_mesh.get_tuple_from_id_simplex(neighbor);
        }
    }

    // BFS on neighboring tets
    while (q_front != q_back) {
        const Tuple tet_tuple = q[q_front++];
        const simplex::Simplex tet(m_mesh, m_pt_top, tet_tuple);

        if (!is_tet_in_topological_offset(tet)) {
            continue;
        }

        tag_accs.at(m_pt_top).scalar_attribute(tet) = m_inside_tag;
        offset_tag_acc.scalar_attribute(tet) = m_offset_tag;
        for (const simplex::Simplex& f : simplex::faces(m_mesh, tet)) {
            tag_accs.at(f.primitive_type()).scalar_attribute(f) = m_inside_tag;
        }

        if (q_back + 4 >= q.size()) {
            q.resize(q.size() * 1.5);
        }

        const auto faces = simplex::faces_single_dimension_tuples(m_mesh, tet, m_pt_face);
        for (Tuple t : faces) {
            if (m_mesh.is_boundary(m_pt_face, t)) {
                continue;
            }
            t = m_mesh.switch_tuple(t, m_pt_top);

            simplex::IdSimplex neighbor = m_mesh.get_id_simplex(t, m_pt_top);
            q[q_back++] = m_mesh.get_tuple_from_id_simplex(neighbor);
        }
    }
}

Eigen::VectorXd TopologicalOffsetGenerator::get_tet_center(const simplex::Simplex& tet)
{
    auto pos_acc = m_mesh.create_const_accessor<double>(m_pos_handle);

    Eigen::VectorXd mid_point;
    mid_point.setZero(pos_acc.dimension());

    const auto vertices = m_mesh.orient_vertices(tet.tuple());
    for (const Tuple& v : vertices) {
        mid_point += pos_acc.const_vector_attribute(v);
    }
    mid_point /= vertices.size();

    return mid_point;
}

double TopologicalOffsetGenerator::get_tet_small_dist(const simplex::Simplex& tet)
{
    auto pos_acc = m_mesh.create_const_accessor<double>(m_pos_handle);

    Eigen::VectorXd mid_point;
    mid_point.setZero(pos_acc.dimension());

    double dv = std::numeric_limits<double>::max();

    const auto vertices = m_mesh.orient_vertices(tet.tuple());
    for (const Tuple& v : vertices) {
        auto p = pos_acc.const_vector_attribute(v);

        dv = std::min(dv, m_bvh->sq_dist(p));
        mid_point += p;
    }
    mid_point /= vertices.size();

    dv = std::sqrt(dv);
    const double d_mid = m_bvh->dist(mid_point);

    if (0.5 * (dv + d_mid) < 0) {
        logger()
            .info("Negative tet distance: {}, dv = {}, d_mid = {}", 0.5 * (dv + d_mid), dv, d_mid);
    }

    return 0.5 * (dv + d_mid);
}

void TopologicalOffsetGenerator::f_to_v_offset_distance()
{
    Mesh& input_mesh = m_bvh->mesh();
    auto f_offset_distance_acc =
        input_mesh.create_accessor<double>(m_bvh->f_offset_distance_handle());
    auto v_offset_distance_acc =
        input_mesh.create_accessor<double>(m_bvh->v_offset_distance_handle());

    for (const Tuple& f_tuple : input_mesh.get_all(input_mesh.top_simplex_type())) {
        const double f_dist = f_offset_distance_acc.const_scalar_attribute(f_tuple);
        const auto vertices = input_mesh.orient_vertices(f_tuple);
        for (const Tuple& v : vertices) {
            const double d = v_offset_distance_acc.const_scalar_attribute(v);
            v_offset_distance_acc.scalar_attribute(v) = std::min(d, f_dist);
        }
    }
}

void TopologicalOffsetGenerator::reset_tet_tags()
{
    Mesh& input_mesh = m_bvh->mesh();

    auto offset_tag_acc = m_mesh.create_accessor<int64_t>(m_offset_tag_handle);

    std::map<PrimitiveType, attribute::Accessor<int64_t>> tag_accs;
    for (const auto& [pt, tag_handle] : m_tag_handles) {
        tag_accs.emplace(pt, m_mesh.create_accessor<int64_t>(tag_handle));
    }

    for (const Tuple& t : m_mesh.get_all(m_pt_top)) {
        const simplex::Simplex tet(m_mesh, m_pt_top, t);
        if (offset_tag_acc.const_scalar_attribute(t) == m_offset_tag) {
            offset_tag_acc.scalar_attribute(t) = m_outside_tag;
            tag_accs.at(m_pt_top).scalar_attribute(t) = m_outside_tag;
            for (const simplex::Simplex& f : simplex::faces(m_mesh, tet)) {
                if (m_mesh.simplex_is_in_child(input_mesh, f)) {
                    continue;
                }
                tag_accs.at(f.primitive_type()).scalar_attribute(f) = m_outside_tag;
            }
        }
        //// DEBUG reset bad tet tag
        // if (tag_accs.at(pt_top).const_scalar_attribute(tet) == 3) {
        //     tag_accs.at(pt_top).scalar_attribute(tet) = m_outside_tag;
        // }
    }
}

} // namespace wmtk::components