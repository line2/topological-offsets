#include "tet_conains_connected_component.hpp"

#include <wmtk/Mesh.hpp>
#include <wmtk/attribute/Accessor.hpp>
#include <wmtk/simplex/IdSimplexCollection.hpp>
#include <wmtk/simplex/faces.hpp>
#include <wmtk/simplex/faces_single_dimension.hpp>
#include <wmtk/simplex/open_star.hpp>
#include <wmtk/simplex/tuples_preserving_primitive_types.hpp>
#include <wmtk/utils/primitive_range.hpp>

namespace wmtk::components::internal::utils {

bool tet_conains_connected_component(
    std::map<PrimitiveType, attribute::MeshAttributeHandle>& tag_handles,
    const int64_t tag_value,
    const simplex::Simplex& s)
{
    Mesh& m = tag_handles.at(PrimitiveType::Vertex).mesh();

    std::map<PrimitiveType, const attribute::Accessor<int64_t>> tag_accs;
    for (const auto& [pt, tag_handle] : tag_handles) {
        tag_accs.emplace(pt, m.create_accessor<int64_t>(tag_handle));
    }

    if (s.primitive_type() == PrimitiveType::Vertex) {
        return tag_accs.at(PrimitiveType::Vertex).const_scalar_attribute(s) == tag_value;
    }

    simplex::SimplexCollection tagged_simplices(m);
    // get all tagged faces of s
    for (const simplex::Simplex& face_simplex : simplex::faces(m, s)) {
        if (tag_accs.at(face_simplex.primitive_type()).const_scalar_attribute(face_simplex) ==
            tag_value) {
            tagged_simplices.add(face_simplex);
        }
    }
    tagged_simplices.sort_and_clean();

    if (tagged_simplices.size() == 0) {
        return false;
    }

    const auto vertices = simplex::faces_single_dimension(m, s, PrimitiveType::Vertex);
    for (const simplex::Simplex& v : vertices) {
        if (tag_accs.at(PrimitiveType::Vertex).const_scalar_attribute(v) != tag_value) {
            continue;
        }

        // find tagged simplices in open star of a vertex
        const auto vertex_tuples_in_tet = simplex::tuples_preserving_primitive_types(
            m,
            v.tuple(),
            v.primitive_type(),
            m.top_simplex_type());

        // add tagged simplices and their faces to a SimplexCollection
        simplex::SimplexCollection sc(m, {v});
        for (const PrimitiveType pt :
             wmtk::utils::primitive_range(PrimitiveType::Edge, m.top_simplex_type())) {
            if (pt == m.top_simplex_type()) {
                continue;
            }
            for (const Tuple& t_coface : vertex_tuples_in_tet) {
                const simplex::Simplex s_cof(m, pt, t_coface);
                if (tag_accs.at(s_cof.primitive_type()).const_scalar_attribute(s_cof) ==
                    tag_value) {
                    sc.add(s_cof);
                    simplex::faces(sc, s_cof, false);
                }
            }
        }
        sc.sort_and_clean();

        // the tet can be tagged if all tagged simplices are included in that SimplexCollection
        if (simplex::SimplexCollection::are_simplex_collections_equal(tagged_simplices, sc)) {
            return true;
        }
    }

    return false;
}

bool tet_conains_face_connected_component(
    std::map<PrimitiveType, attribute::MeshAttributeHandle>& tag_handles,
    const int64_t tag_value,
    const simplex::Simplex& s)
{
    Mesh& m = tag_handles.at(PrimitiveType::Vertex).mesh();

    std::map<PrimitiveType, const attribute::Accessor<int64_t>> tag_accs;
    for (const auto& [pt, tag_handle] : tag_handles) {
        tag_accs.emplace(pt, m.create_accessor<int64_t>(tag_handle));
    }

    if (s.primitive_type() == PrimitiveType::Vertex) {
        return tag_accs.at(PrimitiveType::Vertex).const_scalar_attribute(s) == tag_value;
    }

    simplex::SimplexCollection tagged_simplices(m);
    // get all tagged faces of s
    for (const simplex::Simplex& face_simplex : simplex::faces(m, s)) {
        if (tag_accs.at(face_simplex.primitive_type()).const_scalar_attribute(face_simplex) ==
            tag_value) {
            tagged_simplices.add(face_simplex);
        }
    }
    tagged_simplices.sort_and_clean();

    if (tagged_simplices.size() == 0) {
        return false;
    }

    const auto vertices = simplex::faces_single_dimension(m, s, PrimitiveType::Vertex);
    for (const simplex::Simplex& v : vertices) {
        if (tag_accs.at(PrimitiveType::Vertex).const_scalar_attribute(v) != tag_value) {
            continue;
        }

        // find tagged simplices in open star of a vertex
        const auto vertex_tuples_in_tet = simplex::tuples_preserving_primitive_types(
            m,
            v.tuple(),
            v.primitive_type(),
            m.top_simplex_type());

        // add tagged simplices and their faces to a SimplexCollection
        simplex::SimplexCollection sc(m);

        PrimitiveType pt_face =
            get_primitive_type_from_id(get_primitive_type_id(s.primitive_type()) - 1);

        for (const Tuple& t_coface : vertex_tuples_in_tet) {
            const simplex::Simplex s_cof(m, pt_face, t_coface);
            if (tag_accs.at(s_cof.primitive_type()).const_scalar_attribute(s_cof) == tag_value) {
                sc.add(s_cof);
                simplex::faces(sc, s_cof, false);
            }
        }

        sc.sort_and_clean();

        // the tet can be tagged if all tagged simplices are included in that SimplexCollection
        if (simplex::SimplexCollection::are_simplex_collections_equal(tagged_simplices, sc)) {
            return true;
        }
    }

    return false;
}

} // namespace wmtk::components::internal::utils