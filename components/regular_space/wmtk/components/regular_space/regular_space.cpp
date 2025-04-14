#include "regular_space.hpp"

#include <wmtk/Mesh.hpp>
#include <wmtk/components/utils/get_attributes.hpp>
#include <wmtk/io/Cache.hpp>
#include <wmtk/utils/Logger.hpp>
#include <wmtk/utils/primitive_range.hpp>

#include "internal/RegularSpace.hpp"

namespace wmtk::components {

// auto gather_attributes(io::Cache& cache, Mesh& mesh, const RegularSpaceOptions& options)
//{
//    // collect labels that were there already before this component
//    std::vector<attribute::MeshAttributeHandle> original_attributes;
//
//    std::vector<attribute::MeshAttributeHandle> label_attributes;
//    for (const PrimitiveType& ptype : wmtk::utils::primitive_below(mesh.top_simplex_type())) {
//        std::string attr_name;
//        switch (ptype) {
//        case PrimitiveType::Vertex: {
//            attr_name = "vertex_label";
//            break;
//        }
//        case PrimitiveType::Edge: {
//            attr_name = "edge_label";
//            break;
//        }
//        case PrimitiveType::Triangle: {
//            attr_name = "face_label";
//            break;
//        }
//        case PrimitiveType::Tetrahedron: {
//            attr_name = "tetrahedron_label";
//            break;
//        }
//        default: log_and_throw_error("Unknown primitive type: {}", ptype);
//        }
//
//        if (options.attributes.find(attr_name) == options.attributes.end()) {
//            // no attribute was given for that primitive type
//            auto h = mesh.register_attribute<int64_t>(attr_name + "_regular_space", ptype, 1);
//            label_attributes.emplace_back(h);
//        } else {
//            auto h = mesh.get_attribute_handle<int64_t>(options.attributes.at(attr_name), ptype);
//            label_attributes.emplace_back(h);
//            original_attributes.emplace_back(h);
//        }
//    }
//
//    auto pass_through_attributes = utils::get_attributes(cache, mesh, options.pass_through);
//
//    return std::make_tuple(original_attributes, label_attributes, pass_through_attributes);
//}

void regular_space(Mesh& mesh, const RegularSpaceOptions& options)
{
    using namespace internal;

    // auto [original_attributes, label_attributes, pass_through_attributes] =
    //     gather_attributes(cache, mesh, options);
    //
    //// clean up attributes
    //{
    //    std::vector<attribute::MeshAttributeHandle> keeps = pass_through_attributes;
    //    keeps.insert(keeps.end(), original_attributes.begin(), original_attributes.end());
    //    mesh.clear_attributes(keeps);
    //}
    //
    // std::tie(original_attributes, label_attributes, pass_through_attributes) =
    //    gather_attributes(cache, mesh, options);

    std::vector<attribute::MeshAttributeHandle> tag_attr_vec;
    for (const PrimitiveType& ptype : wmtk::utils::primitive_below(mesh.top_simplex_type())) {
        tag_attr_vec.emplace_back(options.tag_attributes.at(ptype));
    }

    RegularSpace rs(mesh, tag_attr_vec, options.value, options.pass_through_attributes);
    rs.regularize_tags(options.generate_simplicial_embedding);

    // clean up attributes
    {
        std::vector<attribute::MeshAttributeHandle> keeps = options.pass_through_attributes;
        keeps.insert(keeps.end(), tag_attr_vec.begin(), tag_attr_vec.end());
        mesh.clear_attributes(keeps);
    }

    mesh.consolidate();
}

} // namespace wmtk::components