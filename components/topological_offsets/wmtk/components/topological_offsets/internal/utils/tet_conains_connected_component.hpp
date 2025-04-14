#pragma once

#include <wmtk/Mesh.hpp>

namespace wmtk::components::internal::utils {

/**
 * @brief The tagged region represents a single connected component.
 *
 * This check can be performed by iterating through all tagged vertices. If any of the vertices
 * contains all tagged simplices in its open star and or in the faces of the open star, the tagged
 * region is a single connected component. The same is true for any lower-dimensional simplex.
 */
bool tet_conains_connected_component(
    std::map<PrimitiveType, attribute::MeshAttributeHandle>& tag_handles,
    const int64_t tag_value,
    const simplex::Simplex& s);

bool tet_conains_face_connected_component(
    std::map<PrimitiveType, attribute::MeshAttributeHandle>& tag_handles,
    const int64_t tag_value,
    const simplex::Simplex& s);

} // namespace wmtk::components::internal::utils