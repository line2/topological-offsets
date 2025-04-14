#pragma once

#include <wmtk/Mesh.hpp>

namespace wmtk::components {

struct RegularSpaceOptions
{
    std::map<PrimitiveType, attribute::MeshAttributeHandle> tag_attributes;
    int64_t value;
    std::vector<attribute::MeshAttributeHandle> pass_through_attributes;
    bool generate_simplicial_embedding = true;
};

} // namespace wmtk::components