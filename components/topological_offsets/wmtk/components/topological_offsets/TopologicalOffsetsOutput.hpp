#pragma once

#include <nlohmann/json.hpp>
#include <optional>
#include <wmtk/Mesh.hpp>
#include <wmtk/attribute/MeshAttributeHandle.hpp>

namespace wmtk::components {

struct TopologicalOffsetsOutput
{
public:
    std::shared_ptr<Mesh> offset_mesh;
    attribute::MeshAttributeHandle offset_position_handle;

    /**
     * This mesh represents the volume enclosed by the volume and the input.
     */
    std::shared_ptr<Mesh> offset_volume_mesh;

    nlohmann::json report;
};

} // namespace wmtk::components
