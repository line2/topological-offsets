#pragma once

#include <wmtk/Mesh.hpp>

#include "TopologicalOffsetsOptions.hpp"
#include "TopologicalOffsetsOutput.hpp"

namespace wmtk::components {

TopologicalOffsetsOutput topological_offsets(TopologicalOffsetsOptions options);

} // namespace wmtk::components
