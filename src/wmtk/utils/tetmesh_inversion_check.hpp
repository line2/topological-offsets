#pragma once

#include <wmtk/Mesh.hpp>

namespace wmtk::utils {

int64_t tetmesh_inversion_check(const Mesh& m, const std::string& positions);
}