#pragma once

#include <filesystem>
#include <wmtk/Mesh.hpp>

namespace wmtk::utils {

std::shared_ptr<Mesh> igl_input(const std::filesystem::path& file_name);

} // namespace wmtk::utils