#include "igl_input.hpp"

#include <igl/readSTL.h>
#include <fstream>
#include <iostream>
#include <wmtk/TriMesh.hpp>
#include <wmtk/utils/Logger.hpp>
#include <wmtk/utils/mesh_utils.hpp>

namespace wmtk::utils {

std::shared_ptr<Mesh> igl_input(const std::filesystem::path& file_name)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    if (file_name.extension() == ".stl") {
        Eigen::MatrixXd N;

        std::ifstream ifs(file_name);

        const bool success = igl::readSTL(ifs, V, F, N);
        if (!success) {
            log_and_throw_error("Cannot read file `{}`.", file_name.string());
        }
        ifs.close();

        auto m_ptr = std::make_shared<TriMesh>();
        TriMesh& m = *m_ptr;

        m.initialize(F.cast<int64_t>());
        mesh_utils::set_matrix_attribute(V, "vertices", PrimitiveType::Vertex, m);

        return m_ptr;
    }

    log_and_throw_error("Unknown file extension {} in igl_input.", file_name.extension().string());
}

} // namespace wmtk::utils