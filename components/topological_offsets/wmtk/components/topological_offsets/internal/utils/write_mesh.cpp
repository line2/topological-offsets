#include "write_mesh.hpp"

#include <wmtk/Mesh.hpp>
#include <wmtk/io/ParaviewWriter.hpp>
#include <wmtk/utils/Logger.hpp>

namespace wmtk::components::utils {
void write_mesh(
    const attribute::MeshAttributeHandle& pos_handle,
    const std::string& name,
    const bool intermediate_output)
{
    if (!intermediate_output) {
        return;
    }

    const Mesh& mesh = pos_handle.mesh();


    logger().info("Write mesh '{}'", name);

    const bool edge = mesh.top_simplex_type() == PrimitiveType::Edge;
    const bool tri = mesh.top_simplex_type() == PrimitiveType::Triangle;
    const bool tet = mesh.top_simplex_type() == PrimitiveType::Tetrahedron;

    wmtk::io::ParaviewWriter
        writer(name, mesh.get_attribute_name(pos_handle.as<double>()), mesh, false, edge, tri, tet);
    mesh.serialize(writer);
}

void write_mesh(
    const attribute::MeshAttributeHandle& pos_handle,
    const std::string& name,
    const PrimitiveType& pt,
    const bool intermediate_output)
{
    if (!intermediate_output) {
        return;
    }

    const Mesh& mesh = pos_handle.mesh();

    logger().info("Write mesh '{}'", name);

    const bool vert = pt == PrimitiveType::Vertex;
    const bool edge = pt == PrimitiveType::Edge;
    const bool tri = pt == PrimitiveType::Triangle;
    const bool tet = pt == PrimitiveType::Tetrahedron;

    wmtk::io::ParaviewWriter
        writer(name, mesh.get_attribute_name(pos_handle.as<double>()), mesh, vert, edge, tri, tet);
    mesh.serialize(writer);
}
} // namespace wmtk::components::utils