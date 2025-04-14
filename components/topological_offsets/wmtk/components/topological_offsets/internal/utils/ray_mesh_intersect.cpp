#include "ray_mesh_intersect.hpp"

#include <igl/ray_mesh_intersect.h>
#include <wmtk/simplex/faces_single_dimension.hpp>
#include <wmtk/simplex/top_dimension_cofaces.hpp>
#include <wmtk/utils/Logger.hpp>

bool wmtk::components::internal::utils::ray_mesh_intersect(
    const attribute::MeshAttributeHandle& position_handle,
    const Eigen::Vector3d& source,
    const Eigen::Vector3d& dir,
    igl::Hit& hit)
{
    const Mesh& mesh = position_handle.mesh();

    constexpr PrimitiveType PV = PrimitiveType::Vertex;
    constexpr PrimitiveType PE = PrimitiveType::Edge;

    const attribute::Accessor<double> accessor =
        mesh.create_const_accessor(position_handle.as<double>());

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    if (mesh.top_simplex_type() == PrimitiveType::Triangle) {
        int64_t count = 0;
        int64_t index = 0;
        assert(accessor.dimension() == 3);

        const std::vector<Tuple>& face_tuples = mesh.get_all(PrimitiveType::Triangle);

        V.resize(3 * face_tuples.size(), accessor.dimension());
        F.resize(face_tuples.size(), 3);

        for (const Tuple& f : face_tuples) {
            auto p0 = accessor.const_vector_attribute(f);
            auto p1 = accessor.const_vector_attribute(mesh.switch_tuple(f, PV));
            auto p2 = accessor.const_vector_attribute(mesh.switch_tuples(f, {PE, PV}));

            F.row(index) = Eigen::Vector3i(count, count + 1, count + 2);
            V.row(3 * index) = p0;
            V.row(3 * index + 1) = p1;
            V.row(3 * index + 2) = p2;

            count += 3;
            ++index;
        }
    } else if (mesh.top_simplex_type() == PrimitiveType::Edge) {
        int64_t count = 0;
        int64_t index = 0;

        const std::vector<Tuple>& edge_tuples = mesh.get_all(PrimitiveType::Edge);

        V.resize(2 * edge_tuples.size(), accessor.dimension());
        F.resize(edge_tuples.size(), 2);

        for (const Tuple& e : edge_tuples) {
            auto p0 = accessor.const_vector_attribute(e);
            auto p1 = accessor.const_vector_attribute(mesh.switch_tuple(e, PV));

            F.row(index) = Eigen::Vector2i(count, count + 1);
            V.row(2 * index) = p0;
            V.row(2 * index + 1) = p1;

            count += 2;
            ++index;
        }
    } else {
        log_and_throw_error("bvh_from_mesh works only for tri/edges meshes");
    }

    return igl::ray_mesh_intersect(source, dir, V, F, hit);
}
