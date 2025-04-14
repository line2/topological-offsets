#include "tetmesh_inversion_check.hpp"

#include <wmtk/TetMesh.hpp>
#include <wmtk/simplex/faces_single_dimension.hpp>
#include <wmtk/utils/Logger.hpp>

#include "predicates.h"

int64_t wmtk::utils::tetmesh_inversion_check(const Mesh& m, const std::string& positions)
{
    if (m.top_simplex_type() != PrimitiveType::Tetrahedron) {
        logger().warn("This is not a TetMesh, cannot check for inversion.");
        return -1;
    }

    int64_t inversion_counter = 0;

    auto h = m.get_attribute_handle<double>(positions, PrimitiveType::Vertex);
    auto acc = m.create_const_accessor<double>(h);
    const auto tets = m.get_all(PrimitiveType::Tetrahedron);
    for (Tuple t : tets) {
        const auto vs = simplex::faces_single_dimension_tuples(
            m,
            simplex::Simplex::tetrahedron(m, t),
            PrimitiveType::Vertex);
        assert(vs.size() == 4);
        Eigen::Vector3d p0 = acc.const_vector_attribute(vs[0]);
        Eigen::Vector3d p1 = acc.const_vector_attribute(vs[1]);
        Eigen::Vector3d p2 = acc.const_vector_attribute(vs[2]);
        Eigen::Vector3d p3 = acc.const_vector_attribute(vs[3]);
        if (m.is_ccw(t)) {
            if (orient3d(p3.data(), p0.data(), p1.data(), p2.data()) < 0) {
                inversion_counter++;
            }
        } else {
            if (orient3d(p3.data(), p0.data(), p2.data(), p1.data()) < 0) {
                inversion_counter++;
            }
        }
    }

    if (inversion_counter > 0) {
        logger().warn("{} of {} tets are inverted.", inversion_counter, tets.size());
    }

    return inversion_counter;
}
