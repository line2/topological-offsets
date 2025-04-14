#include "apply_distance_function.hpp"

#include <tinyexpr.h>
#include <wmtk/utils/Logger.hpp>
#include <wmtk/utils/bbox_from_mesh.hpp>

namespace wmtk::components::internal::utils {

void apply_distance_function(
    attribute::MeshAttributeHandle& pos_handle,
    attribute::MeshAttributeHandle& distance_handle,
    double distance,
    const std::string& expression)
{
    Mesh& m = pos_handle.mesh();
    if (&m != &distance_handle.mesh()) {
        log_and_throw_error("pos_handle and distance_handle must be from the same mesh.");
    }
    if (pos_handle.primitive_type() != PrimitiveType::Vertex ||
        distance_handle.primitive_type() != PrimitiveType::Vertex) {
        log_and_throw_error("pos_handle and distance_handle must be vertex attributes.");
    }

    const auto pos_acc = m.create_const_accessor<double>(pos_handle);
    auto dist_acc = m.create_accessor<double>(distance_handle);

    const Eigen::MatrixXd bbox = wmtk::utils::bbox_from_mesh(pos_handle);

    double x, y;
    double ux, uy;

    int err;
    logger().info("Evaluating: {}", expression);
    te_expr* n;

    if (pos_handle.dimension() == 2) {
        te_variable vars[] = {{"x", &x}, {"y", &y}, {"ux", &ux}, {"uy", &uy}, {"d", &distance}};

        n = te_compile(expression.data(), vars, 5, &err);

        if (!n) {
            /* Show the user where the error is at. */
            std::string ind(err - 1, ' ');
            logger().info("Error near here\n{:s}\n{:s}{:s}", expression, ind, "^");
            log_and_throw_error("apply_distance_function failed to evaluate expression.");
        }

        for (const Tuple& t : m.get_all(PrimitiveType::Vertex)) {
            const auto p = pos_acc.const_vector_attribute(t);
            x = p(0);
            y = p(1);
            ux = (x - bbox(0, 0)) / (bbox(1, 0) - bbox(0, 0));
            uy = (y - bbox(0, 1)) / (bbox(1, 1) - bbox(0, 1));

            const double r = te_eval(n);
            dist_acc.scalar_attribute(t) = r;
        }
    } else if (pos_handle.dimension() == 3) {
        double z, uz;
        te_variable vars[] = {
            {"x", &x},
            {"y", &y},
            {"z", &z},
            {"ux", &ux},
            {"uy", &uy},
            {"uz", &uz},
            {"d", &distance}};

        n = te_compile(expression.data(), vars, 7, &err);

        if (!n) {
            /* Show the user where the error is at. */
            std::string ind(err - 1, ' ');
            logger().info("Error near here\n{:s}\n{:s}{:s}", expression, ind, "^");
            log_and_throw_error("apply_distance_function failed to evaluate expression.");
        }

        for (const Tuple& t : m.get_all(PrimitiveType::Vertex)) {
            const auto p = pos_acc.const_vector_attribute(t);
            x = p(0);
            y = p(1);
            z = p(2);
            ux = (x - bbox(0, 0)) / (bbox(1, 0) - bbox(0, 0));
            uy = (y - bbox(0, 1)) / (bbox(1, 1) - bbox(0, 1));
            uz = (z - bbox(0, 2)) / (bbox(1, 2) - bbox(0, 2));

            const double r = te_eval(n);
            dist_acc.scalar_attribute(t) = r;
        }
    }

    te_free(n);
}
} // namespace wmtk::components::internal::utils