#pragma once

#include <Eigen/Dense>
#include <wmtk/Mesh.hpp>

#include "OffsetBvh.hpp"

namespace wmtk::components::internal::utils {

/**
 * A wrapper for an Eigen::Matrix4d that represents quadrics on triangles.
 */
class Quadrics
{
public:
    /**
     * @brief Construct a quadric from the plane equation ax+by+cz+d=0.
     */
    Quadrics(const double a, const double b, const double c, const double d);

    /**
     * @brief Construct a quadric from a point on the plane and the normal.
     */
    Quadrics(const Eigen::Vector3d& p0, const Eigen::Vector3d& n);

    /**
     * @brief Load quadric from a face attribute.
     */
    Quadrics(const attribute::MeshAttributeHandle& quadric_handle, const simplex::Simplex& s);

    Quadrics(
        const OffsetBvh& bvh,
        const attribute::MeshAttributeHandle& pos_handle,
        const attribute::MeshAttributeHandle& tet_tag_handle,
        const attribute::MeshAttributeHandle& bvh_mesh_position_handle,
        const simplex::Simplex& s);

    Quadrics(
        const SimpleBVH::BVH& bvh,
        const attribute::MeshAttributeHandle& pos_handle,
        const attribute::MeshAttributeHandle& tet_tag_handle,
        const attribute::MeshAttributeHandle& bvh_mesh_position_handle,
        const simplex::Simplex& s,
        const double offset_distance);

    Quadrics(const Eigen::Ref<const Eigen::Matrix4d> matrix);

    Quadrics(const Quadrics& other);

    Quadrics operator+=(const Quadrics& other);

    Quadrics operator+(const Quadrics& other);

    Quadrics operator*(const double scalar);

    Quadrics operator*=(const double scalar);

    /**
     * @brief Solve quadrics with the SVD pseudo inverse.
     */
    Eigen::Vector3d solve(const Eigen::Ref<const Eigen::Vector3d> p);

    double squared_distance_to_optimum(const Eigen::Ref<const Eigen::Vector3d> p);

    /**
     * @brief Returns the quadrics in a vectorized form suitable for mesh attributes.
     */
    Eigen::VectorXd vectorized();

    /**
     * @brief Compute the distance from the center of gravity of the given simplex to the quadrics
     * optimum.
     */
    double squared_distance_to_optimum(
        const attribute::MeshAttributeHandle& position_handle,
        const simplex::Simplex& s);

    /**
     * @brief Retrieve the quadric of a vertex and convert it into a matrix.
     */
    static Eigen::Matrix4d attribute_to_matrix(
        const attribute::MeshAttributeHandle& quadric_handle,
        const Tuple& t);

    static void matrix_to_attribute(
        attribute::MeshAttributeHandle& quadric_handle,
        const Tuple& t,
        const Eigen::Ref<Eigen::Matrix4d> mat);

    /**
     * @brief Get the squared distance from the simplex assuming that its vertices were moved to
     * their quadrics position.
     *
     * If the simplex is a vertex, it just returns the squared distance to its quadrics position.
     */
    static double squared_distance_of_projected_simplex_to_optimum(
        const attribute::MeshAttributeHandle& position_handle,
        const attribute::MeshAttributeHandle& quadric_handle,
        const simplex::Simplex& s);

private:
    Eigen::Matrix4d m_matrix; // the actual quadrics
};

} // namespace wmtk::components::internal::utils
