#pragma once

#include <map>
#include <wmtk/Mesh.hpp>
#include "utils/OffsetBvh.hpp"

namespace wmtk::components {

/**
 * Generates the topological offset for a given surface.
 *
 * The surface and offset are described as tags on the faces.
 *
 * If the surface is a boundary representation, the inside can be tagged to avoid the generation of
 * the offset in there.
 *
 * Order of execution:
 * - add tag attributes for all primtive types
 * - set input tag
 * - set offset tag
 * - call regularize()
 * - call marching()
 */
class TopologicalOffsetGenerator
{
public:
    TopologicalOffsetGenerator(
        Mesh& mesh,
        attribute::MeshAttributeHandle& pos_handle,
        const double offset_distance,
        bool debug_print = false);

    /**
     * @brief Add everything in that attribute that is tagged with `m_outside_tag` to the input.
     */
    void add_tags(attribute::MeshAttributeHandle& tag_handle);

    /**
     * @brief Tag every simplex of a substructure as input.
     */
    void add_input_substructure(Mesh& child);

    void regularize(bool use_simplicial_embedding);

    /**
     * @brief Regularization can mess with other tags. Here, we transfer the outside tag from the
     * internal one to another.
     */
    void transfer_outside(
        const std::map<PrimitiveType, attribute::MeshAttributeHandle>& tag_handles);

    void register_bvh(const std::shared_ptr<internal::utils::OffsetBvh>& bvh);

    /**
     * @brief Tag tets as inside if they are within the offset distance of the input mesh.
     *
     * Tets on the boundary are never tagged as inside!
     *
     * This function basically changes the offset from an infinitessimal one to a finite offset.
     */
    void tag_tets_as_inside(
        bool keep_topology = true,
        bool use_warm_start = true,
        bool adapt_offset_distance = true);

    void adapt_offset_distance();

    void expand_topological_offset_input();

    void expand_topological_offset_input_for_marching();

    void expand_finite_offset();

    void expand_topological_offset_post();

    void marching(const SimpleBVH::BVH& bvh, double offset_distance);

    void marching_sampling(const SimpleBVH::BVH& bvh, double offset_distance);

    void marching();

    /**
     * @brief Tag all offset top dimension simplices that are in the offset region.
     *
     * The offset region is found by checking for simplices where every incident vertex has either
     * the input or the offset tag.
     *
     * @return the handle holding the offset tag
     */
    attribute::MeshAttributeHandle tag_offset_tets();

    void tag_offset_triangles(
        const attribute::MeshAttributeHandle& tri_tag_handle,
        const int64_t offset_tag_value) const;

    std::shared_ptr<Mesh> generate_substructure_from_offset_tag();

private:
    bool tet_touches_boundary(const simplex::Simplex& tet);

    /**
     * @brief Tag all tets that are (conservatively) fully within the offset distance
     */
    void update_within_od_conservative(
        const attribute::MeshAttributeHandle& within_od_handle,
        const Eigen::MatrixXd& bbox);

    /**
     * @brief Tag all tets that intersect the offset region at least partially.
     *
     * The tagged region is larger than the offset itself.
     */
    void update_within_od_aggressive(
        const attribute::MeshAttributeHandle& within_od_handle,
        const Eigen::MatrixXd& bbox);

    /**
     * @brief Tag all tets with all vertices within offset distance
     */
    void update_within_od_pointwise(const attribute::MeshAttributeHandle& within_od_handle);

    /**
     * @brief Grow the offset region without changing its topology.
     */
    void tag_topological_offset_tets(
        const attribute::MeshAttributeHandle& within_od_handle,
        const attribute::MeshAttributeHandle& tet_touches_boundary_handle,
        bool ignore_boundary = false);

    void tag_topological_offset_tets_through_face(
        const attribute::MeshAttributeHandle& within_od_handle,
        const attribute::MeshAttributeHandle& tet_touches_boundary_handle);

    Eigen::VectorXd get_tet_center(const simplex::Simplex& tet);

    double get_tet_small_dist(const simplex::Simplex& tet);

    void f_to_v_offset_distance();

    /**
     * @brief Reset offset region tag.
     */
    void reset_tet_tags();

private:
    Mesh& m_mesh;
    PrimitiveType m_pt_top;
    PrimitiveType m_pt_face;

    const double m_offset_distance = 0;

    bool m_debug_print = false;
    int64_t m_debug_print_counter = 0;

    attribute::MeshAttributeHandle m_pos_handle;

    std::shared_ptr<internal::utils::OffsetBvh> m_bvh;

    // Inside tags. The offset is around this tag.
    std::map<PrimitiveType, attribute::MeshAttributeHandle> m_tag_handles;
    attribute::MeshAttributeHandle m_offset_tag_handle;
    const int64_t m_inside_tag = 1;
    const int64_t m_outside_tag = 0;
    const int64_t m_offset_tag = 2;
};


} // namespace wmtk::components