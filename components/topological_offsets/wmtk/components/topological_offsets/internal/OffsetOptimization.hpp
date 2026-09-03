#pragma once

#include <SimpleBVH/BVH.hpp>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <wmtk/Mesh.hpp>
#include <wmtk/function/LocalNeighborsSumFunction.hpp>
#include <wmtk/function/PerSimplexFunction.hpp>
#include <wmtk/invariants/Invariant.hpp>
#include <wmtk/invariants/InvariantCollection.hpp>
#include <wmtk/operations/attribute_new/CollapseNewAttributeStrategy.hpp>
#include <wmtk/operations/attribute_update/AttributeTransferStrategy.hpp>

#include "utils/OffsetBvh.hpp"

namespace wmtk::components::internal {

/**
 * Optimize the offset surface.
 *
 * The input and offset mesh and all substructures that should be preserved musst be registered as
 * child meshes.
 */
class OffsetOptimization
{
public:
    OffsetOptimization(Mesh& mesh, attribute::MeshAttributeHandle& pos_handle);

    void set_offset_mesh(
        const std::shared_ptr<Mesh> offset_mesh,
        const std::shared_ptr<Mesh> input_mesh,
        attribute::MeshAttributeHandle offset_pos_handle,
        attribute::MeshAttributeHandle input_pos_handle);

    void set_input_bvh(const std::shared_ptr<utils::OffsetBvh>& input_bvh);

    void set_offset_distance(const double offset_distance);

    void set_normal_deviation(const double min_normal_deviation, const double max_normal_deviation);

    void set_edge_length_constraints(const double min_edge_length, const double max_edge_length);

    void set_tetrahedron_tag_handle(const attribute::MeshAttributeHandle& tet_tag_handle);

    void init_embedding_optimization();
    void init_offset_optimization();

    void optimize_embedding(const int64_t n_iterations);
    void optimize_offset(const int64_t n_iterations);

    /**
     * Smooth offset and embedding.
     */
    void smooth_all(const int64_t n_iterations);

    std::vector<double> get_error_metrics();

    nlohmann::json& metrics_json();

    void set_debug_prints(bool embedding, bool offset, bool smooth);

private:
    /**
     * @brief Update region of interest for tet mesh.
     */
    void roi_update();

    /**
     * @brief Greedy-color the offset vertices so that any two vertices whose
     * parent embedding vertices are equal or share an embedding tet get
     * different colors. Writes the colors into m_offset_color_handle; this
     * makes parallel offset smoothing safe (the smoothing lambda temporarily
     * moves the parent embedding vertex and checks inversions in its
     * tet-neighborhood).
     */
    void color_offset_vertices_by_embedding_adjacency();

    std::tuple<Eigen::Vector3d, double, double> compute_offset_normal(const Eigen::Vector3d& p);
    std::tuple<Eigen::Vector2d, double, double> compute_offset_normal(const Eigen::Vector2d& p);

    /**
     * @brief Check that each offset triangle has exactly one tagged tet as coface.
     */
    void check_tet_tag_validity();

private:
    bool m_debug_print_embedding = false;
    bool m_debug_print_offset = false;
    bool m_debug_print_smooth_all = false;

    nlohmann::json m_metrics_report;

    Mesh& m_mesh;
    std::shared_ptr<Mesh> m_input_mesh;
    std::shared_ptr<Mesh> m_offset_mesh;

    attribute::MeshAttributeHandle m_pos_handle;
    attribute::MeshAttributeHandle m_input_pos_handle;
    attribute::MeshAttributeHandle m_offset_pos_handle;

    attribute::MeshAttributeHandle m_tet_tag_handle;

    std::shared_ptr<utils::OffsetBvh> m_input_bvh;

    double m_offset_distance = 0;
    double m_min_normal_deviation = 0;
    double m_max_normal_deviation = std::numeric_limits<double>::max();
    double m_min_edge_length = 0;
    double m_max_edge_length = std::numeric_limits<double>::max();

    ////////////////////////////////
    // attributes for optimization
    attribute::MeshAttributeHandle m_embedding_edge_length_attribute;
    attribute::MeshAttributeHandle m_embedding_target_edge_length_attribute;
    attribute::MeshAttributeHandle m_embedding_amips_attribute;
    attribute::MeshAttributeHandle m_embedding_roi_attribute;
    attribute::MeshAttributeHandle m_offset_normal_samples_attribute;
    attribute::MeshAttributeHandle m_offset_edge_length_attribute;
    attribute::MeshAttributeHandle m_offset_target_edge_length_attribute;
    attribute::MeshAttributeHandle m_offset_face_normal_deviation_attribute;
    attribute::MeshAttributeHandle m_offset_mean_ratio_metric_attribute;
    attribute::MeshAttributeHandle m_offset_face_distance_attribute;
    attribute::MeshAttributeHandle m_offset_face_target_distance_attribute;
    attribute::MeshAttributeHandle m_offset_vertex_distance_attribute;
    attribute::MeshAttributeHandle m_offset_vertex_target_distance_attribute;
    attribute::MeshAttributeHandle m_offset_vertex_converged_attribute;
    attribute::MeshAttributeHandle m_offset_face_converged_attribute;

    // vertex color attributes for the coloring-based parallel scheduler
    attribute::MeshAttributeHandle m_embedding_color_handle;
    attribute::MeshAttributeHandle m_offset_color_handle;
    attribute::MeshAttributeHandle m_embedding_vertex_color_tmp;


    std::shared_ptr<operations::SingleAttributeTransferStrategy<double, double>>
        m_embedding_edge_length_transfer;

    std::shared_ptr<operations::SingleAttributeTransferStrategy<double, double>>
        m_embedding_target_edge_length_transfer;

    std::shared_ptr<operations::SingleAttributeTransferStrategy<double, double>>
        m_embedding_amips_transfer;

    std::shared_ptr<operations::SingleAttributeTransferStrategy<double, double>>
        m_offset_position_to_embedding_transfer; // offset pos -> embedding pos

    std::shared_ptr<operations::SingleAttributeTransferStrategy<double, double>>
        m_offset_point_to_face_transfer; // offset pos -> normal samples, normal deviation, mean
                                         // ratio metric

    std::shared_ptr<operations::SingleAttributeTransferStrategy<double, double>>
        m_offset_point_to_face_transfer_wout_convergence; // offset pos -> normal samples, normal
                                                          // deviation, mean ratio metric

    std::shared_ptr<operations::SingleAttributeTransferStrategy<double, double>>
        m_offset_edge_length_transfer; // pos -> edge length

    // std::shared_ptr<operations::SingleAttributeTransferStrategy<double, double>>
    //    m_offset_face_normal_deviation_transfer; // pos -> normal deviation
    //
    // std::shared_ptr<operations::SingleAttributeTransferStrategy<double, double>>
    //    m_offset_mean_ratio_metric_transfer; // pos -> mrm

    std::shared_ptr<operations::SingleAttributeTransferStrategy<double, double>>
        m_offset_target_edge_length_transfer; // edge quadric distance -> target edge length

    // priorities
    std::function<double(const simplex::Simplex&)> m_embedding_prio_long_edges_first;
    std::function<double(const simplex::Simplex&)> m_embedding_prio_short_edges_first;

    std::function<double(const simplex::Simplex&)> m_offset_prio_long_edges_first;
    std::function<double(const simplex::Simplex&)> m_offset_prio_short_edges_first;

    // invariants
    std::shared_ptr<wmtk::invariants::InvariantCollection> m_link_conditions;
    std::shared_ptr<wmtk::invariants::Invariant> m_embedding_inversion_invariant;
    std::shared_ptr<wmtk::invariants::InvariantCollection> m_embedding_no_child_invariants;
    std::shared_ptr<wmtk::invariants::Invariant> m_embedding_separate_substructures_invariant;
    std::shared_ptr<wmtk::invariants::Invariant> m_embedding_todo_larger;
    std::shared_ptr<wmtk::invariants::Invariant> m_embedding_todo_smaller;
    std::shared_ptr<wmtk::invariants::Invariant> m_embedding_roi_invariant;

    std::shared_ptr<wmtk::invariants::Invariant> m_embedding_valence_3_invariant;
    std::shared_ptr<wmtk::invariants::Invariant> m_embedding_valence_4_invariant;
    std::shared_ptr<wmtk::invariants::Invariant> m_embedding_swap44_invariant;
    std::shared_ptr<wmtk::invariants::Invariant> m_embedding_swap44_2_invariant;
    std::shared_ptr<wmtk::invariants::Invariant> m_embedding_swap32_invariant;
    std::shared_ptr<wmtk::invariants::Invariant> m_embedding_swap23_invariant;

    std::shared_ptr<wmtk::invariants::Invariant> m_offset_todo_larger;
    std::shared_ptr<wmtk::invariants::Invariant> m_offset_todo_smaller;
    std::shared_ptr<wmtk::invariants::Invariant> m_offset_swap_invariant;
    std::shared_ptr<wmtk::invariants::Invariant> m_normal_deviation_invariant;
    std::shared_ptr<wmtk::invariants::Invariant> m_offset_convergence_invariant;


    // energies
    std::shared_ptr<function::PerSimplexFunction> m_embedding_amips;
    std::shared_ptr<function::LocalNeighborsSumFunction> m_amips_energy;


    // attribute strategies
    std::shared_ptr<operations::CollapseNewAttributeStrategy<double>>
        m_embedding_pos_collapse_strat;
    std::shared_ptr<operations::CollapseNewAttributeStrategy<double>>
        m_embedding_pos_collapse_strat_1;
    std::shared_ptr<operations::CollapseNewAttributeStrategy<double>>
        m_embedding_pos_collapse_strat_2;

    std::shared_ptr<operations::CollapseNewAttributeStrategy<double>> m_offset_pos_collapse_strat_1;
    std::shared_ptr<operations::CollapseNewAttributeStrategy<double>> m_offset_pos_collapse_strat_2;
    std::shared_ptr<operations::CollapseNewAttributeStrategy<int64_t>>
        m_offset_vertex_converged_collapse_strat_1;
    std::shared_ptr<operations::CollapseNewAttributeStrategy<int64_t>>
        m_offset_vertex_converged_collapse_strat_2;

    // offset smoothing function
    std::function<bool(Mesh&, const simplex::Simplex& s)> m_offset_smoothing_function;
};

} // namespace wmtk::components::internal
