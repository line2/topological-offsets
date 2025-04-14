#pragma once

#include <optional>
#include <wmtk/Mesh.hpp>
#include <wmtk/attribute/MeshAttributeHandle.hpp>

namespace wmtk::components {

struct TopologicalOffsetsOptions
{
public:
    /**
     * This attribute holds the positions of the embedding mesh. It also contains a pointer to the
     * embedding mesh itself.
     */
    attribute::MeshAttributeHandle embedding_position_handle;
    /**
     * The mesh to which the offset should be generated.
     */
    std::shared_ptr<Mesh> input_mesh;
    /**
     * Every simplex in the "inside" mesh is tagged as input.
     * For example, if you want to generate the offset of a tet mesh, the inside mesh are the tets
     * and the input mesh the boundary (triangle mesh) of the tet mesh.
     */
    std::shared_ptr<Mesh> inside_mesh;

    /**
     * The (absolute) offset distance.
     */
    double distance;

    /**
     * The min/max normal deviation affects the "roundness" of the offset in convex regions. The
     * smaller the max_normal deviation, the smoother convex regions will be. However, this also
     * increases the mesh complexity.
     */
    double min_normal_deviation = 2;
    double max_normal_deviation = 15;

    /**
     * The min/max edge length are not hard constraints but hints to the optimization when to split
     * or collapse an edge. Too long edges are split, too short ones are collapsed.
     *
     * If no min edge length is specified, it is computed based on the `min_normal_deviation`.
     *
     * If no max edge length is specified, it is set to infinity. As there is a sizing field that
     * creates a smooth transition from short to long edges, the max edge length is usually
     * unnecessary to specify.
     */
    std::optional<double> min_edge_length;
    std::optional<double> max_edge_length;

    /**
     * The maximum average relative distance allowed for a mesh to be converged.
     */
    double convergence_max = 0.5;

    /**
     * @brief Offset distance function.
     *
     * It is possible to initialize the sizing field by defining a volumetric distance function,
     * e.g., x*y*z. Available variables: x,y,z, and d as distance.
     */
    std::string distance_function;

    /**
     * The number of iterations for the entire optimization.
     */
    int64_t passes = 10;
    /**
     * The number of iterations for the embedding optimization within a single pass of the entire
     * optimization.
     */
    int64_t embedding_passes = 2;
    /**
     * The number of iterations for the offset optimization within a single pass of the entire
     * optimization.
     */
    int64_t offset_passes = 2;

    /**
     * Compute a finite offset. By default, this component computes the topological offset instead.
     */
    bool finite_offset = false;
    /**
     * Multiply the distance and min/max edge length by the AABB of the input mesh.
     */
    bool relative_distance_and_length = true;
    /**
     * Print intermediate output of the optimization.
     */
    bool intermediate_output = false;
    /**
     * Use simplicial embedding. This option only exists for testing and should always be true, at
     * least if a topological offset is desired.
     */
    bool use_simplicial_embedding = true;
    /**
     * For the topological offset, add as much tets as possible to the input that are within the
     * offset distance. The input topology will be preserved, so using the warm start just improves
     * performance but does not affect the correctness of the topology.
     *
     * This option is ignored if `finite_offset` is true.
     */
    bool use_warm_start = true;

    bool expand_topo = true;
    /**
     * In case of the warm start, locally adapt the offset distance to avoid colliding and
     * potentially degenerating offsets.
     */
    bool adapt_offset_distance = true;
    /**
     * This option restricts the min edge length to the average edge length of the input. This can
     * be advantageous for very small offset distances, where the computed min edge length would
     * become extremely small.
     */
    bool restrict_min_edge_length_to_input_avg = false;

    bool DEBUG_print_embedding = false;
    bool DEBUG_print_offset = false;
    bool DEBUG_print_smooth = false;
};

} // namespace wmtk::components
