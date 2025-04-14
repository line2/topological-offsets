set(COMPONENT_NAME topological_offsets)
add_component(${COMPONENT_NAME})

if(NOT ${WMTK_ENABLE_COMPONENT_${COMPONENT_NAME}})
    return()
endif()


set(SRC_FILES
    topological_offsets.hpp
    topological_offsets.cpp
    
    TopologicalOffsetsOptions.hpp
    TopologicalOffsetsOutput.hpp

    internal/OffsetOptimization.hpp
    internal/OffsetOptimization.cpp
    internal/TopologicalOffsetGenerator.hpp
    internal/TopologicalOffsetGenerator.cpp

    internal/utils/compute_tet_amips.hpp
    internal/utils/compute_tet_amips.cpp
    internal/utils/compute_edge_length.hpp
    internal/utils/compute_edge_length.cpp
    internal/utils/Quadrics.hpp
    internal/utils/Quadrics.cpp
    internal/utils/normal_angle.hpp
    internal/utils/mean_ratio_metric.hpp
    internal/utils/is_tet_in_offset.hpp
    internal/utils/is_tet_in_offset.cpp
    internal/utils/ConvergenceRate.hpp
    internal/utils/ConvergenceRate.cpp
    internal/utils/write_mesh.hpp
    internal/utils/write_mesh.cpp
    internal/utils/bvh_from_mesh.hpp
    internal/utils/bvh_from_mesh.cpp
    internal/utils/enclosing_sphere.hpp
    internal/utils/enclosing_sphere.cpp
    internal/utils/tet_conains_connected_component.hpp
    internal/utils/tet_conains_connected_component.cpp
    internal/utils/ray_mesh_intersect.hpp
    internal/utils/ray_mesh_intersect.cpp
    internal/utils/bvh_to_mesh_index_map.hpp
    internal/utils/bvh_to_mesh_index_map.cpp
    internal/utils/OffsetBvh.hpp
    internal/utils/OffsetBvh.cpp
    internal/utils/apply_distance_function.hpp
    internal/utils/apply_distance_function.cpp
    
    internal/invariants/OffsetSwapInvariant.hpp
    internal/invariants/OffsetSwapInvariant.cpp
    internal/invariants/LongestEdgeInvariant.hpp
    internal/invariants/LongestEdgeInvariant.cpp
    internal/invariants/NormalDeviationAfterInvariant.hpp
    internal/invariants/NormalDeviationAfterInvariant.cpp
    internal/invariants/MeanRatioAfterInvariant.hpp
    internal/invariants/MeanRatioAfterInvariant.cpp
    internal/invariants/OffsetCollapseBeforeInvariant.hpp
    internal/invariants/OffsetCollapseBeforeInvariant.cpp
    internal/invariants/AmipsBeforeInvariant.hpp
    internal/invariants/AmipsBeforeInvariant.cpp
    internal/invariants/ConvergenceInvariant.hpp
    internal/invariants/ConvergenceInvariant.cpp
    internal/invariants/RoiInvariant.hpp
    internal/invariants/RoiInvariant.cpp
    internal/invariants/NoChildSimplexInvariant.hpp
    internal/invariants/NoChildSimplexInvariant.cpp
    
    internal/functions/AmipsOffsetOptimization.hpp
    internal/functions/AmipsOffsetOptimization.cpp
    internal/functions/AmipsOptimization.hpp
    internal/functions/AmipsOptimization.cpp
    internal/functions/AmipsOffsetProblem.hpp
    internal/functions/AmipsOffsetProblem.cpp
    internal/functions/AmipsProblem.hpp
    internal/functions/AmipsProblem.cpp
    )

include(libigl)
include(tinyexpr)

    
target_sources(wmtk_${COMPONENT_NAME} PRIVATE ${SRC_FILES})

target_link_libraries(wmtk_${COMPONENT_NAME} PRIVATE 
    wmtk::regular_space
    wmtk::marching
    wmtk::multimesh_from_tag
)

target_link_libraries(wmtk_${COMPONENT_NAME} PUBLIC 
    igl::core
    tinyexpr::tinyexpr
)