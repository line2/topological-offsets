set(COMPONENT_NAME regular_space)
add_component(${COMPONENT_NAME})

if(NOT ${WMTK_ENABLE_COMPONENT_${COMPONENT_NAME}})
    return()
endif()

set(SRC_FILES
    RegularSpaceOptions.hpp
    internal/RegularSpace.hpp
    internal/RegularSpace.cpp
	regular_space.hpp
    regular_space.cpp
    )


target_sources(wmtk_${COMPONENT_NAME} PRIVATE ${SRC_FILES})