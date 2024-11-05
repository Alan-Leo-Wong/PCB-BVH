add_library(PCB-BVH INTERFACE)
find_package(Threads)
if (Threads_FOUND)
    # Link with the threading library of the system, which may
    # be required by standard header <thread> on some systems
    target_link_libraries(PCB-BVH INTERFACE Threads::Threads)
endif()

target_include_directories(PCB-BVH INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/PCB-BVH/src>
    $<INSTALL_INTERFACE:include>)

set_target_properties(PCB-BVH PROPERTIES CXX_STANDARD 20)

install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION include/PCB-BVH/v2
    FILES_MATCHING PATTERN "*.h"
    PATTERN "c_api" EXCLUDE)

if (BVH_BUILD_C_API)
    add_subdirectory(c_api)
endif()
