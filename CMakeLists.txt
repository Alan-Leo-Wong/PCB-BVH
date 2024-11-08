option(BVH_BUILD_C_API "Builds the C API library wrapper" OFF)
option(BVH_C_API_STATIC_LINK_STDLIB "Link the C API library statically against the standard C++ library (only supported by clang/gcc)" OFF)
option(BVH_C_API_UNSAFE_CASTS "Enable unsafe casts in C API" OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}PCB-BVH/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}PCB-BVH/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}PCB-BVH/lib)

add_subdirectory(src/bvh/v2)

if (PROJECT_IS_TOP_LEVEL)
    include(CTest)
    if (BUILD_TESTING)
        add_subdirectory(test)
    endif()

    include(cmake/Install.cmake)
endif()
