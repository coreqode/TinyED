# Required by CMake
cmake_minimum_required (VERSION 3.9)

project(TinyED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_LIST_DIR}/cmake")

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}")

if (NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 REQUIRED)
endif()

message("Building TinyED in " ${CMAKE_BUILD_TYPE} " mode")

file(GLOB_RECURSE TINYED_HEADERS ${PROJECT_SOURCE_DIR}/include/*.h)
add_library(TinyED INTERFACE
        include/helpers.h)
      add_custom_target(TinyED-Headers SOURCES ${TINYED_HEADERS}) # IDE shows TinyED headers in one place only
target_include_directories(TinyED INTERFACE ${PROJECT_SOURCE_DIR}/include )
target_compile_options(TinyED INTERFACE $<$<CXX_COMPILER_ID:MSVC>:/bigobj>)
target_link_libraries(TinyED INTERFACE Eigen3::Eigen)

add_executable(TinyEDMain examples/main.cpp)
target_link_libraries(TinyEDMain PRIVATE TinyED)
