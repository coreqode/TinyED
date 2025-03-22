# Required by CMake
cmake_minimum_required (VERSION 3.5)

# Declare your project name, languages, version
project(
    tinyed
    LANGUAGES CXX
    VERSION 1.0
)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_LIST_DIR}/cmake")

set(CMAKE_CXX_STANDARD 17)
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}")

add_subdirectory(extern)
add_subdirectory(src)
add_subdirectory(examples)
