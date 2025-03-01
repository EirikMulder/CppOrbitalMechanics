cmake_minimum_required(VERSION 3.21)

# Set C++ Standard to c++20
set(CMAKE_CXX_STANDARD 20)

# Export compilation commands to clangd can identify added libraries.
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# Project Definition
project(CppOrbital VERSION 0.1)

add_subdirectory(matplotplusplus)

# Executables
add_executable(CppOrbital orbital_propagation.cpp)
target_link_libraries(CppOrbital PUBLIC matplot)
