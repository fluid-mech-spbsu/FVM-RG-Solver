cmake_minimum_required(VERSION 3.14)

set(PROJECT_MAIN "${PROJECT_NAME}_main")
project(${PROJECT_MAIN} LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_executable("Sods-solver" sods-solver-main.cpp)
target_link_libraries("Sods-solver" FVM-RG-Solver)

