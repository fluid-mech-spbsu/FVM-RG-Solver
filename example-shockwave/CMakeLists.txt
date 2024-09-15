cmake_minimum_required(VERSION 3.14)

set(PROJECT_MAIN "${PROJECT_NAME}_main")
project(${PROJECT_MAIN} LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_executable("Shock-wave-ch4-solver" shockwave-ch4-solver-main.cpp)
add_executable("Shock-wave-ar-solver" shockwave-ar-solver-main.cpp)

target_link_libraries("Shock-wave-ch4-solver" FVM-RG-Solver)
target_link_libraries("Shock-wave-ar-solver" FVM-RG-Solver)