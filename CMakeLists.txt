cmake_minimum_required(VERSION 3.15)

project(PT LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIED true)

add_executable(PT main.cpp)
add_executable(p3dif p3_ppm_comparer.cpp)