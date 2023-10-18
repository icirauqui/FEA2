cmake_minimum_required(VERSION 3.5)

project(nlo LANGUAGES CXX)

set(CMAKE_BUILD_TYPE Release)

add_library(nlo
    levenberg_marquardt.cpp
    levenberg_marquardt.hpp
)

target_link_libraries(nlo PRIVATE 
    fea
)













