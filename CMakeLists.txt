cmake_minimum_required(VERSION 3.5)

project(fem_analysis LANGUAGES CXX)

add_executable(fem_analysis
    src/data/data.cpp
    src/data/data.hpp
    src/data/aux.cpp
    src/data/aux.hpp
    src/fea/fea2.cpp
    src/fea/fea2.hpp
    src/fea/fea.cpp
    src/fea/fea.hpp
    src/main.cpp
)

target_link_libraries(fem_analysis PRIVATE)