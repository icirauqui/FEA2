cmake_minimum_required(VERSION 3.5)

project(fem_analysis LANGUAGES CXX)

# Compile in debug mode
set(CMAKE_BUILD_TYPE Debug)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/fea)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/dataset)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/nlo)


add_executable(fem_analysis
    src/main.cpp
)

target_link_libraries(fem_analysis PRIVATE 
    fea
    dataset
    nlo
)