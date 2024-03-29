cmake_minimum_required(VERSION 3.5)

project(fea LANGUAGES CXX)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_VERBOSE_MAKEFILE OFF)


find_package(VTK REQUIRED)
set(VTKFIG_INCLUDE_DIRS "@CONF_INCLUDE_DIRS@;@VTK_INCLUDE_DIRS@;@VTK_PREFIX_PATH@/include/vtk-@VTK_MAJOR_VERSION@.@VTK_MINOR_VERSION@")

find_package(PCL 1.12 REQUIRED)

INCLUDE_DIRECTORIES(
    ${PCL_INCLUDE_DIRS}
    ${VTKFIG_INCLUDE_DIRS}
)

link_directories(
    ${PCL_LIBRARY_DIRS}
    ${PCL_LIBRARY_DIRS}/..
)

add_definitions(
    ${PCL_DEFINITIONS}
)




add_library(fea
    pos.cpp
    pos.hpp
    fem.cpp
    fem.hpp
    fea.cpp
    fea.hpp
)

target_include_directories(fea 
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(fea
  ${VTK_LIBRARIES}
  ${PCL_LIBRARIES}
)