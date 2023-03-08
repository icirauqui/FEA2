cmake_minimum_required(VERSION 3.5)

project(fe LANGUAGES CXX)

find_package(VTK REQUIRED)
set(VTKFIG_INCLUDE_DIRS "@CONF_INCLUDE_DIRS@;@VTK_INCLUDE_DIRS@;@VTK_PREFIX_PATH@/include/vtk-@VTK_MAJOR_VERSION@.@VTK_MINOR_VERSION@")

find_package(PCL REQUIRED)

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


add_library(fe
    fe.cpp
    fe.hpp
)


target_link_libraries(fe
  ${VTK_LIBRARIES}
)