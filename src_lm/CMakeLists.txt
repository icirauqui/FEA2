cmake_minimum_required(VERSION 3.5)

project(fem_analysis LANGUAGES CXX)

find_package(MRPT 2.0 COMPONENTS poses gui OPTIONAL_COMPONENTS vision)

# Compile in debug mode
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_VERBOSE_MAKEFILE OFF)




add_executable(fem_analysis
    src/main.cpp
)

target_link_libraries(fem_analysis
${MRPT_LIBRARIES})