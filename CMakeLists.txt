cmake_minimum_required(VERSION 3.25.2)

project(Spacecraft_maneuver)

set(CMAKE_CXX_STANDARD 17)


add_subdirectory(tests)

find_package(Eigen3 REQUIRED)

target_include_directories(TESTs PRIVATE ${EIGEN3_INCLUDE_DIRS})

target_include_directories(TESTs PRIVATE src)

