cmake_minimum_required(VERSION 2.8)
PROJECT(Homework-2)

# Require MPI for this project:
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
add_compile_options(-g)

include_directories(.)
file(GLOB SOURCES "hw2.cpp")

add_executable(hw2 ${SOURCES})

target_link_libraries(hw2 ${MPI_C_LIBRARIES})
