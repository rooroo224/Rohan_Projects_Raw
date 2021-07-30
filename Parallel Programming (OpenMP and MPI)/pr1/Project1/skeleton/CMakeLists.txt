cmake_minimum_required(VERSION 2.8)
PROJECT(2d_Unsteady)

add_compile_options(-g -fast)

set(VTK_DIR "/usr/lib64/vtk")
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
include_directories(.)
file(GLOB SOURCES "*.cpp")

add_executable(2d_Unsteady ${SOURCES})

target_link_libraries(2d_Unsteady ${VTK_LIBRARIES})
