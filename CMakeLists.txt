cmake_minimum_required(VERSION 3.16)
project(MP)

set(CMAKE_CXX_STANDARD 17)
set(EXECUTABLE_OUTPUT_PATH ".")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
find_package(Eigen3 3.3 REQUIRED)
find_package(OpenCV REQUIRED)



include_directories(${OpenCV_INCLUDE_DIRS})
include_directories( EIGEN3_INCLUDE_DIR )


add_executable(MP main.cpp
)

set(PACKAGE_STRING MP)
target_link_libraries(MP ${OpenCV_LIBS} )
target_link_libraries (MP Eigen3::Eigen)

