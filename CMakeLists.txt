cmake_minimum_required(VERSION 3.16)
project(MP)

set(CMAKE_CXX_STANDARD 17)
set(EXECUTABLE_OUTPUT_PATH ".")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
find_package(Eigen3 3.3 REQUIRED)
find_package(OpenCV REQUIRED)
find_package(OpenMP)
find_package(Boost REQUIRED COMPONENTS filesystem)


include_directories(${OpenCV_INCLUDE_DIRS})
include_directories( EIGEN3_INCLUDE_DIR )
include_directories(${Boost_INCLUDE_DIRS})

add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})
add_executable(MP main.cpp
)

set(PACKAGE_STRING MP)
target_link_libraries(MP ${OpenCV_LIBS} )
target_link_libraries (MP Eigen3::Eigen)
target_link_libraries(MP OpenMP::OpenMP_CXX)
target_link_libraries(MP  Boost::filesystem)