cmake_minimum_required (VERSION 3.27)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_COMPILE g++)

include_directories("include")
find_package(Boost COMPONENTS bimap)
include_directories(${Boost_INCLUDE_DIR})
project(fdm3d)
add_library(eigen_target INTERFACE)
add_subdirectory("src")
add_subdirectory(thirdparty/pybind11)
add_subdirectory(thirdparty/eigen)
add_subdirectory(python)
# add_executable(main main.cpp)
# target_link_libraries(main fdm3dsolverlib)


find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
# find_package(pybind11 REQUIRED HINTS "${Python3_SITELIB}")