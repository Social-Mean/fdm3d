cmake_minimum_required(VERSION 3.27)

pybind11_add_module(fdm3d fdm3dsolver_wrapper.cpp)
target_link_libraries(fdm3d PUBLIC fdm3dsolverlib Eigen3::Eigen eigen_target)
