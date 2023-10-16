add_library(fdm3dsolverlib fdm3dsolver.cpp)
set_target_properties(fdm3dsolverlib PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(fdm3dsolverlib ${Boost_LIBRARIES} Eigen3::Eigen)