# CMakeList.txt: 顶层 CMake 项目文件，在此处执行全局配置
# 并包含子项目。
#
cmake_minimum_required (VERSION 3.8)

include_directories("include")

project ("fdm3d")

# 包含子项目。
add_subdirectory ("src")
add_executable(main main.cpp)
target_link_libraries(main fdm3dsolverlib)
