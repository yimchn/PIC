cmake_minimum_required(VERSION 3.20)

set(PROJECT_NAME ELF)
# set(EMF EMF)
project(${PROJECT_NAME})

set(CMAKE_CXX_FLAGS "-fopenmp")

find_package(Eigen3 CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)

# 添加源文件目录
include_directories(${PROJECT_SOURCE_DIR}/include)
aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_DIR)

# ELF PIC
# add_executable(ELF ${SRC_DIR} example/ELF.cpp)
# target_link_libraries(ELF PRIVATE Eigen3::Eigen spdlog::spdlog spdlog::spdlog_header_only)

# 旋转电磁场
add_executable(FDTD ${SRC_DIR} example/EMF.cpp)
target_link_libraries(FDTD PRIVATE Eigen3::Eigen spdlog::spdlog spdlog::spdlog_header_only)

# 测试程序
add_executable(TEST test/test.cpp)
target_link_libraries(TEST PRIVATE Eigen3::Eigen spdlog::spdlog spdlog::spdlog_header_only)