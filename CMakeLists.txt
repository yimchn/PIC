cmake_minimum_required(VERSION 3.20)

set(PROJECT_NAME ELF)
project(${PROJECT_NAME})

# set(CMAKE_CXX_FLAGS "-fopenmp")

# 添加源文件目录
include_directories(${PROJECT_SOURCE_DIR}/include)
aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_DIR)

# ELF PIC
add_executable(ELF ${SRC_DIR} example/ELF.cpp)

# 旋转电磁场
add_executable(EMF ${SRC_DIR} example/EMF.cpp)

# 测试程序
add_executable(TEST test/test.cpp)