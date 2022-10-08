cmake_minimum_required(VERSION 3.20)

set(PROJECT_NAME ELF)
project(${PROJECT_NAME})

if (UNIX)
    message(STATUS "[UNIX-like os]")
    set(CMAKE_CXX_FLAGS "-fopenmp")
elseif (WIN32)
    message(STATUS "[Windows]")
    find_package(OpenMP REQUIRED)
    if (OPENMP_FOUND)
        message("OPENMP FOUND")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif ()
endif ()


find_package(Eigen3 CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)

# 添加源文件目录
include_directories(${PROJECT_SOURCE_DIR}/include)
aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_DIR)

# ELF PIC
# add_executable(ELF ${SRC_DIR} example/ELF.cpp)
# target_link_libraries(ELF PRIVATE Eigen3::Eigen spdlog::spdlog spdlog::spdlog_header_only)

# 旋转电磁场
add_executable(FDTD ${SRC_DIR} example/FDTD.cpp)
target_link_libraries(FDTD PRIVATE Eigen3::Eigen)
target_link_libraries(FDTD PRIVATE spdlog::spdlog spdlog::spdlog_header_only)

# 测试程序
add_executable(TEST test/decorator.cpp)
target_link_libraries(TEST PRIVATE Eigen3::Eigen)
target_link_libraries(TEST PRIVATE spdlog::spdlog spdlog::spdlog_header_only)