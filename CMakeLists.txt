cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(gRKPM LANGUAGES CXX CUDA)

enable_language(CUDA)
include(CheckLanguage)
check_language(CUDA)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14") # 添加c++11标准支持

# SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb -fopenmp")  # debug模式下 gdb相关选项
# SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall -fopenmp") # release模式下 gdb相关选项
SET(CMAKE_BUILD_TYPE "Debug")
# SET(CMAKE_BUILD_TYPE "Release")

if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    message("Debug mode")
    add_definitions(-DDEBUG) # 相当于代码中 #define DEBUG
else()
    message("Release mode")
endif ()

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin") # 可执行文件输出目录

INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/include") # 头文件包含目录
INCLUDE_DIRECTORIES("~/lib_local/seacas/include") # 头文件包含目录


AUX_SOURCE_DIRECTORY(src DIR_SRCS) # 添加源代码文件夹, 自动扫描所有文件

if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD 14)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

# message("${CUDA_RUNTIME_LIBRARY}")
link_directories(~/lib_local/seacas/lib)
add_executable(gRKPM MACOSX_BUNDLE ${DIR_SRCS})
target_link_libraries(gRKPM ${CUDA_RUNTIME_LIBRARY})
target_link_libraries(${PROJECT_NAME} libexodus.so.2)