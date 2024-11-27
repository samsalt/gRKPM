cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(gRKPM LANGUAGES CXX CUDA)

enable_language(CUDA)
include(CheckLanguage)
check_language(CUDA)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14") 

# SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb -fopenmp")  
# SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall -fopenmp") 
SET(CMAKE_BUILD_TYPE "Debug")
# SET(CMAKE_BUILD_TYPE "Release")

if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    message("Debug mode")
    add_definitions(-DDEBUG) 
else()
    message("Release mode")
endif ()

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin")

INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/include") 
INCLUDE_DIRECTORIES("~/seacas/include") 

AUX_SOURCE_DIRECTORY(src DIR_SRCS) 

if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD 14)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

# message("${CUDA_RUNTIME_LIBRARY}")
# in order to compile, it needs the seacas lib
link_directories(~/seacas/lib)  
add_executable(gRKPM MACOSX_BUNDLE ${DIR_SRCS})
target_link_libraries(gRKPM ${CUDA_RUNTIME_LIBRARY})
target_link_libraries(${PROJECT_NAME} libexodus.so.2)