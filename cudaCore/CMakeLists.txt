cmake_minimum_required(VERSION 3.10)

project(CudaCompute)

find_package(CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CUDA_NVCC_FLAGS -std=c++17)

include_directories(${PROJECT_SOURCE_DIR}$)

file(GLOB CUDASOURCES ./*.cu)

cuda_add_library(core STATIC
        ${CUDASOURCES})