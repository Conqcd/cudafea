cmake_minimum_required(VERSION 3.10)

project(CudaCompute)

find_package(CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CUDA_NVCC_FLAGS -std=c++11)


file(GLOB CUDASOURCES ./*.cu)

cuda_add_library(core SHARED
        ${CUDASOURCES})