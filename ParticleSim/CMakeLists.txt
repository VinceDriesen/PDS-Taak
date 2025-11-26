cmake_minimum_required(VERSION 3.18)

if(NOT CMAKE_CUDA_COMPILER)
    set(CMAKE_CUDA_COMPILER nvcc)
endif()

project(ParticleSim LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES native)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# Output executables in bin/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include_directories(include)

file(GLOB SRC src/*.cpp src/*.cu)

find_package(CUDAToolkit REQUIRED)

add_executable(particlesim ${SRC})

target_link_libraries(particlesim PRIVATE CUDA::cudart)

set_target_properties(particlesim PROPERTIES CUDA_SEPARABLE_COMPILATION ON)