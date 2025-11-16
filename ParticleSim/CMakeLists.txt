cmake_minimum_required(VERSION 3.18)
project(ParticleSim LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# Output executables in bin/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include_directories(include)

file(GLOB SRC src/*.cpp)

add_executable(particlesim ${SRC})
# target_link_libraries(particlesim PRIVATE simcuda)
