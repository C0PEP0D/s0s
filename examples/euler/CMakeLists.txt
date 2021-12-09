cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(solve)

# Build target
add_executable(solve main.cpp func.h)
# Link external Libraries
target_link_libraries(solve tbb)
# Include
target_include_directories(solve PRIVATE .)
# Link s0s
target_include_directories(solve PRIVATE ../../include)
target_include_directories(solve PRIVATE ../../../../thirdparty/eigen)
