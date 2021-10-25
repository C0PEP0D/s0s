cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(s0s)

if (APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
endif()

# Create C++ library
add_library(s0s include/s0s/solver.h)
target_include_directories(s0s PUBLIC include)
set_target_properties(s0s PROPERTIES LINKER_LANGUAGE CXX)
