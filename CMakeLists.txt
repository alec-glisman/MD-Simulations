cmake_minimum_required(VERSION 3.13)
project(MD_Simulations)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_COMPILER /usr/local/bin/g++-8)
set(CMAKE_C_COMPILER /usr/local/bin/gcc-8)
# enable openmp using -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp")

set(Boost_USE_STATIC_LIBS ON) # only find static libs
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost 1.68.0)

if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(MD_Simulations main.cpp Simulation.cpp Simulation.h Utilities.cpp Utilities.h ProgressBar.hpp)
    target_link_libraries(MD_Simulations ${Boost_LIBRARIES})
endif ()

if (NOT Boost_FOUND)
    message(FATAL_ERROR "Could not find boost!")
endif ()

# Compiler Optimization
set(CMAKE_CXX_FLAGS "-O3 -mavx2")  # Better for my MacBook
#set(CMAKE_CXX_FLAGS "-O3 -mavx2 -march=native")  # Better for distributed cluster computing


# For profiling code
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
#set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -O0")
