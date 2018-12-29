cmake_minimum_required(VERSION 3.13)
project(MD_Simulations)

set(CMAKE_CXX_STANDARD 17)

add_executable(MD_Simulations main.cpp Simulation_Class.cpp Simulation_Class.h)