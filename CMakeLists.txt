cmake_minimum_required(VERSION 3.13)
project(MD_Simulations)

set(CMAKE_CXX_STANDARD 17)

# add_executable(MD_Simulations main.cpp Simulation_Class.cpp Simulation_Class.h)


set(Boost_USE_STATIC_LIBS ON) # only find static libs
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost 1.68.0)
# COMPONENTS date_time filesystem system range

if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(MD_Simulations main.cpp Simulation_Class.cpp Simulation_Class.h)
    target_link_libraries(MD_Simulations ${Boost_LIBRARIES})
endif ()

if (NOT Boost_FOUND)
    message(FATAL_ERROR "Could not find boost!")
endif ()