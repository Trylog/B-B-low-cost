cmake_minimum_required(VERSION 3.26)
project(TSP_BB_2)

set(CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc")
set(CMAKE_EXE_LINKER_FLAGS "-static")

add_executable(TSP_BB_2 main.cpp
        GraphMatrix.cpp)

