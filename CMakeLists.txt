cmake_minimum_required(VERSION 3.16)
project(Progetto_Geometria_Computazionale)

set(CMAKE_CXX_STANDARD 20)

include_directories(Progetto_Geometria_Computazionale libs/boost/)
include_directories(Progetto_Geometria_Computazionale libs/cgal/)

add_executable(Progetto_Geometria_Computazionale main.cpp)
