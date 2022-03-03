cmake_minimum_required(VERSION 3.16)
project(Progetto_Geometria_Computazionale)

set(CMAKE_CXX_STANDARD 20)

include_directories(Progetto_Geometria_Computazionale libs/boost_1_77_0/)
include_directories(Progetto_Geometria_Computazionale libs/CGAL-5.4/include)

add_executable(Progetto_Geometria_Computazionale main.cpp)
