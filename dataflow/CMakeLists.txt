cmake_minimum_required(VERSION 3.30)
project(STD_DEGRID C)

set(CMAKE_C_STANDARD 11)

add_executable(STD_DEGRID main.c top.c)
target_link_libraries(STD_DEGRID m)
