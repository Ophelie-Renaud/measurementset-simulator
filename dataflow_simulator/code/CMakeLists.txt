cmake_minimum_required(VERSION 3.30)

# Add definition for relative path into project
#add_definitions( -DPROJECT_ROOT_PATH="${CMAKE_CURRENT_SOURCE_DIR}")

project(STD_DEGRID C)

set(CMAKE_C_STANDARD 11)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/lib/cmake_modules" ${CMAKE_MODULE_PATH})
find_package(FFTW REQUIRED)

if(NOT FFTW_FOUND)
    MESSAGE(FATAL_ERROR "FFTW not found !")
endif()

add_executable(STD_DEGRID main.c top.c)
target_link_libraries(STD_DEGRID m fftw3 fftw3f)
