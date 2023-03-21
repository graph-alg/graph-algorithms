cmake_minimum_required(VERSION 3.10)

project(scnu)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(HOME_DIR ${PROJECT_SOURCE_DIR})

set(LIBRARY_OUTPUT_PATH ${HOME_DIR}/lib)

add_subdirectory(src)
add_subdirectory(test)





