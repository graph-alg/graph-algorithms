project(thread)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

include_directories(${HOME_DIR}/include/)

file(GLOB SRC_LIST src/*.cpp)

add_library(thread SHARED ${SRC_LIST})