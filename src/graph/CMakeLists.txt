project(graph)

include_directories(${HOME_DIR}/include/)

file(GLOB SRC_LIST src/*.cpp)

add_library(graph SHARED ${SRC_LIST})