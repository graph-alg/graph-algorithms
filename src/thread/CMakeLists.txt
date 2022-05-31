project(thread)

include_directories(${HOME_DIR}/include/)

file(GLOB SRC_LIST src/*.cpp)

add_library(thread SHARED ${SRC_LIST})