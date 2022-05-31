project(random)

include_directories(${HOME_DIR}/include/)

file(GLOB SRC_LIST src/*.cpp)

add_library(random SHARED ${SRC_LIST})