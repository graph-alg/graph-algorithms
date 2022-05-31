project(logger)

include_directories(${HOME_DIR}/include/)

file(GLOB SRC_LIST src/*.cpp)

add_library(logger SHARED ${SRC_LIST})