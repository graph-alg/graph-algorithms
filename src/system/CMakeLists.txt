project(system)

include_directories(${HOME_DIR}/include)

file(GLOB SRC_LIST src/*.cpp)

add_library(system SHARED ${SRC_LIST})