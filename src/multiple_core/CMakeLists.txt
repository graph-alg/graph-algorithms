project(multiple_core)

include_directories(${HOME_DIR}/include)

file(GLOB SRC_LIST src/*.cpp)

link_directories(${HOME_DIR}/lib /usr/lib64)

add_library(multiple_core SHARED ${SRC_LIST})