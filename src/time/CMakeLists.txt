include_directories(${HOME_DIR}/include)

file(GLOB SRC_LIST src/*.cpp)

add_library(time SHARED ${SRC_LIST})