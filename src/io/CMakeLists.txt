project(io)

include_directories(${HOME_DIR}/include/)

link_directories(${HOME_DIR}/lib)

link_libraries(graph)

file(GLOB SRC_LIST src/*.cpp)

add_library(io SHARED ${SRC_LIST})