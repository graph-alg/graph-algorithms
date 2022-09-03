project(io)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

include_directories(${HOME_DIR}/include/)

link_directories(${HOME_DIR}/lib)

link_libraries(graph)

file(GLOB SRC_LIST src/*.cpp)

add_library(io SHARED ${SRC_LIST})