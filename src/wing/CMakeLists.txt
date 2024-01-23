project(wing)

include_directories(../../include/)

file(GLOB SRC_LIST src/*.cpp)

set(LIBRARY_OUTPUT_PATH ../../lib)

add_library(wing SHARED ${SRC_LIST} src/bipartite_edge_index.cpp ../../include/wing/bipartite_edge_index.h)