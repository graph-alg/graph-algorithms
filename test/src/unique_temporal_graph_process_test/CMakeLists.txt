project(unique_temporal_graph_process_test)

include_directories(${HOME_DIR}/include ${HOME_DIR}/test/include)

link_directories(${HOME_DIR}/lib /usr/lib64)

link_libraries(graph io string pthread thread stdc++fs)

file(GLOB SRC_LIST src/*.cpp)

add_executable(unique_temporal_graph_process_test ${SRC_LIST})