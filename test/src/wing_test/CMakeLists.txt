project(wing_test)

include_directories(${HOME_DIR}/include ${HOME_DIR}/test/include)

link_directories(${HOME_DIR}/lib)

link_libraries(graph io logger pthread random stdc++fs string thread time wing system)

add_executable(wing_decomposition_test src/wing_decomposition_test.cpp)

add_executable(wing_insertion_test src/wing_insertion_test.cpp)

add_executable(wing_removal_test src/wing_removal_test.cpp)

add_executable(wing_insertion_removal_test src/wing_insertion_removal_test.cpp)

add_executable(wing_insertion_compare_test src/wing_insertion_compare_test.cpp)

add_executable(wing_removal_compare_test src/wing_removal_compare_test.cpp)

add_executable(wing_insertion_size_test src/wing_insertion_size_test.cpp)

add_executable(wing_removal_size_test src/wing_removal_size_test.cpp)

add_executable(wing_insertion_thread_test src/wing_insertion_thread_test.cpp)

add_executable(wing_removal_thread_test src/wing_removal_thread_test.cpp)

add_executable(wing_shuffle_non_temporal_bipartite_edges src/wing_shuffle_non_temporal_bipartite_edges.cpp)

