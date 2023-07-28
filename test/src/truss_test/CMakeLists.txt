project(truss_test)

include_directories(${HOME_DIR}/include ${HOME_DIR}/test/include)

link_directories(${HOME_DIR}/lib /usr/lib64)

link_libraries(time graph io logger string random pthread thread stdc++fs system truss)

add_executable(truss_decomposition_test src/truss_decomposition_test.cpp)

add_executable(truss_insertion_test src/truss_insertion_test.cpp)

add_executable(truss_removal_test src/truss_removal_test.cpp)

add_executable(truss_insertion_removal_test src/truss_insertion_removal_test.cpp)

add_executable(truss_insertion_compare_test src/truss_insertion_compare_test.cpp)

add_executable(truss_removal_compare_test src/truss_removal_compare_test.cpp)

add_executable(truss_insertion_removal_compare_test src/truss_insertion_removal_compare_test.cpp)

add_executable(truss_insertion_size_test src/truss_insertion_size_test.cpp)

add_executable(truss_removal_size_test src/truss_removal_size_test.cpp)

add_executable(truss_insertion_removal_size_test src/truss_insertion_removal_size_test.cpp)

add_executable(truss_insertion_thread_test src/truss_insertion_thread_test.cpp)

add_executable(truss_removal_thread_test src/truss_removal_thread_test.cpp)

add_executable(truss_insertion_removal_thread_test src/truss_insertion_removal_thread_test.cpp)

add_executable(truss_shuffle_non_temporal_edges src/truss_shuffle_non_temporal_edges.cpp)

add_executable(truss_split_non_temporal_edges src/truss_split_non_temporal_edges.cpp)