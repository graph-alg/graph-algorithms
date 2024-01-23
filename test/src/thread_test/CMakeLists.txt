project(thread_test)

include_directories(/usr/local/include ${HOME_DIR}/include ${HOME_DIR}/test/include)

link_directories(${HOME_DIR}/lib)

link_libraries(thread pthread random)

file(GLOB SRC_LIST ./src/*.cpp)

add_executable(thread_test ${SRC_LIST})

