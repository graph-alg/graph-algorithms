project(container_test)

include_directories(${HOME_DIR}/include ${HOME_DIR}/test/include)

file(GLOB SRC_LIST src/*.cpp)

include_directories(${HOME_DIR}/include ${HOME_DIR}/test/include)

link_directories(${HOME_DIR}/lib /usr/lib64)

link_libraries(time logger string random pthread thread stdc++fs)

add_executable(container_test ${SRC_LIST})