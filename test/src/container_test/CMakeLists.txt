project(container_test)

include_directories(${HOME_DIR}/include ${HOME_DIR}/test/include)

file(GLOB SRC_LIST src/*.cpp)

add_executable(container_test ${SRC_LIST})