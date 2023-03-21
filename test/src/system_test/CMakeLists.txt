project(system_test)

include_directories(${HOME_DIR}/include ${HOME_DIR}/test/include)

link_directories(${HOME_DIR}/lib /usr/lib64)

link_libraries(system logger string random)

add_executable(system_test src/system_test.cpp)

