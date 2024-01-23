project(logger_test)

include_directories(/usr/local/include ${HOME_DIR}/include)

link_directories(${HOME_DIR}/lib)

link_libraries(logger)

file(GLOB SRC_LIST src/*.cpp)

add_executable(logger_test ${SRC_LIST})
