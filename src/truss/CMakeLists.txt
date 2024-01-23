project(truss)

include_directories(${HOME_DIR}/include)

file(GLOB SRC_LIST src/*.cpp)

set(LIBRARY_OUTPUT_PATH ${HOME_DIR}/lib)

link_directories(${HOME_DIR}/lib /usr/lib64)

link_libraries(time)

add_library(truss SHARED ${SRC_LIST} src/quasi_truss_maintenance.cpp ../../include/truss/quasi_truss_maintenance.h)